#!/usr/bin/env python3
"""
Calibrate a machine and its msearch build for the paid search: how many
reference CPU-years one instance-hour buys, and how many workers to run.
Linux or macOS, Python 3.8+, standard library only.

usage: machine_cal.py [--bin DIR/msearch] [--out cal.json] [--ks 1,8,16]
                      [--no-bench] [--refs refs.json] [--quiet]

1. Detects the CPU: its ISA (/proc/cpuinfo flags, or sysctl on macOS), the
   physical (on macOS the performance) cores, the logical CPUs and SMT.
2. Exactness gate: bench/quick.txt and bench/full.txt must report every
   instance ok (squares and hashes) and the node totals of the build's
   search path ("carry512": AVX-512BW, or "matrix": the portable path);
   msearch writes its path into its first record. Any mismatch exits 1.
3. Reference workload, the four fixed sums of the ISA study
   (research/ideas.md, "Machine calibration"): T1 plain 12 6 3 2 1 0 1 / 900,
   T2 plain with --r1-stride 4, 13 7 4 3 1 1 / 1900, T3 d-first --d-stride
   16, 12 6 3 2 1 1 / 988, T4 d-first --d-stride 64, 12 6 3 2 1 1 / 1200.
   Their nodes (and the squares of T1, T2) must equal the path's reference.
4. Runs each sum as K concurrent copies, for K = 1, the physical (or
   performance) cores and the logical CPUs: the throughput of the instance
   is K x the reference CPU (the CPU of the sum on the reference machine,
   the fast x86 build) / the mean wall time of a copy, in reference
   CPU-seconds per instance-second, for the plain search (T1 + T2) and the
   d-first search (T3 + T4).
5. Writes cal.json: the ISA, the path, the plain and d-first speeds at the
   best K (reference CPU-hours per instance-hour, and reference CPU-years
   per instance-hour), the recommended worker count (that K), the
   per-process speeds (reference CPU per process CPU second: what
   scheduler.py run --machine scales the records' CPU by) and the
   measured CPU per process against the scheduler's time law, alone and
   under load. scheduler.py forecast --machine cal.json --instance-hours H
   turns it into E and P(>=1 magic square) for H paid instance-hours.

Run it under the conditions of the search (nothing else heavy running).
About 1.5 minutes on a fast x86 instance (each K pass ~25 s of wall); the
portable path is ~4-12x slower per sum.
"""
import argparse
import hashlib
import json
import os
import platform
import re
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
YEAR_H = 8766.0

# The reference workload. ref_cpu: the CPU seconds of the sum ("cpu" of its
# record, the min of two paired runs) of the fast build (-march=native,
# carry512, gcc 13, Sapphire Rapids 2.1 GHz KVM, the development machine:
# research/ideas.md "Machine calibration"); law: the scheduler's shipped
# engine-4 time law for the whole sum at the model's N' and L_raw (plain for
# T1, T2; d-first for T3, T4), which "est_time" (T2-T4: the record's
# estimate of the whole sum's CPU from its sample) or "cpu" (T1) is
# compared to.
TESTS = [
    {"name": "T1", "mode": "plain", "args": ["--sums", "900", "12", "6", "3", "2", "1", "0", "1"],
     "type": "sum", "ref_cpu": 5.40, "law": 5.824, "squares": 18},
    {"name": "T2", "mode": "plain",
     "args": ["--r1-stride", "4", "--sums", "1900", "13", "7", "4", "3", "1", "1"],
     "type": "csum", "ref_cpu": 5.16, "law": 19.887, "squares": 18},
    {"name": "T3", "mode": "dfirst",
     "args": ["--diag-first", "--diag-first-min-n", "0", "--d-stride", "16",
              "--sums", "988", "12", "6", "3", "2", "1", "1"],
     "type": "dsum", "ref_cpu": 6.13, "law": 138.767, "pairs": 0},
    {"name": "T4", "mode": "dfirst",
     "args": ["--diag-first", "--diag-first-min-n", "0", "--d-stride", "64",
              "--sums", "1200", "12", "6", "3", "2", "1", "1"],
     "type": "dsum", "ref_cpu": 8.48, "law": 827.111, "pairs": 0},
]
# the nodes per search path: bench totals and T1-T4 (the search is exact and
# deterministic per path; the ISA changes its speed only)
REFS = {
    "carry512": {"quick": 1770779, "full": 14958507, "prod": 50375738,
                 "T1": 19762185, "T2": 21868831, "T3": 11012428, "T4": 9445350},
    "matrix": {"quick": 2804467, "full": 35176246, "prod": 171112789,
               "T1": 77187256, "T2": 83156821, "T3": 48075085, "T4": 56363950},
}
# a tiny sum that writes a record (bench/quick.txt: 2 squares), for the path
PROBE = ["--sums", "391", "10", "6", "3", "1", "0", "1"]
# the AVX-512 extensions the fast path uses beyond AVX-512BW (without them:
# the "cascadelake" class, ~1.05x the CPU; research/ideas.md)
FAST_EXT = ("avx512vpopcntdq", "avx512vbmi", "gfni", "avx512bitalg")
# canonical ISA names (as msearch's "isa") from /proc/cpuinfo flags
CPUINFO_FLAGS = {"avx2": "avx2", "bmi2": "bmi2", "avx512f": "avx512f", "avx512bw": "avx512bw",
                 "avx512_vpopcntdq": "avx512vpopcntdq", "avx512vbmi": "avx512vbmi",
                 "gfni": "gfni", "avx512_bitalg": "avx512bitalg"}
# ... and from macOS sysctl (hw.optional.* = 1, machdep.cpu.*features words)
SYSCTL_OPTIONAL = {"hw.optional.avx2_0": "avx2", "hw.optional.bmi2": "bmi2",
                   "hw.optional.avx512f": "avx512f", "hw.optional.avx512bw": "avx512bw",
                   "hw.optional.arm64": "arm64"}
SYSCTL_WORDS = {"AVX2": "avx2", "BMI2": "bmi2", "AVX512F": "avx512f", "AVX512BW": "avx512bw",
                "AVX512VPOPCNTDQ": "avx512vpopcntdq", "VPOPCNTDQ": "avx512vpopcntdq",
                "AVX512VBMI": "avx512vbmi", "GFNI": "gfni", "AVX512BITALG": "avx512bitalg",
                "BITALG": "avx512bitalg"}


# ---------------------------------------------------------------- detection

def parse_cpuinfo(text):
    """{"arch", "vendor", "model", "isa": sorted canonical names, "logical",
    "physical" (distinct (physical id, core id), or None), "smt_flag"} from
    the text of /proc/cpuinfo (x86 or ARM)"""
    procs = [b for b in re.split(r"\n\s*\n", text) if b.strip()]
    flags, cores = set(), set()
    vendor = model = None
    logical = 0
    arm = False
    siblings = cpu_cores = None
    for b in procs:
        kv = {}
        for line in b.splitlines():
            if ":" in line:
                k, v = line.split(":", 1)
                kv[k.strip().lower()] = v.strip()
        if "processor" not in kv:
            # (ARM: a trailing block with the hardware name)
            model = model or kv.get("hardware") or kv.get("model name")
            continue
        logical += 1
        vendor = vendor or kv.get("vendor_id") or (
            "ARM implementer " + kv["cpu implementer"] if "cpu implementer" in kv else None)
        model = model or kv.get("model name")
        if "flags" in kv:
            flags |= set(kv["flags"].split())
        if "features" in kv:
            arm = True
            flags |= set(kv["features"].split())
        if "physical id" in kv and "core id" in kv:
            cores.add((kv["physical id"], kv["core id"]))
        if "siblings" in kv and "cpu cores" in kv:
            siblings, cpu_cores = int(kv["siblings"]), int(kv["cpu cores"])
    isa = {CPUINFO_FLAGS[f] for f in flags if f in CPUINFO_FLAGS}
    if arm or "asimd" in flags:
        isa.add("arm64")
        arm = True
    physical = len(cores) or None
    return {"arch": "arm64" if arm else "x86_64", "vendor": vendor, "model": model,
            "isa": sorted(isa), "logical": logical, "physical": physical,
            "smt_siblings": (siblings, cpu_cores) if siblings else None}


def parse_lscpu(text):
    """{"threads_per_core", "cores_per_socket", "sockets"} from lscpu output"""
    out = {}
    keys = {"thread(s) per core": "threads_per_core", "core(s) per socket": "cores_per_socket",
            "socket(s)": "sockets", "core(s) per cluster": "cores_per_socket"}
    for line in text.splitlines():
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        k = k.strip().lower()
        if k in keys and keys[k] not in out:
            try:
                out[keys[k]] = int(v.strip())
            except ValueError:
                pass
    return out


def parse_sysctl(text):
    """{"arch", "vendor", "model", "isa", "logical", "physical", "perf"}
    from `sysctl -a` style text ("key: value" lines) on macOS"""
    kv = {}
    for line in text.splitlines():
        m = re.match(r"^([\w.]+)\s*[:=]\s*(.*)$", line.strip())
        if m:
            kv[m.group(1)] = m.group(2).strip()
    isa = set()
    for k, name in SYSCTL_OPTIONAL.items():
        if kv.get(k, "0").strip() == "1":
            isa.add(name)
    for k in ("machdep.cpu.features", "machdep.cpu.leaf7_features",
              "machdep.cpu.extfeatures"):
        for w in kv.get(k, "").split():
            if w.upper() in SYSCTL_WORDS:
                isa.add(SYSCTL_WORDS[w.upper()])
    arm = "arm64" in isa or kv.get("hw.machine", "").startswith("arm64")
    if arm:
        isa.add("arm64")

    def num(k):
        try:
            return int(kv[k])
        except (KeyError, ValueError):
            return None
    return {"arch": "arm64" if arm else "x86_64",
            "vendor": kv.get("machdep.cpu.vendor") or ("Apple" if arm else None),
            "model": kv.get("machdep.cpu.brand_string"), "isa": sorted(isa),
            "logical": num("hw.logicalcpu") or num("hw.ncpu"),
            "physical": num("hw.physicalcpu"), "perf": num("hw.perflevel0.physicalcpu")}


def isa_class(isa):
    """the class of a CPU (or build) for this search from its ISA names"""
    s = set(isa)
    if "arm64" in s:
        return "arm64 (portable matrix path, ~4-5x slower per core)"
    if "avx512bw" not in s:
        return ("avx2 (no AVX-512BW: matrix path, ~5-12x slower per core)" if "avx2" in s
                else "x86-64 without AVX2 (matrix path)")
    missing = [f for f in FAST_EXT if f not in s]
    if missing:
        return "avx512 cascadelake class (no " + "/".join(missing) + ": ~1.05x the CPU)"
    return "avx512 full (Ice Lake / Sapphire Rapids / Zen 4 class: the fast path)"


def run_text(cmd):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return ""


def cgroup_cpus():
    """the CPU quota of this cgroup (v2 cpu.max or v1 cfs), or None"""
    try:
        q, p = open("/sys/fs/cgroup/cpu.max").read().split()[:2]
        return None if q == "max" else int(q) / int(p)
    except (OSError, ValueError):
        pass
    try:
        q = int(open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read())
        p = int(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read())
        return None if q <= 0 else q / p
    except (OSError, ValueError):
        return None


def detect_cpu():
    if sys.platform == "darwin":
        info = parse_sysctl(run_text(["sysctl", "-a"]))
    else:
        try:
            text = open("/proc/cpuinfo").read()
        except OSError:
            text = ""
        info = parse_cpuinfo(text)
        ls = parse_lscpu(run_text(["lscpu"]))
        if not info["physical"] and ls.get("cores_per_socket"):
            info["physical"] = ls["cores_per_socket"] * ls.get("sockets", 1)
        if ls.get("threads_per_core"):
            info["threads_per_core"] = ls["threads_per_core"]
        info["perf"] = None
    logical = info.get("logical") or os.cpu_count() or 1
    try:
        avail = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        avail = logical
    info["logical"] = logical
    info["available"] = avail
    info["cgroup_cpus"] = cgroup_cpus()
    phys = info.get("perf") or info.get("physical") or logical
    if info.get("threads_per_core", 1) > 1 and not info.get("physical"):
        phys = max(1, logical // info["threads_per_core"])
    info["cores"] = min(phys, avail)
    info["smt"] = logical > (info.get("physical") or logical) or info.get("threads_per_core", 1) > 1
    info["class"] = isa_class(info["isa"])
    return info


# ------------------------------------------------------------------ running

def find_bin(name, msearch):
    p = os.path.join(os.path.dirname(os.path.abspath(msearch)), name)
    return p if os.path.exists(p) else None


def default_msearch():
    for p in (os.path.join(ROOT, "build-laptop", "msearch"), os.path.join(ROOT, "bin", "msearch"),
              os.path.join(ROOT, "cmake-build-release", "src", "c", "msearch")):
        if os.path.exists(p):
            return p
    return os.path.join(ROOT, "bin", "msearch")


def read_records(path):
    recs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("{"):
                recs.append(json.loads(line))
    return recs


def build_of(msearch, tmp):
    """(path, isa) from the first record of a tiny msearch run"""
    out = os.path.join(tmp, "probe.jsonl")
    subprocess.run([msearch, *PROBE, "--out", out], check=True, stdout=subprocess.DEVNULL)
    recs = read_records(out)
    if not recs or "path" not in recs[0]:
        raise SystemExit(f"{msearch}: no \"path\" in its first record (an msearch older than "
                         f"the machine calibration?)")
    sums = [r for r in recs if r.get("type") == "sum"]
    if not sums or sums[0].get("squares") != 2:
        raise SystemExit(f"EXACTNESS FAILED: the probe sum gave {sums}")
    return recs[0]["path"], recs[0].get("isa", [])


def parse_bench(text):
    """(total nodes, all ok, instances) from bench output"""
    total, ok, n = None, True, 0
    for line in text.splitlines():
        if line.startswith("TOTAL"):
            f = line.split()
            total = int(f[1])
            ok = ok and f[2] == "ok"
        elif re.match(r"^\d+x\d+ ", line):
            n += 1
            if " ok " not in line:
                ok = False
        elif "FAIL" in line:
            ok = False
    return total, ok, n


def exactness(bench, path, refs, files=("quick", "full"), say=print):
    """bench each file: every instance ok and the path's node total;
    returns {file: {...}}; raises SystemExit(1) on a mismatch"""
    res = {}
    ref = refs.get(path)
    if ref is None:
        raise SystemExit(f"EXACTNESS FAILED: no reference node counts for search path {path!r}")
    for f in files:
        t0 = time.time()
        p = subprocess.run([bench, os.path.join(ROOT, "bench", f + ".txt")], capture_output=True,
                           text=True)
        total, ok, n = parse_bench(p.stdout)
        good = p.returncode == 0 and ok and total == ref[f]
        res[f] = {"nodes": total, "expected": ref[f], "instances": n, "ok": good,
                  "wall": round(time.time() - t0, 3)}
        say(f"  bench {f}: {n} instances {'ok' if ok else 'FAILED'}, nodes {total} "
            f"(reference {ref[f]} on {path}) -> {'ok' if good else 'MISMATCH'}")
        if not good:
            raise SystemExit(f"EXACTNESS FAILED: bench/{f}.txt on {path}: {res[f]}\n{p.stdout}")
    return res


def check_test(t, rec, path, refs):
    """the nodes (and squares or pairs) of a reference run; returns an error
    string or None"""
    want = refs[path][t["name"]]
    if rec is None:
        return f"{t['name']}: no {t['type']} record"
    if rec.get("nodes") != want:
        return f"{t['name']}: nodes {rec.get('nodes')} != reference {want} ({path})"
    if "squares" in t and rec.get("squares") != t["squares"]:
        return f"{t['name']}: squares {rec.get('squares')} != {t['squares']}"
    if "pairs" in t and rec.get("pairs") != t["pairs"]:
        return f"{t['name']}: pairs {rec.get('pairs')} != {t['pairs']}"
    return None


def run_batch(msearch, t, K, tmp):
    """K concurrent copies of test t: per copy (wall, process CPU, record)"""
    procs = {}
    for i in range(K):
        out = os.path.join(tmp, f"{t['name']}_K{K}_{i}.jsonl")
        if os.path.exists(out):
            os.remove(out)
        p = subprocess.Popen([msearch, *t["args"], "--out", out], stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
        procs[p.pid] = (p, out, time.monotonic())
    res = []
    while len(res) < K:
        pid, status, ru = os.wait4(-1, 0)
        if pid not in procs:
            continue
        p, out, t0 = procs[pid]
        wall = time.monotonic() - t0
        p.returncode = os.waitstatus_to_exitcode(status) if hasattr(
            os, "waitstatus_to_exitcode") else (status >> 8)
        recs = read_records(out) if os.path.exists(out) else []
        rec = next((r for r in recs if r.get("type") == t["type"]), None)
        res.append({"wall": wall, "cpu": ru.ru_utime + ru.ru_stime, "rc": p.returncode,
                    "rec": rec, "first": recs[0] if recs else None})
    return res


def band(rows, names, K):
    ref = sum(t["ref_cpu"] for t in TESTS if t["name"] in names)
    wall = sum(rows[n]["wall"] for n in names)
    cpu = sum(rows[n]["cpu"] for n in names)
    return K * ref / wall, ref / cpu


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--bin", default=None, help="msearch to calibrate (bench beside it)")
    ap.add_argument("--out", default="cal.json")
    ap.add_argument("--ks", default=None,
                    help="concurrent copies to measure, e.g. 1,8,16 (default: 1, the physical "
                         "or performance cores, the logical CPUs)")
    ap.add_argument("--no-bench", action="store_true",
                    help="skip the bench exactness gate (the T1-T4 node checks still run)")
    ap.add_argument("--refs", default=None,
                    help="JSON overriding reference node counts, e.g. {\"carry512\": {\"T1\": 1}}")
    ap.add_argument("--tests", default="T1,T2,T3,T4", help=argparse.SUPPRESS)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()
    say = (lambda *x, **k: None) if a.quiet else (lambda *x, **k: print(*x, **k, flush=True))
    msearch = os.path.abspath(a.bin or default_msearch())
    if not os.path.exists(msearch):
        sys.exit(f"{msearch} not found: build first (scripts/laptop/build.sh or ./build.sh)")
    refs = json.loads(json.dumps(REFS))
    if a.refs:
        with open(a.refs) as f:
            for path, d in json.load(f).items():
                refs.setdefault(path, {}).update(d)
    tests = [t for t in TESTS if t["name"] in a.tests.split(",")]
    t_start = time.time()
    cpu = detect_cpu()
    warnings = []
    say(f"CPU: {cpu.get('model')} ({cpu['arch']}); {cpu['logical']} logical CPUs "
        f"({cpu['available']} available), {cpu.get('physical')} physical cores"
        + (f", {cpu['perf']} performance cores" if cpu.get("perf") else "")
        + (", SMT" if cpu["smt"] else "")
        + (f", cgroup quota {cpu['cgroup_cpus']:g} CPUs" if cpu.get("cgroup_cpus") else ""))
    say(f"ISA: {' '.join(cpu['isa']) or '-'}; class: {cpu['class']}")
    try:
        load0 = os.getloadavg()[0]
    except (AttributeError, OSError):
        load0 = None
    if load0 is not None and load0 > 0.5:
        warnings.append(f"load average {load0:.2f} before the run: other work on this machine "
                        f"makes the speeds below read low (rerun on an idle machine)")
    with tempfile.TemporaryDirectory(prefix="machine_cal_") as tmp:
        path, bisa = build_of(msearch, tmp)
        say(f"build: {msearch}: search path {path}; ISA at compile time: {' '.join(bisa) or '-'}")
        if path != "carry512":
            warnings.append(f"the build's search path is {path!r}, not 'carry512' (AVX-512BW, "
                            f"the fast path): each sum costs ~4-12x the CPU of the reference "
                            f"build" + ("; this CPU has AVX-512BW: rebuild with -march=native"
                                        if "avx512bw" in cpu["isa"] else ""))
        miss = [f for f in FAST_EXT if f not in bisa]
        if path == "carry512" and miss:
            warnings.append(f"the build lacks {', '.join(miss)} (cascadelake class: ~1.05x the "
                            f"CPU)" + ("; the CPU has them: rebuild with -march=native"
                                       if all(f in cpu["isa"] for f in miss) else ""))
        if cpu["isa"] and "arm64" not in cpu["isa"]:
            lost = [f for f in bisa if f not in cpu["isa"] and f != "arm64"]
            if lost:
                warnings.append(f"the build uses {', '.join(lost)}, which /proc/cpuinfo does not "
                                f"list: it may crash on this CPU")
        exact = {}
        if not a.no_bench:
            bench = find_bin("bench", msearch)
            if bench is None:
                sys.exit(f"no bench beside {msearch} (or pass --no-bench)")
            say("exactness gate (bench quick, full):")
            exact = exactness(bench, path, refs, say=say)
        if a.ks:
            ks = sorted({int(k) for k in a.ks.split(",")})
        else:
            ks = sorted({1, cpu["cores"], cpu["available"]})
        runs = []
        for K in ks:
            rows = {}
            for t in tests:
                res = run_batch(msearch, t, K, tmp)
                for r in res:
                    err = (f"{t['name']}: msearch exit {r['rc']}" if r["rc"]
                           else check_test(t, r["rec"], path, refs))
                    if err:
                        raise SystemExit(f"EXACTNESS FAILED: {err}")
                est = [r["rec"].get("est_time", r["rec"].get("cpu")) if t["type"] != "sum"
                       else r["rec"]["cpu"] for r in res]
                rows[t["name"]] = {
                    "wall": sum(r["wall"] for r in res) / K,
                    "wall_max": max(r["wall"] for r in res),
                    "cpu": sum(r["cpu"] for r in res) / K,
                    "rec_cpu": sum(r["rec"]["cpu"] for r in res) / K,
                    "est_time": sum(est) / K,
                    "nodes": res[0]["rec"]["nodes"],
                }
                rows[t["name"]]["law_ratio"] = rows[t["name"]]["est_time"] / t["law"]
                rows[t["name"]]["speed"] = K * t["ref_cpu"] / rows[t["name"]]["wall"]
                say(f"  K={K:<3d} {t['name']}: wall {rows[t['name']]['wall']:7.2f} s (max "
                    f"{rows[t['name']]['wall_max']:.2f}), CPU/process {rows[t['name']]['cpu']:7.2f}"
                    f" s (reference {t['ref_cpu']:.2f}), {rows[t['name']]['speed']:.3f} "
                    f"reference CPU-s per s, est/law {rows[t['name']]['law_ratio']:.3f}")
            slow = [n for n, r in rows.items() if r["cpu"] > 0 and r["wall"] / r["cpu"] > 1.05]
            if K <= cpu["cores"] and slow:
                warnings.append(f"K={K}: wall exceeds process CPU by more than 5% on "
                                f"{', '.join(slow)}: the machine was busy, so these speeds read low")
            run = {"K": K, "tests": rows}
            names = {m: [t["name"] for t in tests if t["mode"] == m] for m in ("plain", "dfirst")}
            for m, ns in names.items():
                if ns:
                    run[m + "_speed"], run[m + "_per_process"] = band(rows, ns, K)
            sp = [run.get("plain_speed"), run.get("dfirst_speed")]
            sp = [s for s in sp if s]
            run["score"] = (sp[0] * sp[-1]) ** 0.5 if sp else 0.0
            runs.append(run)
            say(f"K={K}: plain {run.get('plain_speed', 0):.3f}, d-first "
                f"{run.get('dfirst_speed', 0):.3f} reference CPU-hours per instance-hour")
    best = max(runs, key=lambda r: r["score"])
    k1 = next((r for r in runs if r["K"] == 1), None)
    if best["K"] > cpu["cores"] and cpu["smt"]:
        gain = best["score"] / max(r["score"] for r in runs if r["K"] <= cpu["cores"])
        warnings.append(f"SMT: {best['K']} workers on {cpu['cores']} cores give {gain:.2f}x the "
                        f"throughput of one per core, not {best['K'] / cpu['cores']:.1f}x: a vCPU "
                        f"is not a core (price the instance by the speeds here)")
    load = {}
    if k1 is not None and best is not k1:
        load = {n: best["tests"][n]["cpu"] / k1["tests"][n]["cpu"] for n in best["tests"]}
    ps = best.get("plain_speed")
    ds = best.get("dfirst_speed")
    cal = {
        "version": 1,
        "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "host": platform.node(), "system": platform.platform(), "machine": platform.machine(),
        "cpu": cpu,
        "build": {"msearch": msearch, "sha256": hashlib.sha256(open(msearch, "rb").read()).hexdigest(),
                  "path": path, "isa": bisa},
        "exactness": exact,
        "reference": {t["name"]: {"args": t["args"], "ref_cpu": t["ref_cpu"], "law": t["law"],
                                  "mode": t["mode"]} for t in tests},
        "runs": runs,
        "workers": best["K"],
        # reference CPU-hours per instance-hour (= reference CPU-years per
        # instance-year) at the recommended worker count
        "plain_speed": ps, "dfirst_speed": ds,
        "plain_ref_cpu_years_per_instance_hour": ps / YEAR_H if ps else None,
        "dfirst_ref_cpu_years_per_instance_hour": ds / YEAR_H if ds else None,
        # reference CPU per process CPU second at that load
        "per_process": {"plain": best.get("plain_per_process"),
                        "dfirst": best.get("dfirst_per_process")},
        # the process's CPU (T1) or whole-sum estimate (T2-T4) / the
        # scheduler's shipped law, alone (K = 1) and at the recommended load
        "law_ratio": {"K1": {n: r["law_ratio"] for n, r in k1["tests"].items()} if k1 else None,
                      "best": {n: r["law_ratio"] for n, r in best["tests"].items()}},
        "load_factor": load,
        "warnings": warnings,
        "wall_seconds": round(time.time() - t_start, 1),
    }
    with open(a.out, "w") as f:
        json.dump(cal, f, indent=1)
    say("")
    say(f"recommended workers: {best['K']}")
    if ps:
        say(f"plain search:   {ps:.3f} reference CPU-hours per instance-hour "
            f"({ps / best['K']:.3f} per worker; {ps / YEAR_H:.3g} reference CPU-years per "
            f"instance-hour)")
    if ds:
        say(f"d-first search: {ds:.3f} reference CPU-hours per instance-hour "
            f"({ds / best['K']:.3f} per worker; {ds / YEAR_H:.3g} reference CPU-years per "
            f"instance-hour)")
    if ps and ds:
        say(f"1 reference CPU-year takes ~{YEAR_H / ps:.0f} (plain) to {YEAR_H / ds:.0f} "
            f"(d-first) instance-hours here")
    for w in warnings:
        say("\n*** WARNING: " + w + " ***", file=sys.stderr)
    say(f"wrote {a.out} ({cal['wall_seconds']:.0f} s)")


if __name__ == "__main__":
    main()
