#!/usr/bin/env python3
"""
What leaves the laptop: the scrubbed forms of run.py's meta files,
build_info.txt, units.tsv and progress.log that scripts/laptop/digest.py
and collect.py write. Python 3.8+, standard library only.

* Kept: the machine's architecture, its chip model and performance cores,
  workers, python major.minor, the compiler's family and major version,
  the build flags, the search kernels, the bench/quick.txt times, the
  repo commit of the build, the plan's name and sha256, CPU and wall
  durations.
* Left out: host names, OS and kernel versions, user names and absolute
  paths (a path keeps its last part), time zones and wall-clock times.
  Every time is in seconds from the first run's start (the earliest time
  found in DIR). Only the fields and lines named here are kept: a meta
  field, a build_info.txt line or a progress.log line of any other kind
  is left out (a progress.log line is kept only if all its words are
  run.py's own).
"""
import calendar
import datetime
import json
import os
import platform
import re
import subprocess
import sys
import time

META = re.compile(r"meta_.*\.json$")
# progress.log's time stamps: run.py's local time, or (scrubbed) the time
# from the first run's start
STAMP = re.compile(r"\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\] (.*)$")
REL_STAMP = re.compile(r"\[\+(\d+):(\d\d):(\d\d)\] (.*)$")
# absolute paths (/x, ~/x, ~user/x; also -I/x), and dates
PATH = re.compile(r"(?<![\w.~])~[\w.-]*/[^\s'\"(),;\]]*|(?<![\w.~])~[A-Za-z_][\w.-]*"
                  r"|(?<![\w.~/])/[^\s'\"(),;\]]+")
OPT_PATH = re.compile(r"(?<=-[A-Za-z])/(?:[^\s'\"(),;\]/]+/)+[^\s'\"(),;\]]*")
DATE = re.compile(r"\d{4}-\d\d-\d\d(?:[T ]\d\d:\d\d(?::\d\d(?:\.\d+)?)?(?:\s?[+-]\d\d:?\d\d|Z)?)?")
STATUS = ("done", "incomplete", "cut", "cut-kept")
# the words of run.py's progress lines (all of its versions)
WORDS = set("""
units done this run CPU h in reference x the per core squares with an SP diagonal magic
last min hour stage band so far pairs predicted at f_rho finished stopping finishing nothing
new running unit failed exit it will be retried next stopped by wall clock limit after s
removed runs again machine slept or was busy kept as is not complete times now seems to have
for a meanwhile hours reached them were stall rerun cut U already to workers its
""".split())
NUM = re.compile(r"\d+(?:\.\d+)?")


def last_part(p):
    """a path's last part; a home directory (or one level under /) as a
    whole: <path>"""
    parts = [x for x in p.split("/") if x]
    if p.startswith("~"):
        if len(parts) <= 1:
            return "<path>"
    elif p.startswith("/") and len(parts) <= 2:
        return "<path>"
    return parts[-1] if parts else "<path>"


def scrub_text(s):
    """s without absolute paths (each keeps its last part) and dates"""
    s = PATH.sub(lambda m: last_part(m.group(0)), s)
    s = OPT_PATH.sub(lambda m: last_part(m.group(0)), s)
    return DATE.sub("<date>", s)


# --------------------------------------------------------------------------
# the machine (run.py's meta files; digest.py's made_on)


def sysctl(name):
    try:
        return subprocess.run(["sysctl", "-n", name], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, universal_newlines=True,
                              timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def chip():
    """the CPU's model (e.g. "Apple M1 Pro"), or None"""
    s = ""
    if sys.platform == "darwin":
        s = sysctl("machdep.cpu.brand_string")
    elif sys.platform.startswith("linux"):
        try:
            with open("/proc/cpuinfo") as f:
                for line in f:
                    if line.startswith("model name"):
                        s = line.split(":", 1)[1]
                        break
        except OSError:
            s = ""
    s = re.sub(r"\s+", " ", s).strip()
    return s[:80] if re.match(r"[\w ().@,+-]+$", s) else None


def perf_cores():
    """the performance cores (macOS), else the CPUs"""
    if sys.platform == "darwin":
        for name in ("hw.perflevel0.physicalcpu", "hw.physicalcpu"):
            try:
                return int(sysctl(name))
            except ValueError:
                pass
    return os.cpu_count()


def python_mm(v=None):
    """"3.9" of "3.9.6" """
    m = re.match(r"(\d+)\.(\d+)", str(v if v is not None else platform.python_version()))
    return f"{m.group(1)}.{m.group(2)}" if m else None


def machine_info():
    return {"machine": platform.machine(), "chip": chip(), "perf_cores": perf_cores(),
            "python": python_mm()}


# --------------------------------------------------------------------------
# build_info.txt (scripts/laptop/build.sh)


def compiler_family(s):
    """"Apple clang 15", "clang 18", "gcc 13" (the major version) or None"""
    m = re.search(r"(?:^|\()((?:Apple clang|clang|gcc) \d+)\)?$", s)
    if m:       # (already scrubbed)
        return m.group(1)
    m = re.search(r"Apple (?:LLVM|clang) version (\d+)", s)
    if m:
        return f"Apple clang {m.group(1)}"
    m = re.search(r"clang version (\d+)", s)
    if m:
        return f"clang {m.group(1)}"
    if re.search(r"gcc|GCC", s):
        v = re.findall(r"\b(\d+)\.\d+(?:\.\d+)?\b", s)
        if v:
            return f"gcc {v[-1]}"
    return None


def scrub_compiler(s):
    """build.sh's compiler line: "CC (version line)" (or an older
    build.sh's "version line") as "CC (family major)", CC's paths cut to
    their last part"""
    fam = compiler_family(s)
    m = re.match(r"(.*?) \((.*)\)$", s)
    cc = None
    if m and "version" not in m.group(1) and compiler_family(m.group(1)) is None:
        toks = [last_part(x) if "/" in x else x for x in m.group(1).split()]
        cc = " ".join(x for x in toks if re.match(r"[\w.+=,-]+$", x)) or None
    if cc:
        return f"{cc} ({fam or 'unknown'})"
    return fam or "unknown"


def scrub_build(text):
    """build_info.txt with only these lines: the compiler (family and major
    version), flags, the architecture, the repo commit, the search
    kernels, the checks and the bench/quick.txt times"""
    out = []
    for line in str(text).splitlines():
        line = line.rstrip()
        if line.startswith("compiler: "):
            out.append("compiler: " + scrub_compiler(line[len("compiler: "):]))
        elif re.match(r"[0-9a-f]{40}$", line) or line == "all checks passed":
            out.append(line)
        elif re.match(r"arch: [\w]+$", line):
            out.append(line)
        elif line.startswith(("flags: ", "search kernels: ", "bench/quick.txt search time: ",
                              "time per core: ")):
            out.append(scrub_text(line))
    return "".join(x + "\n" for x in out)


def scrub_bench(text):
    """bench's output (check_quick.txt): its table lines only"""
    keep = re.compile(r"(instance |\d+x\d+ P=|TOTAL |\d+ of \d+ instances)")
    return "".join(scrub_text(x.rstrip()) + "\n" for x in str(text).splitlines() if keep.match(x))


# --------------------------------------------------------------------------
# times: seconds from the first run's start


def parse_started(s):
    """(unix time, UTC offset s, local naive datetime) of run.py's
    "started" (%Y-%m-%dT%H:%M:%S%z), or None"""
    try:
        t = datetime.datetime.strptime(str(s), "%Y-%m-%dT%H:%M:%S%z")
    except ValueError:
        return None
    return t.timestamp(), t.utcoffset().total_seconds(), t.replace(tzinfo=None)


class Runs:
    """the runs of DIR (its meta files) and the time from the first one's
    start. A scrubbed DIR (collect.py's tarballs unpacked: meta files with
    started_s) is already relative"""

    def __init__(self, out_dir):
        self.metas = []         # [(start, meta, local start, UTC offset)] by start
        self.relative = False
        names = sorted(x for x in os.listdir(out_dir) if META.match(x))
        for name in names:
            try:
                with open(os.path.join(out_dir, name)) as f:
                    m = json.load(f)
            except (OSError, ValueError):
                raise ValueError(f"{os.path.join(out_dir, name)}: not a meta file run.py wrote")
            if not isinstance(m, dict):
                raise ValueError(f"{os.path.join(out_dir, name)}: not a meta file run.py wrote")
            if "started_s" in m:
                self.relative = True
                self.metas.append((float(m["started_s"]), m, None, None))
                continue
            p = parse_started(m.get("started"))
            if p is None:
                n = re.match(r"meta_(\d+)\.json$", name)
                p = (float(n.group(1)), None, None) if n else (None, None, None)
            self.metas.append((p[0], m, p[2], p[1]))
        self.metas.sort(key=lambda x: (x[0] is None, x[0] or 0.0))
        self.base = None

    def unix(self, local):
        """the time of a progress.log stamp (run.py's local time): by the
        UTC offset of the run then (its meta file's "started"), else this
        machine's time zone. (A run across a change of daylight saving
        time is an hour off after it.)"""
        dt = datetime.datetime.strptime(local, "%Y-%m-%d %H:%M:%S")
        off = None
        for _, _, start, o in self.metas:
            if o is None:
                continue
            if off is None or (start is not None and start <= dt):
                off = o
        if off is None:
            return time.mktime(dt.timetuple())
        return calendar.timegm(dt.timetuple()) - off

    def set_base(self, times):
        """the first run's start: the earliest of the runs' starts and the
        times given"""
        c = [t for t, _, _, _ in self.metas if t is not None] + [t for t in times if t is not None]
        self.base = 0.0 if self.relative else (min(c) if c else None)

    def rel(self, t):
        """t in seconds from the first run's start (None if unknown)"""
        if t is None or self.base is None:
            return None
        return int(round(t - self.base))

    def rel_now(self):
        """now in seconds from the first run's start (None if unknown: no
        start known, or a scrubbed DIR)"""
        return None if self.relative else self.rel(time.time())


def read_tsv(out_dir):
    """units.tsv's rows as (i, pred_time, cpu, wall, end, status), end
    absolute (run.py's) or relative (a scrubbed units.tsv's end_s), and
    whether they are relative; rows of any other form are left out"""
    rows = []
    rel = False
    try:
        with open(os.path.join(out_dir, "units.tsv")) as f:
            head = f.readline().rstrip("\n").split("\t")
            rel = "end_s" in head
            for line in f:
                x = line.rstrip("\n").split("\t")
                if len(x) != len(head):
                    continue
                r = dict(zip(head, x))
                try:
                    v = (int(r["i"]), float(r["pred_time"]), float(r["cpu"]), float(r["wall"]),
                         int(r["end_s" if rel else "end"]), r["status"])
                except (KeyError, ValueError):
                    continue
                if v[5] in STATUS:
                    rows.append(v)
    except OSError:
        pass
    return rows, rel


def clean_message(msg):
    """a progress line's message if it is one of run.py's (scrubbed), else
    None"""
    if msg.startswith("*** MAGIC SQUARE FOUND *** "):
        try:
            r = json.loads(msg[len("*** MAGIC SQUARE FOUND *** "):])
        except ValueError:
            return None
        return "*** MAGIC SQUARE FOUND *** " + json.dumps(r) if isinstance(r, dict) else None
    m = re.match(r"plan (\S+) \(([0-9a-f]{12})\): (.*)$", msg)
    head = ""
    if m:
        name = last_part(m.group(1)) if "/" in m.group(1) else m.group(1)
        if not re.match(r"[\w.+-]+$", name):
            name = "<plan>"
        head = f"plan {name} ({m.group(2)}): "
        msg = m.group(3)
    words = re.findall(r"[A-Za-z_]+", NUM.sub(" ", msg))
    if not words or any(w not in WORDS for w in words):
        return None
    if re.search(r"[^\w\s.,;:()=/*+%-]", msg):
        return None
    return head + msg


def read_progress(out_dir, runs):
    """progress.log's lines as [(unix time or seconds from the start, message)]
    (run.py's lines only) and the count of lines left out"""
    out = []
    left = 0
    try:
        with open(os.path.join(out_dir, "progress.log"), errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                m = STAMP.match(line)
                r = REL_STAMP.match(line)
                t = msg = None
                try:
                    if m and not runs.relative:
                        t, msg = runs.unix(m.group(1)), m.group(2)
                    elif r and runs.relative:
                        t = int(r.group(1)) * 3600 + int(r.group(2)) * 60 + int(r.group(3))
                        msg = r.group(4)
                except (ValueError, OverflowError):
                    msg = None
                msg = clean_message(msg) if msg is not None else None
                if msg is None:
                    left += 1
                    continue
                out.append((t, msg))
    except OSError:
        pass
    return out, left


def stamp(s):
    """[+H:MM:SS] of s seconds from the start ([+?] if unknown)"""
    if s is None:
        return "[+?]"
    s = max(0, int(s))
    return f"[+{s // 3600}:{s // 60 % 60:02d}:{s % 60:02d}]"


def scrub_meta(m, k, started_s, rootdir=None):
    """a meta file of run.py with only what the analysis needs: the run's
    number (k, by start) and start in seconds from the first, the plan's
    name and sha256, workers, the machine (architecture, chip, performance
    cores), python major.minor, the binary's path from rootdir (or its
    last part) and the build (scrubbed)"""
    out = {"run": k, "started_s": started_s}
    for key in ("plan", "plan_sha256", "workers", "machine", "chip", "perf_cores"):
        if key in m:
            v = m[key]
            out[key] = scrub_text(v) if isinstance(v, str) else v
    if "python" in m:
        out["python"] = python_mm(m["python"])
    b = m.get("bin")
    if isinstance(b, str):
        if os.path.isabs(b) and rootdir:
            r = os.path.relpath(b, rootdir)
            b = last_part(b) if r.startswith("..") else r
        out["bin"] = last_part(b) if os.path.isabs(b) or b.startswith("~") else scrub_text(b)
    if isinstance(m.get("build"), str):
        out["build"] = scrub_build(m["build"])
    return out


def scrub_dir(out_dir, rootdir=None):
    """the scrubbed meta files (by start), units.tsv rows (end in seconds
    from the first run's start) and progress lines ([+H:MM:SS] message) of
    run.py's DIR, and the count of progress lines left out"""
    runs = Runs(out_dir)
    rows, rel = read_tsv(out_dir)
    if rel and not runs.metas:
        runs.relative = True
    lines, left = read_progress(out_dir, runs)
    runs.set_base(([] if rel else [e - w for _, _, _, w, e, _ in rows])
                  + ([] if runs.relative else [t for t, _ in lines[:1]]))
    metas = [scrub_meta(m, k + 1, runs.rel(t), rootdir)
             for k, (t, m, _, _) in enumerate(runs.metas)]
    rows = [(i, p, c, w, e if rel else runs.rel(e), s) for i, p, c, w, e, s in rows]
    lines = [stamp(runs.rel(t)) + " " + msg for t, msg in lines]
    return {"metas": metas, "tsv": rows, "progress": lines, "progress_left_out": left,
            "runs": runs}


def tsv_text(rows):
    return "i\tpred_time\tcpu\twall\tend_s\tstatus\n" + "".join(
        f"{i}\t{p:.3f}\t{c:.3f}\t{w:.1f}\t{e}\t{s}\n" for i, p, c, w, e, s in rows)
