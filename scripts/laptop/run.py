#!/usr/bin/env python3
"""
Run a static plan (scripts/laptop/make_plan.py) on any machine, in the plan's
order, with K parallel msearch processes. Python 3.8+, standard library only.

usage: run.py PLAN(.jsonl or .jsonl.xz) [--bin BUILD/msearch] [--out DIR] [--workers K]
              [--hours H] [--max-units U] [--rerun-cut]

* Resumable: a unit whose output has a "done" record is skipped; a unit that
  was interrupted is run again from the start (its output is replaced).
* Ctrl-C, kill (TERM), closing the terminal (HUP) or --hours stops
  launching, stops the running units, removes their partial output, leaves
  everything finished in place and exits. Run it again to continue.
* Outputs: DIR/U000123.jsonl per unit, DIR/meta_<time>.json (machine, build,
  plan hash), DIR/progress.log. A magic square (best_score >= 14) is
  announced and copied to DIR/MAGIC.txt. U<i> is line i of the plan: a
  finished U<i> whose P or sums are not line i's (another plan's output)
  is refused. Plans that share their first lines (research/stage1's plan
  begins with all of plan-20261010) share those outputs.
* Progress: the units' measured CPU against their predicted reference CPU
  (the fast x86 build's; the ratio also carries the machine's load and the
  predictions' own error, e.g. the fast build itself ran at 1.32x its
  predictions on a busy shared machine), the same over the last 30 minutes
  with the reference CPU-hours done per hour there (throttling and other
  load show in these), and for a stage-1 plan the pairs found in its band
  units (research/stage1.md; scripts/decide.py fits them). DIR/units.tsv
  has a line per finished unit (and per run cut by the wall clock, below):
  i, pred_time, CPU and wall seconds, end time, status.
* The plan's --time-limit is wall-clock time, which runs on while the
  machine sleeps or the unit waits for a core. A unit stopped by it whose
  wall time exceeds its CPU by more than a quarter of the limit lost that
  time to a sleep or other load: it is not kept but run again from its
  start (at most 3 times per run; then it is kept as it is, and later runs
  skip it). Throttling slows CPU and wall time alike, so it opens no such
  gap. A unit stopped by the limit with CPU close to its wall time really
  needed that long and is kept, as are units incomplete for any other
  reason. Finished units do not depend on the speed: the search is exact.
  --rerun-cut also runs again the finished units cut that way (by an older
  run.py, or kept after 3 tries); without it they are counted and named.
"""
import argparse
import hashlib
import json
import lzma
import os
import platform
import signal
import subprocess
try:
    import resource
except ImportError:
    resource = None
import sys
import time


def default_workers():
    if sys.platform == "darwin":  # the performance cores only
        try:
            out = subprocess.run(["sysctl", "-n", "hw.perflevel0.physicalcpu"],
                                 capture_output=True, text=True, check=True).stdout
            return max(1, int(out.strip()))
        except (OSError, ValueError, subprocess.CalledProcessError):
            pass
    return max(1, (os.cpu_count() or 2) - 1)


def done_file(path):
    """the parsed records of a finished unit file, or None"""
    try:
        with open(path) as f:
            recs = [json.loads(line) for line in f if line.strip()]
    except (OSError, ValueError):
        return None
    return recs if any(r.get("type") == "done" for r in recs) else None


# reruns per run of a unit stopped by the wall-clock limit after a stall
RETRIES = 3
# a stall: the wall time ahead of the CPU by this part of the limit
STALL = 0.25
# seconds: the window of the progress line's recent speed
WINDOW = 1800


def time_limit(u):
    """the unit's --time-limit (msearch's wall-clock seconds), or None"""
    args = u.get("args", [])
    for k, x in enumerate(args[:-1]):
        if x == "--time-limit":
            try:
                return float(args[k + 1])
            except ValueError:
                return None
    return None


def records_cpu(recs):
    """the process CPU seconds of a unit's records (sum, csum and dsum
    records' cpu, else their time)"""
    return sum(r.get("cpu", r.get("time", 0.0)) for r in recs
               if r.get("type") in ("sum", "csum", "dsum"))


def cut_by_clock(u, recs, proc_cpu=0.0):
    """(CPU, wall) seconds if the unit stopped at its wall-clock limit
    (--time-limit) with its wall time ahead of its CPU by more than STALL
    of the limit: the machine slept or was busy for that long, and the unit
    would have gone further without it. Else None. The CPU is the larger
    of the records' (which leave out a small enumeration share) and
    proc_cpu, the process's own (from the OS, when run.py ran it)"""
    lim = time_limit(u)
    done = [r for r in recs if r.get("type") == "done"]
    if not lim or not done or done[-1].get("complete", 1):
        return None
    wall = float(done[-1].get("time", 0.0))
    cpu = max(records_cpu(recs), proc_cpu)
    if wall >= 0.99 * lim and wall - cpu > STALL * lim:
        return cpu, wall
    return None


def children_cpu():
    """CPU seconds of the child processes waited for so far"""
    if resource is None:
        return 0.0
    ru = resource.getrusage(resource.RUSAGE_CHILDREN)
    return ru.ru_utime + ru.ru_stime


def unit_matches(u, recs):
    """the records' P and sums are those of plan line u"""
    P = list(u["P"])
    while P and P[-1] == 0:
        P.pop()
    for r in recs:
        if r.get("type") in ("sum", "square"):
            Q = list(r.get("P", []))
            while Q and Q[-1] == 0:
                Q.pop()
            if Q != P or not u["lo"] <= r.get("S", -1) <= u["hi"]:
                return False
    return True


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(os.path.dirname(here))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("plan")
    ap.add_argument("--bin", default=os.path.join(root, "build-laptop", "msearch"))
    ap.add_argument("--out", default=os.path.join(root, "laptop-runs"))
    ap.add_argument("--workers", type=int, default=default_workers())
    ap.add_argument("--hours", type=float, default=0, help="stop after this long (0 = never)")
    ap.add_argument("--max-units", type=int, default=0)
    ap.add_argument("--rerun-cut", action="store_true",
                    help="also run again finished units stopped by the wall-clock limit "
                         "after a stall (an older run.py's, or kept after 3 tries)")
    a = ap.parse_args()

    with open(a.plan, "rb") as f:
        raw = f.read()
    if a.plan.endswith(".xz"):
        raw = lzma.decompress(raw)
    units = [json.loads(line) for line in raw.decode().splitlines() if line.strip()]
    if a.max_units:
        units = units[:a.max_units]
    os.makedirs(a.out, exist_ok=True)
    # one run.py per output directory
    try:
        import fcntl
        lock = open(os.path.join(a.out, ".lock"), "w")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            sys.exit(f"another run.py is running on {a.out}")
    except ImportError:
        pass
    if not os.path.exists(a.bin):
        sys.exit(f"{a.bin} not found: run scripts/laptop/build.sh first")
    meta = {"plan": os.path.basename(a.plan), "plan_sha256": hashlib.sha256(raw).hexdigest(),
            "bin": a.bin, "workers": a.workers, "machine": platform.machine(),
            "system": platform.platform(), "processor": platform.processor(),
            "python": platform.python_version(), "started": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    build_info = os.path.join(os.path.dirname(a.bin), "build_info.txt")
    if os.path.exists(build_info):
        meta["build"] = open(build_info).read()
    with open(os.path.join(a.out, f"meta_{int(time.time())}.json"), "w") as f:
        json.dump(meta, f, indent=1)
    log = open(os.path.join(a.out, "progress.log"), "a")
    tsv_path = os.path.join(a.out, "units.tsv")
    tsv = open(tsv_path, "a")
    if os.path.getsize(tsv_path) == 0:
        tsv.write("i\tpred_time\tcpu\twall\tend\tstatus\n")
        tsv.flush()

    def unit_row(u, cpu, wall, status):
        # status: done, incomplete (kept), cut (stopped by the wall clock
        # after a stall: removed and run again), cut-kept (after RETRIES)
        tsv.write(f"{u['i']}\t{u.get('pred_time', 0.0):.3f}\t{cpu:.3f}\t{wall:.1f}\t"
                  f"{int(time.time())}\t{status}\n")
        tsv.flush()

    def say(msg):
        line = time.strftime("[%Y-%m-%d %H:%M:%S] ") + msg
        try:
            print(line, flush=True)
        except OSError:     # (the terminal is gone: the log still has it)
            pass
        log.write(line + "\n")
        log.flush()

    # partial outputs of an earlier run (U<i>.jsonl.part.<pid>): removed, so
    # that a unit is never written by two processes under one name
    for name in os.listdir(a.out):
        if ".jsonl.part" in name:
            try:
                os.remove(os.path.join(a.out, name))
            except OSError:
                pass
    todo = []
    n_done = 0
    old_cut = []    # finished units stopped by the wall clock after a stall
    for u in units:
        path = os.path.join(a.out, f"U{u['i']:06d}.jsonl")
        recs = done_file(path)
        if recs:
            if not unit_matches(u, recs):
                sys.exit(f"{path} is not the output of line {u['i']} of {a.plan} (another "
                         f"plan's?): use another --out")
            if cut_by_clock(u, recs):
                old_cut.append(u["i"])
                if a.rerun_cut:
                    todo.append(u)
                    continue
            n_done += 1
        else:
            todo.append(u)
    say(f"plan {meta['plan']} ({meta['plan_sha256'][:12]}): {len(units)} units, "
        f"{n_done} already done, {len(todo)} to run with {a.workers} workers")
    if old_cut:
        names = " ".join(f"U{i:06d}" for i in old_cut[:10]) + (" ..." if len(old_cut) > 10 else "")
        say(f"{len(old_cut)} finished units were stopped by the wall-clock limit after a stall "
            f"({names}): " + ("running them again" if a.rerun_cut else
                              "kept; --rerun-cut runs them again"))

    stop = {"now": False}

    def on_signal(signum, frame):
        if stop["now"]:
            return
        stop["now"] = True
        say("stopping: finishing nothing new, stopping the running units")
    signal.signal(signal.SIGINT, on_signal)
    signal.signal(signal.SIGTERM, on_signal)
    if hasattr(signal, "SIGHUP"):
        # (closing the terminal: stop the units too, they run in their own
        # sessions and would otherwise go on as orphans)
        signal.signal(signal.SIGHUP, on_signal)

    t0 = time.time()
    running = {}  # Popen -> (unit, tmp path, output path, launch time)
    stats = {"units": 0, "cpu": 0.0, "ref": 0.0, "squares": 0, "sp": 0, "magic": 0}
    recent = []   # (end time, CPU, pred_time) of the units of the last WINDOW s
    tries = {}    # unit i -> times stopped by the wall clock after a stall
    # (a stage-1 plan: pairs in its band units, all finished ones so far)
    band = {"units": 0, "pairs": 0, "pred": 0.0}
    stage1 = any("stage1" in u for u in units)

    def count_band(u, recs):
        if u.get("stage1"):
            band["units"] += 1
            band["pred"] += u.get("pred_pairs", 0.0)
            band["pairs"] += sum(r.get("sp_count", 0) for r in recs if r.get("type") == "square")
    def magic(r):
        msg = json.dumps(r)
        say("*** MAGIC SQUARE FOUND *** " + msg)
        with open(os.path.join(a.out, "MAGIC.txt"), "a") as f:
            f.write(msg + "\n")
    rerun = set(old_cut) if a.rerun_cut else set()
    for u in units:
        if stage1 and u.get("stage1") and u["i"] not in rerun:
            recs = done_file(os.path.join(a.out, f"U{u['i']:06d}.jsonl"))
            if recs:
                count_band(u, recs)
    last = 0.0
    clocks = (time.time(), time.monotonic())
    while (todo and not stop["now"]) or running:
        # (the wall clock ahead of the monotonic one, which stops while the
        # machine sleeps: log it; best effort, the units' records decide)
        now = (time.time(), time.monotonic())
        gap = (now[0] - clocks[0]) - (now[1] - clocks[1])
        clocks = now
        if gap > 120:
            say(f"the machine seems to have slept for {gap / 60:.0f} min: a unit stopped by its "
                f"wall-clock limit meanwhile is run again")
        if a.hours and time.time() - t0 > a.hours * 3600 and not stop["now"]:
            stop["now"] = True
            say(f"--hours {a.hours} reached: stopping")
        while todo and not stop["now"] and len(running) < a.workers:
            u = todo.pop(0)
            out = os.path.join(a.out, f"U{u['i']:06d}.jsonl")
            tmp = f"{out}.part.{os.getpid()}"
            cmd = [a.bin, *u["args"][:-len(u["P"])], "--out", tmp, *map(str, u["P"])]
            p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                 start_new_session=True)
            running[p] = (u, tmp, out, time.time())
        if stop["now"]:
            for p in list(running):
                try:
                    os.killpg(p.pid, signal.SIGTERM)
                except OSError:
                    pass
        time.sleep(0.05)
        for p in list(running):
            c0 = children_cpu()
            if p.poll() is None:
                continue
            proc_cpu = children_cpu() - c0
            u, tmp, out, t_launch = running.pop(p)
            recs = done_file(tmp) if p.returncode == 0 else None
            if recs is None:
                if os.path.exists(tmp):
                    os.remove(tmp)
                if not stop["now"]:
                    say(f"unit {u['i']} failed (exit {p.returncode}); it will be retried next run")
                continue
            cpu = records_cpu(recs)
            wall = time.time() - t_launch
            cut = cut_by_clock(u, recs, proc_cpu)
            if cut:
                tries[u["i"]] = tries.get(u["i"], 0) + 1
                if tries[u["i"]] <= RETRIES:
                    # (its output is never kept: run again from the start;
                    # a magic square in it is announced all the same)
                    for r in recs:
                        if r.get("type") == "square" and r.get("best_score", 0) >= 14:
                            magic(r)
                    os.remove(tmp)
                    unit_row(u, cpu, wall, "cut")
                    if stop["now"]:
                        say(f"unit {u['i']}: stopped by the wall-clock limit after {cut[0]:.0f} "
                            f"CPU-s in {cut[1]:.0f} s: removed, it runs again next run")
                    else:
                        say(f"unit {u['i']}: stopped by the wall-clock limit after {cut[0]:.0f} "
                            f"CPU-s in {cut[1]:.0f} s (the machine slept or was busy): running "
                            f"it again")
                        todo.insert(0, u)
                    continue
                say(f"unit {u['i']}: stopped by the wall-clock limit {tries[u['i']]} times in "
                    f"this run (now after {cut[0]:.0f} CPU-s in {cut[1]:.0f} s): kept as it is, "
                    f"not complete")
            os.replace(tmp, out)
            done = [r for r in recs if r.get("type") == "done"][-1]
            unit_row(u, cpu, wall, "cut-kept" if cut else
                     "done" if done.get("complete", 1) else "incomplete")
            stats["units"] += 1
            stats["ref"] += u.get("pred_time", 0.0)
            stats["cpu"] += cpu
            recent.append((time.time(), cpu, u.get("pred_time", 0.0)))
            count_band(u, recs)
            for r in recs:
                if r.get("type") == "square":
                    stats["squares"] += 1
                    stats["sp"] += r.get("sp_count", 0) > 0
                    if r.get("best_score", 0) >= 14:
                        stats["magic"] += 1
                        magic(r)
        if time.time() - last > 60:
            last = time.time()
            recent = [x for x in recent if x[0] > last - WINDOW]
            say(progress(n_done, len(units), stats, t0, band if stage1 else None, recent))
    say("finished: " + progress(n_done, len(units), stats, t0, band if stage1 else None))


def progress(n_done, n, stats, t0, band, recent=None):
    now = time.time()
    el = (now - t0) / 3600
    ref = stats["ref"] / 3600
    msg = (f"{n_done + stats['units']}/{n} units done; this run: {stats['units']} units, "
           f"{stats['cpu'] / 3600:.2f} CPU-h in {el:.2f} h = {ref:.2f} reference CPU-h"
           + (f" ({stats['cpu'] / stats['ref']:.2f}x the reference per core)" if stats["ref"] else "")
           + f", {stats['squares']} squares, {stats['sp']} with an SP diagonal, "
           f"{stats['magic']} magic")
    if recent is not None and now - t0 >= 60:
        # (the units finished in the last WINDOW seconds: throttling, sleep
        # and other load show here)
        span = max(min(WINDOW, now - t0), 1.0)
        rc = sum(c for _, c, _ in recent)
        rr = sum(r for _, _, r in recent)
        msg += (f"; last {span / 60:.0f} min: {len(recent)} units, "
                f"{rr / span:.2f} reference CPU-h per hour"
                + (f" ({rc / rr:.2f}x the reference per core)" if rr else ""))
    if band is not None:
        msg += (f"; stage-1 band units so far: {band['units']}, pairs {band['pairs']} "
                f"(predicted {band['pred']:.1f} at f_rho = 1)")
    return msg


if __name__ == "__main__":
    main()
