#!/usr/bin/env python3
"""
Run a static plan (scripts/laptop/make_plan.py) on any machine, in the plan's
order, with K parallel msearch processes. Python 3.8+, standard library only.

usage: run.py PLAN(.jsonl or .jsonl.xz) [--bin BUILD/msearch] [--out DIR] [--workers K]
              [--hours H] [--max-units U]

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
  predictions on a busy shared machine), and for a stage-1 plan the pairs
  found in its band units (research/stage1.md; scripts/decide.py fits them).
"""
import argparse
import hashlib
import json
import lzma
import os
import platform
import signal
import subprocess
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
    for u in units:
        path = os.path.join(a.out, f"U{u['i']:06d}.jsonl")
        recs = done_file(path)
        if recs:
            if not unit_matches(u, recs):
                sys.exit(f"{path} is not the output of line {u['i']} of {a.plan} (another "
                         f"plan's?): use another --out")
            n_done += 1
        else:
            todo.append(u)
    say(f"plan {meta['plan']} ({meta['plan_sha256'][:12]}): {len(units)} units, "
        f"{n_done} already done, {len(todo)} to run with {a.workers} workers")

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
    running = {}  # Popen -> (unit, tmp path)
    stats = {"units": 0, "cpu": 0.0, "ref": 0.0, "squares": 0, "sp": 0, "magic": 0}
    # (a stage-1 plan: pairs in its band units, all finished ones so far)
    band = {"units": 0, "pairs": 0, "pred": 0.0}
    stage1 = any("stage1" in u for u in units)

    def count_band(u, recs):
        if u.get("stage1"):
            band["units"] += 1
            band["pred"] += u.get("pred_pairs", 0.0)
            band["pairs"] += sum(r.get("sp_count", 0) for r in recs if r.get("type") == "square")
    for u in units:
        if stage1 and u.get("stage1"):
            recs = done_file(os.path.join(a.out, f"U{u['i']:06d}.jsonl"))
            if recs:
                count_band(u, recs)
    last = 0.0
    while (todo and not stop["now"]) or running:
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
            running[p] = (u, tmp, out)
        if stop["now"]:
            for p in list(running):
                try:
                    os.killpg(p.pid, signal.SIGTERM)
                except OSError:
                    pass
        time.sleep(0.05)
        for p in list(running):
            if p.poll() is None:
                continue
            u, tmp, out = running.pop(p)
            recs = done_file(tmp) if p.returncode == 0 else None
            if recs is None:
                if os.path.exists(tmp):
                    os.remove(tmp)
                if not stop["now"]:
                    say(f"unit {u['i']} failed (exit {p.returncode}); it will be retried next run")
                continue
            os.replace(tmp, out)
            stats["units"] += 1
            stats["ref"] += u.get("pred_time", 0.0)
            count_band(u, recs)
            for r in recs:
                if r.get("type") in ("sum", "csum", "dsum"):
                    stats["cpu"] += r.get("cpu", r.get("time", 0.0))
                elif r.get("type") == "square":
                    stats["squares"] += 1
                    stats["sp"] += r.get("sp_count", 0) > 0
                    if r.get("best_score", 0) >= 14:
                        stats["magic"] += 1
                        msg = json.dumps(r)
                        say("*** MAGIC SQUARE FOUND *** " + msg)
                        with open(os.path.join(a.out, "MAGIC.txt"), "a") as f:
                            f.write(msg + "\n")
        if time.time() - last > 60:
            last = time.time()
            say(progress(n_done, len(units), stats, t0, band if stage1 else None))
    say("finished: " + progress(n_done, len(units), stats, t0, band if stage1 else None))


def progress(n_done, n, stats, t0, band):
    el = (time.time() - t0) / 3600
    ref = stats["ref"] / 3600
    msg = (f"{n_done + stats['units']}/{n} units done; this run: {stats['units']} units, "
           f"{stats['cpu'] / 3600:.2f} CPU-h in {el:.2f} h = {ref:.2f} reference CPU-h"
           + (f" ({stats['cpu'] / stats['ref']:.2f}x the reference per core)" if stats["ref"] else "")
           + f", {stats['squares']} squares, {stats['sp']} with an SP diagonal, "
           f"{stats['magic']} magic")
    if band is not None:
        msg += (f"; stage-1 band units so far: {band['units']}, pairs {band['pairs']} "
                f"(predicted {band['pred']:.1f} at f_rho = 1)")
    return msg


if __name__ == "__main__":
    main()
