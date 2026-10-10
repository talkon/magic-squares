#!/usr/bin/env python3
"""
Run a static plan (scripts/laptop/make_plan.py) on any machine, in the plan's
order, with K parallel msearch processes. Python 3.8+, standard library only.

usage: run.py PLAN(.jsonl or .jsonl.xz) [--bin BUILD/msearch] [--out DIR] [--workers K]
              [--hours H] [--max-units U]

* Resumable: a unit whose output has a "done" record is skipped; a unit that
  was interrupted is run again from the start (its output is replaced).
* Ctrl-C (or --hours) stops launching, stops the running units and leaves
  everything finished in place. Run it again to continue.
* Outputs: DIR/U000123.jsonl per unit, DIR/meta.json (machine, build, plan
  hash), DIR/progress.log. A magic square (best_score >= 14) is announced and
  copied to DIR/MAGIC.txt.
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
        print(line, flush=True)
        log.write(line + "\n")
        log.flush()

    todo = []
    n_done = 0
    for u in units:
        if done_file(os.path.join(a.out, f"U{u['i']:06d}.jsonl")):
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

    t0 = time.time()
    running = {}  # Popen -> (unit, tmp path)
    stats = {"units": 0, "cpu": 0.0, "squares": 0, "sp": 0, "magic": 0}
    last = 0.0
    while (todo or running):
        if a.hours and time.time() - t0 > a.hours * 3600 and not stop["now"]:
            stop["now"] = True
            say(f"--hours {a.hours} reached: stopping")
        while todo and not stop["now"] and len(running) < a.workers:
            u = todo.pop(0)
            out = os.path.join(a.out, f"U{u['i']:06d}.jsonl")
            tmp = out + ".part"
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
            for r in recs:
                if r.get("type") == "sum":
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
            el = (time.time() - t0) / 3600
            say(f"{n_done + stats['units']}/{len(units)} units done; this run: "
                f"{stats['units']} units, {stats['cpu'] / 3600:.2f} CPU-h in {el:.2f} h, "
                f"{stats['squares']} squares, {stats['sp']} with an SP diagonal, "
                f"{stats['magic']} magic")
    el = (time.time() - t0) / 3600
    say(f"finished: {n_done + stats['units']}/{len(units)} units done; this run "
        f"{stats['units']} units, {stats['cpu'] / 3600:.2f} CPU-h in {el:.2f} h, "
        f"{stats['squares']} squares, {stats['sp']} with an SP diagonal, {stats['magic']} magic")


if __name__ == "__main__":
    main()
