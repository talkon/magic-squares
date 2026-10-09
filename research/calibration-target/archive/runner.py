#!/usr/bin/env python3
"""Runner of the pre-registered calibration search (run.sh feeds it the
commands, one per line:  ID pred=SECONDS limit=SECONDS :: COMMAND...).

* at most --jobs (2) processes at a time, each niced (the commands start
  with nice -n 19; the runner also lowers its own children's priority);
* every process has RLIMIT_CPU = its limit (SIGXCPU then SIGKILL);
* hard budget: before each launch, CPU spent (ledger) + the CPU limits of
  the running processes + this run's limit <= --cap-hours (8); the limit is
  shrunk to fit, and the run is skipped if that leaves less than 1.25 x its
  prediction (reported as "budget");
* resumable: a run whose output already has a "done" record is skipped; a
  run the ledger records as killed (CPU limit) or skipped for budget is not
  retried (--retry-killed to retry); an output without "done" and without a
  ledger entry (the runner was interrupted) is moved aside to
  OUT.partial.N and the run starts again (msearch --out appends);
* runs/ledger.jsonl: one line per finished, killed or interrupted process
  (id, status, cpu, wall, rc); CPU is the child's rusage (user + sys);
* runs/.lock: one runner at a time.
"""
import argparse
import fcntl
import json
import os
import resource
import shlex
import signal
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def parse(lines):
    runs = []
    for ln in lines:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        head, cmd = ln.split("::", 1)
        f = head.split()
        rid = f[0]
        kv = dict(x.split("=", 1) for x in f[1:])
        argv = shlex.split(cmd)
        out = argv[argv.index("--out") + 1]
        runs.append({"id": rid, "pred": float(kv["pred"]), "limit": float(kv["limit"]), "argv": argv,
                     "out": out})
    return runs


def has_done(path):
    if not os.path.exists(path):
        return False
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b'{"type":"done"'):
                return True
    return False


def read_ledger(path):
    led = []
    if os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        led.append(json.loads(line))
                    except ValueError:
                        pass
    return led


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs-dir", default=os.path.join(HERE, "runs"))
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--cap-hours", type=float, default=8.0)
    ap.add_argument("--retry-killed", action="store_true")
    ap.add_argument("--dry-run", action="store_true", help="print what would run, run nothing")
    ap.add_argument("commands", nargs="?", default="-")
    a = ap.parse_args()
    a.jobs = max(1, min(a.jobs, 2))      # the machine is shared: never more than 2
    os.makedirs(a.runs_dir, exist_ok=True)
    lock = open(os.path.join(a.runs_dir, ".lock"), "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        sys.exit("another runner holds runs/.lock")
    src = sys.stdin if a.commands == "-" else open(a.commands)
    runs = parse(src.read().splitlines())
    ledger_path = os.path.join(a.runs_dir, "ledger.jsonl")
    ledger = read_ledger(ledger_path)
    cap = a.cap_hours * 3600.0
    spent = sum(float(e.get("cpu", 0.0)) for e in ledger)
    final = {}
    for e in ledger:
        if e.get("status") in ("killed", "budget", "done", "failed"):
            final[e["id"]] = e["status"]

    def log(msg):
        print(time.strftime("[%H:%M:%S] ") + msg, file=sys.stderr, flush=True)

    def record(e):
        with open(ledger_path, "a") as f:
            f.write(json.dumps(e) + "\n")

    todo = []
    for r in runs:
        if has_done(r["out"]):
            continue
        st = final.get(r["id"])
        if st in ("killed", "failed") and not a.retry_killed:
            continue
        if st == "budget":
            continue
        todo.append(r)
    log(f"{len(runs)} runs, {len(todo)} to do; CPU spent so far {spent / 3600:.3f} h of "
        f"{a.cap_hours} h; predicted for the rest {sum(r['pred'] for r in todo) / 3600:.3f} h")
    if a.dry_run:
        for r in todo:
            print(r["id"], f"pred {r['pred']:.0f} s limit {r['limit']:.0f} s:", shlex.join(r["argv"]))
        return
    running = {}   # pid -> (run, limit, t0)
    stop = {"flag": False}

    def on_signal(sig, frame):
        stop["flag"] = True
        for pid in list(running):
            try:
                os.kill(pid, signal.SIGTERM)
            except OSError:
                pass

    signal.signal(signal.SIGINT, on_signal)
    signal.signal(signal.SIGTERM, on_signal)

    def reap(block=True):
        nonlocal spent
        try:
            pid, status, ru = os.wait4(-1, 0 if block else os.WNOHANG)
        except ChildProcessError:
            return False
        if pid == 0:
            return False
        r, limit, t0 = running.pop(pid)
        cpu = ru.ru_utime + ru.ru_stime
        spent += cpu
        sig = os.WTERMSIG(status) if os.WIFSIGNALED(status) else 0
        rc = os.WEXITSTATUS(status) if os.WIFEXITED(status) else None
        if has_done(r["out"]):
            st = "done"
        elif stop["flag"]:
            st = "interrupted"
        elif sig in (signal.SIGXCPU, signal.SIGKILL) or cpu >= limit - 1:
            st = "killed"
        else:
            st = "failed"
        e = {"id": r["id"], "status": st, "cpu": cpu, "wall": time.time() - t0, "rc": rc, "signal": sig,
             "limit": limit, "pred": r["pred"], "t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
        record(e)
        log(f"{r['id']} {st}: cpu {cpu:.1f} s (pred {r['pred']:.1f}), spent {spent / 3600:.3f} h")
        return True

    for r in todo:
        if stop["flag"]:
            break
        while len(running) >= a.jobs:
            reap()
        if stop["flag"]:
            break
        reserved = sum(lim for _, lim, _ in running.values())
        room = cap - spent - reserved
        need = min(r["limit"], 1.25 * r["pred"])
        limit = min(r["limit"], room)
        if limit < need:
            # wait for the running processes first: their unused limits come back
            while running and limit < need:
                reap()
                reserved = sum(lim for _, lim, _ in running.values())
                limit = min(r["limit"], cap - spent - reserved)
            if limit < need:
                record({"id": r["id"], "status": "budget", "cpu": 0.0, "pred": r["pred"],
                        "room": cap - spent, "t": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
                log(f"{r['id']} skipped: budget ({(cap - spent) / 3600:.3f} h left)")
                continue
        if os.path.exists(r["out"]):
            i = 1
            while os.path.exists(f"{r['out']}.partial.{i}"):
                i += 1
            os.rename(r["out"], f"{r['out']}.partial.{i}")
        lim = int(limit)

        def pre(lim=lim):
            os.setsid()
            os.nice(19)
            resource.setrlimit(resource.RLIMIT_CPU, (lim, lim + 10))

        errp = r["out"][:-len(".jsonl")] + ".err" if r["out"].endswith(".jsonl") else r["out"] + ".err"
        with open(errp, "ab") as err:
            pid = os.fork()
            if pid == 0:
                try:
                    pre()
                    os.dup2(err.fileno(), 2)
                    devnull = os.open(os.devnull, os.O_RDWR)
                    os.dup2(devnull, 0)
                    os.dup2(devnull, 1)
                    os.execvp(r["argv"][0], r["argv"])
                finally:
                    os._exit(127)
        running[pid] = (r, float(lim), time.time())
        log(f"{r['id']} started (pred {r['pred']:.0f} s, CPU limit {lim} s): {shlex.join(r['argv'][3:])}")
    while running:
        reap()
    log(f"finished; CPU spent {spent / 3600:.3f} h")


if __name__ == "__main__":
    main()
