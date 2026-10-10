#!/usr/bin/env python3
"""
Write a static, pre-registered plan for a search on a machine without the
scheduler (e.g. a laptop): the units scheduler v2 would launch first, in its
order, with the model's predictions frozen next to them.

usage: make_plan.py --state DIR --hours H --out plan.jsonl
           [--time-factor F] [--node-factor F] [--dfirst off|auto]

Each line: {"i", "P", "lo", "hi", "mode", "args" (msearch arguments without
the binary and --out), "pred_time" (CPU-s of the fast x86 build),
"pred_squares", "pred_magic", "score"}. --time-limit and --node-limit are
multiplied by --time-factor and --node-factor, because a slower machine
(e.g. an ARM laptop: the same search with the portable kernels, ~2.4x the
fast build's time per core on x86 without AVX-512; the node factor was for
the matrix path's extra nodes, which builds without AVX-512 no longer take)
would otherwise stop units early. Prints the sha256 of the plan for the pre-registration.
"""
import argparse
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import scheduler as sch  # noqa: E402


def scheduler_args(state, dfirst, unit_time):
    """the scheduler's own argument namespace for `emit` (its defaults)"""
    got = {}
    sch.v2_emit = lambda a: got.setdefault("args", a)
    argv = sys.argv
    sys.argv = ["scheduler.py", "--state", state, "emit", "--units", "1",
                "--dfirst", dfirst, "--unit-time", str(unit_time)]
    try:
        sch.main()
    finally:
        sys.argv = argv
    return got["args"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--hours", type=float, default=40.0,
                    help="predicted CPU-hours of the fast x86 build to plan")
    ap.add_argument("--out", required=True)
    ap.add_argument("--time-factor", type=float, default=12.0)
    ap.add_argument("--node-factor", type=float, default=4.0)
    ap.add_argument("--dfirst", default="off")
    ap.add_argument("--unit-time", type=float, default=120)
    a = ap.parse_args()

    args = scheduler_args(a.state, a.dfirst, a.unit_time)
    s = sch.SchedulerV2(args, quiet=True)
    lines = []
    for i, u in enumerate(s.simulate(max_hours=a.hours)):
        cmd = s.command(u, "OUT", check=False)[1:]
        k = cmd.index("--out")
        del cmd[k:k + 2]
        k = cmd.index("--time-limit")
        cmd[k + 1] = f"{float(cmd[k + 1]) * a.time_factor:.0f}"
        k = cmd.index("--node-limit")
        cmd[k + 1] = str(int(int(cmd[k + 1]) * a.node_factor))
        lines.append(json.dumps({
            "i": i, "P": list(u.P), "lo": u.lo, "hi": u.hi, "mode": u.mode,
            "args": cmd, "pred_time": round(u.time, 3),
            "pred_squares": u.squares, "pred_magic": u.magic, "score": u.score},
            sort_keys=True))
    data = ("\n".join(lines) + "\n").encode()
    with open(a.out, "wb") as f:
        f.write(data)
    tot = sum(json.loads(x)["pred_time"] for x in lines) / 3600
    mag = sum(json.loads(x)["pred_magic"] for x in lines)
    print(f"{len(lines)} units, {tot:.1f} predicted CPU-hours (fast x86 build), "
          f"E[magic] {mag:.4f}")
    print("sha256", hashlib.sha256(data).hexdigest())


if __name__ == "__main__":
    main()
