#!/usr/bin/env python3
"""
Write a static, pre-registered plan for a search on a machine without the
scheduler (e.g. a laptop or rented machines): the units scheduler v2 would
launch first, in its order, with the model's predictions frozen next to them.

usage: make_plan.py --state DIR --hours H --out plan.jsonl
           [--time-factor F] [--node-factor F] [--dfirst off|auto]
       make_plan.py --state DIR --stage1 [HOURS[:LO:HI]] --out plan.jsonl
           [--prereg PREREGISTERED.txt] [--decision-table FILE]

Each line: {"i", "P", "lo", "hi", "mode", "args" (msearch arguments without
the binary and --out), "pred_time" (CPU-s of the fast x86 build),
"pred_squares", "pred_magic", "score"}. --time-limit and --node-limit are
multiplied by --time-factor and --node-factor, because a slower machine
would otherwise stop units early. An ARM laptop runs the same search with
the portable kernels: on x86 without AVX-512 they take 2.5-2.6x the fast
build's time per core on plain sums with clang -march=x86-64-v3 and
3.2-3.7x with 128-bit vectors (NEON's width), so an M1 is expected at
~2.5-3.5x (scripts/laptop/build.sh measures it). The scheduler's time
limit is max(2 x unit time, 3 x predicted) of the fast build; the default
factor 12 keeps 3.4-4.8x of that margin at 2.5-3.5x (a factor of ~6 would
do at 3.5x) and still covers a core 6-10x slower (an efficiency core, a
throttled laptop); the cost is only that a runaway unit runs longer. The
node factor was for the matrix path's extra nodes, which builds without
AVX-512 no longer take; the time limit binds first. Prints the sha256 of
the plan for the pre-registration.

--stage1 (research/stage1.md): the plan of scheduler.py run --stage1
HOURS:LO:HI (default 60:3000:6000; --dfirst auto unless given): the first
HOURS of reference CPU (--hours defaults to HOURS), sums with N' in [LO, HI)
plain. Each line also holds the model's predictions of its records at the
shipped-or-learned calibration of the state and f_rho = 1
(scheduler.unit_predictions): pred_S, pred_P (traversals), pred_SP (SP
traversals before the class factors), pred_pairs ((square, SP traversal)
pairs, scripts/decide.py's model), pred_pairs_dfirst (the d loop's, not in
decide.py), lNp (ln N' of its first and last sum) and stage1 (the unit is
in the stage-1 band). The totals, the sha256 of the plan and of the
decision table go into --prereg (default: next to the plan), to be
committed before the first unit runs.
"""
import argparse
import datetime
import hashlib
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import scheduler as sch  # noqa: E402


def scheduler_args(state, dfirst, unit_time, stage1=None):
    """the scheduler's own argument namespace for `emit` (its defaults)"""
    got = {}
    emit = sch.v2_emit
    sch.v2_emit = lambda a: got.setdefault("args", a)
    argv = sys.argv
    sys.argv = (["scheduler.py", "--state", state, "emit", "--units", "1",
                 "--dfirst", dfirst, "--unit-time", str(unit_time)]
                + (["--stage1", stage1] if stage1 else []))
    try:
        sch.main()
    finally:
        sys.argv = argv
        sch.v2_emit = emit
    return got["args"]


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def git_head():
    try:
        return subprocess.run(["git", "-C", HERE, "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--hours", type=float, default=None,
                    help="predicted CPU-hours of the fast x86 build to plan (default 40, "
                         "with --stage1 its HOURS)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--time-factor", type=float, default=12.0)
    ap.add_argument("--node-factor", type=float, default=4.0)
    ap.add_argument("--dfirst", default=None, help="default off (auto with --stage1)")
    ap.add_argument("--unit-time", type=float, default=120)
    ap.add_argument("--stage1", nargs="?", const=":".join(f"{v:g}" for v in sch.STAGE1_DEFAULT),
                    default=None, metavar="HOURS[:LO:HI]")
    ap.add_argument("--prereg", default=None,
                    help="--stage1: the pre-registration file (default PLAN.PREREGISTERED.txt)")
    ap.add_argument("--decision-table", default=os.path.join(
        os.path.dirname(os.path.dirname(HERE)), "research", "stage1", "decision-table.json"),
        help="--stage1: the decision table frozen with the plan (its sha256 is registered)")
    a = ap.parse_args()

    spec = sch.parse_stage1(a.stage1) if a.stage1 else None
    dfirst = a.dfirst or ("auto" if spec else "off")
    hours = a.hours if a.hours is not None else (spec[0] if spec else 40.0)
    args = scheduler_args(a.state, dfirst, a.unit_time, a.stage1)
    s = sch.SchedulerV2(args, quiet=True)
    if spec and s.stage1.spent > 0:
        sys.exit(f"the state's stage 1 has run {s.stage1.spent / 3600:.2f} h: a static stage-1 "
                 f"plan starts from a state without it")
    lines = []
    band_pairs = []   # (--stage1) per unit: the pairs of its sums in the band
    band = (math.log(spec[1]), math.log(spec[2])) if spec else None
    for i, u in enumerate(s.simulate(max_hours=hours)):
        cmd = s.command(u, "OUT", check=False)[1:]
        k = cmd.index("--out")
        del cmd[k:k + 2]
        k = cmd.index("--time-limit")
        cmd[k + 1] = f"{float(cmd[k + 1]) * a.time_factor:.0f}"
        k = cmd.index("--node-limit")
        cmd[k + 1] = str(int(int(cmd[k + 1]) * a.node_factor))
        rec = {"i": i, "P": list(u.P), "lo": u.lo, "hi": u.hi, "mode": u.mode,
               "args": cmd, "pred_time": round(u.time, 3),
               "pred_squares": u.squares, "pred_magic": u.magic, "score": u.score}
        if spec:
            # (s.simulate set the scorer's stage-1 band for this unit)
            p = sch.unit_predictions(s.scorer, u, band=band)
            # (not in the plan's lines, whose format is frozen)
            band_pairs.append(p.pop("pred_pairs_band"))
            p["pred_squares"] = u.squares
            rec.update(p)
            rec["stage1"] = bool(band[0] <= p["lNp"][0] < band[1])
        lines.append(json.dumps(rec, sort_keys=True))
    data = ("\n".join(lines) + "\n").encode()
    with open(a.out, "wb") as f:
        f.write(data)
    recs = [json.loads(x) for x in lines]
    tot = sum(r["pred_time"] for r in recs) / 3600
    mag = sum(r["pred_magic"] for r in recs)
    print(f"{len(lines)} units, {tot:.1f} predicted CPU-hours (fast x86 build), "
          f"E[magic] {mag:.4f}")
    print("sha256", sha256(data))
    if not spec:
        return

    # the pre-registration: totals per N' band and mode
    def band_name(r):
        return sch.NBAND_NAMES[sum(r["lNp"][0] >= math.log(e) for e in sch.NBAND_EDGES)]

    keys = ("pred_time", "pred_squares", "pred_S", "pred_P", "pred_SP", "pred_pairs",
            "pred_pairs_dfirst", "pred_magic")
    by = {}
    for r in recs:
        g = by.setdefault(f"{band_name(r)} {r['mode']}", dict.fromkeys(keys, 0.0) | {"units": 0})
        g["units"] += 1
        for k in keys:
            g[k] += r[k]
    totals = {k: sum(r[k] for r in recs) for k in keys} | {"units": len(recs)}
    inb = [r for r in recs if r["stage1"]]
    tb = {k: sum(r[k] for r in inb) for k in keys} | {"units": len(inb)}
    # reference CPU-hours of the plan at which the band's pairs reach 45:
    # counting the band units only (the units whose first sum is in the
    # band, what decide.py fits without a state), and counting every sum in
    # the band (also the later sums of units that start below LO, what
    # decide.py fits with the state)
    cum = cums = h45 = h45s = 0.0
    h = 0.0
    for r, bp in zip(recs, band_pairs):
        h += r["pred_time"] / 3600
        cums += bp
        if cums >= 45 and not h45s:
            h45s = h
        if r["stage1"]:
            cum += r["pred_pairs"]
            if cum >= 45 and not h45:
                h45 = h
    across = [(r, bp) for r, bp in zip(recs, band_pairs) if bp > 0 and not r["stage1"]]
    nb = [math.exp(x) for r in inb for x in r["lNp"]]
    try:
        with open(a.decision_table, "rb") as f:
            dt = f.read()
        dt_sha = sha256(dt)
    except OSError:
        dt_sha = None
    pre = a.prereg or a.out + ".PREREGISTERED.txt"
    fmt = lambda g: (f"{g['units']:6d} units {g['pred_time'] / 3600:7.2f} h  squares "
                     f"{g['pred_squares']:10.0f}  S {g['pred_S']:9.0f}  P {g['pred_P']:8.0f}  "
                     f"SP {g['pred_SP']:7.1f}  pairs {g['pred_pairs']:7.2f}  "
                     f"(d-first loops {g['pred_pairs_dfirst']:.2f})  magic {g['pred_magic']:.4f}")
    txt = [f"Pre-registration of a stage-1 plan ({datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}, "
           f"before any of its units ran)", "",
           f"plan: {os.path.basename(a.out)}", f"sha256 of the plan: {sha256(data)}",
           f"decision table: {a.decision_table} sha256 {dt_sha}",
           f"code: {git_head()} (scripts/scheduler.py, scripts/amodel.py, scripts/decide.py)",
           f"model: scheduler v2, --stage1 {':'.join(f'{v:g}' for v in spec)}, --dfirst {dfirst}, "
           f"--unit-time {a.unit_time:g}, the state's learned calibration; predictions at f_rho = 1 "
           f"(pairs = sp_kappa r_S r_P 720 e^lpSP per square, scheduler.unit_predictions)",
           f"state: {os.path.abspath(a.state)} ({len(s.summary.files)} unit files, "
           f"{s.summary.totals['squares']} squares)", "",
           "totals:", "  all      " + fmt(totals),
           f"  stage-1 band (N' {spec[1]:g}-{spec[2]:g}) " + fmt(tb), "by N' band and mode:"]
    txt += [f"  {k:13} " + fmt(g) for k, g in sorted(by.items())]
    txt += ["", f"N' of the band units' sums: {min(nb, default=0):.0f}-{max(nb, default=0):.0f}",
            f"pairs predicted in the band's sums: {sum(band_pairs):.2f} ({tb['pred_pairs']:.2f} in "
            f"the band units, {sum(bp for _, bp in across):.2f} in the later sums of "
            f"{len(across)} units that start below N' {spec[1]:g})",
            f"45 pairs predicted in the stage-1 band after {h45:.1f} reference CPU-hours of the "
            f"plan counting the band units only, after {h45s:.1f} counting every sum in the band "
            f"(0 = not within it)",
            "", "analysis plan (scripts/decide.py STATE --plan PLAN): observed / predicted squares, "
                "S and P traversals, pairs per unit run (a prefix of the plan, in order), pooled and "
                "by N' band, with Poisson intervals; the f_rho posterior from the band's plain "
                "pairs (prior: the calibration search's posterior, median 1.06, f_rho^2 ln sd "
                "0.47); the predictive E and P(>=1) at 1 / 3 / 10 CPU-years; the decision table's "
                "rows, evaluated, never acted on by the tool."]
    with open(pre, "w") as f:
        f.write("\n".join(txt) + "\n")
    print(f"pre-registration: {pre}")
    print("\n".join(txt[10:]))


if __name__ == "__main__":
    main()
