#!/usr/bin/env python3
"""
Write the digest of the outputs of scripts/laptop/run.py: one small file
with everything scripts/decide.py needs, in place of the unit files.
Python 3.8+, standard library only.

usage: digest.py [DIR] [--host H] [--out FILE] [--plan PLAN]
                 (DIR: run.py's --out, default laptop-runs)

* Writes research/laptop/digest-<host>.json.xz: one cumulative file per
  host, replaced each time (written to a temporary file, then renamed).
  ~1-2 MB for the whole plan. The unit files stay in DIR. It prints the
  units it covers and the git command that sends the file.
* Contents (JSON, "version" 1): the plan's name and sha256 (from the meta
  files, or --plan for a DIR without them; a DIR with the outputs of two
  plans is refused), the host, when it was made, the meta files (machine,
  build_info text), run.py's speed (units.tsv's rows: i, cpu, wall, end,
  status; totals; the event lines of progress.log). Per finished unit (one
  with a "done" record), as columns: i, complete, sums, first and last S,
  the done record's min_sum, last_sum (less the first and last S) and
  time, squares, their S, P and SP traversals, squares with an SP
  traversal, best score, CPU, truncated sums, n, the file's sha256 and
  size (P is the plan line's). The per-sum rows of decide.py's state mode
  (it rebuilds the records' counts from them): per sum the S step and its
  squares, and for the sums with squares their S, P and SP traversals.
  Every notable square in full (an SP traversal, best_score >= 12, or
  magic) and the lines of MAGIC.txt.
* Counted as scripts/decide.py counts a unit file (decide.count_records):
  a square counts once its sum record is there; the squares of a sum not
  finished count nowhere (but a notable one is kept).
* decide.py reads it: decide.py --plan PLAN --digest FILE [FILE ...], with
  or without the scheduler state; its report is the one of DIR but for
  its first line.
* On the whole plan (2.34M sums, 0.84 GB of synthetic unit files): 1.8 MB,
  ~40 s and 0.3 GB of memory on an x86 core.
"""
import argparse
import hashlib
import json
import lzma
import os
import re
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
from collect import host_name  # noqa: E402
from decide import count_records  # noqa: E402

FORMAT = "magic-squares laptop digest"
VERSION = 1
UNIT = re.compile(r"U(\d+)\.jsonl$")
META = re.compile(r"meta_.*\.json$")
# a notable square: an SP traversal, or a best score of 12 or more (14: magic)
NOTABLE_BEST = 12
MAGIC_BEST = 14
# progress.log: run.py's periodic lines; the others are events
PERIODIC = " units done; this run: "
MAX_EVENTS = 2000

# the per-unit columns (decide.py reads them by these names)
UNIT_COLS = ("i", "complete", "sums", "S_first", "S_last", "min_sum", "last_sum", "done_time",
             "squares", "s", "p", "sp", "sq_sp", "best", "cpu", "truncated", "n", "sha256", "size")


def off(x, base):
    """x - base where both are numbers, else x"""
    if isinstance(x, (int, float)) and isinstance(base, (int, float)):
        return x - base
    return x


def read_unit(path):
    """(the unit's row, its per-sum rows, its notable squares, other record
    types) or None when it has no "done" record (not finished)"""
    with open(path, "rb") as f:
        data = f.read()
    recs = []
    for line in data.split(b"\n"):
        if not line.strip():
            continue
        try:
            r = json.loads(line)
        except ValueError:
            continue
        if isinstance(r, dict):
            recs.append(r)
    done = [r for r in recs if r.get("type") == "done"]
    if not done:
        return None
    o = count_records(recs)
    # per sum, as count_records pairs them: the squares that precede a sum
    # record with its S
    pend = []
    sums = []           # [S, squares, truncated, square records, s, p, sp]
    notable = []
    other = set()
    ns = set()
    Ps = set()
    sq_sp = 0
    for r in recs:
        t = r.get("type")
        if t in ("square", "sum", "done"):
            ns.add(r.get("n"))
            Ps.add(tuple(r.get("P", ())))
        if t == "square":
            pend.append(r)
            if (r.get("sp_count", 0) > 0 or r.get("best_score", 0) >= NOTABLE_BEST):
                notable.append(r)
        elif t == "sum":
            mine = [q for q in pend if q["S"] == r["S"]]
            pend = [q for q in pend if q["S"] != r["S"]]
            row = [int(r["S"]), int(r.get("squares", 0)), int(bool(r.get("truncated"))),
                   len(mine), 0, 0, 0]
            for q in mine:
                row[4] += q["s_count"]
                row[5] += q["p_count"]
                row[6] += q["sp_count"]
                sq_sp += q["sp_count"] > 0
            sums.append(row)
        elif t != "done":
            other.add(str(t))
    d = done[-1]
    row = {"complete": int(o.get("complete", 0)), "sums": len(sums),
           "S_first": sums[0][0] if sums else None, "S_last": sums[-1][0] if sums else None,
           # (the done record's min_sum and last_sum less the first and last
           # S: mostly 0)
           "min_sum": off(d.get("min_sum"), sums[0][0] if sums else None),
           "last_sum": off(d.get("last_sum"), sums[-1][0] if sums else None),
           "done_time": d.get("time"),
           "squares": o["squares"], "s": o["S"], "p": o["P"], "sp": o["pairs"], "sq_sp": sq_sp,
           "best": o.get("best", 0), "cpu": o["cpu"],
           "truncated": sum(x[2] for x in sums), "n": sorted(ns, key=str)[0] if ns else None,
           "sha256": hashlib.sha256(data).hexdigest(), "size": len(data)}
    if len(ns) > 1 or len(Ps) > 1:
        other.add("several n or P")
    return row, sums, notable, other


def speed(out_dir):
    """run.py's speed: units.tsv's rows as columns with totals, and the
    event lines of progress.log (with the periodic line before each)"""
    cols = {k: [] for k in ("i", "pred_time", "cpu", "wall", "end", "status")}
    try:
        with open(os.path.join(out_dir, "units.tsv")) as f:
            head = f.readline().rstrip("\n").split("\t")
            for line in f:
                x = line.rstrip("\n").split("\t")
                if len(x) != len(head):
                    continue
                r = dict(zip(head, x))
                try:
                    v = (int(r["i"]), float(r["pred_time"]), float(r["cpu"]), float(r["wall"]),
                         int(r["end"]), r["status"])
                except (KeyError, ValueError):
                    continue
                for k, x in zip(cols, v):
                    cols[k].append(x)
    except OSError:
        pass
    st = {}
    for s in cols["status"]:
        st[s] = st.get(s, 0) + 1
    tot = {"rows": len(cols["i"]), "cpu": sum(cols["cpu"]), "wall": sum(cols["wall"]),
           "pred_time": sum(cols["pred_time"]), "status": st}
    # (end times as steps from the previous row: they pack better)
    end = cols.pop("end")
    # (pred_time is the plan's: in the totals only)
    cols.pop("pred_time")
    cols["end_step"] = [e - (end[k - 1] if k else 0) for k, e in enumerate(end)]
    events = []
    try:
        with open(os.path.join(out_dir, "progress.log"), errors="replace") as f:
            last = None
            for line in f:
                line = line.rstrip("\n")
                if PERIODIC in line:
                    last = line
                    continue
                if last is not None:
                    events.append(last)
                    last = None
                events.append(line)
            if last is not None:
                events.append(last)
    except OSError:
        pass
    return {"units_tsv": cols, "totals": tot, "progress_events": events[-MAX_EVENTS:],
            "progress_events_left_out": max(0, len(events) - MAX_EVENTS)}


def plan_sha256(plan):
    with open(plan, "rb") as f:
        data = f.read()
    if plan.endswith(".xz"):
        data = lzma.decompress(data)
    return hashlib.sha256(data).hexdigest()


def build(out_dir, host, plan=None):
    """the digest of out_dir (a dict); plan: the plan file, for a
    directory without run.py's meta files"""
    metas = []
    builds = []
    for name in sorted(os.listdir(out_dir)):
        if not META.match(name):
            continue
        try:
            with open(os.path.join(out_dir, name)) as f:
                m = json.load(f)
        except (OSError, ValueError):
            sys.exit(f"{os.path.join(out_dir, name)}: not a meta file run.py wrote")
        m = dict(m, file=name)
        if "build" in m:
            if m["build"] not in builds:
                builds.append(m["build"])
            m["build"] = builds.index(m["build"])
        metas.append(m)
    plans = sorted({(m.get("plan"), m.get("plan_sha256")) for m in metas})
    if plan:
        given = (os.path.basename(plan), plan_sha256(plan))
        if plans and given[1] not in {p[1] for p in plans}:
            sys.exit(f"{out_dir} holds the outputs of {plans[0][0]} ({str(plans[0][1])[:12]}), "
                     f"not of --plan {plan} ({given[1][:12]})")
        plans = plans or [given]
    if not plans:
        sys.exit(f"{out_dir}: no meta_*.json (not run.py's --out?); --plan PLAN names its plan")
    if len({p[1] for p in plans}) > 1:
        sys.exit(f"{out_dir} holds the outputs of {len(plans)} plans: "
                 + ", ".join(f"{p[0]} ({str(p[1])[:12]})" for p in plans)
                 + ": a digest is of one plan's outputs")
    units = sorted(((int(m.group(1)), name) for name in os.listdir(out_dir)
                    for m in [UNIT.match(name)] if m))
    cols = {k: [] for k in UNIT_COLS}
    per = {"dS": [], "q": [], "truncated": [], "nrec": [], "s": [], "p": [], "sp": []}
    notable = []
    unfinished = []
    other = {}
    k = 0               # (the global index of a sum)
    for i, name in units:
        got = read_unit(os.path.join(out_dir, name))
        if got is None:
            unfinished.append(name)
            continue
        row, sums, nb, oth = got
        if oth:
            other[name] = sorted(oth)
        row["i"] = i
        for c in UNIT_COLS:
            cols[c].append(row[c])
        prev = None
        for S, q, tr, nrec, s, p, sp in sums:
            per["dS"].append(0 if prev is None else S - prev)
            prev = S
            per["q"].append(q)
            if tr:
                per["truncated"].append(k)
            if nrec != q:
                per["nrec"].append([k, nrec])
            if nrec:
                per["s"].append(s)
                per["p"].append(p)
                per["sp"].append(sp)
            k += 1
        for r in nb:
            notable.append(dict(r, i=i))
    magic = []
    try:
        with open(os.path.join(out_dir, "MAGIC.txt"), errors="replace") as f:
            magic = [x.rstrip("\n") for x in f if x.strip()]
    except OSError:
        pass
    for r in notable:
        r["magic"] = int(r.get("best_score", 0) >= MAGIC_BEST)
    d = {"format": FORMAT, "version": VERSION, "plan": plans[0][0], "plan_sha256": plans[0][1],
         "host": host, "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
         "created_unix": time.time(), "meta": metas, "builds": builds, "speed": speed(out_dir),
         "units": cols,
         "persum": dict(per, note="per sum of the units in order: dS (S minus the previous "
                        "sum's S, 0 at a unit's first sum, whose S is S_first), q (its record's "
                        "squares); truncated: the sums' indices; nrec: [index, square records] "
                        "where that is not q; s, p, sp: the S, P and SP traversals of the "
                        "squares of each sum with square records"),
         "notable": notable, "magic": magic, "unfinished": unfinished, "other_records": other}
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("dir", nargs="?", default=os.path.join(ROOT, "laptop-runs"))
    ap.add_argument("--host", default=None,
                    help="the host's name in the file name (default this machine's)")
    ap.add_argument("--out", default=None,
                    help="the digest (default research/laptop/digest-<host>.json.xz)")
    ap.add_argument("--plan", default=None,
                    help="the plan run (checked against the meta files; needed only for a "
                         "directory without them)")
    a = ap.parse_args()
    if not os.path.isdir(a.dir):
        sys.exit(f"{a.dir}: no such directory (run.py's --out)")
    host = a.host or host_name()
    out = a.out or os.path.join(ROOT, "research", "laptop", f"digest-{host}.json.xz")
    t0 = time.time()
    d = build(a.dir, host, a.plan)
    if d["other_records"]:
        names = sorted(d["other_records"])
        sys.exit(f"{len(names)} unit files hold records other than square, sum and done "
                 f"({names[0]}: {', '.join(d['other_records'][names[0]])}): the digest is of "
                 f"plain units (use collect.sh)")
    raw = json.dumps(d, separators=(",", ":")).encode()
    data = lzma.compress(raw, preset=9 | lzma.PRESET_EXTREME)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    tmp = out + ".tmp"
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, out)
    u = d["units"]
    nu = len(u["i"])
    ranges = []
    for i in u["i"]:
        if ranges and i == ranges[-1][1] + 1:
            ranges[-1][1] = i
        else:
            ranges.append([i, i])
    span = ", ".join(f"{x}-{y}" if y > x else f"{x}" for x, y in ranges[:8]) + (
        ", ..." if len(ranges) > 8 else "")
    rel = os.path.relpath(out, ROOT)
    if rel.startswith(".."):
        rel = out
    size = f"{len(data) / 1e6:.2f} MB" if len(data) >= 1e5 else f"{len(data) / 1e3:.1f} kB"
    print(f"{rel}: {size}, plan {d['plan']} ({str(d['plan_sha256'])[:12]}), "
          f"made in {time.time() - t0:.0f} s")
    print(f"units: {nu} finished ({nu - sum(u['complete'])} not complete): {span}; "
          f"{sum(u['sums'])} sums, {sum(u['squares'])} squares, {sum(u['sq_sp'])} with an SP "
          f"traversal, {sum(1 for r in d['notable'] if r['magic'])} magic")
    if d["unfinished"]:
        print(f"left out: {len(d['unfinished'])} unit files without a done record "
              f"({' '.join(d['unfinished'][:5])}{' ...' if len(d['unfinished']) > 5 else ''})")
    print(f"to send it: git add {rel} && git commit -m 'laptop digest' && git push")


if __name__ == "__main__":
    main()
