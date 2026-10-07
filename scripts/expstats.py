#!/usr/bin/env python3
"""Summarize msearch JSON-lines output by number of vectors N.

usage: expstats.py [--by P|N] [--bins B1,B2,...] FILE...

For each bin of N (nvecs_raw), prints: sums searched, CPU time, squares,
squares per sum, squares per CPU-hour, and the traversal probabilities of
the squares found (p_S, p_P, p_SP per traversal), plus the implied
"magic mass" per CPU-hour = squares/h * 5400 * p_SP^2 (p_SP from the
observed rates, falling back to rho * p_S * p_P with rho = 2 when no SP
traversal was seen in the bin).
"""
import argparse
import json
import math
import sys


def load(files):
    sums, squares = [], []
    for path in files:
        with open(path) as f:
            for line in f:
                if not line.endswith("\n"):
                    break
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r["type"] == "sum":
                    sums.append(r)
                elif r["type"] == "square":
                    squares.append(r)
    return sums, squares


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--by", default="N", choices=["N", "P", "PN"])
    ap.add_argument("--bins", default="400,700,1000,1400,2000,2800,4000,5600,8000,11000,16000")
    ap.add_argument("files", nargs="+")
    a = ap.parse_args()
    bins = [int(x) for x in a.bins.split(",")]
    sums, squares = load(a.files)

    def nbin(N):
        for i, b in enumerate(bins):
            if N < b:
                return i
        return len(bins)

    def nlabel(i):
        lo = bins[i - 1] if i > 0 else 0
        hi = bins[i] if i < len(bins) else "inf"
        return f"{lo}-{hi}"

    def key(r):
        P = "_".join(map(str, r["P"]))
        if a.by == "N":
            return nbin(r["nvecs_raw"])
        if a.by == "P":
            return P
        return (P, nbin(r["nvecs_raw"]))

    # per (P, S): N for squares
    N_of = {(tuple(r["P"]), r["S"]): r["nvecs_raw"] for r in sums}
    acc = {}
    for r in sums:
        d = acc.setdefault(key(r), dict(sums=0, time=0.0, sq=0, nodes=0, s=0, p=0, sp=0, nsq=0,
                                        sptype=0))
        d["sums"] += 1
        d["time"] += r["time"] + r["setup_time"]
        d["sq"] += r["squares"]
        d["nodes"] += r["nodes"]
    for q in squares:
        N = N_of.get((tuple(q["P"]), q["S"]))
        if N is None:
            continue
        r = {"P": q["P"], "nvecs_raw": N}
        d = acc.get(key(r))
        if d is None:
            continue
        d["nsq"] += 1
        d["s"] += q["s_count"]
        d["p"] += q["p_count"]
        d["sp"] += q["sp_count"]
        d["sptype"] += q["sp_count"] > 0

    def label(k):
        if a.by == "N":
            return nlabel(k)
        if a.by == "P":
            return k
        return f"{k[0]} {nlabel(k[1])}"

    print(f"{'bin':24} {'sums':>5} {'cpu-h':>7} {'sq':>6} {'sq/sum':>7} {'sq/h':>8} {'s/sum(s)':>8} "
          f"{'p_S':>8} {'p_P':>8} {'p_SP':>8} {'#SP':>4} {'mass/h':>9}")
    tot_mass = 0.0
    for k in sorted(acc, key=lambda k: (str(k[0]) if isinstance(k, tuple) else str(k) if a.by == "P" else k,
                                        k[1] if isinstance(k, tuple) else 0)):
        d = acc[k]
        h = d["time"] / 3600
        T = 720 * d["nsq"]
        pS = d["s"] / T if T else 0
        pP = d["p"] / T if T else 0
        pSP = d["sp"] / T if T else 0
        if d["sp"] == 0:
            pSP = 2.0 * pS * pP
        mass = (d["sq"] / h if h else 0) * 5400 * pSP ** 2
        tot_mass += mass * h
        print(f"{label(k):24} {d['sums']:5} {h:7.3f} {d['sq']:6} {d['sq'] / max(d['sums'], 1):7.2f} "
              f"{d['sq'] / h if h else 0:8.1f} {d['time'] / max(d['sums'], 1):8.2f} "
              f"{pS:8.2e} {pP:8.2e} {pSP:8.2e} {d['sp']:4} {mass:9.2e}")
    print(f"total: {len(sums)} sums, {sum(r['time'] + r['setup_time'] for r in sums) / 3600:.2f} CPU-h, "
          f"{len(squares)} squares, expected magic so far {tot_mass:.2e}")


if __name__ == "__main__":
    main()
