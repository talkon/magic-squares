#!/usr/bin/env python3
"""
Candidate pool of scheduler v2: every exponent assignment P over the first
`npr` primes (not only non-increasing exponents) with

    S0 = 6 P^(1/6) in [s0_min, s0_max]    (the AM-GM bound on the sums)
    tau(P) >= tau_min, k_min <= k <= k_max distinct primes,
    at most ones_max exponents equal to 1, every exponent <= emax,
    assignment ratio r = (P / P_sorted)^(1/6) <= ratio_max

where P_sorted puts the same exponent multiset in non-increasing order on
2, 3, 5, ... (the smallest P with those exponents). 55-70% of the expected
magic squares are in assignments that are not sorted (e.g. 13 7 4 3 0 0 1 1,
r = 1.145; research/existence.md 3.2, 6).

Generation is exact: first the sorted multisets (DFS with tau and log P
pruning), then for each multiset every placement on the primes with branch
and bound on the excess log P - log P_sorted, whose exact lower bound is
"the remaining exponents sorted onto the next primes" (rearrangement
inequality), so the work is proportional to the output (~10 us per P).

    python3 scripts/pool.py [--ratio R] [--s0-max X] ...   prints the size
"""
import argparse
import math
import time

import numpy as np

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
LP = [math.log(p) for p in PRIMES]

DEFAULTS = dict(npr=10, s0_max=6000.0, s0_min=0.0, tau_min=1000, k_min=4, k_max=10,
                ratio_max=1.2, ones_max=6, emax=31)


def multisets(npr, lmax, tau_lo, k_lo, k_hi, emax, ones_max):
    """non-increasing exponent vectors a_1 >= ... >= a_k > 0 (on the first k
    primes, so log P is the smallest for the multiset) with log P <= lmax,
    k_lo <= k <= k_hi, tau >= tau_lo, at most ones_max exponents 1"""
    out = []
    cur = []

    def rec(i, prev, lp, t, ones):
        k = len(cur)
        if k >= k_lo and t >= tau_lo:
            out.append((tuple(cur), lp, t))
        if i >= npr or k >= k_hi:
            return
        a = 1
        while a <= prev and lp + a * LP[i] <= lmax + 1e-12:
            if a == 1 and ones >= ones_max:
                a += 1
                continue
            cur.append(a)
            rec(i + 1, a, lp + a * LP[i], t * (a + 1), ones + (a == 1))
            cur.pop()
            a += 1

    rec(0, emax, 0.0, 1, 0)
    return out


def assignments(ms, npr, budget, lmin, lmax, lsorted):
    """every placement of multiset ms (non-increasing) on primes 0..npr-1
    (a prime may stay unused) with log P - lsorted <= budget and log P in
    [lmin, lmax]; returns (exponent tuple, log P) pairs"""
    vals = sorted(set(ms), reverse=True)
    cnt = [ms.count(v) for v in vals]
    k = len(ms)
    res = []
    cur = [0] * npr

    def lb(i, rem):
        return sum(a * LP[i + j] for j, a in enumerate(rem))

    def rec(i, left, lp):
        if left == 0:
            if lp >= lmin - 1e-12:
                res.append((tuple(cur[:i]), lp))
            return
        if npr - i < left:
            return
        rem = [v for v, c in zip(vals, cnt) for _ in range(c)]
        for j, v in enumerate(vals):
            if cnt[j] == 0:
                continue
            cnt[j] -= 1
            rem2 = [x for x, c in zip(vals, cnt) for _ in range(c)]
            l2 = lp + v * LP[i]
            low = l2 + lb(i + 1, rem2)
            if low - lsorted <= budget + 1e-12 and low <= lmax + 1e-12:
                cur[i] = v
                rec(i + 1, left - 1, l2)
                cur[i] = 0
            cnt[j] += 1
        # prime i unused
        low = lp + lb(i + 1, rem)
        if npr - i - 1 >= left and low - lsorted <= budget + 1e-12 and low <= lmax + 1e-12:
            rec(i + 1, left, lp)

    rec(0, k, 0.0)
    return res


def gen_pool(npr=10, s0_max=6000.0, s0_min=0.0, tau_min=1000, k_min=4, k_max=10, ratio_max=1.2,
             ones_max=6, emax=31):
    """the pool as a dict of arrays: exps uint8 [M, npr], S0 float32,
    ratio float32, tau int32, k int8 (rows in lexicographic order of exps)"""
    if npr > len(PRIMES):
        raise ValueError(f"at most {len(PRIMES)} primes")
    lmin = 6 * math.log(max(s0_min, 1e-9) / 6) if s0_min > 0 else -1.0
    lmax = 6 * math.log(s0_max / 6)
    budget = 6 * math.log(ratio_max)
    rows, lps, taus, ks, lss = [], [], [], [], []
    for ms, ls, t in multisets(npr, lmax, tau_min, k_min, k_max, emax, ones_max):
        if ls + budget < lmin:
            continue
        for P, lp in assignments(ms, npr, budget, lmin, lmax, ls):
            rows.append(P + (0,) * (npr - len(P)))
            lps.append(lp)
            taus.append(t)
            ks.append(len(ms))
            lss.append(ls)
    exps = np.array(rows, dtype=np.uint8).reshape(-1, npr)
    lps = np.array(lps)
    order = np.lexsort(exps.T[::-1]) if len(exps) else np.zeros(0, int)
    return {"exps": exps[order], "S0": (6 * np.exp(lps / 6)).astype(np.float32)[order],
            "ratio": np.exp((lps - np.array(lss)) / 6).astype(np.float32)[order],
            "tau": np.array(taus, dtype=np.int32)[order], "k": np.array(ks, dtype=np.int8)[order]}


def brute_force(npr, s0_max, s0_min, tau_min, k_min, k_max, ratio_max, ones_max=99, emax=99):
    """the same set by plain enumeration of all exponent tuples (for tests;
    only for small configurations)"""
    lmin = 6 * math.log(s0_min / 6) if s0_min > 0 else -1.0
    lmax = 6 * math.log(s0_max / 6)
    out = set()
    cur = [0] * npr

    def rec(i, lp):
        if i == npr:
            ms = sorted((a for a in cur if a), reverse=True)
            k = len(ms)
            if not (k_min <= k <= k_max) or lp < lmin - 1e-12:
                return
            if sum(1 for a in ms if a == 1) > ones_max or (ms and ms[0] > emax):
                return
            if math.prod(a + 1 for a in ms) < tau_min:
                return
            ls = sum(a * LP[j] for j, a in enumerate(ms))
            if lp - ls > 6 * math.log(ratio_max) + 1e-12:
                return
            out.add(tuple(cur))
            return
        a = 0
        while lp + a * LP[i] <= lmax + 1e-12:
            cur[i] = a
            rec(i + 1, lp + a * LP[i])
            a += 1
        cur[i] = 0

    rec(0, 0.0)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primes", type=int, default=DEFAULTS["npr"])
    ap.add_argument("--s0-max", type=float, default=DEFAULTS["s0_max"])
    ap.add_argument("--s0-min", type=float, default=DEFAULTS["s0_min"])
    ap.add_argument("--tau-min", type=int, default=DEFAULTS["tau_min"])
    ap.add_argument("--ratio", type=float, default=DEFAULTS["ratio_max"])
    a = ap.parse_args()
    t0 = time.time()
    p = gen_pool(a.primes, a.s0_max, a.s0_min, a.tau_min, ratio_max=a.ratio)
    r = p["ratio"]
    print(f"{len(r)} P ({(r <= 1.0 + 1e-6).sum()} sorted, {(r <= 1.1).sum()} with ratio <= 1.1) "
          f"in {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
