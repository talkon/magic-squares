#!/usr/bin/env python3
"""
The calibration ladder: for every semi-magic square found, how many of each
rarer event on the way to a magic square it "should count for" under several
predictors, against the exact counts recomputed from its grid. Totals rung
by rung, overall and binned by S/S_min, tau(P), N, number of primes, data
source, P source (the first source with a square of that P) and, with the
first search's summary, whether the (P, S) lies inside the range the first
search covered; with Poisson and bootstrap intervals. Also the sub-events
(some coordinates pinned, on one diagonal or both; see below), the
heuristic's variants (its model uncertainty) and E[magic] with the
heuristic's pair factor corrected by the sub-events.

usage:
    python3 scripts/calibrate.py [options] [LABEL=]PATH ...

PATH: msearch JSON-lines output (files, or directories searched recursively
for *.jsonl); default <state>/units. Squares are deduplicated by hash (the
first LABEL that has a square gets it). N comes from the "sum" (or "csum",
"dsum") records in the same files, S_min from the scheduler's pinfo cache (or
bin/enumerate).

Calibration streams and d-first sums (msearch --diag-first): a "csquare"
record (the stream's plain search of every k-th first row) is a square with
weight k (its "weight" field), so that weighted totals estimate the whole
sum's; every sum, rung and sub-event total, its Poisson interval (on the
effective count, see summarize) and its bootstrap use the weights
(--unweighted: weight 1). "dsquare" records ((square, SP diagonal) pairs of
a d-first d loop) are not a sample of the semi-magic squares (each has an SP
diagonal by construction), so they stay out of the ladder; they are listed
in their own section (pairs, est. pairs = pairs x d_stride, partner flags,
best_score) for the rungs above SP.

examples:
    python3 scripts/calibrate.py                       # data/sched/units
    python3 scripts/calibrate.py --jobs 4 seed=runs/seed fresh=data/sched/units
    python3 scripts/calibrate.py --pinfo old/pinfo_6.json --heuristic-cond PATH

options:
    --state DIR      scheduler state (default data/sched): pinfo_6.json (read
                     only), model_6.json (the regression; built-in
                     DEFAULT_MODEL if absent), legacy.json (for --eb)
    --pinfo FILE     more pinfo_6.json caches to read (repeatable)
    --legacy FILE    the first search's stats_short.txt for the per-P check
                     (default stats/stats_short.txt; --legacy none to skip)
    --out DIR        outputs and caches (default data/calibrate)
    --jobs K         processes for the heuristic (default min(4, cores))
    --no-heuristic, --heuristic-cond, --no-sub-events, --eb, --model FILE (repeatable),
    --plugin NAME=FILE[:FUNC[:k=v,...]], --pred NAME=FILE.jsonl[:FIELD],
    --exclude GLOB (repeatable), --hashes FILE (only these squares), --boot B,
    --no-per-square, --w-samples N

outputs (in --out): ladder.md (all tables), results.json (all numbers),
per_square.jsonl (observed and every predictor, per square); caches:
w_table_*.npz, heuristic_<version>.jsonl, pinfo_cache.json.

Events (as square_diag_stats in src/c/square.c):
    traversals (the 720 permutations): S (sum S, incl. SP), P (product P,
        incl. SP), SP (both), S_only, P_only
    partner pairs {sigma, sigma o tau} (tau one of the 15 fixed-point-free
        involutions; 5400 unordered pairs), by exclusive class: 0+0 S+0 P+0
        S+S S+P P+P SP+0 SP+S SP+P SP+SP (scores 0,2,3,4,5,6,7,9,10,14 with
        S = 2, P = 3, SP = 7 per diagonal); inclusive: sum+sum = S+S + SP+S
        + SP+SP (both diagonals sum S), prod+prod = P+P + SP+P + SP+SP
    squares by best pair: best=k, best>=k (k a pair score); magic = SP+SP
    sub-events (square_heuristic.components(sub_events=True)): the sum X and
        d of the k exponent coordinates (cyclic windows of primes) pinned at
        their targets on one diagonal, or on both diagonals of a partner
        pair; the full sets are the rungs. Reported per d: O/E for each
        Edgeworth variant, the X-Y coupling (O/E with X over O/E of the same
        exponents alone), the pair factor given the counts (O / own-count
        null) against the heuristic's (kappa), overall and by the exact
        congruence class of each tau (forbidden / no extra congruence /
        extra congruences), with a bootstrap over squares. beta in
        O/null = kappa^beta is fitted on the exponent windows (d >= 2).

Predictors:
    regression  the scheduler's model: p_S = Model.p_s, p_P = Model.p_p at
                the square's (P, S, S_min); q_SP = rho p_S p_P; traversal
                classes i.i.d. over the 720 traversals, so E[#X+Y pairs] =
                5400 (2) q_X q_Y and P(magic) ~ 5400 q_SP^2 as in the scheduler.
    null        own-count null: the square's own counts of S_only, P_only and
                SP traversals placed on a uniformly random subset of the 720
                (E[#X+Y] = 15/719 n_X n_Y, or 15/719 C(n_X, 2)). It reproduces
                the traversal rungs by construction and tests the pair rungs
                given the counts.
    heuristic   square_heuristic.predict(grid, S, P): from the square's own
                entries (exact moments and lattices, local CLT with Edgeworth
                factors, the exact covariance and mod-3 congruences of each
                pair of diagonals); never looks at the realized counts. No
                parameter is fitted to counts, but its Edgeworth truncation
                was chosen on the P+P rung and the 12-cell product event.
                ~0.25 s per square with the sub-events, cached.
    heur-mixed, heur-marg, heur-gauss, heur-nolat: its variants (bivariate
                factors not within a diagonal's exponents; univariate only;
                plain Gaussian; without the pair congruences), from the same
                components: the model uncertainty.
    heuristic-cond (--heuristic-cond) the same, one rung at a time from the
                realized S and P counts (SP from them, pairs from the realized
                classes), never using the realized SP or pair counts.
    perP        i.i.d. traversals at each P's own pooled observed rates
                (in-sample; used to split the regression's misses into
                rates x spread x pairs-given-counts; its SP+SP is biased up).
    eb (--eb)   regression x the scheduler's per-P empirical Bayes factors
                (posterior means, from these squares + legacy.json; in-sample).

Approximations, explicitly:
  * i.i.d. predictors (regression, perP, eb, plugins returning rates) assume
    the traversal classes are independent across the 720 traversals; pair
    expectations are then exact by linearity. Their best-pair probabilities
    use W(x, y) = P(a random y-set Y and a disjoint random x-set X of
    traversals have a partner pair inside Y or between Y and X) on the exact
    partner graph: exact for y <= 2 and (x, y) = (0, 3), Monte Carlo
    otherwise (s.e. <= 0.0025), Poisson approximation outside the table
    (only in cells of negligible multinomial mass, reported); the class
    counts are mixed over Multinomial(720; q) truncated at a <= 40, b <= 25,
    c <= 6 (captured mass reported).
  * the null's best-pair probabilities use the same W table.
  * the heuristic's best-pair probabilities use Poisson clumping.
  * intervals: "90%" is the exact (Garwood) Poisson interval of the observed
    count divided by E, treating E as exact (Wilson-Hilferty above 300
    events); the bootstrap over squares includes the clustering of pairs
    within a square (it matters for SP+0: 15 pairs per SP traversal).
  * first search (--legacy): only per-P totals were kept (no grids), so the
    check is at the P level: each P's own rates, and the regression at each
    P's midpoint sum (the scheduler's own approximation for legacy totals);
    the 9 SP+S and 1 SP+P squares come from the report (table 11 / 7.1).
    The pooled top rung leaves out the squares here that lie inside the
    first search's range (already in its totals) and is also given with the
    per-P i.i.d. expectation scaled to the own-count null (the ratio of the
    two on these squares).
  * sub-events: bootstrap intervals over squares; beta is a one-parameter
    summary of the trend (weighted least squares through the origin in
    log kappa); the corrected E[magic] applies kappa_c^(beta - 1) to the
    correlation part kappa_c of each square's magic pair factor (keeping
    the exact congruences), an extrapolation from d <= 5 pinned
    coordinates per diagonal to 2 (k + 1).
"""
import argparse
import ast
import fnmatch
import hashlib
import importlib.util
import json
import math
import multiprocessing
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import scheduler as sch  # noqa: E402
import square_heuristic as SH  # noqa: E402

PRIMES = SH.PRIMES

# ---------------------------------------------------------------------------
# geometry of the 6x6 traversals and partner pairs (as in src/c/square.c)

N = 6
NT = SH.NT            # 720 traversals
NPAIRS = SH.NPAIRS    # 5400 unordered partner pairs
CELLIDX = np.arange(N)[None, :] * N + SH.PERMS     # (720, 6) flat cells of traversal i
PAIRS = np.array(sorted({(min(i, j), max(i, j)) for i in range(NT) for j in SH.PIDX[i].tolist()}),
                 dtype=np.int64)
assert len(PAIRS) == NPAIRS
ADJ = np.zeros((NT, NT), dtype=bool)
ADJ[PAIRS[:, 0], PAIRS[:, 1]] = ADJ[PAIRS[:, 1], PAIRS[:, 0]] = True
DEG = SH.NTAU         # 15 partners per traversal
PAIR_FRAC = DEG / (NT - 1)   # P(two random distinct traversals are partners)

# traversal classes: 0 none, 1 S_only, 2 P_only, 3 SP; scores as square.c
CLASS_SCORE = np.array([0, 2, 3, 7])
PAIR_CODES = SH.PAIR_NAMES
PAIR_KEYS = SH.PAIR_KEYS
LEVELS = SH.SCORES
LEVEL_NAME = {0: "none", 2: "S", 3: "P", 4: "S+S", 5: "S+P", 6: "P+P", 7: "SP",
              9: "SP+S", 10: "SP+P", 14: "SP+SP"}
TRAV_KEYS = ["S", "P", "SP", "S_only", "P_only"]
INCL_KEYS = ["sum+sum", "prod+prod"]
BEST_EQ = [f"best={k}" for k in LEVELS]
BEST_GE = [f"best>={k}" for k in LEVELS]
ALL_KEYS = TRAV_KEYS + PAIR_KEYS + INCL_KEYS + BEST_EQ + BEST_GE + ["magic"]
HEURISTIC_KEYS = TRAV_KEYS + PAIR_KEYS + BEST_EQ + BEST_GE + ["none"]   # what square_heuristic returns

# the ladder, in (roughly) decreasing frequency per square
LADDER = [("S", "S traversals"), ("P", "P traversals"), ("S+S", "S+S pairs"),
          ("S+P", "S+P pairs"), ("P+P", "P+P pairs"),
          ("sum+sum", "pairs with both sums S"), ("prod+prod", "pairs with both products P"),
          ("SP", "SP traversals"), ("SP+0", "SP+0 pairs"), ("SP+S", "SP+S pairs"),
          ("SP+P", "SP+P pairs"), ("best>=9", "squares with best >= SP+S"),
          ("SP+SP", "SP+SP pairs (magic)")]
LADDER_KEYS = [k for k, _ in LADDER]


def add_derived(ev):
    """fill the inclusive pair events, magic, S_only/P_only and best=/best>=
    from each other where missing"""
    if "sum+sum" not in ev and all(k in ev for k in ("S+S", "SP+S", "SP+SP")):
        ev["sum+sum"] = ev["S+S"] + ev["SP+S"] + ev["SP+SP"]
    if "prod+prod" not in ev and all(k in ev for k in ("P+P", "SP+P", "SP+SP")):
        ev["prod+prod"] = ev["P+P"] + ev["SP+P"] + ev["SP+SP"]
    if "magic" not in ev and "SP+SP" in ev:
        ev["magic"] = ev["SP+SP"]
    if "S_only" not in ev and "S" in ev and "SP" in ev:
        ev["S_only"] = ev["S"] - ev["SP"]
    if "P_only" not in ev and "P" in ev and "SP" in ev:
        ev["P_only"] = ev["P"] - ev["SP"]
    if all(k in ev for k in BEST_EQ) and not all(k in ev for k in BEST_GE):
        acc = 0.0
        for k in reversed(LEVELS):
            acc = acc + ev[f"best={k}"]
            ev[f"best>={k}"] = acc
    if all(k in ev for k in BEST_GE) and not all(k in ev for k in BEST_EQ):
        for i, k in enumerate(LEVELS):
            nxt = ev[f"best>={LEVELS[i + 1]}"] if i + 1 < len(LEVELS) else 0.0
            ev[f"best={k}"] = ev[f"best>={k}"] - nxt
    return ev


# ---------------------------------------------------------------------------
# W(x, y) on the exact partner graph (for best-pair probabilities)

W_XMAX, W_YMAX = 64, 24
W_TAB = None
W_SE = None
W_USE = {"calls": 0, "approx": 0}


def _w_approx(x, y):
    return 1.0 - np.exp(-PAIR_FRAC * (y * (y - 1) / 2.0 + x * y))


def build_w_table(cache_dir, n_small=200000, n_large=50000, y_small=6, seed=20261007, chunk=10000):
    """W[x, y] for x <= W_XMAX, y <= W_YMAX: Monte Carlo on the exact graph
    (n_small samples for y <= y_small, n_large above), exact for y <= 2 and
    for (0, 3); cached in cache_dir"""
    xmax, ymax = W_XMAX, W_YMAX
    path = os.path.join(cache_dir, f"w_table_{xmax}_{ymax}_{n_small}_{n_large}_{seed}.npz")
    if os.path.exists(path):
        d = np.load(path)
        return d["W"], d["se"]
    t0 = time.time()
    print(f"[building the W table (once, cached in {path})]", file=sys.stderr)
    rng = np.random.default_rng(seed)
    W = np.zeros((xmax + 1, ymax + 1))
    se = np.zeros_like(W)
    for x in range(xmax + 1):          # y = 1: the node's 15 partners among 719
        q = 1.0
        for i in range(x):
            q *= (NT - 1 - DEG - i) / (NT - 1 - i)
        W[x, 1] = 1.0 - q
    within = np.zeros(ymax + 1)
    hist = np.zeros((ymax + 1, xmax + 1))
    nsamp = np.zeros(ymax + 1)
    done = 0
    while done < n_small:
        m = min(chunk, n_small - done)
        nodes = np.argsort(rng.random((m, NT)), axis=1)[:, :ymax + xmax]
        for y in range(2, ymax + 1):
            if y > y_small and done >= n_large:
                continue
            Y, X = nodes[:, :y], nodes[:, y:y + xmax]
            win = ADJ[Y[:, :, None], Y[:, None, :]].any(axis=(1, 2))
            cross = ADJ[Y[:, :, None], X[:, None, :]].any(axis=1)          # (m, xmax)
            first = np.where(cross.any(axis=1), cross.argmax(axis=1), xmax)
            within[y] += win.sum()
            hist[y] += np.bincount(first[~win], minlength=xmax + 1)
            nsamp[y] += m
        done += m
    for y in range(2, ymax + 1):
        cum = np.concatenate([[0.0], np.cumsum(hist[y])[:xmax]])
        W[:, y] = (within[y] + cum) / nsamp[y]
        se[:, y] = np.sqrt(W[:, y] * (1 - W[:, y]) / nsamp[y])
    # y = 2 exactly: Y = {u, v}; if u ~ v the event is certain, else X must
    # avoid N(u) u N(v), of size 30 - cn(u, v) (common neighbours in the
    # Cayley graph of S_6 generated by the 15 involutions)
    cn = np.zeros(NT, dtype=np.int64)
    for t1 in SH.PARTNERS.tolist():
        for t2 in SH.PARTNERS.tolist():
            cn[SH.PERM_INDEX[tuple(t1[t2[i]] for i in range(N))]] += 1
    tset = {SH.PERM_INDEX[tuple(t)] for t in SH.PARTNERS.tolist()}
    ident = SH.PERM_INDEX[tuple(range(N))]
    for x in range(xmax + 1):
        acc = DEG / (NT - 1)
        for gi in range(NT):
            if gi == ident or gi in tset:
                continue
            mnb = 2 * DEG - cn[gi]
            if x <= NT - 2 - mnb:
                avoid = math.exp(math.lgamma(NT - 1 - mnb) - math.lgamma(NT - 1 - mnb - x) -
                                 math.lgamma(NT - 1) + math.lgamma(NT - 1 - x))
            else:
                avoid = 0.0
            acc += (1.0 - avoid) / (NT - 1)
        W[x, 2], se[x, 2] = acc, 0.0
    # (0, 3) exactly: the graph has no triangles
    assert all(cn[i] == 0 for i in tset)
    W[0, 3] = 3 * DEG / (NT - 1) - 3 * DEG * (DEG - 1) / ((NT - 1) * (NT - 2))
    se[0, 3] = 0.0
    os.makedirs(cache_dir, exist_ok=True)
    np.savez(path, W=W, se=se)
    print(f"[W table built in {time.time() - t0:.0f}s]", file=sys.stderr)
    return W, se


def w_func(x, y):
    """vectorized W(x, y): the table, exact for y = 1 at any x, Poisson
    approximation outside the table"""
    x, y = np.broadcast_arrays(np.asarray(x, dtype=np.int64), np.asarray(y, dtype=np.int64))
    out = np.zeros(x.shape)
    inside = (x <= W_XMAX) & (y <= W_YMAX)
    out[inside] = W_TAB[x[inside], y[inside]]
    big = ~inside & (y >= 2)
    out[big] = _w_approx(x[big], y[big])
    one = ~inside & (y == 1)
    if one.any():
        out[one] = 1.0 - np.exp(np.array([
            math.lgamma(NT - DEG) - math.lgamma(NT - DEG - xx) - math.lgamma(NT) + math.lgamma(NT - xx)
            if xx <= NT - DEG - 1 else -np.inf for xx in x[one]]))
    W_USE["calls"] += int(x.size)
    W_USE["approx"] += int(big.sum())
    return out


def best_ge_from_counts(a, b, c):
    """P(best_score >= k | a S_only, b P_only, c SP traversals at uniformly
    random positions) for k in LEVELS; returns an array (10, ...)"""
    a, b, c = np.broadcast_arrays(*(np.asarray(v, dtype=np.int64) for v in (a, b, c)))
    hasc = (c >= 1).astype(float)
    noc = 1.0 - hasc
    out = np.zeros((len(LEVELS),) + a.shape)
    out[0] = 1.0
    out[1] = ((a + b + c) >= 1).astype(float)                  # >= S
    out[2] = np.where((b + c) >= 1, 1.0, w_func(0, a))         # >= P
    out[3] = hasc + noc * w_func(0, a + b)                     # >= S+S
    out[4] = hasc + noc * w_func(a, b)                         # >= S+P
    out[5] = hasc + noc * w_func(0, b)                         # >= P+P
    out[6] = hasc                                              # >= SP
    out[7] = w_func(a + b, c)                                  # >= SP+S
    out[8] = w_func(b, c)                                      # >= SP+P
    out[9] = w_func(0, c)                                      # SP+SP
    return out


def iid_events(qS, qP, qSP, amax=40, bmax=25, cmax=6, chunk=256):
    """all events for i.i.d. traversal classes with per-square probabilities
    (q_S_only, q_P_only, q_SP): pairs exact, best via the multinomial mixture"""
    qS, qP, qSP = (np.asarray(v, dtype=float) for v in (qS, qP, qSP))
    q0 = 1.0 - qS - qP - qSP
    assert (q0 > 0).all() and (qS >= 0).all() and (qP >= 0).all() and (qSP >= 0).all()
    q = [q0, qS, qP, qSP]
    ev = {"S_only": NT * qS, "P_only": NT * qP, "SP": NT * qSP}
    ev["S"] = ev["S_only"] + ev["SP"]
    ev["P"] = ev["P_only"] + ev["SP"]
    for (x, y), name in PAIR_CODES.items():
        ev[name] = NPAIRS * q[x] * q[y] * (1.0 if x == y else 2.0)
    A, B, C = (v.ravel() for v in np.meshgrid(np.arange(amax + 1), np.arange(bmax + 1),
                                              np.arange(cmax + 1), indexing="ij"))
    T = best_ge_from_counts(A, B, C)
    lg = np.array([math.lgamma(v + 1) for v in range(NT + 1)])
    lbase = lg[NT] - lg[A] - lg[B] - lg[C] - lg[NT - A - B - C]
    cnts = [NT - A - B - C, A, B, C]
    approx_cell = ((A + B) > W_YMAX) | ((C >= 2) & ((A + B) > W_XMAX))
    n = len(qS)
    ge = np.zeros((len(LEVELS), n))
    mass, approx_mass = np.zeros(n), np.zeros(n)
    with np.errstate(divide="ignore"):
        lq = [np.log(v) for v in q]
    for lo in range(0, n, chunk):
        sl = slice(lo, lo + chunk)
        lp = np.repeat(lbase[None, :], len(qS[sl]), axis=0)
        with np.errstate(invalid="ignore"):
            for cnt, l in zip(cnts, lq):
                lp += np.where(cnt[None, :] == 0, 0.0, cnt[None, :] * l[sl, None])
        pm = np.exp(lp)
        mass[sl] = pm.sum(1)
        approx_mass[sl] = pm[:, approx_cell].sum(1)
        ge[:, sl] = (pm @ T.T).T
    for i, k in enumerate(LEVELS):
        ev[f"best>={k}"] = ge[i]
    add_derived(ev)
    ev["_mass_captured_min"] = float(mass.min()) if n else 1.0
    ev["_approx_mass_max"] = float(approx_mass.max()) if n else 0.0
    return ev


def null_events(a, b, c):
    """own-count null: the square's own class counts on a random subset"""
    a, b, c = (np.asarray(v, dtype=float) for v in (a, b, c))
    nn = [NT - a - b - c, a, b, c]
    ev = {"S_only": a, "P_only": b, "SP": c, "S": a + c, "P": b + c}
    for (x, y), name in PAIR_CODES.items():
        ev[name] = PAIR_FRAC * (nn[x] * (nn[x] - 1) / 2.0 if x == y else nn[x] * nn[y])
    ge = best_ge_from_counts(a.astype(np.int64), b.astype(np.int64), c.astype(np.int64))
    for i, k in enumerate(LEVELS):
        ev[f"best>={k}"] = ge[i]
    return add_derived(ev)


# ---------------------------------------------------------------------------
# statistics helpers (no scipy)


def _phi(z):
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def gammainc(a, x):
    """regularized lower incomplete gamma P(a, x) (Wilson-Hilferty above a = 300)"""
    if x <= 0:
        return 0.0
    if a <= 0:
        return 1.0
    if a > 300:
        v = 1.0 / (9.0 * a)
        return _phi(((x / a) ** (1.0 / 3.0) - (1.0 - v)) / math.sqrt(v))
    lg = math.lgamma(a)
    if x < a + 1:
        s = term = 1.0 / a
        ap = a
        for _ in range(10000):
            ap += 1
            term *= x / ap
            s += term
            if abs(term) < abs(s) * 1e-15:
                break
        return min(1.0, s * math.exp(-x + a * math.log(x) - lg))
    b, c, d = x + 1 - a, 1e300, 1 / (x + 1 - a)
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        h *= d * c
        if abs(d * c - 1) < 1e-15:
            break
    return max(0.0, 1.0 - math.exp(-x + a * math.log(x) - lg) * h)


def gamma_quantile(a, p):
    lo, hi = 0.0, max(10.0, a * 10 + 50)
    for _ in range(80):
        mid = (lo + hi) / 2
        if gammainc(a, mid) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def poisson_ci(k, conf=0.90):
    """exact (Garwood) interval for a Poisson mean given k observed"""
    al = (1 - conf) / 2
    return (0.0 if k <= 0 else gamma_quantile(k, al)), gamma_quantile(k + 1, 1 - al)


def poisson_tail(k, mu):
    """(P(X >= k), P(X <= k)) for X ~ Poisson(mu)"""
    if mu <= 0:
        return (1.0 if k <= 0 else 0.0), 1.0
    return (1.0 if k <= 0 else gammainc(k, mu)), 1.0 - gammainc(k + 1, mu)


# ---------------------------------------------------------------------------
# S_min and N (read-only use of pinfo caches; new values in our own cache)


class PInfoLite:
    def __init__(self, paths, cache_dir):
        self.paths = [p for p in paths if p and os.path.exists(p)]
        self.cache_path = os.path.join(cache_dir, "pinfo_cache.json")
        self.cache = {}
        if os.path.exists(self.cache_path):
            with open(self.cache_path) as f:
                self.cache = json.load(f)
        self.texts = None
        self.dirty = False

    def _index(self):
        """key -> (file, offset of its value), one regex pass (pinfo files can
        be >100 MB; only the entries used are parsed)"""
        if self.texts is None:
            self.texts, self.offsets = [], {}
            for p in self.paths:
                with open(p) as f:
                    t = f.read()
                k = len(self.texts)
                self.texts.append(t)
                for m in re.finditer(r'"(6:[0-9_]+)": \{"smin"', t):
                    self.offsets.setdefault(m.group(1), (k, m.end() - len('{"smin"')))
        return self.offsets

    def _lookup(self, P, counts=False):
        key = "6:" + "_".join(map(str, P))
        if key in self.cache and (not counts or "counts" in self.cache[key]):
            return self.cache[key]
        ent = None
        off = self._index().get(key)
        if off is not None:
            ent, _ = json.JSONDecoder().raw_decode(self.texts[off[0]], off[1])
        if ent is None:
            out = subprocess.run([os.path.join(ROOT, "bin", "enumerate"), "--vec-size", "6",
                                  "--print-min-sum", *map(str, P)],
                                 capture_output=True, text=True, check=True).stdout.split()
            ent = {"smin": int(out[0]), "counts": {}, "counted_to": 0}
        self.cache[key] = {"smin": ent["smin"]}
        if counts:
            self.cache[key].update({"counts": dict(ent.get("counts", {})),
                                    "counted_to": ent.get("counted_to", 0)})
        self.dirty = True
        return self.cache[key]

    def smin(self, P):
        return self._lookup(P)["smin"]

    def count(self, P, S):
        """number of vectors with sum S (before reduction, = nvecs_raw)"""
        e = self._lookup(P, counts=True)
        if str(S) in e["counts"]:
            return e["counts"][str(S)]
        if S <= e["counted_to"]:
            return 0
        tmp = self.cache_path + f".{os.getpid()}.counts"
        subprocess.run([os.path.join(ROOT, "bin", "enumerate"), "--vec-size", "6", "--reduce", "none",
                        "--counts", "--min-sum", str(S), "--max-sum", str(S), "--file", tmp, *map(str, P)],
                       check=True, stderr=subprocess.DEVNULL)
        c = 0
        with open(tmp) as f:
            for line in f:
                s, v = line.split()
                if int(s) == S:
                    c = int(v)
        os.remove(tmp)
        e["counts"][str(S)] = c
        self.dirty = True
        return c

    def save(self):
        if self.dirty:
            with open(self.cache_path + ".tmp", "w") as f:
                json.dump(self.cache, f)
            os.replace(self.cache_path + ".tmp", self.cache_path)
            self.dirty = False


# ---------------------------------------------------------------------------
# loading square records


def iter_files(path):
    if os.path.isdir(path):
        for root, dirs, files in os.walk(path):
            dirs.sort()
            for f in sorted(files):
                if f.endswith(".jsonl"):
                    yield os.path.join(root, f)
    elif os.path.exists(path):
        yield path


def load_records(specs, excludes, weighted=True):
    """specs: [(label, path)]. Returns (unique 6x6 squares with a grid, in
    first-seen order, each with "weight": 1 for a "square", the stream's
    stride for a "csquare" (1 with weighted=False); {(P, S): nvecs_raw} from
    the sum / csum / dsum records; stats; the "dsquare" records, each with
    "d_stride" from its file's dsum record)"""
    squares, sums = {}, {}
    dsq, dstride = [], {}
    stats = {"files": 0, "records": 0, "dups": 0, "not_n6": 0, "no_grid": 0, "per_label": {},
             "csquares": 0, "dsquares": 0}
    for label, path in specs:
        for f in iter_files(path):
            if any(fnmatch.fnmatch(f, ex) for ex in excludes):
                continue
            stats["files"] += 1
            with open(f) as fh:
                for line in fh:
                    if not line.endswith("\n"):
                        break   # partially written line of a running unit
                    if 'square"' not in line and 'sum"' not in line:
                        continue    # quick filter: other record types (square, csquare,
                        # dsquare; sum, csum, dsum pass)
                    try:
                        r = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if not isinstance(r, dict):
                        continue
                    t = r.get("type")
                    if t in ("sum", "csum", "dsum"):
                        if r.get("n") == 6 and "nvecs_raw" in r:
                            sums[(sch.norm_p(r["P"]), r["S"])] = r["nvecs_raw"]
                        if t == "dsum" and r.get("n") == 6:
                            dstride[(f, sch.norm_p(r["P"]), r["S"])] = int(r.get("d_stride", 1))
                        continue
                    if t == "dsquare":
                        if r.get("n") == 6 and "grid" in r:
                            r["P"] = sch.norm_p(r["P"])
                            r["source"], r["file"] = label, f
                            dsq.append(r)
                            stats["dsquares"] += 1
                        continue
                    if t not in ("square", "csquare"):
                        continue
                    stats["records"] += 1
                    if r.get("n") != 6:
                        stats["not_n6"] += 1
                        continue
                    if "grid" not in r:
                        stats["no_grid"] += 1
                        continue
                    if r["hash"] in squares:
                        stats["dups"] += 1
                        continue
                    r["P"] = sch.norm_p(r["P"])
                    r["source"] = label
                    r["file"] = f
                    r["kind"] = "calib" if t == "csquare" else "plain"
                    r["weight"] = float(r.get("weight", 1)) if (t == "csquare" and weighted) else 1.0
                    stats["csquares"] += t == "csquare"
                    squares[r["hash"]] = r
                    stats["per_label"][label] = stats["per_label"].get(label, 0) + 1
    for r in dsq:
        r["d_stride"] = dstride.get((r["file"], r["P"], r["S"]), 1)
    return list(squares.values()), sums, stats, dsq


# ---------------------------------------------------------------------------
# exact observation (vectorized square_diag_stats)


def observe(squares, chunk=400):
    """exact event counts per square (dict of arrays), best_score per square,
    and the hashes of records that are not valid semi-magic squares"""
    n = len(squares)
    cls_all = np.zeros((n, NT), dtype=np.int8)
    pair_cnt = np.zeros((n, 16), dtype=np.int64)
    best = np.zeros(n, dtype=np.int64)
    bad = []
    for lo in range(0, n, chunk):
        sq = squares[lo:lo + chunk]
        m = len(sq)
        G = np.array([np.array(r["grid"], dtype=np.int64).reshape(-1) for r in sq])  # (m, 36)
        S = np.array([r["S"] for r in sq], dtype=np.int64)
        E = np.zeros((m, len(PRIMES)), dtype=np.int64)
        for i, r in enumerate(sq):
            E[i, :len(r["P"])] = r["P"]
        F = np.zeros((m, 36, len(PRIMES)), dtype=np.int64)
        x = G.copy()
        for j, p in enumerate(PRIMES):
            while True:
                d = (x % p == 0)
                if not d.any():
                    break
                F[:, :, j] += d
                x = np.where(d, x // p, x)
        ok = (x == 1).all(axis=1)
        g3, f4 = G.reshape(m, N, N), F.reshape(m, N, N, len(PRIMES))
        ok &= (g3.sum(2) == S[:, None]).all(1) & (g3.sum(1) == S[:, None]).all(1)
        ok &= (f4.sum(2) == E[:, None, :]).all((1, 2)) & (f4.sum(1) == E[:, None, :]).all((1, 2))
        ok &= np.array([len(set(row.tolist())) == 36 for row in G])
        bad += [sq[i]["hash"] for i in np.nonzero(~ok)[0]]
        isS = G[:, CELLIDX].sum(-1) == S[:, None]
        isP = (F[:, CELLIDX, :].sum(2) == E[:, None, :]).all(-1)
        cls = isS.astype(np.int8) + 2 * isP.astype(np.int8)
        cls_all[lo:lo + m] = cls
        c1 = cls[:, PAIRS[:, 0]].astype(np.int64)
        c2 = cls[:, PAIRS[:, 1]].astype(np.int64)
        code = np.minimum(c1, c2) * 4 + np.maximum(c1, c2)
        pair_cnt[lo:lo + m] = np.bincount((code + 16 * np.arange(m)[:, None]).ravel(),
                                          minlength=16 * m).reshape(m, 16)
        best[lo:lo + m] = (CLASS_SCORE[c1] + CLASS_SCORE[c2]).max(axis=1)
    ev = {"S_only": (cls_all == 1).sum(1).astype(float), "P_only": (cls_all == 2).sum(1).astype(float),
          "SP": (cls_all == 3).sum(1).astype(float)}
    ev["S"] = ev["S_only"] + ev["SP"]
    ev["P"] = ev["P_only"] + ev["SP"]
    for (x, y), name in PAIR_CODES.items():
        ev[name] = pair_cnt[:, x * 4 + y].astype(float)
    for k in LEVELS:
        ev[f"best={k}"] = (best == k).astype(float)
        ev[f"best>={k}"] = (best >= k).astype(float)
    add_derived(ev)
    return ev, best, bad


# ---------------------------------------------------------------------------
# data container


class SquareData:
    def __init__(self, squares, sums, pinfo):
        self.squares = squares
        self.pinfo = pinfo
        self.P = [r["P"] for r in squares]
        self.S = np.array([r["S"] for r in squares], dtype=np.int64)
        self.smin = np.array([pinfo.smin(P) for P in self.P], dtype=np.int64)
        self.N = np.array([sums.get((P, S), -1) for P, S in zip(self.P, self.S.tolist())], dtype=np.int64)
        self.n_from_pinfo = 0
        for i in np.nonzero(self.N < 0)[0]:
            self.N[i] = pinfo.count(self.P[i], int(self.S[i]))
            self.n_from_pinfo += 1
        self.tau = np.array([sch.tau(P) for P in self.P], dtype=np.int64)
        self.k = np.array([sum(1 for e in P if e) for P in self.P], dtype=np.int64)
        self.ratio = self.S / self.smin
        self.source = np.array([r["source"] for r in squares])
        self.w = np.array([float(r.get("weight", 1.0)) for r in squares])
        order = {s: i for i, s in enumerate(dict.fromkeys(self.source.tolist()))}
        first = {}
        for P, s in zip(self.P, self.source.tolist()):
            if P not in first or order[s] < order[first[P]]:
                first[P] = s
        self.psource = np.array([first[P] for P in self.P])
        self.searched = None     # set by set_legacy_range
        self.obs = None

    def set_legacy_range(self, legacy_rows):
        """searched[i]: the first search covered this (P, S), i.e. it searched P
        up to max_S >= S (stats_short.txt; max_S >= 99999 = all sums)"""
        maxs = {}
        for r in legacy_rows:
            maxs[r["P"]] = max(maxs.get(r["P"], 0), r["maxS"])
        self.searched = np.array([S <= maxs.get(P, -1) or maxs.get(P, -1) >= 99999
                                  for P, S in zip(self.P, self.S.tolist())])


# ---------------------------------------------------------------------------
# predictors: predict(data) -> {event: array over squares}


class RegressionPredictor:
    def __init__(self, model_path=None, name="regression"):
        self.name = name
        if model_path and os.path.exists(model_path):
            with open(model_path) as f:
                self.coef = json.load(f)
            self.src = model_path
        else:
            self.coef = dict(sch.DEFAULT_MODEL)
            self.src = "scheduler.DEFAULT_MODEL"
        self.note = (f"scheduler model ({self.src}): p_S, p_P per square, q_SP = rho p_S p_P "
                     f"(rho = {self.coef['rho']:.4g}), i.i.d. traversals")

    def model(self):
        sch.CLAMP.update(sch.WIDE)
        return sch.Model(dict(self.coef))   # sets the model's own clamps

    def probs(self, data):
        m = self.model()
        ps = np.array([m.p_s(P, S, smin) for P, S, smin in zip(data.P, data.S, data.smin)])
        pp = np.array([m.p_p(P, S, smin) for P, S, smin in zip(data.P, data.S, data.smin)])
        return ps, pp, self.coef["rho"]

    def predict(self, data):
        ps, pp, rho = self.probs(data)
        qsp = np.minimum(rho * ps * pp, np.minimum(ps, pp))
        ev = iid_events(ps - qsp, pp - qsp, qsp)
        ev["_pS"], ev["_pP"] = ps, pp
        return ev


class EBPredictor(RegressionPredictor):
    """regression x per-P empirical Bayes factors of scheduler.Scorer
    (posterior means, no exploration bonus), from the loaded squares and the
    legacy per-P totals; in-sample by construction"""

    def __init__(self, model_path=None, legacy_json=None, name="eb"):
        super().__init__(model_path, name)
        self.legacy_json = legacy_json
        self.note = ("regression x per-P factors f_S, f_P (scheduler.Scorer, posterior mean, from these "
                     "squares" + (" + legacy totals" if legacy_json else "") + "; in-sample)")

    def predict(self, data):
        m = self.model()
        res = sch.Results(6)
        for r, P, o_s, o_p in zip(data.squares, data.P, data.obs["S"], data.obs["P"]):
            res.squares[r["hash"]] = {"P": P, "S": r["S"], "s_count": int(o_s), "p_count": int(o_p)}
        if self.legacy_json and os.path.exists(self.legacy_json):
            with open(self.legacy_json) as f:
                ours = set(data.P)
                for v in json.load(f):
                    P = sch.norm_p(v["P"])
                    if P in ours and v.get("sols", 0) > 0:
                        res.legacy[P] = v
        scorer = sch.Scorer(m, res, data.pinfo, optimistic=False)
        fac = {P: scorer.factor(P, False) for P in set(data.P)}
        ps0, pp0, rho = self.probs(data)
        ps = np.minimum(1.0, np.array([fac[P][1] for P in data.P]) * ps0)
        pp = np.minimum(1.0, np.array([fac[P][2] for P in data.P]) * pp0)
        qsp = np.minimum(rho * ps * pp, np.minimum(ps, pp))
        return iid_events(ps - qsp, pp - qsp, qsp)


class NullPredictor:
    name = "null"
    note = "own-count null: the square's own S_only/P_only/SP counts on a random subset of the 720"

    def predict(self, data):
        o = data.obs
        return null_events(o["S_only"], o["P_only"], o["SP"])


class PerPPredictor:
    name = "perP"
    note = ("i.i.d. traversals at each P's own pooled observed S_only/P_only/SP rates "
            "(in-sample; what the first search's per-P totals allow)")

    def predict(self, data):
        o = data.obs
        acc = {}
        for i, P in enumerate(data.P):
            a = acc.setdefault(P, [0.0, 0.0, 0.0, 0])
            a[0] += o["S_only"][i]
            a[1] += o["P_only"][i]
            a[2] += o["SP"][i]
            a[3] += 1
        q = np.array([[acc[P][j] / (NT * acc[P][3]) for j in range(3)] for P in data.P])
        return iid_events(q[:, 0], q[:, 1], q[:, 2])


def code_version(path):
    """hash of a module's code without its docstrings and comments (so that
    editing the documentation keeps the cache)"""
    with open(path) as f:
        tree = ast.parse(f.read())
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)) and node.body \
                and isinstance(node.body[0], ast.Expr) and isinstance(node.body[0].value, ast.Constant) \
                and isinstance(node.body[0].value.value, str):
            node.body = node.body[1:] or [ast.Pass()]
    return hashlib.sha1(ast.dump(tree).encode()).hexdigest()[:12]


# heuristic variants computed (cheaply, from the same components) and cached
# for every square: key -> combine() arguments. "moment" is the default
# heuristic; the others measure its model uncertainty.
HEUR_VARIANTS = {
    "moment": dict(mode="moment"),                      # pairwise (Kirkwood) Edgeworth: the default
    "conditional": dict(mode="conditional"),
    "moment:mixed": dict(mode="moment", ew="mixed"),    # no bivariate factor for 2 exponents of a diagonal
    "moment:marg": dict(mode="moment", ew="marg"),      # univariate Edgeworth factors only
    "moment:gauss": dict(mode="moment", ew=None),       # plain Gaussian local CLT
    "moment:nolat": dict(mode="moment", pair_lattice=False),   # without the congruences of a pair
}
VARIANT_PRED = {"moment:mixed": "heur-mixed", "moment:marg": "heur-marg", "moment:gauss": "heur-gauss",
                "moment:nolat": "heur-nolat"}


def _round_sig(x, sig=7):
    if isinstance(x, list):
        return [_round_sig(v, sig) for v in x]
    if isinstance(x, float) and math.isfinite(x) and x != 0.0:
        return float(f"{x:.{sig}g}")
    return x


def _heuristic_work(item):
    h, grid, S, P, sub = item
    try:
        comp = SH.components(grid, S, P, sub_events=sub)
        out = {"hash": h, **{k: SH.combine(comp, **kw) for k, kw in HEUR_VARIANTS.items()}}
        if sub:
            out["sub"] = [{k: _round_sig(v) for k, v in row.items()} for row in comp["sub"]]
        return out
    except Exception as e:  # noqa: BLE001  (reported, the square is left out)
        return {"hash": h, "error": repr(e)}


class HeuristicPredictor:
    """square_heuristic per square: every variant of HEUR_VARIANTS and (with
    sub_events) the sub-event rows of square_heuristic.components, computed
    together and cached in cache_dir by hash and module version"""

    def __init__(self, cache_dir, jobs=1, sub_events=True):
        self.version = code_version(SH.__file__)
        self.cache_path = os.path.join(cache_dir, f"heuristic_{self.version}.jsonl")
        self.jobs = max(1, jobs)
        self.keys = tuple(HEUR_VARIANTS) + (("sub",) if sub_events else ())
        self.sub_events = sub_events
        self.results = None
        self.errors = {}

    def compute(self, data):
        if self.results is not None:
            return self.results
        cache = {}
        if os.path.exists(self.cache_path):
            with open(self.cache_path) as f:
                for line in f:
                    if not line.endswith("\n"):
                        break
                    try:
                        r = json.loads(line)
                    except json.JSONDecodeError:
                        continue
                    if "error" not in r:
                        cache.setdefault(r["hash"], {}).update({m: r[m] for m in r if m != "hash"})
        todo = [(r["hash"], r["grid"], int(r["S"]), list(r["P"]), self.sub_events)
                for r in data.squares if not all(m in cache.get(r["hash"], {}) for m in self.keys)]
        if todo:
            per = 0.3 if self.sub_events else 0.12
            print(f"[heuristic: {len(todo)} squares to compute ({len(data.squares) - len(todo)} cached), "
                  f"{self.jobs} processes, ~{per * len(todo) / self.jobs / 60:.1f} min]", file=sys.stderr)
            t0 = time.time()
            os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
            with open(self.cache_path, "a") as fo:
                if self.jobs > 1:
                    pool = multiprocessing.get_context("fork").Pool(self.jobs)
                    it = pool.imap_unordered(_heuristic_work, todo, chunksize=4)
                else:
                    pool, it = None, map(_heuristic_work, todo)
                try:
                    for i, r in enumerate(it):
                        if "error" in r:
                            self.errors[r["hash"]] = r["error"]
                        else:
                            cache.setdefault(r["hash"], {}).update({m: r[m] for m in r if m != "hash"})
                            fo.write(json.dumps(r) + "\n")
                        if (i + 1) % 500 == 0:
                            fo.flush()
                            print(f"  [heuristic] {i + 1}/{len(todo)} {time.time() - t0:.0f}s", file=sys.stderr)
                finally:
                    if pool:
                        pool.close()
                        pool.join()
            if self.errors:
                print(f"[heuristic failed on {len(self.errors)} squares, e.g. "
                      f"{next(iter(self.errors.items()))}; they get NaN]", file=sys.stderr)
        self.results = cache
        return cache

    def predict(self, data, mode="moment"):
        cache = self.compute(data)
        ev = {k: np.array([cache.get(r["hash"], {}).get(mode, {}).get(k, np.nan) for r in data.squares])
              for k in HEURISTIC_KEYS}
        return add_derived(ev)

    def sub_rows(self, data):
        """per square, its sub-event rows (None if not computed)"""
        cache = self.compute(data)
        return [cache.get(r["hash"], {}).get("sub") for r in data.squares]


class PluginPredictor:
    """FUNC(grid, S, P_exponents, **kw) from a module file, per square: a dict
    of event expectations (keys as above; missing derived keys are filled in)
    or {q_S_only, q_P_only, q_SP} for an i.i.d. traversal model (then every
    event is computed here). Cached per hash in cache_dir."""

    def __init__(self, spec, cache_dir):
        name, rest = spec.split("=", 1)
        parts = rest.split(":")
        self.name, self.path = name, parts[0]
        self.func = parts[1] if len(parts) > 1 and parts[1] else "predict"
        self.kw = {}
        if len(parts) > 2 and parts[2]:
            for kv in parts[2].split(","):
                k, v = kv.split("=", 1)
                try:
                    v = json.loads(v)
                except json.JSONDecodeError:
                    pass
                self.kw[k] = v
        self.cache_path = os.path.join(cache_dir, f"plugin_{name}.json")
        self.note = f"plugin {self.path}:{self.func}({self.kw})"

    def predict(self, data):
        d = os.path.dirname(os.path.abspath(self.path))
        if d not in sys.path:
            sys.path.insert(0, d)
        spec = importlib.util.spec_from_file_location(f"plugin_{self.name}", self.path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        fn = getattr(mod, self.func)
        with open(self.path, "rb") as f:
            tag = hashlib.sha1(f.read() + json.dumps([self.func, self.kw], sort_keys=True).encode()).hexdigest()
        cache = {}
        if os.path.exists(self.cache_path):
            with open(self.cache_path) as f:
                c = json.load(f)
            if c.get("tag") == tag:
                cache = c["rows"]
        rows = []
        for r in data.squares:
            if r["hash"] not in cache:
                out = fn(r["grid"], r["S"], list(r["P"]), **self.kw)
                cache[r["hash"]] = {k: float(v) for k, v in out.items()
                                    if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool)}
            rows.append(cache[r["hash"]])
        with open(self.cache_path, "w") as f:
            json.dump({"tag": tag, "rows": cache}, f)
        if rows and all(k in rows[0] for k in ("q_S_only", "q_P_only", "q_SP")):
            return iid_events(*[np.array([row[k] for row in rows]) for k in ("q_S_only", "q_P_only", "q_SP")])
        keys = set().union(*[set(r) for r in rows]) if rows else set()
        return add_derived({k: np.array([row.get(k, np.nan) for row in rows]) for k in keys})


class PrecomputedPredictor:
    """per-square predictions computed elsewhere: JSON lines
    {"hash": ..., FIELD: {event: value}}"""

    def __init__(self, spec):
        name, rest = spec.split("=", 1)
        path, _, field = rest.partition(":")
        self.name, self.field = name, field or name
        self.rows = {}
        with open(path) as f:
            for line in f:
                if not line.endswith("\n"):
                    break
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if self.field in r and "hash" in r:
                    self.rows[r["hash"]] = r[self.field]
        self.note = f"precomputed {path} [{self.field}] ({len(self.rows)} squares)"

    def covers(self, h):
        return h in self.rows

    def predict(self, data):
        keys = set().union(*[set(self.rows[r["hash"]]) for r in data.squares]) if data.squares else set()
        ev = {}
        for k in keys:
            vals = [self.rows[r["hash"]].get(k, np.nan) for r in data.squares]
            if all(isinstance(v, (int, float)) for v in vals):
                ev[k] = np.array(vals, dtype=float)
        return add_derived(ev)


# ---------------------------------------------------------------------------
# summaries


BINNINGS = {
    "S/S_min": ("ratio", [1.0, 1.1, 1.2, 1.3, 1.5, 2.0, float("inf")]),
    "tau(P)": ("tau", [0, 2000, 4000, 8000, 16000, float("inf")]),
    "N": ("N", [0, 1000, 2000, 4000, 8000, float("inf")]),
    "primes": ("k", [0, 4, 5, 6, 7, float("inf")]),
}


def bin_labels(data):
    out = {}
    for name, (field, edges) in BINNINGS.items():
        def label(i):
            lo, hi = edges[i], edges[i + 1]
            if name == "primes":
                return f"{int(lo)}" if hi == lo + 1 else (f"<= {int(hi) - 1}" if lo == 0 else f">= {int(lo)}")
            f = (lambda z: f"{z:g}") if name == "S/S_min" else (lambda z: f"{int(z)}")
            return f"[{f(lo)}, {f(hi)})" if hi != float("inf") else f">= {f(lo)}"

        v = getattr(data, field).astype(float)
        idx = np.clip(np.searchsorted(edges, v, side="right") - 1, 0, len(edges) - 2)
        out[name] = (np.array([label(i) for i in idx]), {label(i): i for i in range(len(edges) - 1)})
    srcs = list(dict.fromkeys(data.source.tolist()))
    out["source"] = (data.source, {s: i for i, s in enumerate(srcs)})
    # P source: the first source (in the order given) with a square of this P,
    # so that "seed P" means the P, whichever run found the square
    psrc = np.array([f"P of {s}" for s in data.psource])
    out["P source"] = (psrc, {f"P of {s}": i for i, s in enumerate(srcs)})
    if data.searched is not None:
        out["first search"] = (np.where(data.searched, "searched", "new"), {"searched": 0, "new": 1})
    return out


# per-square weights (SquareData.w) of the summaries; None = all 1
SQW = None


def summarize(obs, preds, mask, events, conf=0.90):
    """per event: observed total and, per predictor, E, O/E, Poisson interval
    and tail probabilities. Squares where a predictor is NaN are left out of
    its comparison (O is then recomputed on the rest; "n_missing").
    With weights (SQW): O and E are weighted totals, and the interval and
    tails are those of the effective count O / s against E / s, where
    s = sum w^2 E / sum w E is the scale of one event (Var(sum w o) =
    sum w^2 E under the predictor), so a csquare of weight k counts as ~k
    squares but carries the uncertainty of one."""
    res = {}
    wt = np.ones(len(mask)) if SQW is None else SQW
    for e in events:
        if e not in obs:
            continue
        O = float((wt * obs[e])[mask].sum())
        row = {"obs": O, "pred": {}}
        for name, ev in preds.items():
            if e not in ev:
                continue
            v = np.asarray(ev[e], dtype=float)
            ok = mask & np.isfinite(v)
            if not ok.any():
                continue
            Op = float((wt * obs[e])[ok].sum())
            Ev = float((wt * v)[ok].sum())
            sc = float((wt * wt * v)[ok].sum() / Ev) if (SQW is not None and Ev > 0) else 1.0
            sc = sc if sc > 0 else 1.0
            lo, hi = poisson_ci(int(round(Op / sc)), conf)
            ge, le = poisson_tail(int(round(Op / sc)), Ev / sc)
            p = {"E": Ev, "ratio": (Op / Ev) if Ev > 0 else None,
                 "ci": [lo * sc / Ev, hi * sc / Ev] if Ev > 0 else None, "p_ge": ge, "p_le": le}
            if sc != 1.0:
                p["event_scale"] = sc
            if ok.sum() < mask.sum():
                p.update({"obs": Op, "n_missing": int(mask.sum() - ok.sum())})
            row["pred"][name] = p
        res[e] = row
    return res


def bootstrap(obs, preds, events, B=1000, conf=0.90, seed=7):
    """bootstrap over squares of O/E (squares resampled with replacement)"""
    n = len(next(iter(obs.values())))
    rng = np.random.default_rng(seed)
    Wt = rng.multinomial(n, np.full(n, 1.0 / n), size=B).astype(float)
    if SQW is not None:
        Wt = Wt * SQW[None, :]
    out = {}
    for e in events:
        if e not in obs:
            continue
        Ob = Wt @ obs[e]
        out[e] = {}
        for name, ev in preds.items():
            if e not in ev or not np.all(np.isfinite(ev[e])):
                continue
            Eb = Wt @ np.asarray(ev[e], dtype=float)
            with np.errstate(divide="ignore", invalid="ignore"):
                r = Ob / Eb
            r = r[np.isfinite(r)]
            if len(r):
                out[e][name] = [float(np.quantile(r, (1 - conf) / 2)), float(np.quantile(r, 1 - (1 - conf) / 2))]
    return out


# ---------------------------------------------------------------------------
# sub-events: the local-CLT machinery and the pair factor below the rungs

SUB_V = SH.SUB_VARIANTS                       # gauss, marg, pair, mixed
SUB_Q1 = ["O1"] + [f"E1_{v}" for v in SUB_V] + ["CO", "CE", "YO", "YE"]
SUB_Q2 = ["O2"] + [f"E2_{v}" for v in SUB_V] + ["E2n", "NULL", "NK", "NKN"]
# by congruence class c of the tau (SH.CONG_CLASSES): observed pairs, null,
# null x kappa (with / without the pair congruences), E, and E of independent
# pairs at the heuristic's own rate
SUB_QC = [f"{q}_c{c}" for c in range(3) for q in ("O2", "NULL", "NK", "NKN", "E2", "IID", "NT")]


def sub_event_arrays(sub_rows):
    """per square and per (kind, d) the sums over its windows: kind "Y" =
    d exponent coordinates, "XY" = the sum X and d - 1 exponents (on one
    diagonal; on both diagonals of a partner pair for the *2 quantities).
    NULL = own-count null of the pair count (15/719 C(n1, 2)), NK / NKN =
    NULL x the heuristic's pair factor kappa = e2 / (5400 (e1 / 720)^2) with /
    without the pair congruences; CO, CE / YO, YE = O, E of an XY window and of
    the Y window of the same exponents (the X-Y coupling). Squares without
    rows get NaN."""
    keys = sorted({("XY" if 0 in r["I"] else "Y", len(r["I"])) for rows in sub_rows if rows for r in rows})
    kidx = {k: j for j, k in enumerate(keys)}
    n = len(sub_rows)
    A = {q: np.zeros((n, len(keys))) for q in SUB_Q1 + SUB_Q2 + SUB_QC}
    for i, rows in enumerate(sub_rows):
        if not rows:
            for q in A:
                A[q][i] = np.nan
            continue
        ywin = {}
        for r in rows:
            if 0 not in r["I"]:
                ywin[tuple(r["I"])] = r
        for r in rows:
            vals = r["e1"] + r["e2"] + [r["e2n"]]
            if not all(isinstance(v, (int, float)) and math.isfinite(v) for v in vals):
                continue
            j = kidx[("XY" if 0 in r["I"] else "Y", len(r["I"]))]
            n1, n2 = r["n1"], r["n2"]
            A["O1"][i, j] += n1
            A["O2"][i, j] += n2
            for v, e1, e2 in zip(SUB_V, r["e1"], r["e2"]):
                A[f"E1_{v}"][i, j] += e1
                A[f"E2_{v}"][i, j] += e2
            A["E2n"][i, j] += r["e2n"]
            null = PAIR_FRAC * n1 * (n1 - 1) / 2.0
            p1 = r["e1"][2] / NT
            A["NULL"][i, j] += null
            if p1 > 0:
                A["NK"][i, j] += null * r["e2"][2] / (NPAIRS * p1 * p1)
                A["NKN"][i, j] += null * r["e2n"] / (NPAIRS * p1 * p1)
            if "ntc" in r:
                tnull = n1 * (n1 - 1) / (2.0 * (NT - 1))       # own-count null per tau (360 pairs)
                for c in range(3):
                    nt = r["ntc"][c]
                    A[f"NT_c{c}"][i, j] += nt
                    A[f"O2_c{c}"][i, j] += r["n2c"][c]
                    A[f"NULL_c{c}"][i, j] += nt * tnull
                    A[f"E2_c{c}"][i, j] += r["e2c"][c]
                    A[f"IID_c{c}"][i, j] += nt * NT / 2 * p1 * p1
                    if p1 > 0:
                        A[f"NK_c{c}"][i, j] += tnull * r["e2c"][c] / (NT / 2 * p1 * p1)
                        A[f"NKN_c{c}"][i, j] += tnull * r["e2nc"][c] / (NT / 2 * p1 * p1)
            if 0 in r["I"]:
                w = tuple(x for x in r["I"] if x != 0)
                y = ywin.get(w) if w else {"n1": NT, "e1": [NT] * 4}
                if y is not None and math.isfinite(y["e1"][2]):
                    A["CO"][i, j] += n1
                    A["CE"][i, j] += r["e1"][2]
                    A["YO"][i, j] += y["n1"]
                    A["YE"][i, j] += y["e1"][2]
    return keys, A


def sub_event_summary(keys, A, mask, B=1000, seed=11, conf=0.90):
    """totals and ratios of the sub-events over the squares in mask, with a
    bootstrap over squares (the windows of a square are nested and its pairs
    cluster, so Poisson intervals would be too narrow)"""
    ok = mask & np.all(np.isfinite(A["O1"]), axis=1)
    idx = np.nonzero(ok)[0]
    n = len(idx)
    if not n:
        return {"squares": 0, "rows": []}
    rng = np.random.default_rng(seed)
    Wt = rng.multinomial(n, np.full(n, 1.0 / n), size=B).astype(float) if B else np.ones((1, n))
    T = {q: A[q][idx].sum(0) for q in A}
    Tb = {q: Wt @ A[q][idx] for q in A}
    al = (1 - conf) / 2

    def ratio(qn, qd, j):
        """T[qn] / T[qd] in column j and its bootstrap interval"""
        with np.errstate(divide="ignore", invalid="ignore"):
            v = T[qn][j] / T[qd][j] if T[qd][j] > 0 else float("nan")
            b = Tb[qn][:, j] / Tb[qd][:, j]
        b = b[np.isfinite(b)]
        ci = [float(np.quantile(b, al)), float(np.quantile(b, 1 - al))] if len(b) else [float("nan")] * 2
        return float(v), ci

    rows = []
    for j, (kind, d) in enumerate(keys):
        r = {"kind": kind, "d": d, "squares": int((A["O1"][idx, j] + A["E1_pair"][idx, j] > 0).sum()),
             "O1": float(T["O1"][j]), "O2": float(T["O2"][j])}
        for v in SUB_V:
            r[f"E1_{v}"] = float(T[f"E1_{v}"][j])
            r[f"E2_{v}"] = float(T[f"E2_{v}"][j])
        r["oe1"], r["oe1_ci"] = ratio("O1", "E1_pair", j)
        r["oe2"], r["oe2_ci"] = ratio("O2", "E2_pair", j)
        r["null"] = float(T["NULL"][j])
        r["kappa"] = float(T["NK"][j] / T["NULL"][j]) if T["NULL"][j] > 0 else float("nan")
        r["kappa_nocong"] = float(T["NKN"][j] / T["NULL"][j]) if T["NULL"][j] > 0 else float("nan")
        r["o_null"], r["o_null_ci"] = ratio("O2", "NULL", j)
        r["o_nk"], r["o_nk_ci"] = ratio("O2", "NK", j)
        with np.errstate(divide="ignore", invalid="ignore"):
            res = (T["O2"][j] / T["E2_pair"][j]) / (T["O1"][j] / T["E1_pair"][j]) ** 2
            resb = (Tb["O2"][:, j] / Tb["E2_pair"][:, j]) / (Tb["O1"][:, j] / Tb["E1_pair"][:, j]) ** 2
        resb = resb[np.isfinite(resb)]
        r["moment_residual"] = float(res)
        r["moment_residual_ci"] = ([float(np.quantile(resb, al)), float(np.quantile(resb, 1 - al))]
                                   if len(resb) else [float("nan")] * 2)
        r["by_class"] = []
        for c in range(3):
            cc = {"class": SH.CONG_CLASSES[c], "taus": float(T[f"NT_c{c}"][j]), "O2": float(T[f"O2_c{c}"][j]),
                  "null": float(T[f"NULL_c{c}"][j]), "E2": float(T[f"E2_c{c}"][j]),
                  "E2_over_iid": float(T[f"E2_c{c}"][j] / T[f"IID_c{c}"][j]) if T[f"IID_c{c}"][j] > 0 else float("nan"),
                  "kappa": float(T[f"NK_c{c}"][j] / T[f"NULL_c{c}"][j]) if T[f"NULL_c{c}"][j] > 0 else float("nan"),
                  "kappa_nocong": (float(T[f"NKN_c{c}"][j] / T[f"NULL_c{c}"][j]) if T[f"NULL_c{c}"][j] > 0
                                   else float("nan"))}
            cc["o_null"], cc["o_null_ci"] = ratio(f"O2_c{c}", f"NULL_c{c}", j)
            cc["o_nk"], cc["o_nk_ci"] = ratio(f"O2_c{c}", f"NK_c{c}", j)
            r["by_class"].append(cc)
        if kind == "XY":
            with np.errstate(divide="ignore", invalid="ignore"):
                c = (T["CO"][j] / T["CE"][j]) / (T["YO"][j] / T["YE"][j])
                cb = (Tb["CO"][:, j] / Tb["CE"][:, j]) / (Tb["YO"][:, j] / Tb["YE"][:, j])
            cb = cb[np.isfinite(cb)]
            r["coupling"] = float(c)
            r["coupling_ci"] = ([float(np.quantile(cb, al)), float(np.quantile(cb, 1 - al))]
                                if len(cb) else [float("nan")] * 2)
        rows.append(r)
    return {"squares": n, "rows": rows}


def fit_pair_exponent(summ, min_events=20):
    """beta in O/null = kappa^beta, fitted (through the origin, weighted by the
    events) on the exponent-only sub-events on both diagonals: the share of
    the heuristic's pair factor that the data show (beta = 1: all of it,
    beta = 0: none, i.e. pairs independent given the counts)"""
    num = den = 0.0
    for r in summ["rows"]:
        if r["kind"] != "Y" or r["d"] < 2 or r["O2"] < min_events or not r["kappa"] > 1:
            continue
        lk, lo = math.log(r["kappa"]), math.log(r["o_null"])
        num += r["O2"] * lk * lo
        den += r["O2"] * lk * lk
    return num / den if den > 0 else float("nan")


def pair_residual_by_k(summ, min_events=20, fit_from=3):
    """the data's pair factor over the heuristic's, O / (null x kappa), for d
    exponent coordinates pinned on both diagonals, as a function of d:
    measured where there are >= min_events pairs, extrapolated beyond from
    the last measured d with the slope of a log-linear fit in d (weighted by
    events, over the measured d >= fit_from), capped at 1.
    Returns ({d: (r, measured?)}, slope, intercept)."""
    pts = {r["d"]: (r["o_nk"], r["O2"]) for r in summ["rows"]
           if r["kind"] == "Y" and r["O2"] >= min_events and r["o_nk"] > 0}
    use = sorted(d for d in pts if d >= fit_from)
    slope = icpt = float("nan")
    if len(use) >= 2:
        w = np.array([pts[d][1] for d in use], float)
        X = np.vstack([np.array(use, float), np.ones(len(use))]).T
        y = np.log([pts[d][0] for d in use])
        slope, icpt = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * y))
    out = {}
    last = max(pts) if pts else None
    for d in range(1, 12):
        if d in pts:
            out[d] = (float(pts[d][0]), True)
        elif last is not None and math.isfinite(slope) and d > last:
            out[d] = (float(min(1.0, pts[last][0] * math.exp(slope * (d - last)))), False)
    return out, float(slope), float(icpt)


def _ci(v, ci, f="{:.3f}"):
    if v is None or not math.isfinite(v):
        return "-"
    if ci and all(math.isfinite(c) for c in ci):
        return f"{f.format(v)} [{f.format(ci[0])}, {f.format(ci[1])}]"
    return f.format(v)


def md_sub_events(summ):
    rows = summ["rows"]
    out = ["One diagonal (O = traversals with the pinned coordinates at their target; O/E per Edgeworth "
           "variant; the X-Y coupling is O/E of X + the exponents over O/E of the same exponents alone):", "",
           "| pinned | d | squares | O | E (pair) | O/E pair [boot 90%] | mixed | marg | gauss | X-Y coupling [boot] |",
           "|---|---:|---:|---:|---:|---|---:|---:|---:|---|"]
    for r in rows:
        lab = "exponents" if r["kind"] == "Y" else "X + exponents"
        oe = {v: (r["O1"] / r[f"E1_{v}"] if r[f"E1_{v}"] > 0 else float("nan")) for v in SUB_V}
        out.append(f"| {lab} | {r['d']} | {r['squares']} | {fnum(r['O1'])} | {fnum(r['E1_pair'])} | "
                   f"{_ci(r['oe1'], r['oe1_ci'])} | {oe['mixed']:.3f} | {oe['marg']:.3f} | {oe['gauss']:.3f} | "
                   f"{_ci(r.get('coupling'), r.get('coupling_ci')) if r['kind'] == 'XY' else ''} |")
    out += ["", "Both diagonals of a partner pair (O = pairs with the coordinates pinned on both; the pair "
            "factor given the counts is O / null, the heuristic's is kappa; O/(null kappa) is their ratio):", "",
            "| pinned on each | d | O | E (pair) | O/E pair [boot] | mixed | marg | gauss | O/null [boot] | kappa | "
            "kappa w/o pair congruences | O/(null kappa) [boot] | (O2/E2)/(O1/E1)^2 [boot] |",
            "|---|---:|---:|---:|---|---:|---:|---:|---|---:|---:|---|---|"]
    for r in rows:
        if r["E2_pair"] <= 0 and r["O2"] == 0:
            continue
        lab = "exponents" if r["kind"] == "Y" else "X + exponents"
        oe = {v: (r["O2"] / r[f"E2_{v}"] if r[f"E2_{v}"] > 0 else float("nan")) for v in SUB_V}
        out.append(f"| {lab} | {r['d']} | {fnum(r['O2'])} | {fnum(r['E2_pair'])} | {_ci(r['oe2'], r['oe2_ci'])} | "
                   f"{oe['mixed']:.3f} | {oe['marg']:.3f} | {oe['gauss']:.3f} | "
                   f"{_ci(r['o_null'], r['o_null_ci'])} | {r['kappa']:.3f} | {r['kappa_nocong']:.3f} | "
                   f"{_ci(r['o_nk'], r['o_nk_ci'])} | {_ci(r['moment_residual'], r['moment_residual_ci'])} |")
    out += ["", "Both diagonals, by the exact congruence class of each partner involution tau (forbidden: the "
            "target is outside the pair's coset; extra congruences: the pair lattice is finer than the two "
            "diagonals' lattices, and the heuristic boosts these pairs by the index). O/null tests the "
            "heuristic's factor kappa class by class; in the class without extra congruences kappa is only "
            "the correlation between the diagonals (and Edgeworth terms):", "",
            "| pinned on each | d | class | taus | O | null | O/null [boot] | heuristic kappa | E / independent | "
            "O/(null kappa) [boot] |", "|---|---:|---|---:|---:|---:|---|---:|---:|---|"]
    for r in rows:
        if r["O2"] == 0 and r["E2_pair"] < 0.5:
            continue
        lab = "exponents" if r["kind"] == "Y" else "X + exponents"
        for cc in r.get("by_class", []):
            if cc["taus"] == 0:
                continue
            out.append(f"| {lab} | {r['d']} | {cc['class']} | {fnum(cc['taus'])} | {fnum(cc['O2'])} | "
                       f"{fnum(cc['null'])} | {_ci(cc['o_null'], cc['o_null_ci'])} | {cc['kappa']:.3f} | "
                       f"{cc['E2_over_iid']:.3f} | {_ci(cc['o_nk'], cc['o_nk_ci'])} |")
    return "\n".join(out)


def fnum(x, digits=3):
    if x is None:
        return "-"
    if x == 0:
        return "0"
    ax = abs(x)
    if float(x).is_integer() and ax < 1e7:
        return f"{int(x):,}"
    if ax >= 1000:
        return f"{x:,.0f}"
    if ax >= 100:
        return f"{x:.0f}"
    if ax >= 10:
        return f"{x:.1f}"
    if ax >= 0.01:
        return f"{x:.{digits}g}" if ax < 1 else f"{x:.2f}"
    return f"{x:.2e}"


def fratio(p):
    if p is None or p.get("ratio") is None:
        return "-"
    lo, hi = p["ci"]
    r = p["ratio"]
    f = (lambda z: f"{z:.2f}") if r < 10 and hi < 10 else (lambda z: f"{z:.3g}")
    return f"{f(r)} [{f(lo)}, {f(hi)}]"


def md_ladder(summ, pred_names, by_construction=("null",)):
    lines = ["| rung | observed | " + " | ".join(f"{p}: E | O/E [90%]" for p in pred_names) + " |",
             "|---|---:|" + "---:|---|" * len(pred_names)]
    for e, label in LADDER:
        if e not in summ:
            continue
        r = summ[e]
        cells = [label, fnum(r["obs"])]
        for p in pred_names:
            pp = r["pred"].get(p)
            if pp is None:
                cells += ["-", "-"]
            elif p in by_construction and e in TRAV_KEYS:
                cells += ["(= obs)", "-"]
            else:
                cells += [fnum(pp["E"]), fratio(pp)]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def md_boot(boot, summ, pred_names):
    lines = ["| rung | " + " | ".join(f"{p}: O/E, bootstrap 90%" for p in pred_names) + " |",
             "|---|" + "---|" * len(pred_names)]
    for e, label in LADDER:
        if e not in boot:
            continue
        cells = [label]
        for p in pred_names:
            b, pp = boot[e].get(p), summ[e]["pred"].get(p)
            cells.append(fratio({"ratio": pp["ratio"], "ci": b}) if b and pp and pp["ratio"] else "-")
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def md_best(summ, pred_names):
    lines = ["| best pair at least | observed | " + " | ".join(f"{p}: E, O/E [90%]" for p in pred_names) + " |",
             "|---|---:|" + "---|" * len(pred_names)]
    for k in LEVELS[1:]:
        ge = summ.get(f"best>={k}")
        cells = [f"{LEVEL_NAME[k]} ({k})", fnum(ge["obs"])]
        for p in pred_names:
            pp = ge["pred"].get(p)
            cells.append(f"{fnum(pp['E'])}, {fratio(pp)}" if pp else "-")
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


BIN_EVENTS = ["S", "P", "SP", "S+S", "sum+sum", "best>=9", "SP+SP"]


def md_bins(binsumm, nsq, pred_names):
    lines = ["| bin | squares | event | observed | " + " | ".join(f"{p}: E | O/E [90%]" for p in pred_names) + " |",
             "|---|---:|---|---:|" + "---:|---|" * len(pred_names)]
    for b, s in binsumm.items():
        first = True
        for e in BIN_EVENTS:
            r = s.get(e)
            if r is None:
                continue
            cells = [b if first else "", str(nsq[b]) if first else "", e, fnum(r["obs"])]
            first = False
            for p in pred_names:
                pp = r["pred"].get(p)
                if pp is None or (p == "null" and e in TRAV_KEYS):
                    cells += ["-", "-"]
                else:
                    cells += [fnum(pp["E"]), fratio(pp) if e != "SP+SP" else "-"]
            lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# the first search: per-P totals of stats_short.txt


def parse_legacy(path):
    """rows of stats_short.txt (as scheduler.py import-legacy) with the
    hSS = sum C(s,2), hSP = sum s p, hPP = sum C(p,2) columns"""
    lines = open(path).read().split("\n")
    start = next(i for i, l in enumerate(lines) if l.startswith("Statistics for each P, sorted"))
    rows = []
    for l in lines[start + 2:]:
        if not l.startswith("  ") or l.startswith("  P "):
            break
        parts = l.split()
        # P_val num_vecs max max_S count time sols #S #P #SP hSS hSP hPP hM [best]
        for t in (15, 14):
            if len(parts) <= t:
                continue
            try:
                P = tuple(int(x) for x in parts[:len(parts) - t])
            except ValueError:
                continue
            if not P or len(P) > len(sch.PRIMES) or str(sch.p_value(P)) != parts[len(parts) - t]:
                continue
            v = parts[len(parts) - t:]
            rows.append(dict(P=sch.norm_p(P), maxS=int(v[3]), sols=int(v[6]), nS=int(v[7]), nP=int(v[8]),
                             nSP=int(v[9]), hSS=int(v[10]), hSP=int(v[11]), hPP=int(v[12]), hM=int(v[13])))
            break
    return rows


def legacy_section(path, pinfo, reg):
    """P-level checks on the first search's totals: each P's own rates
    (i.i.d. traversals and squares) and the regression at each P's midpoint
    sum; observed best pairs from the report (9 SP+S, 1 SP+P, 0 magic)"""
    rows = [r for r in parse_legacy(path) if r["sols"] > 0]
    tot = {k: sum(r[k] for r in rows) for k in ("sols", "nS", "nP", "nSP", "hSS", "hSP", "hPP", "hM")}
    pi = dict.fromkeys(["hSS", "hPP", "hSP", "SP+S", "SP+P", "SP+SP"], 0.0)
    for r in rows:
        n = r["sols"]
        ps, pp, psp = r["nS"] / NT / n, r["nP"] / NT / n, r["nSP"] / NT / n
        pi["hSS"] += n * NT * (NT - 1) / 2 * ps ** 2
        pi["hPP"] += n * NT * (NT - 1) / 2 * pp ** 2
        pi["hSP"] += n * (NT * psp + NT * (NT - 1) * ps * pp)
        pi["SP+S"] += n * NPAIRS * 2 * psp * (ps - psp)
        pi["SP+P"] += n * NPAIRS * 2 * psp * (pp - psp)
        pi["SP+SP"] += n * NPAIRS * psp ** 2
    m, rho = reg.model(), reg.coef["rho"]
    pr = dict.fromkeys(["S", "P", "SP", "SP+S", "SP+P", "SP+SP", "hSS", "hPP"], 0.0)
    for r in rows:
        smin = pinfo.smin(r["P"])
        S_mid = (smin + min(r["maxS"], 3 * smin)) / 2
        ps, pp = m.p_s(r["P"], S_mid, smin), m.p_p(r["P"], S_mid, smin)
        qsp, n = rho * ps * pp, r["sols"]
        pr["S"] += n * NT * ps
        pr["P"] += n * NT * pp
        pr["SP"] += n * NT * qsp
        pr["SP+S"] += n * NPAIRS * 2 * qsp * (ps - qsp)
        pr["SP+P"] += n * NPAIRS * 2 * qsp * (pp - qsp)
        pr["SP+SP"] += n * NPAIRS * qsp ** 2
        pr["hSS"] += n * NT * (NT - 1) / 2 * ps ** 2
        pr["hPP"] += n * NT * (NT - 1) / 2 * pp ** 2
    rho_legacy = {"own_square": tot["nSP"] / ((tot["hSP"] - tot["nSP"]) / (NT - 1)),
                  "global": tot["nSP"] / (tot["nS"] * tot["nP"] / (NT * tot["sols"]))}
    return {"n_P": len(rows), "totals": tot, "pred_per_P_iid": pi, "pred_regression": pr,
            "rho": rho_legacy, "observed_best": {"SP+S": 9, "SP+P": 1, "SP+SP": 0},
            "report_estimates": {"SP+S": 4.58, "SP+P": 1.46}}


# ---------------------------------------------------------------------------
# main


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="see the module docstring for the events, predictors and approximations")
    ap.add_argument("paths", nargs="*", help="[LABEL=]PATH: msearch output files or directories")
    ap.add_argument("--state", default=os.path.join(ROOT, "data", "sched"),
                    help="scheduler state: pinfo_6.json, model_6.json, legacy.json, units/ (default data/sched)")
    ap.add_argument("--pinfo", action="append", default=[], help="more pinfo_6.json caches (read only)")
    ap.add_argument("--legacy", default=os.path.join(ROOT, "stats", "stats_short.txt"),
                    help="stats_short.txt of the first search ('none' to skip)")
    ap.add_argument("--out", default=os.path.join(ROOT, "data", "calibrate"),
                    help="outputs and caches (default data/calibrate)")
    ap.add_argument("--jobs", type=int, default=min(4, os.cpu_count() or 1),
                    help="processes for the heuristic (~0.1 s per square, cached)")
    ap.add_argument("--exclude", action="append", default=[], help="fnmatch pattern on file paths")
    ap.add_argument("--model", action="append", default=[],
                    help="model json for the regression (repeatable; default <state>/model_6.json, "
                         "else the built-in DEFAULT_MODEL; 'default' = the built-in one)")
    ap.add_argument("--no-heuristic", action="store_true", help="skip square_heuristic")
    ap.add_argument("--heuristic-cond", action="store_true", help="also the heuristic's conditional mode")
    ap.add_argument("--no-sub-events", action="store_true",
                    help="skip the sub-events (pinning some of the coordinates; ~2.5x the heuristic's cost)")
    ap.add_argument("--eb", action="store_true", help="also the regression with per-P empirical Bayes factors")
    ap.add_argument("--plugin", action="append", default=[], help="NAME=FILE[:FUNC[:k=v,...]]")
    ap.add_argument("--pred", action="append", default=[],
                    help="NAME=FILE.jsonl[:FIELD]: precomputed per-square predictions keyed by hash "
                         "(the analysis is restricted to the squares they all cover)")
    ap.add_argument("--hashes", default=None,
                    help="only the squares listed in this file (one hash per line, or JSON lines with "
                         "a 'hash' field, e.g. a per_square.jsonl): to repeat an analysis on a fixed set")
    ap.add_argument("--boot", type=int, default=1000, help="bootstrap replicates (0 = none)")
    ap.add_argument("--no-per-square", action="store_true", help="do not write per_square.jsonl")
    ap.add_argument("--unweighted", action="store_true",
                    help="csquare records with weight 1 instead of their stride")
    ap.add_argument("--w-samples", type=int, default=200000,
                    help="Monte Carlo samples for the partner-graph table (a quarter of them above y = 6)")
    args = ap.parse_args()

    global W_TAB, W_SE, SQW
    t0 = time.time()
    os.makedirs(args.out, exist_ok=True)
    specs = []
    for p in args.paths or [os.path.join(args.state, "units")]:
        if "=" in p and not os.path.exists(p):
            lab, path = p.split("=", 1)
        else:
            ap_ = os.path.abspath(p)
            lab, path = (os.path.relpath(ap_, ROOT) if ap_.startswith(ROOT + os.sep) else p), p
        specs.append((lab, path))
    squares, sums, lstats, dsquares = load_records(specs, args.exclude, weighted=not args.unweighted)
    print(f"loaded {len(squares)} unique 6x6 squares from {lstats['files']} files "
          f"({lstats['records']} records, {lstats['dups']} duplicates)", file=sys.stderr)
    if not squares:
        raise SystemExit("no square records found")

    pre = [PrecomputedPredictor(s) for s in args.pred]
    restricted = None
    only = None
    if args.hashes:
        with open(args.hashes) as f:
            only = {json.loads(line)["hash"] if line.lstrip().startswith("{") else line.strip()
                    for line in f if line.strip()}
    if pre or only is not None:
        keep = [r for r in squares if all(p.covers(r["hash"]) for p in pre) and (only is None or r["hash"] in only)]
        if len(keep) < len(squares):
            by = [p.name for p in pre] + ([args.hashes] if only is not None else [])
            restricted = {"before": len(squares), "after": len(keep), "by": by}
            print(f"restricting to the {len(keep)} of {len(squares)} squares covered by {by}", file=sys.stderr)
            squares = keep
    pinfo = PInfoLite([os.path.join(args.state, "pinfo_6.json")] + args.pinfo, args.out)
    data = SquareData(squares, sums, pinfo)
    pinfo.save()
    SQW = data.w if np.any(data.w != 1.0) else None
    if SQW is not None:
        print(f"weighted: {int((data.w != 1).sum())} csquares of weight {data.w[data.w != 1].min():g}-"
              f"{data.w.max():g} (total weight {data.w.sum():.0f} for {len(squares)} squares)", file=sys.stderr)
    use_legacy = bool(args.legacy and args.legacy != "none" and os.path.exists(args.legacy))
    if use_legacy:
        data.set_legacy_range(parse_legacy(args.legacy))
    obs, best, bad = observe(squares)
    data.obs = obs
    mism = {"s_count": 0, "p_count": 0, "sp_count": 0, "best_score": 0}
    for i, r in enumerate(squares):
        for key, val in (("s_count", obs["S"][i]), ("p_count", obs["P"][i]), ("sp_count", obs["SP"][i]),
                         ("best_score", best[i])):
            if key in r and r[key] != int(val):
                mism[key] += 1
    print(f"recomputed vs records: mismatches {mism}, invalid grids {len(bad)}", file=sys.stderr)
    if bad:
        print("  invalid (not semi-magic) records are kept in the totals: " + ", ".join(bad[:5]), file=sys.stderr)

    W_TAB, W_SE = build_w_table(args.out, n_small=args.w_samples, n_large=max(1, args.w_samples // 4))

    preds, notes = {}, {}
    models = args.model or [os.path.join(args.state, "model_6.json")]
    for i, mp in enumerate(models):
        name = "regression" if i == 0 else f"regression{i + 1}"
        rp = RegressionPredictor(None if mp == "default" else mp, name)
        preds[name], notes[name] = rp.predict(data), rp.note
    reg0 = RegressionPredictor(None if models[0] == "default" else models[0])
    preds["null"], notes["null"] = NullPredictor().predict(data), NullPredictor.note
    hp = None
    if not args.no_heuristic:
        hp = HeuristicPredictor(args.out, args.jobs, sub_events=not args.no_sub_events)
        preds["heuristic"] = hp.predict(data, "moment")
        notes["heuristic"] = (f"square_heuristic.predict (moment mode, version {hp.version}): from each "
                              f"square's entries (moments, lattices, Edgeworth factors); no parameter fitted to "
                              f"counts, but its Edgeworth truncation was chosen on the P+P rung and the 12-cell "
                              f"product event")
        if args.heuristic_cond:
            preds["heuristic-cond"] = hp.predict(data, "conditional")
            notes["heuristic-cond"] = ("square_heuristic.predict (conditional mode): each rung from the "
                                       "realized rung below it (realized S, P counts)")
        for key, name in VARIANT_PRED.items():
            preds[name] = hp.predict(data, key)
            notes[name] = f"heuristic variant {HEUR_VARIANTS[key]} (model uncertainty)"
    preds["perP"], notes["perP"] = PerPPredictor().predict(data), PerPPredictor.note
    if args.eb:
        legacy_json = os.path.join(args.state, "legacy.json")
        eb = EBPredictor(None if models[0] == "default" else models[0],
                         legacy_json if os.path.exists(legacy_json) else None)
        preds["eb"], notes["eb"] = eb.predict(data), eb.note
    for p in pre:
        preds[p.name], notes[p.name] = p.predict(data), p.note
    for spec in args.plugin:
        pl = PluginPredictor(spec, args.out)
        preds[pl.name], notes[pl.name] = pl.predict(data), pl.note
    pinfo.save()

    events = [k for k in ALL_KEYS if k in obs]
    allmask = np.ones(len(squares), bool)
    overall = summarize(obs, preds, allmask, events)
    boot = bootstrap(obs, preds, LADDER_KEYS, B=args.boot) if args.boot else {}
    binlab = bin_labels(data)
    bins_out, nsq_bins = {}, {}
    for bname, (labs, order) in binlab.items():
        bins_out[bname], nsq_bins[bname] = {}, {}
        for b in sorted(set(labs.tolist()), key=lambda z: order.get(z, 99)):
            mask = labs == b
            bins_out[bname][b] = summarize(obs, preds, mask, events)
            nsq_bins[bname][b] = int(mask.sum())

    magic = {name: {"total": float(np.nansum(data.w * ev["SP+SP"])),
                    "by_source": {s: float(np.nansum((data.w * np.asarray(ev["SP+SP"]))[data.source == s]))
                                  for s in dict.fromkeys(data.source.tolist())}}
             for name, ev in preds.items() if "SP+SP" in ev}

    s_, p_, sp_ = data.w * obs["S"], data.w * obs["P"], data.w * obs["SP"]
    groups = {}
    for i, P in enumerate(data.P):
        groups.setdefault(P, []).append(i)
    indep_own = float(np.sum(data.w * (obs["S"] * obs["P"] - obs["SP"]) / (NT - 1)))
    rho_est = {"observed_SP": float(sp_.sum()),
               "own_square": float(sp_.sum() / indep_own) if indep_own > 0 else None,
               "own_square_expected_indep": indep_own,
               "per_P_pooled": float(sp_.sum() / sum(s_[ix].sum() * p_[ix].sum() / (NT * data.w[ix].sum())
                                                     for ix in map(np.array, groups.values()))),
               "model_rates": float(sp_.sum() / np.sum(data.w * preds["regression"]["_pS"] * preds["regression"]["_pP"] * NT)),
               "model_rho": float(reg0.coef["rho"]),
               "ci90_SP_count": list(poisson_ci(int(sp_.sum())))}
    legacy = None
    if use_legacy:
        legacy = legacy_section(args.legacy, pinfo, reg0)
        pinfo.save()
        # the top rung pooled with the first search at the P level, without the
        # squares of this data set that lie inside its searched range (those
        # are already in its totals)
        new = ~data.searched
        top = ("SP+S", "SP+P")
        e_leg = sum(legacy["pred_per_P_iid"][k] for k in top)
        o_leg = sum(legacy["observed_best"][k] for k in top)
        o_new = float(obs["best>=9"][new].sum())
        e_new = float(np.nansum(preds["null"]["best>=9"][new]))
        # per-P i.i.d. rates vs the own-count null on these squares: the bias
        # of the first search's per-P expectation
        e_pp = float(sum(np.nansum(preds["perP"][k]) for k in top))
        e_nl = float(sum(np.nansum(preds["null"][k]) for k in top))
        bias = e_nl / e_pp if e_pp > 0 else float("nan")
        o_tot = o_leg + o_new
        lo, hi = poisson_ci(int(o_tot))
        pooled = {"obs_legacy": o_leg, "E_legacy_perP": e_leg, "obs_new": o_new, "E_new_null": e_new,
                  "perP_vs_null_bias": bias, "obs": o_tot}
        for tag, E in (("", e_leg + e_new), ("_bias_corrected", e_leg * bias + e_new)):
            pooled["E" + tag] = E
            pooled["ratio" + tag] = o_tot / E
            pooled["ci" + tag] = [lo / E, hi / E]
        legacy["top_rung_pooled"] = pooled
        legacy["squares_in_searched_range"] = int(data.searched.sum())

    # sub-events: pin some of the coordinates, on one diagonal or both
    sub, beta = None, float("nan")
    if hp is not None and not args.no_sub_events:
        skeys, SA = sub_event_arrays(hp.sub_rows(data))
        if SQW is not None:
            SA = {q: a * SQW[:, None] for q, a in SA.items()}
        sub = {"all": sub_event_summary(skeys, SA, allmask, B=args.boot)}
        for bname in ("P source", "first search", "source"):
            if bname in binlab:
                labs, order = binlab[bname]
                for b in sorted(set(labs.tolist()), key=lambda z: order.get(z, 99)):
                    sub[f"{bname}: {b}"] = sub_event_summary(skeys, SA, labs == b, B=args.boot)
        beta = fit_pair_exponent(sub["all"])

    # the heuristic's E[magic] with its pair factor corrected by the
    # sub-events. (1) "measured residual": x r(k), the data's pair factor over
    # the heuristic's for the k exponent coordinates of the square's P pinned
    # on both diagonals (measured for d <= 5, extrapolated beyond); the sum
    # coordinate is not corrected (its S+S residual is 1.19, unexplained).
    # (2) "congruences only": the correlation part kappa_c of the pair
    # factor (pairs without the pair congruences over independent pairs at
    # its own SP rate) set to 1, the exact congruences kept.
    magic_corr, pair_corr = {}, {"beta": beta}
    if hp is not None:
        H, Hn = preds["heuristic"], preds["heur-nolat"]
        iid = NPAIRS * (H["SP"] / NT) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            kc = np.where((Hn["SP+SP"] > 0) & (iid > 0), Hn["SP+SP"] / iid, 1.0)
        variants = {}
        if sub:
            rk, slope, icpt = pair_residual_by_k(sub["all"])
            pair_corr["residual_by_d"] = {d: {"r": v[0], "measured": v[1]} for d, v in rk.items()}
            pair_corr["residual_fit"] = {"slope": slope, "intercept": icpt}
            fac = np.array([rk.get(k, (1.0, False))[0] for k in data.k.tolist()])
            variants["heuristic x measured pair residual r(k) (sub-events; extrapolated beyond d = "
                     f"{max([d for d, v in rk.items() if v[1]] or [0])})"] = H["SP+SP"] * fac
        variants["heuristic, pairs given the counts independent except for the congruences"] = H["SP+SP"] / kc
        for lab, v in variants.items():
            magic[lab] = {"total": float(np.nansum(data.w * v)),
                          "by_source": {s: float(np.nansum((data.w * v)[data.source == s]))
                                        for s in dict.fromkeys(data.source.tolist())}}
            magic_corr[lab] = v
        pair_corr["kappa_c_overall"] = float(np.nansum(H["SP+SP"]) / np.nansum(H["SP+SP"] / kc))
        pair_corr["pair_factor_overall"] = float(np.nansum(H["SP+SP"]) / np.nansum(iid))
        pair_corr["congruence_factor_overall"] = float(np.nansum(H["SP+SP"] / kc) / np.nansum(iid))

    # d-first (square, SP diagonal) pairs: the rungs above SP (each pair's
    # SP diagonal d and its 15 partners; partner = an SP partner, i.e. magic)
    dsec = {"records": len(dsquares), "pairs": len(dsquares),
            "est_pairs": float(sum(r["d_stride"] for r in dsquares)),
            "squares": len({r["hash"] for r in dsquares}),
            "by_best_score": {}, "partner_sp": int(sum(1 for r in dsquares if r.get("partner"))),
            "magic": int(sum(1 for r in dsquares if r.get("magic")))}
    for r in dsquares:
        b = LEVEL_NAME.get(int(r.get("best_score", 0)), str(r.get("best_score")))
        dsec["by_best_score"][b] = dsec["by_best_score"].get(b, 0) + 1
    if dsquares:
        g = observe([dict(r, grid=r["grid"]) for r in dsquares])[0]
        dsec["recomputed"] = {k: float(g[k].sum()) for k in ("SP", "SP+0", "SP+S", "SP+P", "SP+SP")}
        dsec["est"] = {k: float((np.array([r["d_stride"] for r in dsquares]) * g[k]).sum())
                       for k in ("SP+S", "SP+P", "SP+SP")}

    if not args.no_per_square:
        with open(os.path.join(args.out, "per_square.jsonl"), "w") as f:
            for i, r in enumerate(squares):
                row = {"hash": r["hash"], "source": r["source"], "kind": r.get("kind", "plain"),
                       "weight": float(data.w[i]), "P_source": str(data.psource[i]),
                       "first_search": (None if data.searched is None else
                                        ("searched" if data.searched[i] else "new")),
                       "P": list(r["P"]), "S": int(r["S"]),
                       "smin": int(data.smin[i]), "N": int(data.N[i]), "tau": int(data.tau[i]),
                       "obs": {k: float(obs[k][i]) for k in TRAV_KEYS + PAIR_KEYS + INCL_KEYS}, "best": int(best[i])}
                for name, ev in preds.items():
                    row[name] = {k: _round_sig(float(np.asarray(ev[k])[i]), 9) for k in
                                 TRAV_KEYS + PAIR_KEYS + INCL_KEYS + BEST_GE if k in ev}
                for lab, v in magic_corr.items():
                    row.setdefault("magic_corrected", {})[lab] = _round_sig(float(v[i]), 9)
                f.write(json.dumps(row) + "\n")

    summary = {
        "data": {"squares": len(squares), "files": lstats["files"], "records": lstats["records"],
                 "duplicates": lstats["dups"], "per_source": {s: int((data.source == s).sum())
                                                              for s in dict.fromkeys(data.source.tolist())},
                 "N_from_sum_records": int(len(squares) - data.n_from_pinfo),
                 "N_from_pinfo_or_enumerate": int(data.n_from_pinfo), "mismatch_vs_records": mism,
                 "invalid_grids": bad, "specs": specs, "restricted": restricted, "excludes": args.exclude,
                 "values_of_P": len(groups), "csquares": lstats["csquares"],
                 "weights": {"total": float(data.w.sum()), "max": float(data.w.max())}},
        "dsquares": dsec,
        "predictors": notes,
        "heuristic_errors": hp.errors if hp else {},
        "w_table": {"xmax": W_XMAX, "ymax": W_YMAX, "max_se": float(W_SE.max()),
                    "fraction_of_W_evaluations_approximated": W_USE["approx"] / max(1, W_USE["calls"])},
        "mass_captured_min": {k: v["_mass_captured_min"] for k, v in preds.items() if "_mass_captured_min" in v},
        "approx_mass_max": {k: v["_approx_mass_max"] for k, v in preds.items() if "_approx_mass_max" in v},
        "rho_estimates": rho_est, "overall": overall, "bootstrap": boot, "bins": bins_out,
        "bin_sizes": nsq_bins, "magic_predicted": magic, "pair_correction": pair_corr, "legacy": legacy,
        "sub_events": sub, "elapsed_s": time.time() - t0,
    }
    with open(os.path.join(args.out, "results.json"), "w") as f:
        json.dump(summary, f, indent=1)

    # markdown
    main_preds = [p for p in ["regression", "null", "heuristic", "heuristic-cond"] if p in preds]
    other = [p for p in preds if p not in main_preds]
    md = ["# Calibration ladder (generated by scripts/calibrate.py)", "",
          f"{len(squares):,} unique squares ({len(groups)} values of P) from {lstats['files']} files; "
          "sources: " + ", ".join(f"{s} {n:,}" for s, n in summary["data"]["per_source"].items()) +
          f". Recomputed counts vs the records: mismatches {mism}, invalid grids {len(bad)}." +
          (f" {lstats['csquares']:,} of them are calibration-stream csquares, weighted by their stride "
           f"(total weight {data.w.sum():,.0f}); intervals use the effective counts." if SQW is not None else ""),
          "",
          f"d-first (square, SP diagonal) pairs (dsquare records, not in the ladder): {dsec['pairs']} "
          f"(est. {dsec['est_pairs']:g} with the d strides) on {dsec['squares']} squares; by best pair "
          f"{dsec['by_best_score'] or '-'}; SP partners (magic) {dsec['partner_sp']}.", "",
          "Predictors: " + "; ".join(f"**{k}**: {v}" for k, v in notes.items()) + ".", "",
          "## The ladder (sums over all squares)", "",
          "O/E with the exact 90% Poisson interval of O (E treated as exact). The null reproduces "
          "the traversal rungs by construction.", "",
          md_ladder(overall, main_preds), ""]
    if other:
        md += ["Other predictors:", "", md_ladder(overall, other), ""]
    if boot:
        md += ["## Bootstrap over squares (90%)", "", md_boot(boot, overall, list(preds)), ""]
    md += ["## Squares by best pair", "", md_best(overall, list(preds)), ""]
    md += ["## The regression's misses: rates x spread x pairs given counts", "",
           "O/E(regression) = [E(perP)/E(regression): each P's own traversal rates vs the model] x "
           "[E(null)/E(perP): spread of the counts across the squares of a P] x [O/E(null)].", "",
           "| event | O | E regression | E perP | E null | O/E reg | rates | spread | pairs given counts [90%] |",
           "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    for e in ["S+S", "S+P", "P+P", "sum+sum", "prod+prod", "SP+0", "SP+S", "SP+P", "best>=4", "best>=9"]:
        r = overall.get(e)
        if not r or not all(k in r["pred"] for k in ("regression", "perP", "null")):
            continue
        Er, Ep, En = (r["pred"][k]["E"] for k in ("regression", "perP", "null"))
        dv = (lambda a, b: f"{a / b:.2f}" if b > 0 else "-")
        md.append(f"| {e} | {fnum(r['obs'])} | {fnum(Er)} | {fnum(Ep)} | {fnum(En)} | {dv(r['obs'], Er)} | "
                  f"{dv(Ep, Er)} | {dv(En, Ep)} | {fratio(r['pred']['null'])} |")
    md += ["", "## SP rate (rho)", "",
           f"{rho_est['observed_SP']:.0f} SP traversals (90% interval {rho_est['ci90_SP_count'][0]:.1f}-"
           f"{rho_est['ci90_SP_count'][1]:.1f}). rho = SP / (SP if sum and product were independent): within "
           f"each square (sum (s p - sp) / 719) {rho_est['own_square']:.2f}, per P pooled "
           f"{rho_est['per_P_pooled']:.2f}, with the model's p_S p_P {rho_est['model_rates']:.2f}; the model "
           f"uses rho = {rho_est['model_rho']:.3f}.", "",
           "## Expected magic squares (sum of E[#SP+SP pairs]) over these squares", "",
           "| predictor | total | " + " | ".join(dict.fromkeys(data.source.tolist())) + " |",
           "|---|---:|" + "---:|" * len(set(data.source.tolist()))]
    for name, mg in magic.items():
        md.append(f"| {name} | {mg['total']:.3g} | " + " | ".join(f"{v:.2g}" for v in mg["by_source"].values()) + " |")
    md.append("")
    if "kappa_c_overall" in pair_corr:
        md += [f"The heuristic's pair factor at the magic rung is {pair_corr['pair_factor_overall']:.2f} over these "
               f"squares (E[magic] / 5400 q_SP^2 at its own SP rate): congruences "
               f"{pair_corr['congruence_factor_overall']:.2f} x correlation {pair_corr['kappa_c_overall']:.2f}."]
        if "residual_by_d" in pair_corr:
            md[-1] += (" Measured pair residual r(d) = O/(null kappa) for d exponents on both diagonals: " +
                       ", ".join(f"{d}: {v['r']:.3f}{'' if v['measured'] else ' (extrapolated)'}"
                                 for d, v in pair_corr["residual_by_d"].items() if d <= 8) + ".")
        md.append("")
    if sub:
        md += ["## Sub-events: pinning some of the coordinates", "",
               "A traversal's coordinates are its sum X and the exponents of the k primes of P. The sub-events pin "
               "d of them at their targets (windows of d consecutive primes, with or without X), on one diagonal or "
               "on both diagonals of a partner pair. The full sets are the rungs (S, P, SP; both sums S, both "
               "products P, magic); the partial ones are 10-10^5x commoner and test the heuristic's local CLT and "
               "its pair factor below the magic rung. Intervals: bootstrap over squares.", "",
               f"Summary of the pair factor given the counts against the heuristic's: O/null = kappa^beta with "
               f"beta = {beta:.2f} (fitted on exponent windows on both diagonals, d >= 2; beta = 1 would confirm "
               "the heuristic's pair factor, 0 would mean pairs independent given the counts). The table by "
               "congruence class separates its two parts.", ""]
        for g, sm in sub.items():
            md += [f"### {g} ({sm['squares']:,} squares)", "", md_sub_events(sm), ""]
    if legacy:
        L, t = legacy, legacy["totals"]
        exp_leg = L["rho"]["own_square"] * rho_est["own_square_expected_indep"]
        md += ["## First search (per-P totals of stats_short.txt; best-pair squares from the report)", "",
               f"{L['n_P']} values of P with squares, {t['sols']:,} squares; traversals S {t['nS']:,}, "
               f"P {t['nP']:,}, SP {t['nSP']} (inclusive). rho within squares {L['rho']['own_square']:.2f}, "
               f"global {L['rho']['global']:.2f}. With the first search's within-square rho the squares above "
               f"would have {exp_leg:.1f} SP traversals; observed {rho_est['observed_SP']:.0f} "
               f"(P(X >= obs) = {poisson_tail(int(rho_est['observed_SP']), exp_leg)[0]:.3f}). Not an independent "
               f"test where the squares lie inside the first search's searched range "
               f"({L.get('squares_in_searched_range', 0):,} of them): those are partly the same squares.", "",
               "| event | observed | per-P i.i.d. (own rates) | O/E | regression (midpoint sum) | O/E | report |",
               "|---|---:|---:|---:|---:|---:|---:|"]
        pi, pr = L["pred_per_P_iid"], L["pred_regression"]
        for e, ob, a, b, rep in [("S traversals", t["nS"], None, pr["S"], None),
                                 ("P traversals", t["nP"], None, pr["P"], None),
                                 ("SP traversals", t["nSP"], None, pr["SP"], None),
                                 ("sum C(s,2)", t["hSS"], pi["hSS"], pr["hSS"], None),
                                 ("sum C(p,2)", t["hPP"], pi["hPP"], pr["hPP"], None),
                                 ("sum s p", t["hSP"], pi["hSP"], None, None),
                                 ("SP+S squares", 9, pi["SP+S"], pr["SP+S"], 4.58),
                                 ("SP+P squares", 1, pi["SP+P"], pr["SP+P"], 1.46),
                                 ("SP+SP squares", 0, pi["SP+SP"], pr["SP+SP"], None)]:
            md.append(f"| {e} | {fnum(ob)} | {fnum(a)} | {f'{ob / a:.2f}' if a else '-'} | {fnum(b)} | "
                      f"{f'{ob / b:.2f}' if b else '-'} | {fnum(rep)} |")
        md.append("")
        tp = L.get("top_rung_pooled")
        if tp:
            md += [f"Top rung (best >= SP+S) pooled with the first search at the P level, leaving out the "
                   f"{L['squares_in_searched_range']:,} squares of this data set inside the first search's searched "
                   f"range (already in its totals): observed {tp['obs_legacy']:.0f} + {tp['obs_new']:.0f}, "
                   f"expected {tp['E_legacy_perP']:.2f} (per-P i.i.d.) + {tp['E_new_null']:.2f} (own-count null) = "
                   f"{tp['E']:.2f}: O/E {fratio(tp)}. On these squares the per-P i.i.d. rates give "
                   f"{1 / tp['perP_vs_null_bias']:.2f}x the own-count null's SP+S + SP+P; corrected for that: "
                   f"E {tp['E_bias_corrected']:.2f}, O/E "
                   f"{fratio({'ratio': tp['ratio_bias_corrected'], 'ci': tp['ci_bias_corrected']})}.", ""]
    for bname in bins_out:
        md += [f"## By {bname}", "", md_bins(bins_out[bname], nsq_bins[bname], main_preds), ""]
    with open(os.path.join(args.out, "ladder.md"), "w") as f:
        f.write("\n".join(md))
    print(f"wrote {args.out}/ladder.md, results.json" + ("" if args.no_per_square else ", per_square.jsonl") +
          f" ({time.time() - t0:.0f}s)", file=sys.stderr)


if __name__ == "__main__":
    main()
