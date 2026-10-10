#!/usr/bin/env python3
"""
Adaptive scheduler for the search over (P, S).

Instead of fixing P and sweeping S, this picks the (P, S) pairs expected to
produce the most *magic* squares per CPU-second, runs bin/msearch on them in
parallel, and keeps learning from the results.

    expected magic squares per second at (P, S)
        = (semi-magic squares per sum at (P, S)) / (CPU-seconds per sum)
        * P(a semi-magic square at (P, S) is magic)

Two models are available (--model).

analytic (the default, "v2"): the existence study's calibrated analytic
model (scripts/amodel.py; research/existence/analytic.md, existence.md 2-4).
For each P it predicts, as a function of S, the number of vectors N, the
semi-magic squares per sum and P(magic | square) = D_n kappa p_pair from a
max-entropy local limit theorem with exact lattice factors and ~3 fitted
constants; it was tested on held-out data up to N = 45k (the regression
below stops growing at N ~ 5.5k). On top of it:
  * the review factors of existence.md 4.1, revised on held-out data
    (squares: x1.23 at N' >= 12k, x0.80 at N' 6-12k, x0.78 below N' 1k,
    x0.90 at x = ln(S / S_min) < 0.1; P(magic): x1.22 at N' >= 3k, x0.85
    at 7+ primes, x0.92 at x < 0.1, pair factor x0.87) as the priors of
    three main-effects quasi-Poisson GLMs (squares, S traversals, P
    traversals) over 216 calibration cells (N' band x primes x assignment
    ratio x x band), refit every --refit-every units by MAP; the dispersion
    of a cell comes from the per-P factors below (counts of one P move
    together), so a few P cannot pin a class factor;
  * per-P gamma factors (a = 15 on squares, i.e. sd 0.26 between P, with
    the squares of one P tempered by their overdispersion 2.5; a = 30 on
    traversals), posterior means (no exploration bonus by default);
  * time per sum ln t = th . [1, ln(N'/4000), ln(L/150), max(0,
    ln(N'/8000)), [labels > 128], N' band offsets] with N' the model's N
    (bias corrected), L its predicted number of distinct entries and the
    step where the label bitsets need a third 64-bit word; refit per msearch
    engine on the process CPU time msearch reports ("cpu"), learning only an
    intercept, the step and the band offsets (the shape stays at the
    prior: a refit dominated by cheap small sums mispredicted the large
    ones by 1.4-2x);
  * d-first units (--dfirst auto, the default): a sum is searched
    diagonal-first (msearch --diag-first, which finds every magic square
    but not the other semi-magic squares) where the measured d-first /
    plain CPU ratio (amodel.DFIRST_RATIO_PRIOR, its level learned online)
    with a plain calibration stream of ~7% (--calib-r1-stride k, whose
    sampled squares enter the class and per-P factors, weighted, and
    whose plain-time estimate the plain law and the ratio) is below 1:
    from N' ~ 3.7k on (engine 4; 4.3k with engine 3). The d-first time
    law (amodel.DFIRST_TIME_PRIOR with the engine's shift, its
    level learned) prices them. A d-first sum longer than 1.5 units is
    split into units of d (--d-range lo:hi); the summary merges their
    "dchunk" records and the planner continues a sum from its first d not
    searched, crediting a part with its share of the sum's CPU and E.
The candidates are a wide pool (scripts/pool.py, --pool wide): every
exponent assignment over 2..29 with S0 = 6 P^(1/6) <= 6000, tau >= 1000,
4-10 primes and (P / P_sorted)^(1/6) <= 1.2 (377,908 P), since most of the
yield is in exponents that are not non-increasing. Each P's model is
evaluated once on 24 sums S0 (1 + u), u = 0.01..1 (`profile`: ~1 CPU-hour
for the pool, 2 niced processes, resumable) and interpolated in S. Units
start at ceil(S0) (no S_min computation) and stop at 2 S0. A lazy greedy
planner keeps all P in a heap keyed by an upper bound of their density and
replaces it by the exact score of the next unit when popped. bin/enumerate
is never called.

regression ("v1", --model regression --pool classic): the number of squares
per sum is a Poisson regression on the number of vectors N (piecewise
linear in log N), how far S is above the smallest possible sum S_min(P),
and the divisor structure of P, fitted to our own runs (one observation per
sum) and to the per-P totals of the first search. The time per sum is a
power law in N. The probability of being magic is estimated from the
traversals of the squares found: if p_S and p_P are the probabilities that a
traversal (a possible diagonal) has the magic sum / product, and rho corrects
for these not being independent, then a square has 5400 (unordered) pairs of
possible diagonals and P(magic) ~= 5400 * (rho * p_S * p_P)^2. These are
multiplied by per-P factors (empirical Bayes, with an optimistic bonus), and
the candidates are gen_candidates (non-increasing exponents on 2, 3, 5, 7
plus 0/1 on at most two of 11..19), with S_min and vector counts from
bin/enumerate.

With both, each P is searched in increasing order of S, the scheduler keeps
a frontier S per P and repeatedly runs the next chunk of sums (a unit: about
--unit-time predicted CPU-seconds, cut where the predicted density falls
below --unit-drop x the unit's best) of the P whose next unit has the
highest predicted magic squares per CPU-second.

All state lives in a directory (default data/sched):
    units/*.jsonl[.gz]  raw msearch output, one file per unit of work
    legacy.json         per-P totals of the first search (`import-legacy`)
    pinfo_6.json        regression: S_min(P) and vector counts (read by
                        analytic for exact S_min where present)
    model_6.json        regression: fitted model (else a built-in default)
    pool_6.npz          analytic: the wide pool (rebuilt if its options change)
    profiles_6/         analytic: model profiles per pool (memmap) + extra.json
    summary_6.json/npz  analytic: everything learned from units/, updated
                        incrementally (rebuilt when missing or the model changes)
    calib_6.json, time_6.json   analytic: learned class factors, time laws
    launched_6.jsonl    analytic: every unit launched, with its predictions
The state is rebuilt from units/ on every start, so runs can be interrupted
and resumed at any time; partially finished units keep their completed sums.

usage:
    scheduler.py import-legacy stats/stats_short.txt   (once; skips the sums
                                          the first search already did)
    scheduler.py profile [--workers K]    analytic: model profiles of the pool
                                          (otherwise done by the first command)
    scheduler.py run [--workers K] [--hours H] [--unit-time T] [--machine cal.json]
    scheduler.py plan [--top K]           show the current ranking
    scheduler.py forecast [--hours H]     predicted squares and magic squares
                                          for the next H CPU-hours (analytic:
                                          composition, uncertainty band;
                                          --sample F for a faster estimate;
                                          --machine cal.json --instance-hours H:
                                          E and P(>=1) after H hours of a
                                          machine, scripts/machine_cal.py)
    scheduler.py fit                      refit the model from all results
    scheduler.py report                   summary of results, calibration
                                          (obs/pred), best squares
    scheduler.py emit [--units K]         print msearch arguments for a static
                                          plan (e.g. SuperCloud, see
                                          submit-sc-plan.sh)
    scheduler.py ingest FILE...           add msearch outputs run elsewhere
    scheduler.py pool [--stats]           analytic: the candidate pool
    scheduler.py compact                  gzip finished unit files
Global options (before the command): --state DIR, --vec-size N. Model and
pool options (after the command): --model, --pool, --pool-ratio,
--pool-s0-max, --pool-primes, --tau-min, --explore, --profile-workers;
d-first (analytic): --dfirst {auto,off,on}, --dfirst-min-n, --calib-frac.
Needs numpy (always for analytic; for `fit` and refits with regression).
"""
import argparse
import collections
import glob
import hashlib
import json
import math
import os
import shutil
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

# numpy work here is small and the machine is shared with the searches: one
# BLAS thread (set before numpy is first imported; profile workers inherit it)
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

# the primes msearch / enumerate accept (ENUM_PRIMES in src/c/enumerate.h;
# len(PRIMES) == ENUM_MAX_PRIMES, checked by test_scheduler.py)
PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
NUM_DIAG_PAIRS = {5: 15 * 120 // 2, 6: 5400, 7: 105 * 5040 // 2}
NUM_TRAVERSALS = {5: 120, 6: 720, 7: 5040}
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# --------------------------------------------------------------------------
# small helpers


def p_value(P):
    v = 1
    for p, e in zip(PRIMES, P):
        v *= p ** e
    return v


def tau(P):
    t = 1
    for e in P:
        t *= e + 1
    return t


def norm_p(P):
    P = list(P)
    while P and P[-1] == 0:
        P.pop()
    return tuple(P)


def p_str(P, sep=" "):
    return sep.join(str(e) for e in P)


def parse_p(s):
    return norm_p(int(x) for x in s.replace("_", " ").replace(",", " ").split())


def log(msg):
    print(time.strftime("[%H:%M:%S] ") + msg, file=sys.stderr, flush=True)


def binary(name, check=True):
    path = os.path.join(ROOT, "bin", name)
    if check and not os.path.exists(path):
        sys.exit(f"{path} not found: run ./build.sh first")
    return path


# --------------------------------------------------------------------------
# per-P information: S_min and number of vectors per sum (cached)


class PInfo:
    def __init__(self, path, n):
        self.path = path
        self.n = n
        self.data = {}
        if os.path.exists(path):
            with open(path) as f:
                self.data = json.load(f)

    def save(self):
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump(self.data, f)
        os.replace(tmp, self.path)

    def key(self, P):
        return f"{self.n}:{p_str(P, '_')}"

    # vector counts computed together with S_min: up to WINDOW * S_min, or
    # until some sum has UNTIL vectors (enough to see where squares appear)
    WINDOW = 1.6
    UNTIL = 1500

    def smin(self, P):
        k = self.key(P)
        if k not in self.data:
            out = subprocess.run([binary("enumerate"), "--vec-size", str(self.n),
                                  "--print-min-sum", "--window", str(self.WINDOW),
                                  "--until", str(self.UNTIL), *map(str, P)],
                                 capture_output=True, text=True, check=True).stdout.split("\n")
            smin = int(out[0])
            counts = {}
            counted_to = int(self.WINDOW * smin)
            for line in out[1:]:
                if line.startswith("# counted to"):
                    counted_to = int(line.split()[-1])
                elif line.strip():
                    S, c = line.split()
                    counts[S] = int(c)
            self.data[k] = {"smin": smin, "counts": counts, "counted_to": counted_to,
                            "window": True}
        return self.data[k]["smin"]

    def counts(self, P, lo, hi):
        """number of vectors (before reduction) for each sum in [lo, hi]"""
        self.smin(P)
        d = self.data[self.key(P)]
        if d["counted_to"] < hi:
            start = max(lo, d["counted_to"] + 1)
            # extend in big steps to avoid many small calls
            end = max(hi, start + 300)
            tmp = f"{self.path}.{os.getpid()}.{abs(hash(P))}.counts"
            subprocess.run([binary("enumerate"), "--vec-size", str(self.n), "--reduce", "none",
                            "--counts", "--min-sum", str(start), "--max-sum", str(end),
                            "--file", tmp, *map(str, P)],
                           check=True, stderr=subprocess.DEVNULL)
            with open(tmp) as f:
                for line in f:
                    s, c = line.split()
                    d["counts"][s] = int(c)
            os.remove(tmp)
            d["counted_to"] = end
        return {S: d["counts"].get(str(S), 0) for S in range(lo, hi + 1)}


# --------------------------------------------------------------------------
# results


def split_range(lo, hi, holes):
    """[lo, hi] without the sums in holes, as a list of ranges"""
    out = []
    for h in sorted(x for x in set(holes) if lo <= x <= hi):
        if h > lo:
            out.append((lo, h - 1))
        lo = h + 1
    if lo <= hi:
        out.append((lo, hi))
    return out


class Results:
    """everything learned from msearch output files"""

    def __init__(self, n):
        self.n = n
        self.sums = {}       # (P, S) -> sum record
        self.squares = {}    # hash -> square record
        # P -> list of (lo, hi): ranges of sums known to be searched (a done
        # record covers [min_sum, last_sum], a sum record its own S). Units
        # of one P can run concurrently (plans run as job arrays) and stop
        # early, so the searched sums need not be a prefix.
        self.covered = {}
        self.legacy = {}     # P -> legacy aggregate (searched up to maxS)
        # d-first sums (msearch --diag-first): (P, S) -> their "dsum"
        # records. They have no "sum" record, and their squares come as
        # (square, diagonal) pairs ("dsquare"), not as the sum's semi-magic
        # squares, so the fits (squares, traversals, time) do not use them
        # until the model is taught to; a sum searched in full ("complete")
        # counts as covered, a part of one (a --d-range unit, a --d-stride
        # sample, a truncated run) does not, even inside a "done" range.
        self.dsums = {}
        self.dmagic = {}     # hash -> "dsquare" record of a magic square

    def add_file(self, path):
        try:
            if path.endswith(".gz"):
                import gzip
                f = gzip.open(path, "rt")
            else:
                f = open(path)
        except OSError:
            return
        # the sums of this file that its "done" records must not cover:
        # d-first sums not searched in full, the plain sums that a
        # --d-range unit left to the unit of its range that starts at d 0,
        # and r1-sampled plain sums ("csum" records of mode "sampled")
        partial = {}
        with f:
            for line in f:
                if not line.endswith("\n"):
                    break  # partially written line of a running unit
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                if r.get("n") != self.n:
                    continue
                P = norm_p(r["P"])
                if r["type"] == "csum" and r.get("mode") == "sampled":
                    partial.setdefault(P, set()).add(r["S"])
                if any(k in r for k in SAMPLED_KEYS):
                    # sampled records (and the done record of a run that
                    # samples its plain sums) are not those of the sum
                    continue
                if r["type"] == "sum":
                    self.sums[(P, r["S"])] = r
                    self.covered.setdefault(P, []).append((r["S"], r["S"]))
                elif r["type"] == "square":
                    r["P"] = P
                    self.squares[r["hash"]] = r
                elif r["type"] == "dsum":
                    self.dsums.setdefault((P, r["S"]), []).append(r)
                    if r.get("complete"):
                        self.covered.setdefault(P, []).append((r["S"], r["S"]))
                    else:
                        partial.setdefault(P, set()).add(r["S"])
                elif r["type"] == "skip":
                    partial.setdefault(P, set()).add(r["S"])
                elif r["type"] == "dsquare":
                    if r.get("magic") or r.get("partner"):
                        r["P"] = P
                        self.dmagic[r["hash"]] = r
                elif r["type"] == "done" and r["last_sum"] >= r["min_sum"]:
                    for lo, hi in split_range(r["min_sum"], r["last_sum"],
                                              partial.get(P, ())):
                        self.covered.setdefault(P, []).append((lo, hi))

    def frontier(self, P, pinfo):
        """next sum to search for P: the first one not covered, from S_min
        (or past the legacy search's range)"""
        f = max(pinfo.smin(P), self.legacy.get(P, {}).get("maxS", 0) + 1)
        for lo, hi in sorted(self.covered.get(P, ())):
            if lo > f:
                break
            f = max(f, hi + 1)
        return f

    def next_covered(self, P, S):
        """the first covered sum > S (None if there is none), where a unit
        starting at S has to stop"""
        nxt = [lo for lo, hi in self.covered.get(P, ()) if lo > S]
        return min(nxt) if nxt else None

    def mark_covered(self, P, lo, hi):
        self.covered.setdefault(P, []).append((lo, hi))

    def by_p(self):
        out = {}
        for (P, S), r in self.sums.items():
            out.setdefault(P, []).append(r)
        for v in out.values():
            v.sort(key=lambda r: r["S"])
        return out


# --------------------------------------------------------------------------
# model

# features are clamped to the range seen when fitting, so that the model is
# not extrapolated far outside its data
WIDE = {"lN": (0.0, 50.0), "x": (0.0, 50.0), "lt": (0.0, 50.0), "k": (0.0, 50.0),
        "lp": (0.0, 50.0)}
CLAMP = dict(WIDE)


def clamp(v, key):
    lo, hi = CLAMP[key]
    return min(max(v, lo), hi)


# below this many vectors the model predicts no squares: none were found in
# ~7000 sums of 36 values of P near S_min (the smallest N with a square was
# 452, Morgenstern's square), and extrapolating the fit there is unreliable
MIN_N_SQUARES = 400

# knots (in log N) of the piecewise linear dependence of the number of
# squares per sum on N: per CPU-second it rises with N, is roughly flat
# around N ~ 1000-2000, then falls
RATE_KNOTS = (math.log(700), math.log(1300), math.log(2000), math.log(3000))


def p_features(P):
    """features of P alone: log number of divisors, number of distinct primes,
    log of the largest prime (clamped to the range seen when fitting)"""
    lt = clamp(math.log(tau(P)), "lt")
    k = clamp(float(sum(1 for e in P if e)), "k")
    lp = clamp(math.log(PRIMES[len(P) - 1]), "lp")
    return lt, k, lp


def sum_features(P, S, N, smin):
    """features for the number of semi-magic squares at (P, S)"""
    lN = clamp(math.log(max(N, 1)), "lN")
    x = clamp(math.log(S / smin), "x")
    lt, k, lp = p_features(P)
    hinges = [max(0.0, lN - kn) for kn in RATE_KNOTS]
    return [1.0, lN, *hinges, x, lt, k, lp, k * hinges[1], k * hinges[2], k * x]


def trav_features(P, S, smin):
    """features for the probability that a traversal is S- or P-magic"""
    lt, k, lp = p_features(P)
    return [1.0, math.log(S), clamp(math.log(S / smin), "x"), lt, k, lp]


def time_features(N):
    lN = math.log(max(N, 1))
    return [1.0, lN, lN * lN]


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def fit_poisson(X, y, exposure, groups=None, weights=None, ridge=1.0, iters=200):
    """Poisson regression with optional grouped observations.

    Row i has features X[i] and exposure[i], with mean mu_i = exposure[i] *
    exp(X[i] . beta). Without groups, y[i] ~ Poisson(mu_i). With groups (an
    int per row), only the total of each group is observed:
    y[g] ~ Poisson(sum of mu_i over rows i in group g) -- e.g. the number of
    squares found by the legacy program over a whole range of sums.

    Each observation can have a weight (its log-likelihood is multiplied by
    it). The first column of X must be the intercept; the others are
    standardized and get a ridge penalty, so constant or collinear columns
    stay harmless. Fitted by Fisher scoring with step halving."""
    import numpy as np
    X = np.asarray(X, float)
    expo = np.asarray(exposure, float)
    if groups is None:
        groups = np.arange(len(X))
    groups = np.asarray(groups)
    y = np.asarray(y, float)
    G = len(y)
    wt = np.ones(G) if weights is None else np.asarray(weights, float)
    mean = X[:, 1:].mean(axis=0)
    sd = X[:, 1:].std(axis=0)
    sd[sd < 1e-9] = np.inf
    Z = np.hstack([X[:, :1], (X[:, 1:] - mean) / sd])
    d = Z.shape[1]
    pen = np.full(d, ridge)
    pen[0] = 1e-9
    off = np.log(expo)

    def evaluate(b):
        mu = np.exp(np.clip(Z @ b + off, -60, 60))
        M = np.bincount(groups, weights=mu, minlength=G)
        Mpos = np.maximum(M, 1e-300)
        ll = float((wt * (y * np.log(Mpos) - M)).sum() - 0.5 * (pen * b * b).sum())
        return ll, mu, Mpos

    b = np.zeros(d)
    b[0] = math.log(max(y.sum(), 0.5) / expo.sum())
    ll, mu, M = evaluate(b)
    for _ in range(iters):
        F = np.zeros((G, d))
        np.add.at(F, groups, mu[:, None] * Z)
        grad = ((wt * (y / M - 1.0))[:, None] * F).sum(axis=0) - pen * b
        info = (F * (wt / M)[:, None]).T @ F + np.diag(pen)
        step = np.linalg.solve(info, grad)
        t = 1.0
        while True:
            ll2, mu2, M2 = evaluate(b + t * step)
            if ll2 >= ll - 1e-9 or t < 1e-6:
                break
            t /= 2
        b, ll, mu, M = b + t * step, ll2, mu2, M2
        if np.abs(t * step).max() < 1e-8:
            break
    beta = np.zeros_like(b)
    beta[1:] = b[1:] / sd
    beta[0] = b[0] - (beta[1:] * mean).sum()
    return [float(v) for v in beta]


def fit_lsq(X, y, ridge=1e-6):
    import numpy as np
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    beta = np.linalg.solve(X.T @ X + ridge * np.eye(X.shape[1]), X.T @ y)
    resid = y - X @ beta
    return [float(b) for b in beta], float(resid.std())


# a unit of work ends where its predicted density (magic squares per
# CPU-second) falls below this fraction of its best (see Scorer.next_unit;
# --unit-drop; 0 = only by --unit-time). In 50-hour forecast simulations,
# 0.5 made the result nearly independent of --unit-time: with 600-s units it
# collected 1.4-2.1x more than without, with 120-s units about the same
# (+18-25% in the first 15 hours, -3% at 50 hours).
UNIT_DROP = 0.5

# version of msearch's search (its "engine" field; records without it are
# version 1): the time model uses only timings of the newest version (v2:
# each engine's law starts from the previous engine's posterior, see
# fit_time_models). 3 = cx/integrated (per-r1 widths, carried bitsets up to
# 512 labels, the pretest); 4 = round 2 of the d-first search (the class
# support in the V_d searches, gated, and the star cover K = 4; the plain
# search unchanged): its d-first law and ratio start from engine 3's with
# amodel.DFIRST_ENGINE_SHIFT / DFIRST_RATIO_ENGINE_SHIFT
ENGINE = 4


def engine_of(r):
    return r.get("engine", 1)


def msearch_env():
    """the environment of the msearch units: ours without bench's SAMPLE_*
    measurement variables (msearch ignores them since cx/integrated, but an
    older binary would turn every unit into an r1-sampled run that never
    covers its sums)"""
    return {k: v for k, v in os.environ.items() if not k.startswith("SAMPLE_")}


# default model (6x6), used until `fit` writes model_6.json in the state
# directory: fitted on msearch runs of 36 values of P from their smallest sum
# (4 minutes each: 7278 sums, 4920 squares) together with the per-P totals of
# the first search (767 values of P, 657410 squares)
DEFAULT_MODEL = {
    "n": 6,
    "squares": [10.2951, 5.81154, 4.10269, -2.73385, 1.67464, -2.95712, -3.28715, -7.38849, 1.82764, -0.692608, 0.23595, -0.616791, -3.36582],
    "min_n": 400,
    "time": [-21.9209, 1.23053, 0.249218],
    "time_sd": 0.468397,
    "time_range": [5.70378, 7.82445],
    "rate_cap": 3.55782,
    "ps": [-0.819939, -1.13315, -1.41967, 0.375084, -0.161966, 0.0577261],
    "pp": [-0.345985, 0.175057, -1.66299, -0.689746, -0.551649, 0.024616],
    "rho": 2.14761,
    "clamp": {"lN": [5.99146, 8.62173], "x": [0.0315903, 1.09861], "lt": [6.79794, 8.77648], "k": [4.0, 6.0], "lp": [1.94591, 2.94444]},
}


class Model:
    def __init__(self, coef):
        self.c = coef
        if "clamp" in coef:
            CLAMP.update({k: tuple(v) for k, v in coef["clamp"].items()})

    @staticmethod
    def load(path, n=6):
        if os.path.exists(path):
            with open(path) as f:
                return Model(json.load(f))
        if n != DEFAULT_MODEL["n"]:
            return None
        return Model(dict(DEFAULT_MODEL))

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.c, f, indent=1)

    def squares(self, P, S, N, smin):
        """expected number of semi-magic squares with sum S"""
        if N < self.c.get("min_n", MIN_N_SQUARES):
            return 0.0
        return math.exp(min(30, dot(self.c["squares"], sum_features(P, S, N, smin))))

    def rate(self, P, S, N, smin):
        """semi-magic squares per CPU-second"""
        r = self.squares(P, S, N, smin) / self.sum_time(N)
        return min(r, self.c.get("rate_cap", r))

    def sum_time(self, N):
        """CPU-seconds to search one sum: quadratic in log N over the fitted
        range, extended as a power law (linearly in log N) beyond it"""
        c = self.c["time"]
        lo, hi = self.c.get("time_range", (0.0, 50.0))
        lN = min(max(math.log(max(N, 1)), lo), hi)
        logt = c[0] + c[1] * lN + c[2] * lN * lN
        if math.log(max(N, 1)) > hi:
            slope = max(3.0, c[1] + 2 * c[2] * hi)
            logt += slope * (math.log(N) - hi)
        return max(math.exp(min(50.0, logt + self.c["time_sd"] ** 2 / 2)), 1e-4)

    def p_s(self, P, S, smin):
        return min(1.0, math.exp(dot(self.c["ps"], trav_features(P, S, smin))))

    def p_p(self, P, S, smin):
        return min(1.0, math.exp(dot(self.c["pp"], trav_features(P, S, smin))))

    def p_magic(self, n, ps, pp):
        return NUM_DIAG_PAIRS[n] * (self.c["rho"] * ps * pp) ** 2


LEGACY_MIN_SQUARES = 0   # use legacy totals of P with this many squares
                         # (keep 0: P that found few squares are informative)
LEGACY_MAX_RATIO = 3.0   # ... and with all sums searched below this * S_min


def legacy_range(lg, smin, max_ratio=LEGACY_MAX_RATIO):
    """the sums whose squares a legacy total counts: [smin, maxS] for P
    searched up to maxS <= max_ratio * smin; for P searched exhaustively, we
    take [smin, max_ratio * smin] (nearly all their squares are there); None
    for long sweeps that stopped above max_ratio * smin"""
    if lg["maxS"] >= 10 ** 9:
        return int(max_ratio * smin)
    if lg["maxS"] <= max_ratio * smin:
        return lg["maxS"]
    return None


def fit_model(results, pinfo, min_n=MIN_N_SQUARES, legacy_min_squares=LEGACY_MIN_SQUARES,
              legacy_max_ratio=LEGACY_MAX_RATIO, min_n_time=300, legacy_weight=0.1):
    """fit all parts of the model; returns (Model, diagnostics).

    Uses every searched sum in results, and, for values of P searched by the
    legacy program (results.legacy, with vector counts available in pinfo),
    the total number of squares and traversals over their range of sums."""
    n = results.n
    T = NUM_TRAVERSALS[n]
    CLAMP.update(WIDE)  # refit from scratch

    # time per sum, from our own runs with the newest version of the search
    # (msearch's "engine"; older timings would overstate the time), or the
    # default model's if there are too few of those. It grows like N^5 in
    # the range that matters; small N are dominated by overheads and not
    # worth modeling.
    engine = max([engine_of(r) for r in results.sums.values()] or [ENGINE])
    Xt, yt = [], []
    for (P, S), r in results.sums.items():
        t = r["time"] + r["setup_time"] + r.get("enum_time", 0)
        if (engine_of(r) == engine and r["nvecs_raw"] >= min_n_time
                and not r["truncated"] and t > 1e-4):
            Xt.append(time_features(r["nvecs_raw"]))
            yt.append(math.log(t))
    if len(yt) >= 20:
        lNs = [x[1] for x in Xt]
        time_range = (min(lNs), max(lNs))
        tcoef, tsd = fit_lsq(Xt, yt)
    else:
        tcoef, tsd = DEFAULT_MODEL["time"], DEFAULT_MODEL["time_sd"]
        time_range = DEFAULT_MODEL["time_range"]

    # squares per sum: our sums (one observation each) ...
    X, expo, groups, y = [], [], [], []
    own_time = 0.0
    for (P, S), r in results.sums.items():
        N = r["nvecs_raw"]
        if N < min_n or r["truncated"]:
            continue
        X.append(sum_features(P, S, N, pinfo.smin(P)))
        expo.append(1.0)
        groups.append(len(y))
        y.append(r["squares"])
        own_time += r["time"] + r["setup_time"]
    own_rows, own_squares = len(y), sum(y)
    # ... and legacy totals per P (one observation per P, over its sums)
    legacy_used = []
    own_p = results.by_p()
    for P, lg in sorted(results.legacy.items()):
        if lg["sols"] < legacy_min_squares or P in own_p:
            continue
        k = pinfo.key(P)
        if k not in pinfo.data or pinfo.data[k]["smin"] == 0:
            continue
        smin = pinfo.smin(P)
        hi = legacy_range(lg, smin)
        if hi is None or pinfo.data[k]["counted_to"] < hi:
            continue
        counts = pinfo.counts(P, smin, hi)
        g = len(y)
        rows = 0
        for S, N in counts.items():
            if N >= min_n:
                X.append(sum_features(P, S, N, smin))
                expo.append(1.0)
                groups.append(g)
                rows += 1
        if rows:
            y.append(lg["sols"])
            legacy_used.append(P)
    if own_squares + sum(y[own_rows:]) < 20:
        raise SystemExit(f"not enough data to fit: {own_squares} squares")
    for key, col in (("lN", 1), ("x", 6), ("lt", 7), ("k", 8), ("lp", 9)):
        CLAMP[key] = (min(x[col] for x in X), max(x[col] for x in X))
    # cap on the predicted rate: a few times the best rate of any P we ran
    # (with the newest search), or the default model's cap
    per_p = {}
    for (P, S), r in results.sums.items():
        if engine_of(r) != engine:
            continue
        a = per_p.setdefault(P, [0, 0.0])
        a[0] += r["squares"]
        a[1] += r["time"] + r["setup_time"]
    rates = [sq / t for sq, t in per_p.values() if sq >= 20 and t > 0]
    rate_cap = 3 * max(rates) if rates else DEFAULT_MODEL["rate_cap"]
    # the legacy totals are dominated by large sums, far from where the
    # scheduler works, so they get a smaller weight than our own data
    wts = [1.0] * own_rows + [legacy_weight] * (len(y) - own_rows)
    sq_coef = fit_poisson(X, y, expo, groups, wts)
    model = Model({"n": n, "squares": sq_coef, "min_n": min_n, "time": tcoef,
                   "time_sd": tsd, "time_range": time_range, "rate_cap": rate_cap})

    # traversal probabilities: our squares (one row each) ...
    Xs, es, gs, ys, yp = [], [], [], [], []
    ysp = 0
    for q in results.squares.values():
        P, S = q["P"], q["S"]
        Xs.append(trav_features(P, S, pinfo.smin(P)))
        es.append(T)
        gs.append(len(ys))
        ys.append(q["s_count"])
        yp.append(q["p_count"])
        ysp += q["sp_count"]
    # ... and legacy totals, spread over the sums of P in proportion to the
    # expected number of squares at each sum
    legacy_sp = 0
    for P in legacy_used:
        lg = results.legacy[P]
        smin = pinfo.smin(P)
        counts = pinfo.counts(P, smin, legacy_range(lg, smin))
        if lg["sols"] == 0:
            continue
        w = {S: model.squares(P, S, N, smin) for S, N in counts.items() if N >= min_n}
        tot = sum(w.values())
        if tot <= 0:
            continue
        g = len(ys)
        for S, wS in w.items():
            Xs.append(trav_features(P, S, smin))
            es.append(T * lg["sols"] * wS / tot)
            gs.append(g)
        ys.append(lg["nS"])
        yp.append(lg["nP"])
        legacy_sp += lg["nSP"]
    ps = fit_poisson(Xs, ys, es, gs)
    pp = fit_poisson(Xs, yp, es, gs)
    expected_sp = sum(e * math.exp(dot(ps, x)) * math.exp(dot(pp, x))
                      for x, e in zip(Xs, es))
    # rho: SP traversals are more common than if S and P were independent (the
    # report measured ~2x); shrink towards 2 when there are few SP traversals
    rho = (ysp + legacy_sp + 2.0 * 5) / (expected_sp + 5)
    model = Model({"n": n, "squares": sq_coef, "min_n": min_n, "time": tcoef,
                   "time_sd": tsd, "time_range": time_range, "rate_cap": rate_cap,
                   "ps": ps, "pp": pp, "rho": rho,
                   "clamp": {k: list(v) for k, v in CLAMP.items()}})
    diag = {"own_sums": own_rows, "own_squares": int(own_squares),
            "own_cpu_seconds": own_time,
            "legacy_P": len(legacy_used), "legacy_squares": int(sum(y[own_rows:])),
            "traversals_S": int(sum(ys)), "traversals_P": int(sum(yp)),
            "traversals_SP": int(ysp + legacy_sp), "expected_SP_if_independent": expected_sp}
    return model, diag


def describe_fit(diag):
    return (f"refit model on {diag['own_sums']} sums / {diag['own_squares']} squares of ours "
            f"and {diag['legacy_P']} legacy values of P / {diag['legacy_squares']} squares")


# --------------------------------------------------------------------------
# per-P empirical Bayes factors and scores

Unit = collections.namedtuple("Unit", "P lo hi score time squares magic")


PRIOR_SQUARES = 3.0     # prior strength of the per-P rate factor, in squares
PRIOR_TRAV = 20.0       # prior strength of the per-P traversal factors
EXPLORE = 1.0           # optimism: posterior mean + EXPLORE * posterior sd


class Scorer:
    def __init__(self, model, results, pinfo, optimistic=True):
        self.optimistic = optimistic  # add the exploration bonus to scores
        self.model = model
        self.results = results
        self.pinfo = pinfo
        self.n = results.n
        T = NUM_TRAVERSALS[self.n]
        # observed vs expected, per P
        acc = {}
        for (P, S), r in results.sums.items():
            a = acc.setdefault(P, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            a[0] += r["squares"]
            a[1] += model.squares(P, S, r["nvecs_raw"], pinfo.smin(P))
        for q in results.squares.values():
            P, S = q["P"], q["S"]
            smin = pinfo.smin(P)
            a = acc.setdefault(P, [0.0] * 6)
            a[2] += q["s_count"]
            a[3] += T * model.p_s(P, S, smin)
            a[4] += q["p_count"]
            a[5] += T * model.p_p(P, S, smin)
        for P, lg in results.legacy.items():
            # the legacy search only gives totals per P: use its traversal
            # counts (its squares are used by the global model)
            if lg["sols"] == 0:
                continue
            smin = pinfo.smin(P)
            S_mid = (smin + min(lg["maxS"], 3 * smin)) / 2
            a = acc.setdefault(P, [0.0] * 6)
            a[2] += lg["nS"]
            a[3] += T * lg["sols"] * model.p_s(P, S_mid, smin)
            a[4] += lg["nP"]
            a[5] += T * lg["sols"] * model.p_p(P, S_mid, smin)
        self.acc = acc

    def factor(self, P, optimistic=True):
        """per-P multipliers (rate, p_S, p_P) from gamma posteriors"""
        a = self.acc.get(P, [0.0] * 6)

        def post(obs, exp, prior):
            # gamma(prior, prior) prior on the multiplier; observed counts are
            # Poisson with mean multiplier * expected
            shape, rate = prior + obs, prior + exp
            mean = shape / rate
            sd = math.sqrt(shape) / rate
            return mean + (EXPLORE * sd if optimistic else 0.0)

        return (post(a[0], a[1], PRIOR_SQUARES), post(a[2], a[3], PRIOR_TRAV),
                post(a[4], a[5], PRIOR_TRAV))

    def score(self, P, S, N):
        """expected magic squares per CPU-second at (P, S), and its parts"""
        smin = self.pinfo.smin(P)
        if N < 2 * self.n:
            return 0.0, 0.0, 0.0
        fr, fs, fp = self.factor(P, self.optimistic)
        rate = fr * self.model.rate(P, S, N, smin)
        pm = self.model.p_magic(self.n, min(1, fs * self.model.p_s(P, S, smin)),
                                min(1, fp * self.model.p_p(P, S, smin)))
        return rate * pm, rate, pm

    def prior(self, P):
        """rough score of P from the vector counts already known (computed
        with S_min): the best predicted yield among the next sums"""
        self.pinfo.smin(P)
        d = self.pinfo.data[self.pinfo.key(P)]
        lo = self.results.frontier(P, self.pinfo)
        hi = min(lo + 400, d["counted_to"])
        best = 0.0
        for S in range(lo, hi + 1, max(1, (hi - lo) // 24)):
            N = d["counts"].get(str(S), 0)
            if N >= 2 * self.n:
                best = max(best, self.score(P, S, N)[0])
        return best

    def next_unit(self, P, unit_time, max_sums=2000, drop=None):
        """the next chunk of sums of P (about unit_time CPU-seconds), with its
        predicted CPU time, semi-magic squares and magic squares; its score
        is the predicted magic squares per CPU-second.

        The chunk also ends where the predicted magic squares per CPU-second
        of the next sum fall below `drop` times the best of the chunk so far:
        past a P's best sums the density falls fast (time grows like N^5.5,
        P(magic) falls like S^-8), and the rest of the stretch should compete
        again with the other values of P rather than ride along (in forecasts,
        units of a fixed 600 s collected ~20% less than ideal at 2000
        CPU-hours, and 3600 s ~55% less)."""
        if drop is None:
            drop = UNIT_DROP
        lo = self.results.frontier(P, self.pinfo)
        stop = self.results.next_covered(P, lo)  # fill a gap, don't redo
        counts = self.pinfo.counts(P, lo, lo + 200)
        total_t, total_sq, total_m, S = 0.0, 0.0, 0.0, lo
        best = 0.0
        hi = lo
        while stop is None or S < stop:
            if S not in counts:
                counts.update(self.pinfo.counts(P, S, S + 200))
            N = counts[S]
            if N >= 2 * self.n:
                t = self.model.sum_time(N)
                m, rate, _ = self.score(P, S, N)
                if total_t > 0 and total_t + t > unit_time:
                    break
                if total_t > 0 and m < drop * best:
                    break
                best = max(best, m)
                total_t += t
                total_sq += rate * t
                total_m += m * t
            hi = S
            S += 1
            if S - lo >= max_sums:
                break
        return Unit(P, lo, hi, total_m / total_t if total_t > 0 else 0.0, total_t,
                    total_sq, total_m)


# --------------------------------------------------------------------------
# candidate values of P


def gen_candidates(n=6, tau_min=800, tau_max=12000, max_big_primes=2):
    """exponent tuples over 2, 3, 5, 7, 11, 13, 17, 19 that look like the
    values of P that worked well before: decreasing exponents for 2, 3, 5, 7,
    and exponent 0 or 1 for at most max_big_primes of 11, 13, 17, 19"""
    out = []
    if n == 6:
        ranges = [range(8, 19), range(3, 10), range(2, 7), range(1, 5)]
    elif n == 5:
        ranges = [range(5, 15), range(2, 8), range(1, 5), range(0, 4)]
    else:
        ranges = [range(8, 19), range(4, 11), range(2, 7), range(1, 5)]
    for e2 in ranges[0]:
        for e3 in ranges[1]:
            if e3 > e2:
                continue
            for e5 in ranges[2]:
                if e5 > e3:
                    continue
                for e7 in ranges[3]:
                    if e7 > e5:
                        continue
                    for mask in range(16):
                        tail = [(mask >> i) & 1 for i in range(4)]
                        if sum(tail) > max_big_primes:
                            continue
                        P = norm_p([e2, e3, e5, e7] + tail)
                        if tau_min <= tau(P) <= tau_max:
                            out.append(P)
    return sorted(set(out))


# --------------------------------------------------------------------------
# running


class Scheduler:
    def __init__(self, args, need_model=True):
        self.args = args
        self.n = args.vec_size
        self.dir = args.state
        os.makedirs(os.path.join(self.dir, "units"), exist_ok=True)
        self.pinfo = PInfo(os.path.join(self.dir, f"pinfo_{self.n}.json"), self.n)
        self.results = Results(self.n)
        self.model_path = os.path.join(self.dir, f"model_{self.n}.json")
        legacy_path = os.path.join(self.dir, "legacy.json")
        if os.path.exists(legacy_path) and not getattr(args, "no_legacy", False):
            with open(legacy_path) as f:
                for r in json.load(f):
                    if r.get("n", 6) == self.n:
                        self.results.legacy[norm_p(r["P"])] = r
        for path in unit_files(os.path.join(self.dir, "units")):
            self.results.add_file(path)
        self.model = Model.load(self.model_path, self.n)
        if self.model is None and need_model:
            sys.exit(f"no model: run `scheduler.py fit` once some results exist in "
                     f"{self.dir}/units (or `ingest` some msearch output)")

    def candidates(self):
        cands = set(gen_candidates(self.n, self.args.tau_min, self.args.tau_max))
        if self.args.only:
            cands = {parse_p(p) for p in self.args.only.split(";")}
        cands |= {P for (P, S) in self.results.sums}
        # values of P whose sums were all searched by the legacy program
        exhausted = {P for P, r in self.results.legacy.items() if r["maxS"] >= 10 ** 9}
        return sorted(cands - exhausted)

    def parallel(self, fn, items):
        with ThreadPoolExecutor(max(1, self.args.workers)) as ex:
            return list(ex.map(fn, items))

    def prepare(self, cands, scorer=None, keep=()):
        """S_min for all candidates, then vector counts for the most promising
        ones; returns the active set of P: the best by prior, the P with
        results, and the P in keep"""
        def known(P):
            d = self.pinfo.data.get(self.pinfo.key(P))
            return d is not None and (d.get("window")
                                      or d["counted_to"] >= int(PInfo.WINDOW * d["smin"]))

        missing = [P for P in cands if not known(P)]
        if missing:
            log(f"computing S_min and vector counts near it for {len(missing)} values of P...")

            def smin(P):
                info = PInfo(self.pinfo.path + ".scratch", self.n)
                info.data = {}
                info.smin(P)
                return info.data

            for d in self.parallel(smin, missing):
                self.pinfo.data.update(d)
            self.pinfo.save()
        scorer = scorer or Scorer(self.model, self.results, self.pinfo)
        ranked = sorted(cands, key=lambda P: -scorer.prior(P))
        active = ranked[:self.args.active]
        active = sorted(set(active) | {P for (P, S) in self.results.sums if P in cands}
                        | set(keep))
        need = [P for P in active
                if self.pinfo.data[self.pinfo.key(P)]["counted_to"]
                < self.results.frontier(P, self.pinfo) + 100]
        if need:
            log(f"computing vector counts for {len(need)} values of P...")

            def counts(P):
                info = PInfo(self.pinfo.path + ".scratch", self.n)
                info.data = {info.key(P): self.pinfo.data[info.key(P)]}
                lo = self.results.frontier(P, self.pinfo)
                info.counts(P, lo, lo + 300)
                return info.data

            for d in self.parallel(counts, need):
                self.pinfo.data.update(d)
            self.pinfo.save()
        return active

    def ranking(self, cands, scorer):
        units = [scorer.next_unit(P, self.args.unit_time) for P in cands]
        return sorted(units, key=lambda u: -u.score)

    def unit_path(self, P, lo, hi):
        seq = len(glob.glob(os.path.join(self.dir, "units", "*.jsonl")))
        name = f"{int(time.time())}_{seq:06d}_{p_str(P, '_')}_{lo}_{hi}.jsonl"
        return os.path.join(self.dir, "units", name)

    def command(self, P, lo, hi, out):
        return [binary("msearch"), "--vec-size", str(self.n), "--min-sum", str(lo),
                "--max-sum", str(hi), "--time-limit", str(self.args.unit_time * 2),
                "--node-limit", str(self.args.node_limit), "--out", out, *map(str, P)]

    def run(self):
        all_cands = self.candidates()
        scorer = Scorer(self.model, self.results, self.pinfo)
        cands = self.prepare(all_cands, scorer)
        deadline = time.time() + self.args.hours * 3600 if self.args.hours else None
        running = {}  # Popen -> (P, path)
        units_done = 0
        seen_squares = set(self.results.squares)
        stopping = False
        cache = {}  # P -> next Unit of P

        def refresh(P):
            cache[P] = scorer.next_unit(P, self.args.unit_time)

        def handle_sigint(sig, frame):
            nonlocal stopping
            stopping = True
            log("stopping: waiting for running units (Ctrl-C again to kill them)")
            signal.signal(signal.SIGINT, signal.default_int_handler)

        signal.signal(signal.SIGINT, handle_sigint)
        log(f"{len(all_cands)} candidate values of P, {len(cands)} active, "
            f"{self.args.workers} workers")
        for P in cands:
            refresh(P)
        try:
            while True:
                # collect finished units
                for proc in list(running):
                    if proc.poll() is None:
                        continue
                    P, path = running.pop(proc)
                    units_done += 1
                    self.results.add_file(path)
                    for h, q in self.results.squares.items():
                        if h in seen_squares:
                            continue
                        seen_squares.add(h)
                        tag = "MAGIC SQUARE" if q["best_score"] >= 14 else "square"
                        if q["best_score"] >= self.args.announce_score:
                            log(f"{tag}: P={p_str(q['P'])} S={q['S']} best={q['best_score']} "
                                f"#S={q['s_count']} #P={q['p_count']} #SP={q['sp_count']} "
                                f"{q['grid']}")
                    scorer = Scorer(self.model, self.results, self.pinfo)
                    if units_done % self.args.refit_every == 0:
                        try:
                            self.model, diag = fit_model(self.results, self.pinfo)
                            self.model.save(self.model_path)
                            log(describe_fit(diag))
                        except SystemExit as e:
                            log(str(e))
                        scorer = Scorer(self.model, self.results, self.pinfo)
                        cands = self.prepare(all_cands, scorer)
                        for Q in cands:
                            refresh(Q)
                    else:
                        refresh(P)
                if stopping or (deadline and time.time() > deadline):
                    if not running:
                        break
                    time.sleep(0.5)
                    continue
                # start new units
                busy = {P for P, _ in running.values()}
                while len(running) < self.args.workers:
                    free = [cache[P] for P in cands if P not in busy and P in cache]
                    if not free:
                        break
                    u = max(free, key=lambda u: u.score)
                    if u.score <= 0:
                        break
                    path = self.unit_path(u.P, u.lo, u.hi)
                    proc = subprocess.Popen(self.command(u.P, u.lo, u.hi, path), env=msearch_env(),
                                            stdout=subprocess.DEVNULL)
                    running[proc] = (u.P, path)
                    busy.add(u.P)
                    log(f"start P={p_str(u.P)} S={u.lo}..{u.hi} (predicted {u.time:.0f}s, "
                        f"{u.squares:.1f} squares, {u.score * 3.15e7:.3g} magic/CPU-year)")
                time.sleep(0.5)
        except KeyboardInterrupt:
            for proc in running:
                proc.terminate()
            for proc, (P, path) in running.items():
                proc.wait()
                self.results.add_file(path)
        self.pinfo.save()
        log("stopped")
        report(self.results, self.pinfo, self.model, top=10)


# --------------------------------------------------------------------------
# reporting


def report(results, pinfo, model, top=20):
    n = results.n
    sums = list(results.sums.values())
    cpu = sum(r["time"] + r["setup_time"] for r in sums)
    sq = list(results.squares.values())
    print(f"\n{len(results.by_p())} values of P, {len(sums)} sums, {cpu / 3600:.2f} CPU-hours, "
          f"{len(sq)} semi-magic squares ({3600 * len(sq) / max(cpu, 1):.0f}/CPU-hour)")
    if results.dsums:
        dcpu = sum(r.get("cpu", r["time"]) for v in results.dsums.values() for r in v)
        # the pairs of the sums searched in full (their first complete
        # record), estimated over every d (dsum_est_pairs)
        full = [next((r for r in v if r.get("complete")), None) for v in results.dsums.values()]
        est = [dsum_est_pairs(r) for r in full if r is not None]
        print(f"d-first (not in the fits): {len(results.dsums)} sums, {dcpu / 3600:.2f} "
              f"CPU-hours, {len(results.dmagic)} magic squares, "
              f"{sum(e for e in est if e is not None):.1f} (square, SP traversal) pairs "
              f"estimated over {sum(e is not None for e in est)} sums searched in full")
        for q in sorted(results.dmagic.values(), key=lambda q: (q["S"], q["hash"])):
            print(f"  MAGIC (d-first) P={p_str(q['P']):18} S={q['S']:5}  {q.get('grid')}")
    types = {}
    for q in sq:
        types[q["best_score"]] = types.get(q["best_score"], 0) + 1
    names = {0: "0", 2: "S", 3: "P", 4: "S+S", 5: "S+P", 6: "P+P", 7: "SP",
             9: "SP+S", 10: "SP+P", 14: "MAGIC (SP+SP)"}
    print("best pair of diagonals: " + ", ".join(
        f"{names.get(k, k)}: {v}" for k, v in sorted(types.items(), reverse=True)))
    print(f"traversals: S {sum(q['s_count'] for q in sq)}, P {sum(q['p_count'] for q in sq)}, "
          f"SP {sum(q['sp_count'] for q in sq)}")
    if model:
        exp_magic = 0.0
        for q in sq:
            smin = pinfo.smin(q["P"])
            exp_magic += model.p_magic(n, model.p_s(q["P"], q["S"], smin),
                                       model.p_p(q["P"], q["S"], smin))
        print(f"expected magic squares among these (model): {exp_magic:.2e}")
    best = sorted(sq, key=lambda q: (-q["best_score"], -q["sp_count"], q["S"]))[:top]
    if best:
        print("\nbest squares:")
    for q in best:
        print(f"  P={p_str(q['P']):18} S={q['S']:5} best={q['best_score']:2} "
              f"#S={q['s_count']} #P={q['p_count']} #SP={q['sp_count']}  {q['grid']}")


# --------------------------------------------------------------------------
# commands


def cmd_fit(args):
    sch = Scheduler.__new__(Scheduler)
    sch.args = args
    sch.n = args.vec_size
    pinfo = PInfo(os.path.join(args.state, f"pinfo_{args.vec_size}.json"), args.vec_size)
    results = Results(args.vec_size)
    for path in unit_files(os.path.join(args.state, "units")):
        results.add_file(path)
    legacy_path = os.path.join(args.state, "legacy.json")
    if os.path.exists(legacy_path) and not args.no_legacy:
        with open(legacy_path) as f:
            for r in json.load(f):
                if r.get("n", 6) == args.vec_size:
                    results.legacy[norm_p(r["P"])] = r
    model, diag = fit_model(results, pinfo)
    model.save(os.path.join(args.state, f"model_{args.vec_size}.json"))
    pinfo.save()
    print(describe_fit(diag))
    print(json.dumps(diag, indent=1))
    print(json.dumps(model.c, indent=1))


def cmd_plan(args):
    sch = Scheduler(args)
    scorer = Scorer(sch.model, sch.results, sch.pinfo)
    cands = sch.prepare(sch.candidates(), scorer)
    units = sch.ranking(cands, scorer)
    print(f"{'P':20} {'S range':>13} {'pred. time':>10} {'squares/h':>9} "
          f"{'P(magic)':>9} {'magic/CPU-year':>14}")
    for u in units[:args.top]:
        h = u.time / 3600
        print(f"{p_str(u.P):20} {u.lo:6}-{u.hi:<6} {u.time:9.0f}s "
              f"{u.squares / max(h, 1e-9):9.0f} {u.magic / max(u.squares, 1e-30):9.2e} "
              f"{u.score * 3.15e7:14.3g}")


def simulate(sch, scorer, all_cands, unit_time, max_units=None, max_hours=None,
             refresh_every=200):
    """greedily simulate the scheduler, assuming every unit takes its
    predicted time and finds its predicted number of squares; yields units.
    Like `run`, the active set of P is refreshed regularly, so new values of
    P come in as the best ones are used up."""
    started = set()
    cands = sch.prepare(all_cands, scorer)
    cache = {P: scorer.next_unit(P, unit_time) for P in cands}
    hours, i = 0.0, 0
    while (max_units is None or i < max_units) and (max_hours is None or hours < max_hours):
        if i and i % refresh_every == 0:
            for P in sch.prepare(all_cands, scorer, keep=started):
                if P not in cache:
                    cache[P] = scorer.next_unit(P, unit_time)
        u = max(cache.values(), key=lambda u: u.score)
        if u.score <= 0:
            break
        yield u
        hours += u.time / 3600
        i += 1
        started.add(u.P)
        sch.results.mark_covered(u.P, u.lo, u.hi)
        cache[u.P] = scorer.next_unit(u.P, unit_time)


def cmd_emit(args):
    """static plan: print msearch arguments for the units the scheduler would
    run (one per line), simulating it with predicted times and squares"""
    sch = Scheduler(args)
    scorer = Scorer(sch.model, sch.results, sch.pinfo)
    cands = sch.candidates()
    units_dir = os.path.join(args.state, "units")
    if os.path.abspath(units_dir).startswith(ROOT + os.sep):
        units_dir = os.path.relpath(units_dir, ROOT)  # portable plan
    for i, u in enumerate(simulate(sch, scorer, cands, args.unit_time, max_units=args.units)):
        out = os.path.join(units_dir, f"plan{i:06d}_{p_str(u.P, '_')}_{u.lo}_{u.hi}.jsonl")
        print(" ".join(sch.command(u.P, u.lo, u.hi, out)[1:]))


def cmd_forecast(args):
    """predicted squares and magic squares for the next CPU-hours of search"""
    if args.hours is None:
        args.hours = 8766.0
    sch = Scheduler(args)
    cands = sch.candidates()
    # predicted numbers of squares, without the exploration bonus
    scorer = Scorer(sch.model, sch.results, sch.pinfo, optimistic=False)
    checkpoints = [args.hours * f for f in (0.01, 0.03, 0.1, 0.3, 1.0)]
    hours = squares = magic = 0.0
    used = set()
    print(f"{'CPU-hours':>10} {'squares':>12} {'magic squares':>14} {'P used':>7}  latest unit")
    for u in simulate(sch, scorer, cands, args.unit_time, max_hours=args.hours):
        hours += u.time / 3600
        squares += u.squares
        magic += u.magic
        used.add(u.P)
        while checkpoints and hours >= checkpoints[0]:
            checkpoints.pop(0)
            print(f"{hours:10.1f} {squares:12.0f} {magic:14.3g} {len(used):7}  "
                  f"P={p_str(u.P)} S={u.lo}..{u.hi}")
    if magic > 0:
        print(f"predicted CPU-years per magic square at this pace: {hours / magic / 8766:.3g}")


def cmd_report(args):
    sch = Scheduler(args, need_model=False)
    report(sch.results, sch.pinfo, sch.model, top=args.top)


def cmd_ingest(args):
    os.makedirs(os.path.join(args.state, "units"), exist_ok=True)
    for path in args.files:
        dst = os.path.join(args.state, "units", "ingested_" + os.path.basename(path))
        if not dst.endswith(".jsonl"):
            dst += ".jsonl"
        shutil.copy(path, dst)
        print(f"{path} -> {dst}")


def cmd_import_legacy(args):
    """import per-P totals from a postprocess.py summary (stats_short.txt):
    sums S_min(P) <= S <= max_S were searched by the old program"""
    rows = []
    lines = open(args.stats).read().split("\n")
    start = next(i for i, l in enumerate(lines)
                 if l.startswith("Statistics for each P, sorted by factorization"))
    for l in lines[start + 2:]:
        if not l.startswith("  ") or l.startswith("  P "):
            break
        parts = l.split()
        for t in (15, 14):  # with / without the "best" column
            if len(parts) <= t:
                continue
            try:
                P = tuple(map(int, parts[:len(parts) - t]))
            except ValueError:
                continue
            if not P or len(P) > len(PRIMES) or max(P) > 63:
                continue
            if str(p_value(P)) != parts[len(parts) - t]:
                continue
            v = parts[len(parts) - t:]
            max_s = int(v[3])
            rows.append({"n": args.vec_size, "P": list(norm_p(P)),
                         "maxS": 10 ** 9 if max_s >= 99999 else max_s,
                         "count": int(v[4]), "time": float(v[5]), "sols": int(v[6]),
                         "nS": int(v[7]), "nP": int(v[8]), "nSP": int(v[9])})
            break
    os.makedirs(args.state, exist_ok=True)
    with open(os.path.join(args.state, "legacy.json"), "w") as f:
        json.dump(rows, f)
    print(f"imported {len(rows)} values of P "
          f"({sum(r['sols'] for r in rows)} squares) from {args.stats}")
    if getattr(args, "model", "regression") == "analytic":
        return  # the analytic model needs no vector counts
    # vector counts over the legacy ranges, so `fit` can use these totals
    pinfo = PInfo(os.path.join(args.state, f"pinfo_{args.vec_size}.json"), args.vec_size)
    todo = [norm_p(r["P"]) for r in rows if r["sols"] >= LEGACY_MIN_SQUARES and r["P"]]
    print(f"computing vector counts over the legacy ranges of {len(todo)} values of P...")

    def work(P):
        info = PInfo(pinfo.path + ".scratch", args.vec_size)
        k = info.key(P)
        info.data = {k: pinfo.data[k]} if k in pinfo.data else {}
        smin = info.smin(P)
        r = next(r for r in rows if norm_p(r["P"]) == P)
        hi = legacy_range(r, smin) if smin else None
        if hi is not None:
            info.counts(P, smin, hi)
        return info.data

    with ThreadPoolExecutor(args.workers) as ex:
        for d in ex.map(work, todo):
            pinfo.data.update(d)
    pinfo.save()


# ==========================================================================
# scheduler v2 (--model analytic, the default): the analytic model of the
# existence study over a wide pool of P, a lazy greedy planner, learning by
# cells, and an incremental summary of the results


def _np():
    import numpy
    return numpy


def _am():
    import amodel
    return amodel


# --------------------------------------------------------------------------
# calibration cells: N' band x k x assignment ratio x x = ln(S / smin_approx)

NBAND_EDGES = (1000, 3000, 6000, 12000, 24000)
XBAND_EDGES = (0.1, 0.25)
NBAND_NAMES = ("<1k", "1-3k", "3-6k", "6-12k", "12-24k", ">=24k")
KBAND_NAMES = ("k<=5", "k=6", "k=7", "k>=8")
RBAND_NAMES = ("sorted", "r<=1.1", "r>1.1")
XBAND_NAMES = ("x<0.1", "x 0.1-0.25", "x>=0.25")
NCELLS = 6 * 4 * 3 * 3
LN12K = math.log(12000)
LN8K = math.log(8000)
# class rates of S and P traversals that the shipped P(magic | square)
# corresponds to (analytic.md 4.3: 0.86 and 0.64 near S_min)
TRAV_BASE = (0.86, 0.64)
# prior strength of the per-P gamma factors: A_SQ = 6 is a between-P sd of
# 1/sqrt(6) = 0.41 in ln(squares obs/pred), as measured by the calibration
# search in the target region (research/calibration-target.md 3.1: sd 0.46
# over its plain sums, 0.40 after stratum means; its GLM refit uses A_SQ 6).
# It was 15 (sd 0.26, the per-unit sd of earlier held-out units), with
# which 15-18% of those sums fell outside their 90% interval. Squares within
# one P are overdispersed too (per-sum Pearson chi2/df 2.5), so their counts
# enter the posterior divided by PHI_SUM. The measured 0.41 is the spread
# beyond Poisson only (analyze.py's between_sd subtracts 1/E), so it also
# holds the per-sum overdispersion that PHI_SUM models: for a one-sum P of
# E ~30 squares the dispersion is 2.5 + 30/6 = 7.5 against the measured
# Pearson 5.4 at 3-6k, a little wide. Taking PHI_SUM out, 0.41^2 - (PHI_SUM -
# 1)/E gives A ~8-9; 6 is the write-up's refit choice (its glm.out), kept
A_SQ, A_TRAV = 6.0, 30.0
PHI_SUM = 2.5
# selection: the forecast's E is no longer discounted (it was x0.8, from
# held-out units v2 ranked highest at 0.7-0.9 of prediction). Within the
# 3-6k band the calibration search found no winner's curse: the top
# quartile by predicted magic per CPU-second came in at 0.96 [0.82, 1.11],
# the slope of ln(obs/pred) on ln density +0.020 +- 0.027 (244 sums); the
# earlier shortfall is the 6-12k squares correction, which the prior now
# carries (calibration-target.md 5). What remains goes into the band as
# lognormal(0, SELECTION_SD), the write-up's selection term.
SELECTION_SD = 0.12
EXPLORE_V2 = 0.0            # optimism: posterior mean + EXPLORE * sd

# d-first units (msearch --diag-first; research/scheduler-v2.md, "d-first
# units"). --dfirst auto: a sum is searched d-first where the measured
# d-first / plain ratio (amodel.DFIRST_RATIO_PRIOR, its level learned
# online) with the calibration stream is below 1, (1 + CALIB_FRAC) r(N') <
# 1, i.e. N' >= ~3.7k with engine 4 (~4.3k with engine 3), and its N' >=
# --dfirst-min-n. The two time laws
# only price the sums (their quotient is not label-free and the plain law
# under-predicts at 5-8k: it put the switch at ~9-10k).
DFIRST_POLICIES = ("auto", "off", "on")
DFIRST_MIN_NP = 2000.0
# the calibration stream of a d-first sum (--calib-r1-stride k, a plain
# search of every k-th first row: semi-magic squares for the models) costs
# about this share of the sum's predicted d-first CPU; k = round(1 /
# (CALIB_FRAC r(N'))) (~18 at N' 6k, ~41 at 25k), by the ratio law, so that
# the plain law's errors do not change the stream's cost. retrospective.md
# 7 assumed 10% (T11: 5% gives E x1.01 at 1-1,000 CPU-years); at the
# frontier's ~10^3 squares per sum it records 20-60 squares per sum.
CALIB_FRAC = 0.07
# ... but the stream (which runs in one piece, after the unit's d loop)
# takes at most DFIRST_STREAM_MAX x --unit-time (k raised to fit; above N'
# ~15k at --unit-time 120 it then costs less than CALIB_FRAC), and the
# first unit of a split sum, which carries it, gets at least
# DFIRST_FIRST_MIN x --unit-time of d loop besides
DFIRST_STREAM_MAX = 1.0
DFIRST_FIRST_MIN = 0.25
# a unit of d (--d-range) writes a "dchunk" checkpoint every 1/DFIRST_CHUNKS
# of its d (at least 8 d): --time-limit acts between chunks, and the parts
# of a killed unit count
DFIRST_CHUNKS = 8
# a d-first sum predicted to take more than DFIRST_SPLIT x --unit-time is
# split into units of d (--d-range lo:hi, each ~--unit-time); the last part
# is extended to the end of the sum rather than leave a tail below
# DFIRST_TAIL x a unit
DFIRST_SPLIT = 1.5
DFIRST_TAIL = 0.25
# a part of a d-first sum ("dsum" with complete 0) enters the d-first time
# law when it searched at least this many d (its CPU scaled to the whole sum
# by the cost profile amodel.DFIRST_COST_PROFILE, with weight its share of
# the d loop, so that a sum split into many parts weighs as one)
DFIRST_TIME_MIN_ND = 32
# ... and its sum has at least this many d (vectors before reduction): the
# law was fitted at N 4-32k, and below ~2k the overheads outside the d loop
# (enumeration, process start) dominate
DFIRST_TIME_MIN_N = 2000
# the calibration stream's squares enter the squares GLM with their
# expectation E_sq / k, weighted down by their sampling dispersion phi_s =
# se_squares^2 / (k est_squares) (squares sharing a first row come
# together): w = min(1, PHI_SUM / (phi_s + PHI_SUM / k)); phi_s is taken as
# CALIB_PHI where the stream sampled fewer than 5 squares
CALIB_PHI = 2.5

# stage 1 (--stage1 HOURS[:LO:HI], research/stage1.md): for the first HOURS
# of reference CPU (the units' predicted CPU-seconds at the scheduler's laws,
# the fast x86 build's clock, counted from the start of stage 1 in the
# state's stage1_<n>.json) the sums with N' in [LO, HI) are searched plain
# (as --dfirst off for them), in the scheduler's own order; everything else
# is unchanged. Plain search records every (square, SP traversal) pair of a
# sum at weight 1 (d-first records them through the star cover, 1 in K d,
# and only where the summary holds d-first pairs), so the SP coupling f_rho
# is learned fastest from plain squares at N' 3-6k, where d-first / plain
# is 0.9-1.0 and almost no E is lost (scripts/decide.py fits it).
STAGE1_DEFAULT = (60.0, 3000.0, 6000.0)


def parse_stage1(spec):
    """--stage1 HOURS[:LO:HI] -> (hours, lo, hi), or None for off / 0 h"""
    if spec is None or str(spec).strip().lower() in ("", "off", "none", "0"):
        return None
    parts = str(spec).split(":")
    if len(parts) not in (1, 3):
        raise ValueError(f"--stage1 {spec}: HOURS or HOURS:LO:HI")
    h = float(parts[0])
    lo, hi = (STAGE1_DEFAULT[1], STAGE1_DEFAULT[2]) if len(parts) == 1 else map(float, parts[1:])
    if h <= 0:
        return None
    if not 0 < lo < hi:
        raise ValueError(f"--stage1 {spec}: need 0 < LO < HI")
    return (h, lo, hi)


class Stage1:
    """the stage-1 clock: the spec (hours, lo, hi) and the reference
    CPU-seconds of the units launched under it (persisted in the state for
    `run`; `plan`, `emit` and `forecast` start from the state's value and
    count their simulated units without saving)"""

    def __init__(self, state, n, spec, scale=1.0):
        self.path = os.path.join(state, f"stage1_{n}.json")
        self.spec = spec
        self.scale = scale        # forecast --sample: the budget is scaled
        self.spent = 0.0
        try:
            with open(self.path) as f:
                d = json.load(f)
            if spec is not None and tuple(d.get("spec", ())) == tuple(spec):
                self.spent = float(d.get("spent", 0.0))
        except (OSError, ValueError):
            pass

    @property
    def limit(self):
        return 0.0 if self.spec is None else self.spec[0] * 3600.0 * self.scale

    @property
    def active(self):
        return self.spec is not None and self.spent < self.limit

    def band(self):
        """(ln LO, ln HI) while stage 1 is on, else None"""
        if not self.active:
            return None
        return (math.log(self.spec[1]), math.log(self.spec[2]))

    def charge(self, seconds):
        """count a unit; True when this ends stage 1"""
        was = self.active
        self.spent += float(seconds)
        return was and not self.active

    def save(self):
        if self.spec is None:
            return
        tmp = self.path + ".tmp"
        with open(tmp, "w") as f:
            json.dump({"spec": list(self.spec), "spent": self.spent,
                       "hours": self.spent / 3600.0}, f)
        os.replace(tmp, self.path)


def ratio_bin(r):
    np = _np()
    r = np.asarray(r, float)
    return np.where(r <= 1.0 + 1e-6, 0, np.where(r <= np.float32(1.1), 1, 2))


def cell_index(lNp, k, rbin, x):
    """calibration cell (0..215) of sums with bias-corrected log N', k
    distinct primes, ratio bin and x = ln(S / smin_approx)"""
    np = _np()
    nb = np.searchsorted(np.log(NBAND_EDGES), lNp, side="right")
    kb = np.clip(np.asarray(k) - 5, 0, 3)
    xb = np.searchsorted(XBAND_EDGES, x, side="right")
    return ((nb * 4 + kb) * 3 + rbin) * 3 + xb


def cell_parts(c):
    c = int(c)
    return c // 36, (c // 9) % 4, (c // 3) % 3, c % 3


# columns of the main-effects design (one row per cell)
EFFECTS = ("intercept", "N'>=3k", "N'<1k", "N' 6-12k", "N' 12-24k", "N'>=24k",
           "k<=5", "k>=7", "k>=8", "r<=1.1", "r>1.1", "x<0.1", "x>=0.25")


def cell_design():
    np = _np()
    Z = np.zeros((NCELLS, len(EFFECTS)))
    for c in range(NCELLS):
        nb, kb, rb, xb = cell_parts(c)
        Z[c] = [1, nb >= 2, nb == 0, nb == 3, nb == 4, nb == 5, kb == 0, kb >= 2, kb == 3,
                rb == 1, rb == 2, xb == 0, xb == 2]
    return Z


def _prior(intercept, sd0, sd, **effects):
    mean = [intercept] + [0.0] * (len(EFFECTS) - 1)
    sds = [sd0] + [sd] * (len(EFFECTS) - 1)
    for name, (m, s) in effects.items():
        i = EFFECTS.index(name)
        mean[i], sds[i] = m, s
    return mean, sds


# priors of the three GLMs (log rate = intercept + effects): squares relative
# to the analytic model x SQ12; S and P traversals per 720 e^{lpS}, e^{lpP}.
# They are the posterior of the calibration search's GLM refit
# (research/calibration-target.md section 5, forecast_update.py "glm" with
# A_SQ 6, archive/results/glm.out and glm_update.json's calib_post; means
# and sds per effect, the correlations dropped): 321 pre-registered sums at
# N' 3-45k fitted on top of the earlier priors, which were: squares N' < 1k
# -0.25 (0.74 on 146 P), x < 0.1 -0.10 (S_min, x < 0.05: 0.65), 6-12k -0.22
# (held-out live units at 0.78, analytic.md 3.4's 0.78 at N 3.4-8.9k); S
# ln 0.86 + 0.04 at N' >= 3k (the review's +0.09 lowered by the live units);
# P ln 0.64 + 0.06 at N' >= 3k, the review factors K7 (e^-0.16) and X01
# (e^-0.08) on P(magic | square) as k >= 7 -0.08 and x < 0.1 -0.04. The
# search's raw ratios against those: squares 1.04 at 3-6k and 1.28 at 6-12k
# (the x0.80 refits to x0.93), S traversals per square 0.939 [0.918, 0.960]
# at 3-6k and 0.973 at 6-12k, P 1.02 and 1.08 (sections 3.1, 3.2). The
# refit puts the S deficit mostly into k <= 5, ratio <= 1.1 and x < 0.1 and
# the P excess into N' >= 3k; these main effects also act below 3k, where
# nothing was measured. At the means m(c) below is e^{2 (0.024 + 0.109)} =
# 1.30 x PAIR at N' 3-6k for k 6 and the middle ratio and x bands (it was
# e^{2 (0.04 + 0.06)} = 1.22).
# effects: intercept, N'>=3k, N'<1k, N' 6-12k, N' 12-24k, N'>=24k, k<=5,
#          k>=7, k>=8, r<=1.1, r>1.1, x<0.1, x>=0.25 (EFFECTS)
GLM_PRIORS = {
    "sq": ([-0.003, -0.002, -0.25, -0.078, 0.058, -0.061, -0.006, 0.033, 0.021, 0.060, -0.017,
            -0.088, 0.060],
           [0.196, 0.194, 0.150, 0.085, 0.190, 0.232, 0.073, 0.112, 0.237, 0.067, 0.117, 0.078,
            0.079]),
    "S": ([-0.159, 0.032, 0.0, -0.018, -0.054, 0.009, -0.041, -0.038, -0.004, -0.048, -0.008,
           -0.067, 0.002],
          [0.073, 0.073, 0.100, 0.050, 0.093, 0.099, 0.043, 0.062, 0.098, 0.039, 0.066, 0.050,
           0.045]),
    "P": ([-0.422, 0.085, 0.0, 0.042, 0.008, -0.002, -0.042, -0.072, -0.003, -0.006, 0.029,
           -0.102, 0.018],
          [0.075, 0.075, 0.100, 0.063, 0.098, 0.100, 0.054, 0.085, 0.100, 0.051, 0.084, 0.067,
           0.056]),
}
# (the earlier priors, for the record: shipped until the round-2 integration)
GLM_PRIORS_V2 = {
    "sq": _prior(0.0, 0.3, 0.25, **{"N'<1k": (-0.25, 0.15), "x<0.1": (-0.1, 0.15),
                                    "N' 6-12k": (-0.22, 0.15)}),
    "S": _prior(math.log(TRAV_BASE[0]), 0.1, 0.1, **{"N'>=3k": (0.04, 0.1)}),
    "P": _prior(math.log(TRAV_BASE[1]), 0.1, 0.1, **{"N'>=3k": (0.06, 0.1), "k>=7": (-0.08, 0.1),
                                                     "x<0.1": (-0.04, 0.1)}),
}


def fit_poisson_map(Z, y, offset, prior_mean, prior_sd, iters=50, phi=None):
    """MAP of a (quasi-)Poisson GLM, log mu = offset + Z beta, with
    independent normal priors on beta (Newton with step halving); phi: the
    dispersion per row (variance = phi mu; default 1), which divides each
    row's log-likelihood; returns (beta, the Laplace covariance)"""
    np = _np()
    Z = np.asarray(Z, float)
    y = np.asarray(y, float)
    off = np.asarray(offset, float)
    m = np.asarray(prior_mean, float)
    lam = 1.0 / np.asarray(prior_sd, float) ** 2
    w = np.ones(len(y)) if phi is None else 1.0 / np.maximum(np.asarray(phi, float), 1e-9)

    def logpost(b):
        eta = off + Z @ b
        return float((w * (y * eta - np.exp(np.minimum(eta, 50)))).sum()
                     - 0.5 * (lam * (b - m) ** 2).sum())

    b = m.copy()
    lp = logpost(b)
    H = Z.T @ ((w * np.exp(np.minimum(off + Z @ b, 50)))[:, None] * Z) + np.diag(lam)
    for _ in range(iters):
        mu = np.exp(np.minimum(off + Z @ b, 50))
        g = Z.T @ (w * (y - mu)) - lam * (b - m)
        H = Z.T @ ((w * mu)[:, None] * Z) + np.diag(lam)
        step = np.linalg.solve(H, g)
        t = 1.0
        while True:
            lp2 = logpost(b + t * step)
            if lp2 >= lp - 1e-10 or t < 1e-8:
                break
            t /= 2
        b, lp = b + t * step, lp2
        if np.abs(t * step).max() < 1e-10:
            break
    mu = np.exp(np.minimum(off + Z @ b, 50))
    H = Z.T @ ((w * mu)[:, None] * Z) + np.diag(lam)
    return b, np.linalg.inv(H)


def dispersion(e, ee, key):
    """quasi-Poisson dispersion of a count with expectation e summed over P
    whose squared per-P expectations sum to ee: per-P gamma factors (between-P
    variance 1/A) and, for squares, PHI_SUM within a P:
    phi = phi_within + (1/A) ee / e"""
    np = _np()
    a, phi0 = (A_SQ, PHI_SUM) if key == "sq" else (A_TRAV, 1.0)
    e = np.asarray(e, float)
    return phi0 + np.asarray(ee, float) / a / np.maximum(e, 1e-12)


class Calibration:
    """learned class factors: per cell, g_sq(c) on squares and m(c) on
    P(magic | square), from three main-effects Poisson GLMs over the cells"""

    def __init__(self, beta=None, cov=None, data=None):
        np = _np()
        self.beta = {k: np.array(beta[k] if beta else GLM_PRIORS[k][0], float) for k in GLM_PRIORS}
        self.cov = {k: (np.array(cov[k], float) if cov else np.diag(np.array(GLM_PRIORS[k][1]) ** 2))
                    for k in GLM_PRIORS}
        self.data = data or {}

    @staticmethod
    def fit(table, ee=None):
        """table: [NCELLS, 9] = sums, O_sq, E_sq, O_S, E_S, O_P, E_P, O_SP, E_SP;
        ee: [NCELLS, 3] sums over P of the squared per-P expectations E_sq,
        E_S, E_P (Summary.cell_ee), for the quasi-Poisson dispersion (None:
        Poisson)"""
        np = _np()
        Z = cell_design()
        beta, cov, data = {}, {}, {}
        for j, (key, (o, e)) in enumerate((("sq", (1, 2)), ("S", (3, 4)), ("P", (5, 6)))):
            use = table[:, e] > 0
            phi = None if ee is None else dispersion(table[use, e], ee[use, j], key)
            b, C = fit_poisson_map(Z[use], table[use, o], np.log(table[use, e]), *GLM_PRIORS[key],
                                   phi=phi)
            beta[key], cov[key] = b, C
            data[key] = [float(table[use, o].sum()), float(table[use, e].sum())]
        return Calibration(beta, cov, data)

    def rates(self, beta=None):
        """per cell: ln g_sq, ln r_S, ln r_P, ln m"""
        am = _am()
        Z = cell_design()
        beta = beta or self.beta
        lg, lS, lP = Z @ beta["sq"], Z @ beta["S"], Z @ beta["P"]
        lm = math.log(am.PAIR) + 2 * (lS + lP - math.log(TRAV_BASE[0]) - math.log(TRAV_BASE[1]))
        return lg, lS, lP, lm

    def draws(self, rng, k):
        """k draws of (ln g_sq, ln m) per cell from the Laplace posterior"""
        out = []
        for _ in range(k):
            b = {key: rng.multivariate_normal(self.beta[key], self.cov[key]) for key in self.beta}
            lg, _, _, lm = self.rates(b)
            out.append((lg, lm))
        return out

    def to_json(self):
        return {"version": _am().AMODEL_VERSION, "effects": list(EFFECTS),
                "beta": {k: [float(v) for v in b] for k, b in self.beta.items()},
                "cov": {k: c.tolist() for k, c in self.cov.items()}, "data": self.data}

    @staticmethod
    def from_json(d):
        return Calibration(d["beta"], d["cov"], d.get("data"))


def fit_time_models(stats, prior=None):
    """TimeModel per msearch engine and mode from sufficient statistics
    {"engine:mode": stats}: engines in increasing order, each starting from
    the posterior of the previous engine (engine 2, the shipped prior's,
    from the shipped prior); returns {key: TimeModel}"""
    am = _am()
    np = _np()
    out = {}
    last = {}
    for key in sorted(stats, key=lambda s: (int(s.split(":")[0]), s)):
        engine, mode = key.split(":")
        engine = int(engine)
        st = {k: (np.array(v) if isinstance(v, list) else v) for k, v in stats[key].items()}
        # (engine 3 on: the previous engine's posterior, shifted by
        # amodel.ENGINE_TIME_SHIFT; the d-first law starts at engine 3 from
        # amodel.DFIRST_TIME_PRIOR)
        first = _prior_engine(mode)
        base = am.time_prior(engine, last.get(mode), mode) if engine >= first else None
        tm = am.TimeModel(engine=engine, mode=mode).fit(st, prior=base)
        tm.engine = engine
        out[key] = tm
        if engine >= first:
            last[mode] = tm
    return out


def current_ratio_stats(ratio, engine=None, with_engine=False):
    """the d-first / plain ratio pairs [n, sum x, sum y, sum y^2] of the
    newest engine that has any (at least ENGINE; an older engine's pairs
    stand in until the current one has some, like the time laws' priors);
    with_engine: (pairs, their engine)"""
    out = (None, None)
    if ratio:
        engines = sorted(int(e) for e in ratio)
        engine = max([ENGINE] + engines) if engine is None else engine
        if str(engine) in ratio:
            out = (ratio[str(engine)], engine)
        else:
            older = [e for e in engines if e < engine]
            if older:
                out = (ratio[str(max(older))], max(older))
    return out if with_engine else out[0]


def current_ratio_level(ratio, engine=None):
    """(the ratio's level a0 for engine, default the newest at least ENGINE,
    and the pairs of the newest engine that has any): the engines' pairs
    chained in increasing order, as fit_time_models chains the time laws:
    each engine's level is the posterior of its own pairs under a prior at
    the previous engine's level shifted to it (amodel.dfirst_ratio_coefs;
    the oldest from its shipped prior), and the newest is shifted to
    engine (amodel.dfirst_ratio_level)"""
    am = _am()
    engs = sorted(int(e) for e in (ratio or {}))
    engine = max([ENGINE] + engs) if engine is None else engine
    level, prev = None, None
    for e in engs:
        if e > engine:
            break
        prior = None if prev is None else level + am.dfirst_ratio_coefs(e)[0] - am.dfirst_ratio_coefs(prev)[0]
        level, prev = am.dfirst_ratio_level(ratio[str(e)], e, e, prior=prior), e
    st = current_ratio_stats(ratio, engine)
    if prev is None:
        return am.dfirst_ratio_level(None, engine), st
    return level + am.dfirst_ratio_coefs(engine)[0] - am.dfirst_ratio_coefs(prev)[0], st


def _prior_engine(mode):
    am = _am()
    return (am.DFIRST_TIME_PRIOR if mode == "dfirst" else am.TIME_PRIOR)["engine"]


def current_time_model(models, engine=None, mode="plain"):
    """the time model of the newest engine (at least ENGINE) for a search
    mode ("plain", "dfirst"): the fit of that engine, or the posterior of
    the newest older engine with a prior (>= 2 plain, >= 3 d-first), or the
    shipped prior (either shifted by amodel.ENGINE_TIME_SHIFT, see
    amodel.time_prior)"""
    am = _am()
    engines = [int(k.split(":")[0]) for k in models if k.endswith(":" + mode)]
    engine = max([ENGINE] + engines) if engine is None else engine
    key = f"{engine}:{mode}"
    if key in models:
        return models[key]
    older = [e for e in engines if _prior_engine(mode) <= e < engine]
    if older:
        return am.time_prior(engine, models[f"{max(older)}:{mode}"], mode)
    return am.time_prior(engine, mode=mode)


# --------------------------------------------------------------------------
# the pool of candidate P (pool.py, or gen_candidates for --pool classic)


class Pool:
    def __init__(self, arrays, kind, params):
        self.exps = arrays["exps"]
        self.S0 = arrays["S0"]
        self.ratio = arrays["ratio"]
        self.tau = arrays["tau"]
        self.k = arrays["k"]
        self.kind = kind
        self.params = params
        self._index = None

    def __len__(self):
        return len(self.exps)

    def key(self):
        import hashlib
        h = hashlib.sha1(self.exps.tobytes())
        h.update(json.dumps([self.kind, self.params], sort_keys=True).encode())
        return h.hexdigest()[:16]

    def index(self):
        """exponent string 'e_e_e' (as p_str(P, '_')) -> row"""
        if self._index is None:
            self._index = {p_str(norm_p(int(v) for v in row), "_"): i
                           for i, row in enumerate(self.exps.tolist())}
        return self._index

    @staticmethod
    def arrays_for(Ps, width=len(PRIMES)):
        """pool arrays for an explicit list of P"""
        np = _np()
        am = _am()
        Ps = [tuple(P) for P in Ps]
        exps = np.zeros((len(Ps), width), np.uint8)
        for i, P in enumerate(Ps):
            exps[i, :len(P)] = P
        return {"exps": exps, "S0": np.array([am.s0(P) for P in Ps], np.float32),
                "ratio": np.array([am.ratio(P) for P in Ps], np.float32),
                "tau": np.array([tau(P) for P in Ps], np.int32),
                "k": np.array([sum(1 for a in P if a) for P in Ps], np.int8)}

    @staticmethod
    def load(args, state):
        """the configured pool (wide: cached in state/pool_6.npz)"""
        np = _np()
        if args.pool == "classic":
            cands = gen_candidates(args.vec_size, args.tau_min, args.tau_max)
            params = {"tau_min": args.tau_min, "tau_max": args.tau_max}
            return Pool(Pool.arrays_for(cands), "classic", params)
        if args.vec_size != 6:
            sys.exit("--pool wide is for 6x6 only (pool.py uses S0 = 6 P^(1/6)): use --pool classic")
        import pool as poolmod
        params = dict(poolmod.DEFAULTS, npr=args.pool_primes, s0_max=args.pool_s0_max,
                      tau_min=args.tau_min, ratio_max=args.pool_ratio)
        if params["npr"] > len(PRIMES):
            sys.exit(f"--pool-primes {params['npr']}: msearch accepts {len(PRIMES)} primes")
        path = os.path.join(state, f"pool_{args.vec_size}.npz")
        meta = json.dumps({"params": params, "generator": poolmod.POOL_VERSION}, sort_keys=True)
        if os.path.exists(path):
            with np.load(path) as z:
                if str(z["meta"]) == meta:
                    return Pool({k: z[k] for k in ("exps", "S0", "ratio", "tau", "k")}, "wide",
                                params)
        t0 = time.time()
        arrays = poolmod.gen_pool(**params)
        if arrays["exps"].shape[1] < len(PRIMES):
            pad = len(PRIMES) - arrays["exps"].shape[1]
            arrays["exps"] = np.pad(arrays["exps"], ((0, 0), (0, pad)))
        log(f"generated the pool: {len(arrays['exps'])} values of P ({time.time() - t0:.1f} s)")
        tmp = path + ".tmp.npz"
        np.savez(tmp, meta=np.array(meta), **arrays)
        os.replace(tmp, path)
        return Pool(arrays, "wide", params)


# --------------------------------------------------------------------------
# profiles of the model on the grid, per P (state/profiles_6/)


class ProfileStore:
    """the analytic model of every pool P on the grid S0 (1 + GRID_U),
    computed once (ProcessPoolExecutor, resumable) and kept as a float16
    memmap [nP, len(FIELDS), len(GRID_U)], with valid flags per point.
    P outside the pool (--only, or P with results) are kept in extra.json."""

    CHUNK = 2000

    def __init__(self, state, pool, n):
        np = _np()
        am = _am()
        self.n = n
        self.pool = pool
        self.dir = os.path.join(state, f"profiles_{n}")
        os.makedirs(self.dir, exist_ok=True)
        self.nf, self.ng = len(am.FIELDS), len(am.GRID_U)
        self.extra_path = os.path.join(self.dir, "extra.json")
        tag = {"version": am.PROFILE_VERSION, "hash": am.profile_hash(),
               "grid": [float(u) for u in am.GRID_U], "fields": list(am.FIELDS)}
        self.extra = {}
        self.extra_dirty = False
        if os.path.exists(self.extra_path):
            with open(self.extra_path) as f:
                d = json.load(f)
            if d.get("tag") == tag:
                self.extra = d["P"]
        self.tag = tag
        self.guards = {g: 0 for g in am.GUARDS}
        self.arr = self.valid = self.done = self.pdir = None
        if pool is None:
            return
        meta = dict(tag, pool=pool.key(), nP=len(pool), params=pool.params, kind=pool.kind)
        # one subdirectory per pool, so that switching pools (e.g. --pool-ratio)
        # never discards a filled store
        self.pdir = os.path.join(self.dir, f"{pool.kind}-{pool.key()}")
        os.makedirs(self.pdir, exist_ok=True)
        mpath = os.path.join(self.pdir, "meta.json")
        old = None
        if os.path.exists(mpath):
            with open(mpath) as f:
                old = json.load(f)
        files = [os.path.join(self.pdir, x) for x in ("arrays.f16", "valid.u1", "done.u1")]
        fresh = (old is None or {k: old.get(k) for k in meta} != meta
                 or not all(os.path.exists(x) for x in files))
        shape = (len(pool), self.nf, self.ng)
        mode = "w+" if fresh else "r+"
        if fresh:
            for x in files:
                if os.path.exists(x):
                    os.remove(x)
        if len(pool):
            self.arr = np.memmap(files[0], np.float16, mode, shape=shape)
            self.valid = np.memmap(files[1], np.uint8, mode, shape=(len(pool), self.ng))
            self.done = np.memmap(files[2], np.uint8, mode, shape=(len(pool),))
        else:
            self.arr = np.zeros(shape, np.float16)
            self.valid = np.zeros((0, self.ng), np.uint8)
            self.done = np.zeros(0, np.uint8)
        if fresh:
            meta["guards"] = self.guards
            self._write_meta(meta)
        else:
            self.guards = old.get("guards", self.guards)
        self.meta = meta

    def _write_meta(self, meta):
        path = os.path.join(self.pdir, "meta.json")
        with open(path + ".tmp", "w") as f:
            json.dump(meta, f)
        os.replace(path + ".tmp", path)

    def save_extra(self):
        if not self.extra_dirty:
            return
        with open(self.extra_path + ".tmp", "w") as f:
            json.dump({"tag": self.tag, "P": self.extra}, f)
        os.replace(self.extra_path + ".tmp", self.extra_path)
        self.extra_dirty = False

    def row_of(self, key):
        if self.pool is None:
            return None
        return self.pool.index().get(key)

    def fill(self, rows, workers=2, quiet=False):
        """compute the missing profiles among pool rows"""
        np = _np()
        am = _am()
        rows = np.asarray(sorted(set(int(r) for r in rows)), np.int64)
        if len(rows) == 0:
            return
        rows = rows[self.done[rows] == 0]
        if len(rows) == 0:
            return
        chunks = [rows[i:i + self.CHUNK] for i in range(0, len(rows), self.CHUNK)]
        t0 = time.time()
        if not quiet:
            log(f"profiling {len(rows)} values of P with the analytic model "
                f"({workers} workers, niced)...")

        def Ps(ch):
            return [norm_p(int(v) for v in self.pool.exps[r]) for r in ch]

        def store(ch, res):
            arr, val, gc = res
            self.arr[ch] = arr
            self.valid[ch] = val
            self.arr.flush()
            self.valid.flush()
            self.done[ch] = 1
            self.done.flush()
            for j, g in enumerate(am.GUARDS):
                self.guards[g] = self.guards.get(g, 0) + int(gc[:, j].sum())
            self.meta["guards"] = self.guards
            self._write_meta(self.meta)

        if workers <= 1 or len(chunks) == 1:
            am.worker_init() if workers <= 1 and len(rows) > 50 else None
            results = (am.profile_rows(Ps(ch), self.n) for ch in chunks)
            for i, (ch, res) in enumerate(zip(chunks, results)):
                store(ch, res)
        else:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor
            ctx = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(workers, mp_context=ctx, initializer=am.worker_init) as ex:
                futs = [ex.submit(am.profile_rows, Ps(ch), self.n) for ch in chunks]
                for i, (ch, fut) in enumerate(zip(chunks, futs)):
                    store(ch, fut.result())
                    if not quiet and (i + 1) % 5 == 0:
                        done = (i + 1) * self.CHUNK
                        el = time.time() - t0
                        log(f"  {min(done, len(rows))}/{len(rows)} P, {el:.0f} s, "
                            f"~{el / done * (len(rows) - done) / 60:.0f} min left")
        self.meta["guards"] = self.guards
        self._write_meta(self.meta)
        if not quiet:
            log(f"profiled {len(rows)} P in {time.time() - t0:.0f} s")

    def fill_extra(self, Ps):
        am = _am()
        todo = [P for P in Ps if p_str(P, "_") not in self.extra]
        if not todo:
            return
        arr, val, gc = am.profile_rows(todo, self.n)
        for P, a, v in zip(todo, arr, val):
            self.extra[p_str(P, "_")] = {"a": a.astype(float).tolist(), "v": v.tolist()}
        self.extra_dirty = True

    def get(self, P, compute=True):
        """(arr float32 [fields, grid], valid bool [grid]) of P, or None"""
        np = _np()
        key = p_str(P, "_")
        r = self.row_of(key)
        if r is not None:
            if not self.done[r]:
                if not compute:
                    return None
                self.fill([r], workers=1, quiet=True)
            return np.asarray(self.arr[r], np.float32), np.asarray(self.valid[r], bool)
        if key not in self.extra:
            if not compute:
                return None
            self.fill_extra([P])
        e = self.extra[key]
        return np.array(e["a"], np.float32), np.array(e["v"], bool)


# --------------------------------------------------------------------------
# summary of all results, updated incrementally (state/summary_6.*)


def _merge(intervals):
    out = []
    for lo, hi in sorted(intervals):
        if out and lo <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def _merge_half(intervals):
    """merge half-open integer intervals [lo, hi) (overlapping or touching)"""
    out = []
    for lo, hi in sorted(intervals):
        if hi <= lo:
            continue
        if out and lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def _covers(iv, nd):
    """the merged half-open intervals iv cover [0, nd)"""
    return nd is not None and bool(iv) and iv[0][0] <= 0 and iv[0][1] >= nd


def _in_cover(cov, S):
    return any(lo <= S <= hi for lo, hi in cov)


def _uncovered_len(iv, lo, hi):
    """the length of [lo, hi) outside the merged half-open intervals iv"""
    n = max(hi - lo, 0)
    for a, b in iv:
        n -= max(0, min(b, hi) - max(a, lo))
    return max(n, 0)


# the star cover of d-first sums (msearch --dfirst-star K, research/ideas.md
# "The star cover of the d loop"): v2 passes K = DFIRST_STAR_K to every
# d-first unit of an even n (0, no star cover, for odd n), and the K and x*
# of a sum's parts to the units that continue it (--dfirst-star-x)
DFIRST_STAR_K = 4


def dfirst_star_k(n):
    """the star cover's K of v2's d-first units for vec size n"""
    return DFIRST_STAR_K if n % 2 == 0 else 0


def dchunk_star(r):
    """the star cover of a "dchunk" or "dsum" record as (x*, K): (None, 0)
    without one (K 0, or a record of engine 3); None for records that
    cannot count towards a sum's coverage: --dfirst-star-only (only the
    star d), or a star cover without star_x (the c2/star prototype's
    chunks)"""
    k = int(r.get("star_k", 0) or 0)
    if r.get("star_only"):
        return None
    if k == 0:
        return (None, 0)
    if r.get("star_x") is None:
        return None
    return (int(r["star_x"]), k)


def dchunk_est_pairs(r, star):
    """the estimated (square, SP traversal) pairs of every d of a d_stride-1
    "dchunk" record's range, under its star cover star = (x*, K): the d
    without x* as found plus K x the pairs of the star d it searched (the
    star d of rank % K == 0, ranked over the whole sum, so the chunks of a
    sum add up to msearch's est_pairs of the whole sum); None with K = -1
    (no star d searched, no estimate)"""
    k = star[1]
    p = float(r.get("pairs", 0))
    if k == 0:
        return p
    if k < 0:
        return None
    ps = float(r.get("pairs_star", 0))
    return p - ps + k * ps


def dfirst_law_star(r, n):
    """does a "dsum" record teach the d-first time law (learnt as run)? Not
    with --dfirst-star-only, and from engine 4 on only with v2's own K
    (dfirst_star_k); engine 3 had no star cover"""
    if r.get("star_only"):
        return False
    k = int(r.get("star_k", 0) or 0)
    return k == (dfirst_star_k(n) if engine_of(r) >= 4 else 0)


def _dgroup_key(star):
    return "-" if star[1] == 0 else f"{star[0]}:{star[1]}"


def dsum_span(r):
    """[d_lo, end) of the d a d_stride-1 "dsum" record went through: its nd
    searched plus the d the star filter skipped (msearch stops between
    chunks under --time-limit, so d_hi can be beyond)"""
    end = r["d_lo"] + r.get("nd", 0) + r.get("nd_star_skipped", 0) + r.get("nd_other_skipped", 0)
    return r["d_lo"], min(end, r.get("nvecs_raw", end))


def open_text(path):
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, "rb")
    return open(path, "rb")


def unit_files(units_dir):
    """unit files (*.jsonl, *.jsonl.gz), a compressed file standing for the
    plain one of the same name"""
    plain = glob.glob(os.path.join(units_dir, "*.jsonl"))
    gz = [p for p in glob.glob(os.path.join(units_dir, "*.jsonl.gz")) if p[:-3] not in plain]
    return sorted(plain + gz)


# keys that mark a sampled msearch record (r1-sampling research builds,
# msearch's "csum" records and the "done" record of an --r1-* run)
SAMPLED_KEYS = ("sample", "stride", "r1_stride", "r1_sample")
# the format of the summary (bump it when update_file / _ingest read the
# records differently, so that an old summary is rebuilt): 2 = d-first
# sums (msearch --diag-first); 3 = distinct d-first sums, notable squares
# deduplicated, the CPU of sampled records; 4 = the parts of d-first sums
# merged (dcov), the d-first time law, the calibration streams in the cells;
# 5 = parts weighted in the d-first law, the streams in the plain law and
# the d-first / plain ratio pairs, the CPU of killed d-first units, SP-type
# squares found d-first; 6 = the star cover: parts merged per (x*, K), the
# d-first pairs and their estimates, the span of a part with skipped d;
# 7 = a d-first / plain ratio pair from the parts of one engine only
SUMMARY_VERSION = 7


def dsum_est_pairs(r):
    """the estimated (square, SP traversal) pairs over all the d of a "dsum"
    record's range (est_pairs; a record without it: pairs), or None where
    the record estimates only some of them: with the star cover (msearch
    --dfirst-star K, star_k = K) the star d are searched every K-th and
    their pairs weighed K x in est_pairs, but with K = -1 none is searched
    (est_pairs covers the other d only) and with --dfirst-star-only only
    they are"""
    if r.get("star_only") or r.get("star_k", 0) < 0:
        return None
    return float(r.get("est_pairs", r.get("pairs", 0)))


# the CPU and time fields of msearch records that a machine record scales
# to reference CPU (scheduler.py run --machine)
MACHINE_TIME_KEYS = ("cpu", "time", "setup_time", "enum_time", "index_time", "est_time",
                     "se_time", "reduce_time", "vd_time", "search_time")
DFIRST_TYPES = ("dsum", "dchunk", "dsquare")


def machine_record(P, S, mid, speed):
    """the record run --machine writes at the top of a unit file: the
    machine's per-process speeds (reference CPU seconds per process CPU
    second, scripts/machine_cal.py's per_process), by which the summary
    scales the CPU of the file's records, so that the time laws learn
    reference CPU from every machine"""
    return {"type": "machine", "n": None, "P": list(P), "S": S, "machine": mid,
            "speed": {"plain": float(speed["plain"]), "dfirst": float(speed["dfirst"])}}


def scale_record(r, speed):
    """a msearch record with its CPU and times in reference CPU seconds
    (speed: a machine record's per-process speeds; d-first records at the
    d-first speed, the others at the plain speed)"""
    if not speed or r.get("type") in ("square", "csquare", "done", "skip"):
        return r
    f = float(speed["dfirst" if r.get("type") in DFIRST_TYPES else "plain"])
    if f == 1.0:
        return r
    r = dict(r)
    for k in MACHINE_TIME_KEYS:
        if isinstance(r.get(k), (int, float)):
            r[k] = r[k] * f
    return r


def sum_cpu(r):
    """CPU seconds of a msearch sum record: the process CPU time where
    msearch reports it ("cpu"), else its wall times"""
    if "cpu" in r:
        return r["cpu"]
    return r["time"] + r["setup_time"] + r.get("enum_time", 0.0)


class Summary:
    """everything the scheduler needs from the msearch output files, without
    keeping the records: per file the bytes parsed; per touched P merged
    coverage, observed counts and per-cell expectations; the cell table;
    time and bias statistics; notable squares in full. Rebuilt from units/
    when missing or when the model version changes."""

    NOTABLE = 7
    FIELDS9 = ("sums", "O_sq", "E_sq", "O_S", "E_S", "O_P", "E_P", "O_SP", "E_SP")

    def __init__(self, state, n):
        self.state = state
        self.n = n
        self.path = os.path.join(state, f"summary_{n}")
        self.reset()

    def reset(self):
        am = _am()
        self.tag = {"version": am.AMODEL_VERSION, "summary": SUMMARY_VERSION,
                    "hash": am.model_hash(),
                    "grid": [float(u) for u in am.GRID_U],
                    "cells": [list(NBAND_EDGES), list(XBAND_EDGES), NCELLS]}
        self.files = {}      # basename (without .gz) -> per-file stats
        self.cover = {}      # 'e_e_e' -> merged [[lo, hi], ...]
        # the d-first sums searched in part: 'e_e_e' -> {"S": {"nd": number
        # of d or None, "g": {star cover "x:K" or "-": {"x", "k", "iv":
        # merged [[lo, hi), ...] of d searched, "pairs", "est"}}}}, from the
        # "dchunk" records (d_stride 1; _dcover); a sum leaves it for cover
        # once one group covers [0, nd)
        self.dcov = {}
        # 'e_e_e' -> {"o": [sq, S, P, SP], "cpu", "nsums", "cells": {c: [9]},
        #  "pcov": merged [[lo, hi]] of the plain sums with a "sum" record,
        #  "calS": {"S": the calibration stream's contribution [cell, 8
        #  values] to the cells and its [squares, est_squares, model
        #  squares], or None}, "dfirst_S": [complete d-first S],
        #  "dT": {"S": [weight, weighted ln whole-sum d-first CPU, ln(N'/4000),
        #  engine] of the parts so far}, "calT": {"S": ln of the stream's plain CPU
        #  estimate}, "rpair": [S paired in self.ratio]}
        self.perP = {}
        self.time = {}       # 'engine:mode' -> {"n", "FF", "Fy", "yy"} (weighted)
        # 'engine' -> [n, sum x, sum y, sum y^2]: per d-first sum searched in
        # full with a calibration stream, y = ln(its d-first CPU / the
        # stream's plain estimate), x = ln(N'/4000) (amodel.dfirst_ratio_level)
        self.ratio = {}
        self.tband = {}      # 'engine:mode:band:W' -> [n, sum of the NT features, sum y]
        self.nbias = {}      # band -> [n, sum, sumsq] of ln(nvecs_raw / N')
        self.lbias = {}      # band -> [n, sum, sumsq] of ln(labels / labels_obs)
        self.notable = []
        self.types = {}
        self.totals = {"sums": 0, "squares": 0, "cpu": 0.0, "trav_S": 0, "trav_P": 0, "trav_SP": 0,
                       "outside": 0, "truncated": 0, "other_n": 0,
                       # d-first sums (msearch --diag-first): searched in
                       # full, parts of one (a --d-range unit, a --d-stride
                       # sample, a truncated run), their CPU seconds, and
                       # the distinct magic squares found
                       "dfirst_sums": 0, "dfirst_partial": 0, "dfirst_cpu": 0.0,
                       "dfirst_magic": 0, "dfirst_truncated": 0,
                       # of the d-first sums searched in full: the (square,
                       # SP traversal) pairs found, their estimate over
                       # every d (the star cover searches every K-th d
                       # through x*) and the sums it covers, the sums with
                       # the star cover; chunks left out of the coverage
                       # (--dfirst-star-only, a star cover without star_x)
                       "dfirst_pairs": 0.0, "dfirst_est_pairs": 0.0, "dfirst_est_sums": 0,
                       "dfirst_star_sums": 0, "dfirst_star_skipped": 0,
                       # CPU seconds of the sampled records ("csum" of an
                       # --r1-* run), in no fit
                       "sampled_cpu": 0.0,
                       # calibration streams (--calib-r1-stride, "csum"
                       # mode calib): sums in the cells, their sampled
                       # squares, CPU, est_squares and the model's squares
                       # (x SQ12, before class factors) there, streams of a
                       # sum already in the cells (not counted again)
                       "calib_sums": 0, "calib_squares": 0, "calib_cpu": 0.0,
                       "calib_est": 0.0, "calib_pred": 0.0, "calib_dup": 0}
        self.dirty = False

    @staticmethod
    def load(state, n):
        s = Summary(state, n)
        js = s.path + ".json"
        if os.path.exists(js):
            with open(js) as f:
                d = json.load(f)
            if d.get("tag") == s.tag:
                for k in ("files", "cover", "dcov", "time", "ratio", "tband", "nbias", "lbias",
                          "notable", "types", "totals"):
                    setattr(s, k, d[k])
                s.perP = d["perP"]
                s.types = {int(k): v for k, v in s.types.items()}
                for st in s.perP.values():
                    st["cells"] = {int(c): v for c, v in st["cells"].items()}
        return s

    def save(self):
        d = {"tag": self.tag, "files": self.files, "cover": self.cover, "dcov": self.dcov,
             "time": self.time, "ratio": self.ratio,
             "tband": self.tband, "nbias": self.nbias, "lbias": self.lbias,
             "notable": self.notable, "types": self.types, "totals": self.totals,
             "perP": self.perP}
        np = _np()
        tmp = self.path + ".json.tmp"
        with open(tmp, "w") as f:
            json.dump(d, f, default=lambda x: x.tolist() if hasattr(x, "tolist") else float(x))
        os.replace(tmp, self.path + ".json")
        tmp = self.path + ".tmp.npz"
        np.savez(tmp, cells=self.cell_table(), tag=np.array(json.dumps(self.tag)))
        os.replace(tmp, self.path + ".npz")
        self.dirty = False

    def cell_table(self):
        np = _np()
        t = np.zeros((NCELLS, 9))
        for st in self.perP.values():
            for c, v in st["cells"].items():
                t[c] += v
        return t

    def cell_ee(self, groups=None):
        """[NCELLS, 3]: per cell, the sum over P of the squared per-P
        expectations E_sq, E_S, E_P (the dispersion of the cell's counts);
        with groups (cell -> group index), the same per group, a P's cells
        in one group summed first"""
        np = _np()
        ng = NCELLS if groups is None else int(max(groups)) + 1
        out = np.zeros((ng, 3))
        for st in self.perP.values():
            acc = {}
            for c, v in st["cells"].items():
                g = c if groups is None else int(groups[c])
                a = acc.setdefault(g, [0.0, 0.0, 0.0])
                a[0] += v[2]
                a[1] += v[4]
                a[2] += v[6]
            for g, a in acc.items():
                out[g] += np.square(a)
        return out

    def update(self, units_dir, store, extra_files=()):
        """parse the new complete lines of every unit file; returns the set
        of P keys that got new records"""
        touched = set()
        for path in list(unit_files(units_dir)) + list(extra_files):
            touched |= self.update_file(path, store)
        store.save_extra()
        return touched

    def update_file(self, path, store):
        name = os.path.basename(path)
        if name.endswith(".gz"):
            name = name[:-3]
        st = self.files.get(name)
        if st is None:
            st = {"off": 0, "sums": 0, "squares": 0, "cpu": 0.0, "done": 0, "complete": 0,
                  "P": None, "last_sum": None}
        if st.get("closed"):
            return set()
        try:
            f = open_text(path)
        except OSError:
            return set()
        with f:
            if path.endswith(".gz"):
                data = f.read()
                if len(data) <= st["off"]:
                    st["closed"] = 1
                    self.files[name] = st
                    return set()
                data = data[st["off"]:]
            else:
                f.seek(st["off"])
                data = f.read()
        end = data.rfind(b"\n") + 1
        if end <= 0:
            return set()
        recs = {}
        # msearch writes a sum's squares before the sum's record: a square
        # counts only once its sum record has arrived (a unit killed in the
        # middle of a sum leaves orphan squares, which the rerun of that sum
        # writes again); squares still waiting are kept in st["pending"]
        pending = st.pop("pending", [])
        # the same for the "csquare" records of a sampled plain search: they
        # precede their "csum"; those of a calibration stream (csum mode
        # "calib") are ingested with it, the others are dropped
        cpend = st.pop("cpending", [])
        for line in data[:end].split(b"\n"):
            if not line:
                continue
            try:
                r = json.loads(line)
            except ValueError:
                continue
            if not isinstance(r, dict) or "type" not in r:
                continue
            if r["type"] == "machine":
                # (run --machine: the records after it ran on that machine)
                st["speed"] = r.get("speed")
                continue
            if r.get("n") != self.n:
                self.totals["other_n"] += 1
                continue
            if st.get("speed"):
                r = scale_record(r, st["speed"])
            # the sums of this file that its "done" record must not cover
            # (fst["holes"], kept across the incremental reads): d-first
            # sums not searched in full, the plain sums that a --d-range
            # unit left to the unit from d 0 ("skip"), r1-sampled plain sums
            t = r["type"]
            if (t == "dsum" and not r.get("complete")) or t == "skip" or (
                    t == "csum" and r.get("mode") == "sampled"):
                h = st.setdefault("holes", {}).setdefault(p_str(norm_p(r["P"]), "_"), [])
                if r["S"] not in h:
                    h.append(r["S"])
            if t == "csquare":
                cpend.append(r)
                continue
            if t == "csum":
                mine = [q for q in cpend if q["S"] == r["S"] and q["P"] == r["P"]]
                cpend = [q for q in cpend if not (q["S"] == r["S"] and q["P"] == r["P"])]
                if r.get("mode") == "calib":
                    # a d-first sum's calibration stream: into the cells
                    # (see _ingest), not a sampled research record
                    self.totals["calib_cpu"] += r.get("cpu", 0.0)
                    recs.setdefault(norm_p(r["P"]), []).append(dict(r, _csq=mine))
                    continue
            if any(k in r for k in SAMPLED_KEYS):
                # a sampled search (every k-th first row, research builds):
                # its counts and times are not those of the sum
                self.totals["sampled"] = self.totals.get("sampled", 0) + 1
                if t == "csum":
                    self.totals["sampled_cpu"] = self.totals.get("sampled_cpu", 0.0) + r.get("cpu", 0.0)
                continue
            if r["type"] == "square":
                pending.append(r)
                continue
            P = norm_p(r["P"])
            if r["type"] == "sum":
                keep = []
                for q in pending:
                    if q["S"] == r["S"] and norm_p(q["P"]) == P:
                        recs.setdefault(P, []).append(q)
                    else:
                        keep.append(q)
                pending = keep
            elif r["type"] == "done" and pending:
                self.totals["orphans"] = self.totals.get("orphans", 0) + len(pending)
                pending = []
            if r["type"] == "done":
                cpend = []
            recs.setdefault(P, []).append(r)
        if pending:
            st["pending"] = pending
        if cpend:
            st["cpending"] = cpend
        st["off"] += end
        for P, rs in recs.items():
            self._ingest(P, rs, st, store)
        self.files[name] = st
        self.dirty = True
        return {p_str(P, "_") for P in recs}

    def _ingest(self, P, rs, fst, store):
        np = _np()
        am = _am()
        key = p_str(P, "_")
        fst["P"] = fst["P"] or key
        cov = []
        sums = [r for r in rs if r["type"] == "sum"]
        sqs = [r for r in rs if r["type"] == "square"]
        holes = fst.get("holes", {}).get(key, ())
        ps = self.perP.setdefault(key, {"o": [0, 0, 0, 0], "cpu": 0.0, "nsums": 0, "cells": {}})
        for r in rs:
            if r["type"] == "sum":
                cov.append((r["S"], r["S"]))
            elif r["type"] == "done":
                fst["done"] = 1
                fst["complete"] = int(r.get("complete", 0))
                fst["last_sum"] = r["last_sum"]
                if r["last_sum"] >= r["min_sum"]:
                    cov.extend(split_range(r["min_sum"], r["last_sum"], holes))
        # d-first sums searched in full: one complete "dsum", or parts (the
        # "dchunk" records of --d-range units, also of killed ones) that
        # together cover every d
        for S in self._dcover(key, ps, rs):
            cov.append((S, S))
        if cov:
            self.cover[key] = _merge([tuple(x) for x in self.cover.get(key, [])] + cov)
        if sums:
            # the plain sums (a calibration stream of the same sum is no
            # longer needed in the cells: the full search replaces it)
            ps["pcov"] = _merge([tuple(x) for x in ps.get("pcov", [])]
                                + [(r["S"], r["S"]) for r in sums])
            dc = self.dcov.get(key)
            for r in sums:
                self._uncalib(ps, r["S"])
                if dc:
                    dc.pop(str(r["S"]), None)
            if dc is not None and not dc:
                self.dcov.pop(key, None)
        for r in sums:
            t = sum_cpu(r)
            ps["cpu"] += t
            ps["nsums"] += 1
            fst["sums"] += 1
            fst["cpu"] += t
            self.totals["sums"] += 1
            self.totals["cpu"] += t
        notable_seen = None
        for q in sqs:
            ps["o"][0] += 1
            ps["o"][1] += q["s_count"]
            ps["o"][2] += q["p_count"]
            ps["o"][3] += q["sp_count"]
            fst["squares"] += 1
            self.totals["squares"] += 1
            self.totals["trav_S"] += q["s_count"]
            self.totals["trav_P"] += q["p_count"]
            self.totals["trav_SP"] += q["sp_count"]
            self.types[q["best_score"]] = self.types.get(q["best_score"], 0) + 1
            if q["best_score"] >= self.NOTABLE:
                # once per square (a sum searched again, or a magic square
                # found d-first before), by its hash where it has one
                if notable_seen is None:
                    notable_seen = {(tuple(x["P"]), x["S"], x["hash"]) for x in self.notable
                                    if "hash" in x}
                kq = (tuple(P), q["S"], q.get("hash"))
                if q.get("hash") is None or kq not in notable_seen:
                    notable_seen.add(kq)
                    self.notable.append(dict(q, P=list(P)))
        self._ingest_dfirst(P, rs, fst)
        # the parts of d-first sums that teach the d-first time law, and the
        # calibration streams
        dsums = [r for r in rs if r["type"] == "dsum" and r.get("d_stride", 1) == 1
                 and not r.get("truncated") and "cpu" in r
                 and r.get("nvecs_raw", 0) >= DFIRST_TIME_MIN_N
                 and r.get("nd", 0) >= min(DFIRST_TIME_MIN_ND, r["nvecs_raw"])
                 and dfirst_law_star(r, self.n)]
        calibs = [r for r in rs if r["type"] == "csum" and "_csq" in r]
        if not sums and not sqs and not dsums and not calibs:
            self._pair_ratio(ps, rs)
            return
        try:
            prof = store.get(P)
        except Exception:
            prof = None
        if prof is None or len(P) > len(am.PRIMES) or not prof[1].any():
            self.totals["outside"] += len(sums)
            for r in calibs:
                self._calib_seen(ps, r, None)
            return
        A, v = prof
        Sg = am.grid_sums(P, self.n)[v]
        k = sum(1 for a in P if a)
        rb = int(ratio_bin(np.float32(am.ratio(P, self.n))))
        smin = am.smin_approx(P, self.n)

        def at(S):
            S = np.asarray(S, float)
            out = {f: np.interp(S, Sg, A[j, v]) for j, f in enumerate(am.FIELDS)}
            out["lNp"] = out["lN"] + am.n_bias(out["lN"])
            out["cell"] = cell_index(out["lNp"], k, rb, np.log(S / smin))
            out["inside"] = (S >= Sg[0] - 1e-9) & (S <= Sg[-1] + 1e-9)
            return out

        cells = ps["cells"]

        def acc(c, j, val):
            row = cells.get(c)
            if row is None:
                row = cells[c] = [0.0] * 9
            row[j] += val

        good = [r for r in sums if not r["truncated"]]
        self.totals["truncated"] += len(sums) - len(good)
        if good:
            a = at([r["S"] for r in good])
            esq = np.exp(a["lEs"]) * np.where(a["lNp"] >= LN12K, am.SQ12, 1.0)
            F = am.time_features(a["lNp"], a["lLraw"], k)
            lab = am.labels_obs(a["lNp"], a["lLraw"], k)
            for i, r in enumerate(good):
                if not a["inside"][i]:
                    self.totals["outside"] += 1
                    continue
                c = int(a["cell"][i])
                acc(c, 0, 1)
                acc(c, 1, r["squares"])
                acc(c, 2, float(esq[i]))
                N = r["nvecs_raw"]
                if N < 300:
                    continue
                band = NBAND_NAMES[int(np.searchsorted(NBAND_EDGES, N, side="right"))]
                for d, val in ((self.nbias, math.log(N) - a["lNp"][i]),
                               (self.lbias, math.log(max(r["labels"], 1)) - math.log(lab[i]))):
                    b = d.setdefault(band, [0, 0.0, 0.0])
                    b[0] += 1
                    b[1] += float(val)
                    b[2] += float(val) ** 2
                t = sum_cpu(r)
                if t < 0.01:
                    continue
                self._time_row(f"{engine_of(r)}:{r.get('mode', 'plain')}", F[i], math.log(t),
                               band, r["labels"])
        if dsums:
            # the whole-sum CPU of a part [d_lo, end) of a d-first sum
            # (dsum_span: its d searched or skipped by the star filter):
            # its overhead (reduction, enumeration share, the d index) plus
            # its d loop scaled by the cost profile along d; the row's
            # weight is the part's share of the d loop (a sum split into m
            # parts weighs 1, not m: its parts share its N' error, and small
            # parts are noisy)
            a = at([r["S"] for r in dsums])
            F = am.time_features(a["lNp"], a["lLraw"], k)
            dT = ps.setdefault("dT", {})
            for i, r in enumerate(dsums):
                if not a["inside"][i]:
                    continue
                raw = r["nvecs_raw"]
                frac = am.dfirst_cost_frac(*dsum_span(r), raw)
                if frac <= 0:
                    continue
                loop = max(r["time"] - r.get("index_time", 0.0), 0.0)
                t = max(r["cpu"] - loop, 0.0) + loop / frac
                if t < 0.01:
                    continue
                band = NBAND_NAMES[int(np.searchsorted(NBAND_EDGES, raw, side="right"))]
                self._time_row(f"{engine_of(r)}:dfirst", F[i], math.log(t), band, r["labels"],
                               w=frac)
                e = dT.setdefault(str(r["S"]), [0.0, 0.0, float(a["lNp"][i]) - math.log(4000),
                                                engine_of(r)])
                # (one engine per ratio pair: parts of another engine than
                # the sum's first part, e.g. older units of another coverage
                # group, stay out of it)
                if e[3] != engine_of(r):
                    continue
                e[0] += frac
                e[1] += frac * math.log(t)
        if calibs:
            a = at([r["S"] for r in calibs])
            F_cal = am.time_features(a["lNp"], a["lLraw"], k)
            T = NUM_TRAVERSALS[self.n]
            for i, r in enumerate(calibs):
                if not a["inside"][i] or r.get("truncated"):
                    self._calib_seen(ps, r, None)
                    continue
                if not self._calib_seen(ps, r, True):
                    continue
                kk = max(int(r.get("r1_stride", 1)), 1)
                esq = math.exp(a["lEs"][i]) * (am.SQ12 if a["lNp"][i] >= LN12K else 1.0)
                o = int(r["squares"])
                est, se = float(r.get("est_squares", o * kk)), float(r.get("se_squares", 0.0))
                phi = se * se / (kk * est) if o >= 5 and est > 0 and se > 0 else CALIB_PHI
                w = min(1.0, PHI_SUM / (phi + PHI_SUM / kk))
                c = int(a["cell"][i])
                acc(c, 0, w)
                self._calib_time(ps, r, a, i, k, F_cal)
                vals = [w * o, w * esq / kk, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                for q in r["_csq"]:
                    vals[2] += w * q["s_count"]
                    vals[3] += w * T * math.exp(a["lpS"][i])
                    vals[4] += w * q["p_count"]
                    vals[5] += w * T * math.exp(a["lpP"][i])
                    vals[6] += w * q["sp_count"]
                    vals[7] += w * T * math.exp(a["lpSP"][i])
                for j, val in enumerate(vals):
                    acc(c, 1 + j, val)
                ps["calS"][str(r["S"])] = [c] + vals + [o, est, esq, w]
                self.totals["calib_sums"] += 1
                self.totals["calib_squares"] += o
                self.totals["calib_est"] += est
                self.totals["calib_pred"] += esq
        if sqs:
            a = at([q["S"] for q in sqs])
            T = NUM_TRAVERSALS[self.n]
            for i, q in enumerate(sqs):
                if not a["inside"][i]:
                    continue
                c = int(a["cell"][i])
                acc(c, 3, q["s_count"])
                acc(c, 4, T * math.exp(a["lpS"][i]))
                acc(c, 5, q["p_count"])
                acc(c, 6, T * math.exp(a["lpP"][i]))
                acc(c, 7, q["sp_count"])
                acc(c, 8, T * math.exp(a["lpSP"][i]))
        self._pair_ratio(ps, rs)

    def _time_row(self, tk, f, y, band, labels, w=1.0):
        """one row (weight w) of the time law tk ("engine:mode")"""
        np = _np()
        nt = len(f)
        ts = self.time.setdefault(tk, {"n": 0, "FF": [[0.0] * nt for _ in range(nt)],
                                       "Fy": [0.0] * nt, "yy": 0.0})
        ts["n"] += w
        ts["FF"] = (np.array(ts["FF"]) + w * np.outer(f, f)).tolist()
        ts["Fy"] = (np.array(ts["Fy"]) + w * f * y).tolist()
        ts["yy"] += w * y * y
        W = (max(labels, 1) + 63) // 64
        tb = self.tband.setdefault(f"{tk}:{band}:{W}", [0.0] * (nt + 2))
        tb[0] += w
        for j in range(nt):
            tb[1 + j] += w * float(f[j])
        tb[nt + 1] += w * y

    def _calib_time(self, ps, r, a, i, k, F):
        """a calibration stream as data on the plain search of its sum: its
        unbiased estimate of the plain CPU (est_time plus the reduction and
        the enumeration share, as a "sum" record's cpu) is a row of the
        plain time law, weighted by its precision (the law learns where no
        plain sum runs), and one half of the sum's d-first / plain pair;
        the sum's N and label bias too (once per sum: the stream runs once)"""
        np = _np()
        am = _am()
        est = float(r.get("est_time", 0.0) or 0.0)
        raw = r.get("nvecs_raw", 0)
        if raw >= 300:
            lab = float(am.labels_obs(a["lNp"][i], a["lLraw"][i], k))
            band = NBAND_NAMES[int(np.searchsorted(NBAND_EDGES, raw, side="right"))]
            for d, val in ((self.nbias, math.log(raw) - a["lNp"][i]),
                           (self.lbias, math.log(max(r["labels"], 1)) - math.log(lab))):
                b = d.setdefault(band, [0, 0.0, 0.0])
                b[0] += 1
                b[1] += float(val)
                b[2] += float(val) ** 2
        if est <= 0 or raw < 300:
            return
        t = est + r.get("reduce_time", 0.0) + r.get("enum_time", 0.0)
        rel = float(r.get("se_time", 0.0) or 0.0) / est
        w = 1.0 / (1.0 + (rel / am.TIME_PRIOR["sd"]) ** 2)
        self._time_row(f"{engine_of(r)}:plain", F[i], math.log(t), band, r["labels"], w=w)
        ps.setdefault("calT", {})[str(r["S"])] = math.log(t)

    def _pair_ratio(self, ps, rs):
        """d-first sums searched in full that have both their d-first CPU
        (dT, weight >= 0.5) and the plain estimate of their stream (calT):
        one pair each of the d-first / plain ratio (self.ratio, by engine)"""
        dT, calT = ps.get("dT", {}), ps.get("calT", {})
        if not dT or not calT:
            return
        fin = set(ps.get("dfirst_S", ()))
        done = ps.setdefault("rpair", [])
        for S in list(dT):
            if int(S) not in fin or S not in calT or int(S) in done:
                continue
            w, wy, x, eng = dT[S]
            if w < 0.5:
                continue
            y = wy / w - calT[S]
            st = self.ratio.setdefault(str(eng), [0.0, 0.0, 0.0, 0.0])
            st[0] += 1
            st[1] += x
            st[2] += y
            st[3] += y * y
            done.append(int(S))
            dT.pop(S)
            calT.pop(S)

    def _calib_seen(self, ps, r, use):
        """record a calibration stream of sum r["S"]; with use, whether its
        squares may enter the cells: not when the sum already has a stream
        there or was searched plain (no double counting)"""
        cal = ps.setdefault("calS", {})
        S = str(r["S"])
        if not use:
            cal.setdefault(S, None)
            return False
        if cal.get(S) is not None or _in_cover(ps.get("pcov", ()), r["S"]):
            self.totals["calib_dup"] += 1
            return False
        return True

    def _uncalib(self, ps, S):
        """a plain search of sum S replaces its calibration stream in the
        cells"""
        row = ps.get("calS", {}).get(str(S))
        if not row:
            return
        c, vals = int(row[0]), row[1:9]
        cell = ps["cells"].get(c)
        if cell is not None:
            for j, val in enumerate(vals):
                cell[1 + j] -= val
            cell[0] -= row[12]
        ps["calS"][str(S)] = None
        T = self.totals
        T["calib_sums"] -= 1
        T["calib_squares"] -= row[9]
        T["calib_est"] -= row[10]
        T["calib_pred"] -= row[11]

    def _dcover(self, key, ps, rs):
        """merge the d-first records of one P into dcov; returns the sums
        newly searched in full (a complete "dsum", or parts covering every
        d of the sum, duplicates and overlaps merged; only d_stride 1).

        The parts of a sum merge only within one star cover (x*, K) of
        msearch --dfirst-star: two parts with different x* could each skip
        one diagonal of a magic square, so the chunks of each (x*, K) are a
        coverage of their own ("g": key "x:K", or "-" without a star cover),
        and the sum is covered once one of them covers every d. The planner
        continues the group with the most of the sum's d-loop CPU
        (dfirst_part), passing its x* and K. --dfirst-star-only chunks
        (only the star d) and star chunks without star_x count towards no
        coverage. Each group also adds up the pairs of its chunks and their
        estimate over every d (dchunk_est_pairs), each chunk weighted by
        the share of its range not yet covered in the group (a duplicate
        adds nothing), so a sum's estimate is that of one pass over its d:
        exact for disjoint chunks and duplicates, the planner's units (it
        plans only the gaps of a group); a chunk that partly overlaps the
        group (a hand-made unit) adds its pairs by that share, as if they
        were spread evenly along its d (a dchunk record has no per-d
        pairs), so the sum's pairs and est_pairs are then approximate."""
        dc = self.dcov.get(key, {})
        done = []
        fin = set(ps.get("dfirst_S", ()))
        cov = self.cover.get(key, ())
        T = self.totals
        for r in rs:
            t = r["type"]
            if t not in ("dchunk", "dsum"):
                continue
            S = r["S"]
            # (records of a sum searched in full already, d-first or plain:
            # a duplicate unit, or the dsum of the unit whose chunks
            # completed the sum, read after them; no dcov entry again)
            if S in fin or _in_cover(cov, S) or S in [x[0] for x in done]:
                if t == "dchunk":
                    T["dfirst_truncated"] += int(r.get("truncated", 0))
                continue
            if t == "dsum" and r.get("complete"):
                e = dsum_est_pairs(r)
                done.append((S, float(r.get("pairs", 0)), e, int(r.get("star_k", 0) or 0)))
                continue
            if r.get("d_stride", 1) != 1:
                continue
            e = dc.setdefault(str(S), {"nd": None, "g": {}})
            if r.get("nvecs_raw"):
                e["nd"] = int(r["nvecs_raw"])
            if t != "dchunk":
                continue
            # (a V_d that hit --node-limit: counted as searched, like a
            # truncated plain sum)
            T["dfirst_truncated"] += int(r.get("truncated", 0))
            star = dchunk_star(r)
            if star is None:
                T["dfirst_star_skipped"] += 1
                continue
            g = e["g"].setdefault(_dgroup_key(star), {"x": star[0], "k": star[1], "iv": [],
                                                       "pairs": 0.0, "est": 0.0})
            lo, hi = int(r["d_lo"]), int(r["d_hi"])
            if hi > lo:
                w = _uncovered_len(g["iv"], lo, hi) / (hi - lo)
                g["pairs"] += w * float(r.get("pairs", 0))
                ce = dchunk_est_pairs(r, star)
                g["est"] = None if (ce is None or g["est"] is None) else g["est"] + w * ce
            g["iv"] = _merge_half(g["iv"] + [[lo, hi]])
            if _covers(g["iv"], e["nd"]):
                done.append((S, g["pairs"], g["est"], g["k"]))
        out = []
        for S, pairs, est, k in done:
            dc.pop(str(S), None)
            fin = ps.setdefault("dfirst_S", [])
            if S not in fin:
                fin.append(S)
                T["dfirst_sums"] += 1
                T["dfirst_pairs"] += pairs
                if est is not None:
                    T["dfirst_est_pairs"] += est
                    T["dfirst_est_sums"] += 1
                if k:
                    T["dfirst_star_sums"] += 1
                out.append(S)
        if dc:
            self.dcov[key] = dc
        else:
            self.dcov.pop(key, None)
        return out

    def _dgroup(self, key, S):
        """the coverage group of d-first sum S of P (key) that the planner
        continues: the one with the most of the sum's d-loop CPU searched
        (the first of equals), as (nd, group) or None"""
        e = self.dcov.get(key, {}).get(str(S))
        if not e:
            return None
        nd, best, bf = e["nd"], None, -1.0
        for g in e["g"].values():
            f = (sum(_am().dfirst_cost_frac(lo, hi, nd) for lo, hi in g["iv"]) if nd
                 else float(sum(hi - lo for lo, hi in g["iv"])))
            if f > bf:
                best, bf = g, f
        return nd, best

    def dfirst_part(self, key, S):
        """(nd, merged d intervals, share of the sum's d-loop CPU searched,
        star cover (x* or None, K) or None) of a d-first sum searched in
        part, by its group that the planner continues (_dgroup), or None"""
        e = self._dgroup(key, S)
        if e is None:
            return None
        nd, g = e
        if g is None:
            return nd, [], 0.0, None
        frac = (sum(_am().dfirst_cost_frac(lo, hi, nd) for lo, hi in g["iv"]) if nd else 0.0)
        return nd, g["iv"], frac, (g["x"], g["k"])

    def _ingest_dfirst(self, P, rs, fst=None):
        """d-first records (msearch --diag-first): a "dsum" per sum (or part
        of one) and a "dsquare" per (square, SP diagonal) pair. They have no
        semi-magic squares of the sum (only those with an SP traversal, each
        once per such traversal), so they stay out of the cells (squares and
        traversal fits; the sum's calibration stream, a "csum" of mode
        calib, goes there) and out of the plain time law (the d-first law
        learns from them); their CPU is counted apart, and their notable
        squares (magic, or best_score >= NOTABLE) join the notable squares
        once. Complete sums are counted in _dcover."""
        fst = fst if fst is not None else {}
        # d-first CPU: each chunk's d loop as it completes (also those of
        # killed units, whose chunks count as searched), the rest of a dsum
        # (index, reduction, enumeration share) with it
        dct = fst.setdefault("dct", {})
        for r in rs:
            if r["type"] == "dchunk":
                t = float(r.get("time", 0.0))
                self.totals["dfirst_cpu"] += t
                kS = str(r["S"])
                dct[kS] = dct.get(kS, 0.0) + t
            elif r["type"] == "dsum":
                if not r.get("complete"):
                    self.totals["dfirst_partial"] += 1
                c = r.get("cpu", r.get("time", 0.0))
                self.totals["dfirst_cpu"] += max(c - dct.pop(str(r["S"]), 0.0), 0.0)
                fst["dcpu"] = fst.get("dcpu", 0.0) + c
            elif r["type"] == "csum" and "_csq" in r:
                fst["dcpu"] = fst.get("dcpu", 0.0) + r.get("cpu", 0.0)
        if not dct:
            fst.pop("dct", None)
        # notable squares found d-first: the magic ones (a "dsquare" with
        # magic or partner, twice each) and those with an SP traversal and
        # best_score >= NOTABLE (SP, SP+S, SP+P: the only record of them in
        # a d-first sum, once per SP traversal), once per (P, S, hash)
        dm = [q for q in rs if q["type"] == "dsquare" and (
            q.get("best_score", 0) >= self.NOTABLE or q.get("magic") or q.get("partner"))]
        cq = [q for r in rs if r["type"] == "csum" for q in r.get("_csq", ())
              if q["best_score"] >= self.NOTABLE]
        if dm or cq:
            known = {(tuple(q["P"]), q["S"], q.get("hash")) for q in self.notable}
            for q, flag in [(q, "dfirst") for q in dm] + [(q, "calib") for q in cq]:
                kq = (tuple(P), q["S"], q["hash"])
                if kq not in known:
                    known.add(kq)
                    self.notable.append(dict(q, P=list(P), **{flag: 1}))
                    if flag == "dfirst" and (q.get("magic") or q.get("partner")):
                        self.totals["dfirst_magic"] += 1

    def time_stats(self):
        np = _np()
        return {k: {"n": v["n"], "FF": np.array(v["FF"]), "Fy": np.array(v["Fy"]), "yy": v["yy"]}
                for k, v in self.time.items()}


# --------------------------------------------------------------------------
# candidates, scorer and planner

# a unit of v2: P, sums lo..hi; mode "plain" or "dfirst" (msearch
# --diag-first --diag-first-min-n 0); a d-first unit with dlo set is one sum
# (lo == hi) searched on the d with index in [dlo, dhi) of nd (dhi None: to
# the end of the sum, nd only predicted), the share frac of the sum's
# d-loop CPU and E; calib: the stride of its calibration stream (0: none);
# time: predicted CPU, score: magic / CPU with the sum's calibration stream
# spread over its parts; star_k, star_x: the star cover of a d-first unit
# (msearch --dfirst-star K, and --dfirst-star-x x* when it continues a sum
# searched in part with that x*; None: msearch chooses)
UnitV2 = collections.namedtuple(
    "UnitV2", Unit._fields + ("a", "nodes", "cells", "mvec", "mode", "dlo", "dhi", "nd", "calib",
                              "frac", "star_k", "star_x"),
    defaults=(None, 0.0, None, None, "plain", None, None, 0, 0, 1.0, 0, None))


class Cands:
    """the values of P the planner chooses from: pool rows (minus those the
    legacy search exhausted) plus extras, or exactly the --only list"""

    def __init__(self, exps, src, pool_ratio=None, n=6):
        np = _np()
        am = _am()
        self.n = n
        self.exps = exps                      # uint8 [nA, len(PRIMES)]
        self.src = np.asarray(src, np.int64)  # pool row, or -1
        lp = exps.astype(float) @ np.log(np.array(PRIMES, float))
        self.S0 = n * np.exp(lp / n)
        self.tau = np.prod(exps.astype(np.int64) + 1, axis=1)
        self.k = (exps > 0).sum(1)
        self.smin = self.S0 * (1 + 10.0 / self.tau)
        if pool_ratio is None:
            pool_ratio = np.array([am.ratio(self.P(a), n) for a in range(len(exps))], np.float32)
        self.ratio = np.asarray(pool_ratio, np.float32)
        self.rbin = ratio_bin(self.ratio)
        self.end = np.floor(2 * self.S0).astype(np.int64)
        self._index = None

    def __len__(self):
        return len(self.exps)

    def P(self, a):
        return norm_p(int(v) for v in self.exps[a])

    def key(self, a):
        return p_str(self.P(a), "_")

    def index(self):
        if self._index is None:
            self._index = {p_str(norm_p(int(v) for v in row), "_"): a
                           for a, row in enumerate(self.exps.tolist())}
        return self._index


class AnalyticScorer:
    """expected magic squares per CPU-second at (P, S) from the profiles:

        lEm = lEs + ln g_sq(c) + ln SQ12 [N' >= 12k] + min(lPm + ln m(c), ln 1e-6)
              + ln f_sq(P) + ln F_m(P),        time from TimeModel(lN', L_raw)

    with c the calibration cell of the sum and f_sq, F_m the per-P factors.
    Each sum is searched plain (time law tm) or d-first (tmd, plus its
    calibration stream, about calib_frac of that), by policy: "off" plain,
    "on" d-first at N' >= dmin, "auto" d-first at N' >= dmin where the
    measured d-first / plain ratio r(N') (amodel.DFIRST_RATIO_PRIOR, level
    lr0) with the stream is below 1 (see _choose)."""

    def __init__(self, cands, store, calib, tm, n=6, tmd=None, policy="auto",
                 calib_frac=CALIB_FRAC, dmin=DFIRST_MIN_NP, lr0=None):
        np = _np()
        self.c = cands
        self.store = store
        self.n = n
        self.tm = tm
        self.tmd = tmd if tmd is not None else _am().time_prior(ENGINE, mode="dfirst")
        if policy not in DFIRST_POLICIES:
            raise ValueError(f"d-first policy {policy}")
        self.policy = policy
        self.calib_frac = float(calib_frac)
        self.ldmin = math.log(max(float(dmin), 1.0))
        # the level of ln(d-first / plain CPU) (amodel.dfirst_log_ratio) of
        # msearch engine self.engine (its slope)
        self.engine = ENGINE
        self.lr0 = _am().dfirst_ratio_coefs(self.engine)[0] if lr0 is None else float(lr0)
        # stage 1: (ln LO, ln HI) of the N' band searched plain, or None
        self.stage1 = None
        self.set_calibration(calib)
        self.lnfsq = np.zeros(len(cands))
        self.lnFm = np.zeros(len(cands))
        self._cache = {}

    def set_calibration(self, calib):
        self.calib = calib
        self.ln_gsq, self.ln_rS, self.ln_rP, self.ln_m = calib.rates()

    @staticmethod
    def _fast(tm):
        np = _np()
        th = [float(v) for v in tm.th]
        # ln t = c0 + th1 lN' + th2 lL + th3 max(0, lN' - ln 8000) + th4 [labels > 128]
        #        + band offset, + sd^2/2 for the mean (am.time_features, written out)
        return (th[0] - th[1] * math.log(4000) - th[2] * math.log(150) + 0.5 * tm.sd ** 2,
                th[1], th[2], th[3], th[4], np.array(th[5:]))

    @property
    def tm(self):
        return self._tm

    @tm.setter
    def tm(self, tm):
        self._tm = tm
        am = _am()
        np = _np()
        self._tc = self._fast(tm)
        # labels_obs > 128  <=>  lL + LAB1 lN' > _w3c - LAB2 (k - 6)
        self._w3c = math.log(128) - am.LAB0 + am.LAB1 * math.log(4000)
        self._tbe = np.log(np.array(am.TIME_BAND_EDGES))

    @property
    def tmd(self):
        return self._tmd

    @tmd.setter
    def tmd(self, tm):
        self._tmd = tm
        self._tcd = self._fast(tm)

    def prof(self, a):
        """per cand: S of the valid grid points, lN, lEs, lPm, L_raw there
        (float64), and the cell offset of its k and ratio bin"""
        hit = self._cache.get(a)
        if hit is not None:
            return hit
        np = _np()
        am = _am()
        src = self.c.src[a]
        if src >= 0:
            A = np.asarray(self.store.arr[src], np.float64)
            v = np.asarray(self.store.valid[src], bool)
        else:
            A, v = self.store.get(self.c.P(a))
            A = A.astype(np.float64)
        kb = min(max(int(self.c.k[a]) - 5, 0), 3)
        out = (self.c.S0[a] * (1 + am.GRID_U[v]), A[0, v].copy(), A[1, v].copy(), A[2, v].copy(),
               A[6, v].copy(), kb * 9 + int(self.c.rbin[a]) * 3, math.log(self.c.smin[a]),
               int(self.c.k[a]))
        if len(self._cache) > 100000:
            self._cache.clear()
        self._cache[a] = out
        return out

    def log_ratio(self, lNp):
        """ln(d-first CPU / plain CPU) of a sum at ln N'"""
        return _am().dfirst_log_ratio(lNp, self.lr0, self.engine)

    def _choose(self, lNp, tp, td):
        """per sum: d-first?, and the calibration stream's predicted CPU
        (calib_frac of the d-first CPU). auto: by the ratio law, label-free
        and measured on the same sums in both modes, not by td against tp
        (two separately fitted laws whose errors do not cancel: see
        amodel.DFIRST_RATIO_PRIOR); tp and td only price the sums"""
        np = _np()
        cf = max(self.calib_frac, 0.0)
        tcal = cf * np.asarray(td, float)
        if self.policy == "off":
            dm = np.zeros(np.shape(lNp), bool)
        elif self.policy == "on":
            dm = lNp >= self.ldmin
        else:
            dm = (lNp >= self.ldmin) & ((1 + cf) * np.exp(self.log_ratio(lNp)) < 1)
        if self.stage1 is not None:
            dm = dm & ~((lNp >= self.stage1[0]) & (lNp < self.stage1[1]))
        return dm, tcal

    def grid_density(self, rows):
        """log magic squares per CPU-second on the grid, [len(rows), grid]
        (-inf at invalid points), with each point's search mode"""
        np = _np()
        am = _am()
        rows = np.asarray(rows, np.int64)
        A = np.zeros((len(rows), len(am.FIELDS), len(am.GRID_U)), np.float32)
        V = np.zeros((len(rows), len(am.GRID_U)), bool)
        src = self.c.src[rows]
        inpool = src >= 0
        if inpool.any():
            idx = np.nonzero(inpool)[0]
            order = np.argsort(src[idx])
            A[idx[order]] = self.store.arr[src[idx][order]]
            V[idx[order]] = self.store.valid[src[idx][order]].astype(bool)
        for i in np.nonzero(~inpool)[0]:
            A[i], V[i] = self.store.get(self.c.P(rows[i]))
        lN, lEs, lPm, lL = A[:, 0], A[:, 1], A[:, 2], A[:, 6]
        lNp = lN + am.n_bias(lN)
        S = self.c.S0[rows, None] * (1 + am.GRID_U[None, :])
        x = np.log(S / self.c.smin[rows, None])
        cell = cell_index(lNp, self.c.k[rows, None], self.c.rbin[rows, None], x)
        lem = (lEs + self.ln_gsq[cell] + math.log(am.SQ12) * (lNp >= LN12K)
               + np.minimum(lPm + self.ln_m[cell], am.LOG_PM_CAP)
               + (self.lnfsq[rows] + self.lnFm[rows])[:, None])
        tp = self.tm.time(lNp, lL, self.c.k[rows, None])
        if self.policy != "off":
            td = self.tmd.time(lNp, lL, self.c.k[rows, None])
            dm, tcal = self._choose(lNp, tp, td)
            tp = np.where(dm, td + tcal, tp)
        return np.where(V, lem - np.log(tp), -np.inf)

    _NB = None

    def eval_modes(self, a, S):
        """per sum: squares, magic squares, plain CPU, d-first CPU (whole
        sum, without its calibration stream), the stream's CPU, d-first?,
        lN', L_raw, cell (as grid_density, at any S, by linear
        interpolation of the profile's logs in S; written for speed)"""
        np = _np()
        am = _am()
        Sg, gN, gEs, gPm, gL, base, lsmin, k = self.prof(a)
        S = np.asarray(S, float)
        if len(Sg) == 0:
            z = np.zeros(len(S))
            t = z + am.SUM_OVERHEAD
            return z, z, t, t, z, np.zeros(len(S), bool), z, z, np.zeros(len(S), int)
        if AnalyticScorer._NB is None:
            AnalyticScorer._NB = (np.log(np.array(NBAND_EDGES, float)), np.array(XBAND_EDGES))
        nbe, xbe = AnalyticScorer._NB
        lN = np.interp(S, Sg, gN)
        lNp = lN + am.NBIAS_SLOPE * np.maximum(0.0, lN - am.NBIAS_KNOT)
        lL = np.interp(S, Sg, gL)
        cell = (np.searchsorted(nbe, lNp, side="right") * 36 + base
                + np.searchsorted(xbe, np.log(S) - lsmin, side="right"))
        lsq = np.interp(S, Sg, gEs) + self.ln_gsq[cell] + self.lnfsq[a]
        big = lNp >= LN12K
        if big.any():
            lsq = lsq + math.log(am.SQ12) * big
        lpm = np.minimum(np.interp(S, Sg, gPm) + self.ln_m[cell], am.LOG_PM_CAP) + self.lnFm[a]
        sq = np.exp(lsq)
        lo, hi = Sg[0] - 1e-9, Sg[-1] + 1e-9
        below = S < lo
        outside = below | (S > hi)
        c0, c1, c2, c3, c4, cb = self._tc
        w3 = lL + am.LAB1 * lNp > self._w3c - am.LAB2 * (k - 6)
        band = np.searchsorted(self._tbe, lNp, side="right")
        hinge = np.maximum(0.0, lNp - LN8K)
        tp = np.exp(c0 + c1 * lNp + c2 * lL + c3 * hinge + c4 * w3 + cb[band]) + am.SUM_OVERHEAD
        if self.policy != "off":
            c0, c1, c2, c3, c4, cb = self._tcd
            td = np.exp(c0 + c1 * lNp + c2 * lL + c3 * hinge + c4 * w3 + cb[band]) + am.SUM_OVERHEAD
            dm, tcal = self._choose(lNp, tp, td)
        else:
            td, tcal, dm = tp, np.zeros(len(S)), np.zeros(len(S), bool)
        if outside.any():
            sq[outside] = 0.0
            tp[below] = td[below] = am.SUM_OVERHEAD
            tcal[below] = 0.0
            dm = dm & ~below
        m = sq * np.exp(lpm)
        return sq, m, tp, td, tcal, dm, lNp, lL, cell

    def eval_sums(self, a, S):
        """per sum: squares, magic squares, CPU seconds in the chosen mode
        (d-first: with the calibration stream), lN', L_raw, cell"""
        np = _np()
        sq, m, tp, td, tcal, dm, lNp, lL, cell = self.eval_modes(a, S)
        return sq, m, np.where(dm, td + tcal, tp), lNp, lL, cell

    def calib_stride(self, lNp, td, unit_time=None):
        """(--calib-r1-stride k, the stream's predicted CPU) for a d-first
        sum at ln N' with d-first CPU td: the stream costs ~ t_plain / k with
        t_plain = td / r(N') (the ratio law, not the plain law, which is
        off by 0.3-1.35x at N' >= 5k), k ~ 1 / (calib_frac r) so ~ calib_frac
        x td; with unit_time, k raised so that it takes at most
        DFIRST_STREAM_MAX x unit_time ((0, 0) without a stream)"""
        if self.calib_frac <= 0:
            return 0, 0.0
        r = math.exp(float(self.log_ratio(lNp)))
        tpr = td / r
        k = max(1, int(round(1 / (self.calib_frac * r))))
        if unit_time:
            k = max(k, int(math.ceil(tpr / (DFIRST_STREAM_MAX * unit_time))))
        return k, tpr / k

    def next_unit(self, a, lo, hi_max, unit_time, drop, detail=False, dpart=None, calib_done=()):
        """the next unit of cand a: sums lo.. (at most hi_max) in the mode
        of sum lo, while the cumulative predicted time stays within
        unit_time (at least one sum) and each sum's magic squares per
        CPU-second stay >= drop x the best of the unit so far (see
        Scorer.next_unit); None if lo > hi_max. A d-first sum predicted to
        take more than DFIRST_SPLIT x unit_time, or searched in part before
        (dpart: S -> (nd or None, merged d intervals, its star cover (x*,
        K) or None)), gets a unit of d (_drange_unit). calib_done: the sums whose calibration stream has
        run (or is planned)."""
        np = _np()
        am = _am()
        if lo > hi_max:
            return None
        dpart = dpart or {}
        # the mode of the unit is that of its first sum inside the model's
        # grid (the sums below it have no predicted squares and go with
        # either mode; if that first sum gets units of d, the sums below
        # are skipped, as Planner._walk skips such gaps)
        Sg = self.prof(a)[0]
        g0 = int(math.ceil(Sg[0] - 1e-9)) if len(Sg) else lo
        first = max(lo, g0) if g0 <= hi_max else lo
        sq0, m0, tp0, td0, tc0, dm0, lN0, _, cell0 = self.eval_modes(a, [float(first)])
        dfirst = bool(dm0[0])
        if first in dpart and self.policy != "off":
            dfirst = True
        if not dfirst:
            # (--dfirst off: a sum searched in part is searched plain, whole)
            dpart = {}
        cal0 = first in calib_done
        if dfirst and (first in dpart
                       or td0[0] + (0 if cal0 else tc0[0]) > DFIRST_SPLIT * unit_time):
            return self._drange_unit(a, first, unit_time, dpart.get(first), cal0, sq0[0], m0[0],
                                     tp0[0], td0[0], tc0[0], lN0[0], cell0[0], detail)
        parts = np.array(sorted(dpart), float) if dpart else None
        cals = np.array(sorted(calib_done), float) if (dfirst and calib_done) else None
        # the calibration stride of a d-first unit (one for all its sums)
        k = 0 if (not dfirst or cal0) else self.calib_stride(lN0[0], td0[0], unit_time)[0]
        acc_t = acc_m = acc_sq = 0.0
        best = 0.0
        hi = lo - 1
        start = lo
        size = 64
        cells, mv = [], []
        while start <= hi_max:
            S = np.arange(start, min(start + size, hi_max + 1), dtype=float)
            sq, m, tp, td, tc, dm, lNp, lL, cell = self.eval_modes(a, S)
            if dfirst:
                # every sum of the unit gets the stream, or none (the
                # squares recorded are the stream's)
                t = td + (tc if k else 0.0)
                same = dm.copy()
                if cals is not None:
                    same &= np.isin(S, cals) == cal0
                sq = sq / k if k else np.zeros_like(sq)
            else:
                t = tp
                same = ~dm
            same |= S < g0
            if parts is not None:
                same &= ~np.isin(S, parts)
            dens = m / t
            ct = np.cumsum(t)
            ct += acc_t
            rm = np.maximum.accumulate(dens)
            np.maximum(rm, best, out=rm)
            # (the sums below the grid go with the first sum inside it)
            ok = ((ct <= unit_time) | (S <= g0)) & same
            ok[1:] &= dens[1:] >= drop * rm[:-1]
            ok[0] &= dens[0] >= drop * best
            if acc_t == 0:
                ok[0] = True
            cut = len(S) if ok.all() else int(ok.argmin())
            if cut:
                acc_t = float(ct[cut - 1])
                acc_m += float(m[:cut].sum())
                acc_sq += float(sq[:cut].sum())
                best = float(rm[cut - 1])
                hi = int(S[cut - 1])
                if detail:
                    cells.append(cell[:cut])
                    mv.append(m[:cut])
            if cut < len(S):
                break
            start = int(S[-1]) + 1
            size = min(size * 4, 2048)
        tot_t = acc_t + am.UNIT_OVERHEAD
        return UnitV2(self.c.P(a), lo, hi, acc_m / tot_t, tot_t, acc_sq, acc_m, a, None,
                      np.concatenate(cells) if detail else None,
                      np.concatenate(mv) if detail else None,
                      "dfirst" if dfirst else "plain", None, None, 0, k, 1.0,
                      dfirst_star_k(self.n) if dfirst else 0, None)

    def _drange_unit(self, a, S, unit_time, part, cal_done, sq, m, tp, td, tc, lNp, cell, detail):
        """a unit of d of the d-first sum S: from the first d not searched
        yet (part: (nd or None, merged d intervals, star cover)) for about
        unit_time (the d loop by the cost profile along d, plus the
        calibration stream), to the end of the gap; the sum's calibration
        stream goes with the first unit that runs. Its star cover is that of
        the parts (their x* and K, so that its chunks merge with theirs),
        else v2's K with msearch's x*"""
        np = _np()
        am = _am()
        nd, iv, star = part if part else (None, [], None)
        if star is None:
            star = (None, dfirst_star_k(self.n))
        # the number of d: known from the records, else predicted (N' is
        # the model's prediction of nvecs_raw)
        ndp = nd if nd else max(int(round(math.exp(float(lNp)))), 1)
        # (the score keeps the model's N', as grid_density and the upper
        # bounds of the lazy heap)
        td_s, tc_s = td, tc
        if nd:
            # the sum's own N for its CPU (the d-first law's slope; E stays
            # the model's)
            lNx = math.log(nd)
            f = math.exp(float(self.tmd.th[1]) * (lNx - float(lNp)))
            td = am.SUM_OVERHEAD + (td - am.SUM_OVERHEAD) * f
            tc *= f
            lNp = lNx
        g0, gap_end = 0, ndp
        for lo_, hi_ in iv:
            if lo_ <= g0:
                g0 = max(g0, hi_)
            else:
                gap_end = lo_
                break
        if g0 >= ndp:
            # more d than predicted (only without a known nd): to the end
            ndp = g0 + 1
            gap_end = ndp
        # the stream (at most DFIRST_STREAM_MAX units) runs in this unit
        # after its d loop: the d loop gets the rest of the unit, at least
        # DFIRST_FIRST_MIN of one
        k, tcal = (0, 0.0) if cal_done else self.calib_stride(lNp, td, unit_time)
        w0 = float(am.dfirst_cost_cum(g0 / ndp))
        unit = unit_time / max(td, 1e-9)
        want = max(unit_time - tcal, DFIRST_FIRST_MIN * unit_time) / max(td, 1e-9)
        dhi = int(math.ceil(ndp * float(am.dfirst_cost_inv(w0 + want))))
        dhi = min(max(dhi, g0 + 1), gap_end)
        if am.dfirst_cost_frac(dhi, gap_end, ndp) < DFIRST_TAIL * unit:
            dhi = gap_end
        frac = am.dfirst_cost_frac(g0, dhi, ndp)
        time_ = frac * td + tcal + am.UNIT_OVERHEAD
        magic = frac * m
        # the stream's cost spread over the sum's parts in the score (as in
        # grid_density), so the parts of a sum rank alike
        score = magic / (frac * (td_s + (0.0 if cal_done else tc_s)) + am.UNIT_OVERHEAD)
        open_end = dhi >= ndp and not nd
        return UnitV2(self.c.P(a), S, S, score, time_, sq / k if k else 0.0, magic, a, None,
                      np.array([cell]) if detail else None, np.array([magic]) if detail else None,
                      "dfirst", g0, None if open_end else dhi, ndp, k, frac, star[1],
                      star[0] if star[1] else None)

    def unit_nodes(self, u):
        """predicted search nodes of the unit's largest (last) sum"""
        am = _am()
        _, _, _, lNp, lL, _ = self.eval_sums(u.a, [float(u.hi)])
        return float(math.exp(am.log_nodes(lNp[0], lL[0], int(self.c.k[u.a]))))

    def limits(self, u, unit_time):
        """msearch --node-limit (per sum) and --time-limit (per unit)"""
        # (at least 1 s: msearch reads --time-limit 0 as none)
        return int(max(2e10, 10 * self.unit_nodes(u))), max(2 * unit_time, 3 * u.time, 1.0)

    def set_factors(self, summary, explore=EXPLORE_V2):
        """per-P gamma factors f_sq and F_m from each P's observed counts and
        its per-cell expectations at the current class rates"""
        self.explore = explore
        self.lnfsq[:] = 0.0
        self.lnFm[:] = 0.0
        for key in summary.perP:
            self.set_factor(summary, key)

    def set_factor(self, summary, key):
        """per-P factors of one P (by key) after new results"""
        np = _np()
        a = self.c.index().get(key)
        st = summary.perP.get(key)
        if a is None or not st or not st["cells"]:
            return
        self._cache.pop(a, None)
        cs = np.fromiter(st["cells"].keys(), int)
        v = np.array(list(st["cells"].values()))
        e_sq = float((v[:, 2] * np.exp(self.ln_gsq[cs])).sum())
        e_S = float((v[:, 4] * np.exp(self.ln_rS[cs])).sum())
        e_P = float((v[:, 6] * np.exp(self.ln_rP[cs])).sum())
        # observed counts over the same sums as the expectations (inside the
        # model's grid, not truncated)
        o_sq, o_S, o_P = float(v[:, 1].sum()), float(v[:, 3].sum()), float(v[:, 5].sum())
        f_sq, f_S, f_P = per_p_factors(o_sq, e_sq, o_S, e_S, o_P, e_P,
                                       getattr(self, "explore", EXPLORE_V2))
        # E[(f_S f_P)^2] under the gamma posteriors, relative to the prior's
        # (so a P without results keeps factor 1)
        F_m = (f_S ** 2 * f_P ** 2 * (1 + 1 / (A_TRAV + o_S)) * (1 + 1 / (A_TRAV + o_P))
               / (1 + 1 / A_TRAV) ** 2)
        self.lnfsq[a] = math.log(f_sq)
        self.lnFm[a] = math.log(F_m)


def per_p_factors(o_sq, e_sq, o_S, e_S, o_P, e_P, explore=0.0):
    """gamma posterior means (+ explore x sd) of the per-P factors on
    squares and on S and P traversals: (a + observed) / (a + expected), the
    squares tempered by their overdispersion (o / PHI_SUM, e / PHI_SUM)"""
    def post(o, e, a):
        mean = (a + o) / (a + e)
        return mean + explore * math.sqrt(a + o) / (a + e)
    return (post(o_sq / PHI_SUM, e_sq / PHI_SUM, A_SQ), post(o_S, e_S, A_TRAV),
            post(o_P, e_P, A_TRAV))


class Planner:
    """lazy greedy over the candidates: a max-heap keyed by an upper bound
    of each P's density (1.1 x the best grid point at or after its
    frontier), replaced by the exact score of its next unit when popped"""

    UB_MARGIN = 1.1

    def __init__(self, cands, scorer, f0, cover, unit_time, drop, max_sums=2000, dpart=None,
                 calib=None):
        np = _np()
        self.c = cands
        self.sc = scorer
        self.f0 = np.asarray(f0, np.int64)     # start: max(S_min bound, legacy maxS + 1)
        self.cover = cover                     # a -> merged [[lo, hi], ...]
        # d-first sums searched in part: a -> {S: (nd or None, merged [[lo, hi), ...] of d)}
        self.dpart = dpart if dpart is not None else {}
        # sums whose calibration stream has run (or is planned): a -> set(S)
        self.calib = calib if calib is not None else {}
        self.unit_time = unit_time
        self.drop = drop
        self.max_sums = max_sums
        self.frontier = self.f0.copy()
        for a in cover:
            self.frontier[a] = self._walk(a, self.f0[a])
        self.heap = []

    GAP_NEGLIGIBLE = 1e-3

    def _walk(self, a, f):
        """the first sum >= f not covered; a gap below the model's first
        valid point (no squares predicted there, e.g. between ceil(S0) and
        the S_min another search started from) counts as covered, and so
        does a gap whose predicted density stays below GAP_NEGLIGIBLE x the
        P's best (else its tiny score would hold back the P's later sums)"""
        cov = self.cover.get(a, ())
        if not cov:
            return f
        Sg = self.sc.prof(a)[0]
        zero_below = Sg[0] if len(Sg) else float("inf")
        for lo, hi in cov:
            if lo > f and lo - 1 >= zero_below and not self._negligible(a, f, lo - 1):
                break
            f = max(f, hi + 1)
        return f

    def _negligible(self, a, lo, hi):
        np = _np()
        best = float(np.exp(self.sc.grid_density([a])[0].max()))
        if not best > 0:
            return True
        S = np.unique(np.linspace(lo, hi, 9).round())
        sq, m, t, _, _, _ = self.sc.eval_sums(a, S)
        return float((m / t).max()) < self.GAP_NEGLIGIBLE * best

    def set_cover(self, a, intervals, dpart=None, calib=None):
        self.cover[a] = [list(x) for x in intervals]
        if dpart is not None:
            self.dpart[a] = dpart
        if calib is not None:
            self.calib[a] = calib
        self.frontier[a] = self._walk(a, self.f0[a])

    def advance(self, a, hi):
        """planned (simulated) coverage up to hi"""
        self.frontier[a] = self._walk(a, max(self.frontier[a], hi + 1))

    def advance_unit(self, u):
        """planned (simulated) coverage after unit u: its sums, or its part
        of a d-first sum (the sum is covered once its parts cover every d)"""
        a = u.a
        if u.calib:
            self.calib.setdefault(a, set()).update(range(u.lo, u.hi + 1))
        if u.dlo is None:
            self.advance(a, u.hi)
            return
        dp = self.dpart.setdefault(a, {})
        nd, iv, star = dp.get(u.lo, (None, [], None))
        end = u.dhi if u.dhi is not None else max(u.nd, nd or 0)
        iv = _merge_half([list(x) for x in iv] + [[u.dlo, end]])
        if u.dhi is None or _covers(iv, nd or u.nd):
            dp.pop(u.lo, None)
            self.cover[a] = _merge([tuple(x) for x in self.cover.get(a, [])] + [(u.lo, u.lo)])
            self.frontier[a] = self._walk(a, self.frontier[a])
        else:
            # (the simulated parts keep the unit's K; x* stays as known)
            dp[u.lo] = (nd, iv, star if star is not None else (u.star_x, u.star_k))

    def hi_max(self, a):
        f = int(self.frontier[a])
        hm = min(f + self.max_sums - 1, int(self.c.end[a]))
        for lo, hi in self.cover.get(a, ()):
            if lo > f:
                hm = min(hm, lo - 1)
                break
        return hm

    def unit(self, a, detail=False):
        return self.sc.next_unit(a, int(self.frontier[a]), self.hi_max(a), self.unit_time,
                                 self.drop, detail, self.dpart.get(a), self.calib.get(a, ()))

    def upper_bounds(self, rows, chunk=20000):
        np = _np()
        am = _am()
        rows = np.asarray(rows, np.int64)
        out = np.zeros(len(rows))
        for i in range(0, len(rows), chunk):
            r = rows[i:i + chunk]
            d = self.sc.grid_density(r)
            uf = self.frontier[r] / self.c.S0[r] - 1
            jf = np.searchsorted(am.GRID_U, uf) - 1
            d = np.where(np.arange(d.shape[1])[None, :] >= jf[:, None], d, -np.inf)
            # (grid_density charges every d-first sum its calibration
            # stream; a sum whose stream has run scores up to 1 +
            # calib_frac higher)
            ub = self.UB_MARGIN * np.exp(d.max(1))
            if self.sc.policy != "off":
                ub *= 1 + max(self.sc.calib_frac, 0.0)
            out[i:i + chunk] = np.where(self.frontier[r] <= self.c.end[r], ub, 0.0)
        return out

    def build(self, busy=()):
        import heapq
        np = _np()
        ub = self.upper_bounds(np.arange(len(self.c)))
        busy = set(busy)
        self.heap = [(-float(u), int(a), None) for a, u in enumerate(ub)
                     if u > 0 and a not in busy]
        heapq.heapify(self.heap)

    def pop(self):
        """the best next unit (its P leaves the heap), or None"""
        import heapq
        while self.heap:
            negk, a, u = heapq.heappop(self.heap)
            if u is not None:
                return u
            u = self.unit(a)
            if u is None or not u.score > 0:
                continue
            if not self.heap or u.score >= -self.heap[0][0]:
                return u
            heapq.heappush(self.heap, (-u.score, a, u))
        return None

    def push(self, a):
        import heapq
        u = self.unit(a)
        if u is not None and u.score > 0:
            heapq.heappush(self.heap, (-u.score, int(a), u))


# --------------------------------------------------------------------------
# the v2 scheduler


def read_legacy(state, n, no_legacy=False):
    out = {}
    path = os.path.join(state, "legacy.json")
    if os.path.exists(path) and not no_legacy:
        with open(path) as f:
            for r in json.load(f):
                if r.get("n", 6) == n:
                    out[norm_p(r["P"])] = r
    return out


def read_exact_smin(state, n):
    """exact S_min of the P already in pinfo (no enumerate call)"""
    path = os.path.join(state, f"pinfo_{n}.json")
    out = {}
    if os.path.exists(path):
        with open(path) as f:
            for k, d in json.load(f).items():
                nn, e = k.split(":")
                if int(nn) == n and d.get("smin"):
                    out[e] = d["smin"]
    return out


def dfirst_args(u):
    """msearch options of a d-first unit: --diag-first with threshold 0 (the
    scheduler chose the mode of every sum of the unit), its d range with
    DFIRST_CHUNKS checkpoints (--d-chunk; msearch's default 256 d would make
    a unit of < 256 d one chunk, without a checkpoint or a --time-limit
    stop inside it), its calibration stream, and its star cover, always
    explicit (--dfirst-star K, so that a unit's cost and records do not
    depend on the binary's default, and --dfirst-star-x x* for a unit that
    continues the parts of a sum searched with that x*)"""
    if u.mode != "dfirst":
        return []
    out = ["--diag-first", "--diag-first-min-n", "0", "--dfirst-star", str(u.star_k)]
    if u.star_k and u.star_x is not None:
        out += ["--dfirst-star-x", str(u.star_x)]
    if u.dlo is not None:
        end = u.dhi if u.dhi is not None else max(u.nd, u.dlo + 1)
        chunk = max(8, -(-(end - u.dlo) // DFIRST_CHUNKS))
        out += ["--d-range", f"{u.dlo}:{'' if u.dhi is None else u.dhi}", "--d-chunk", str(chunk)]
    if u.calib:
        out += ["--calib-r1-stride", str(u.calib)]
    return out


def describe_dunit(u):
    """" d-first ..." for log lines and plans ("" for a plain unit)"""
    if u.mode != "dfirst":
        return ""
    out = " d-first"
    if u.dlo is not None:
        end = "end" if u.dhi is None else u.dhi - 1
        out += f" d {u.dlo}..{end} of {u.nd} ({100 * u.frac:.0f}%)"
    if u.calib:
        out += f" +calib 1/{u.calib}"
    if u.star_x is not None:
        out += f" star x*={u.star_x} K={u.star_k}"
    return out


def unit_tag(u):
    """the suffix of a unit file's name: "" plain, "D" d-first sums,
    "d<lo>-<hi>" a d range (hi "end": to the end of the sum)"""
    if u is None or u.mode != "dfirst":
        return ""
    if u.dlo is None:
        return "D"
    return f"d{u.dlo}-{'end' if u.dhi is None else u.dhi}"


class SchedulerV2:
    def __init__(self, args, quiet=False):
        np = _np()
        self.args = args
        self.n = args.vec_size
        if self.n != 6:
            sys.exit("--model analytic is calibrated for 6x6 only: use --model regression")
        self.dir = args.state
        self.units = os.path.join(self.dir, "units")
        os.makedirs(self.units, exist_ok=True)
        self.quiet = quiet
        self.legacy = read_legacy(self.dir, self.n, getattr(args, "no_legacy", False))
        self.smin_exact = read_exact_smin(self.dir, self.n)
        only = getattr(args, "only", None)
        self.pool = None if only else Pool.load(args, self.dir)
        self.store = ProfileStore(self.dir, self.pool, self.n)
        self.summary = Summary.load(self.dir, self.n)
        if self.pool is not None and len(self.pool):
            # profile the pool P with results first (the summary needs them)
            self._prefill_touched()
        self.summary.update(self.units, self.store)
        # candidates
        if only:
            Ps = [parse_p(p) for p in only.split(";")]
            arr = Pool.arrays_for(Ps)
            self.cands = Cands(arr["exps"], [-1] * len(Ps), arr["ratio"], self.n)
        else:
            keep = np.ones(len(self.pool), bool)
            idx = self.pool.index()
            removed = 0
            for P, lg in self.legacy.items():
                r = idx.get(p_str(P, "_"))
                if r is not None and lg["maxS"] >= math.floor(2 * am_s0(P)):
                    keep[r] = False
                    removed += 1
            self.legacy_removed = removed
            rows = np.nonzero(keep)[0]
            extra = [k for k in set(self.summary.perP) | set(self.summary.cover)
                     if k not in idx and len(k.split("_")) <= len(PRIMES)]
            extra = sorted(parse_p(k) for k in extra)
            ea = Pool.arrays_for(extra)
            exps = np.concatenate([self.pool.exps[rows], ea["exps"]])
            src = np.concatenate([rows, -np.ones(len(extra), np.int64)])
            ratio = np.concatenate([self.pool.ratio[rows], ea["ratio"]])
            self.cands = Cands(exps, src, ratio, self.n)
        self.fill_profiles()
        self.stage1 = Stage1(self.dir, self.n, parse_stage1(getattr(args, "stage1", None)))
        self.refit(save=False)
        self.scorer.stage1 = self.stage1.band()

    def _prefill_touched(self):
        idx = self.pool.index()
        rows = []
        for path in unit_files(self.units):
            name = os.path.basename(path)
            if name.endswith(".gz"):
                name = name[:-3]
            if name in self.summary.files:
                continue
            # the P of a scheduler unit is in its name; others are found when parsed
            parts = name[:-6].split("_")
            for i in range(len(parts)):
                r = idx.get("_".join(parts[i:-2]))
                if r is not None:
                    rows.append(r)
                    break
        if rows:
            self.store.fill(rows, getattr(self.args, "profile_workers", 2), quiet=self.quiet)

    def fill_profiles(self):
        np = _np()
        src = self.cands.src
        rows = src[src >= 0]
        if len(rows) and self.store.done is not None:
            missing = rows[self.store.done[rows] == 0]
            if len(missing):
                self.store.fill(missing, getattr(self.args, "profile_workers", 2), quiet=self.quiet)
        extra = [self.cands.P(a) for a in np.nonzero(src < 0)[0]]
        self.store.fill_extra(extra)
        self.store.save_extra()

    def f0(self):
        np = _np()
        c = self.cands
        f = np.ceil(c.S0 - 1e-9).astype(np.int64)
        idx = c.index() if (self.smin_exact or self.legacy) else {}
        for key, smin in self.smin_exact.items():
            a = idx.get(key)
            if a is not None:
                f[a] = smin
        for P, lg in self.legacy.items():
            a = idx.get(p_str(P, "_"))
            if a is not None:
                f[a] = max(f[a], min(lg["maxS"], 10 ** 9) + 1)
        return f

    def cover(self):
        idx = self.cands.index()
        return {idx[k]: v for k, v in self.summary.cover.items() if k in idx}

    def dpart_of(self, key):
        """the d-first sums of P (key) searched in part: {S: (nd, iv, star)},
        iv and star = (x* or None, K) those of the coverage group the
        planner continues (Summary.dfirst_part)"""
        out = {}
        for S in self.summary.dcov.get(key, {}):
            nd, iv, _, star = self.summary.dfirst_part(key, S)
            out[int(S)] = (nd, [list(x) for x in iv], star)
        return out

    def calib_of(self, key):
        """the sums of P (key) whose calibration stream has run"""
        return {int(S) for S in self.summary.perP.get(key, {}).get("calS", {})}

    def dpart(self):
        idx = self.cands.index()
        return {idx[k]: self.dpart_of(k) for k in self.summary.dcov if k in idx}

    def calib_done(self):
        idx = self.cands.index()
        return {idx[k]: self.calib_of(k) for k, st in self.summary.perP.items()
                if k in idx and st.get("calS")}

    def refit(self, save=True):
        """class factors, time models and per-P factors from the summary"""
        table = self.summary.cell_table()
        self.calib = Calibration.fit(table, self.summary.cell_ee())
        self.time_models = fit_time_models(self.summary.time_stats())
        self.tm = current_time_model(self.time_models)
        self.tmd = current_time_model(self.time_models, mode="dfirst")
        self.lr0, self.ratio_st = current_ratio_level(self.summary.ratio)
        explore = self.args.explore if self.args.explore is not None else EXPLORE_V2
        if not hasattr(self, "scorer"):
            a = self.args
            self.scorer = AnalyticScorer(self.cands, self.store, self.calib, self.tm, self.n,
                                         self.tmd, getattr(a, "dfirst", "auto") or "auto",
                                         getattr(a, "calib_frac", CALIB_FRAC),
                                         getattr(a, "dfirst_min_n", DFIRST_MIN_NP), self.lr0)
        else:
            self.scorer.set_calibration(self.calib)
            self.scorer.tm = self.tm
            self.scorer.tmd = self.tmd
            self.scorer.lr0 = self.lr0
        self.scorer.set_factors(self.summary, explore)
        if save:
            self.save()

    def save(self):
        with open(os.path.join(self.dir, f"calib_{self.n}.json"), "w") as f:
            json.dump(self.calib.to_json(), f, indent=1)
        with open(os.path.join(self.dir, f"time_{self.n}.json"), "w") as f:
            json.dump({k: v.to_json() for k, v in self.time_models.items()}
                      | {"current": self.tm.to_json(), "current_dfirst": self.tmd.to_json(),
                         "dfirst_ratio": {"a0": self.lr0, "a1": _am().dfirst_ratio_coefs(ENGINE)[1],
                                          "pairs": self.ratio_st}},
                      f, indent=1)
        self.summary.save()

    def planner(self, busy=()):
        p = Planner(self.cands, self.scorer, self.f0(), self.cover(), self.args.unit_time,
                    self.args.unit_drop, dpart=self.dpart(), calib=self.calib_done())
        p.build(busy)
        return p

    def command(self, u, out, check=True):
        nl, tl = self.scorer.limits(u, self.args.unit_time)
        # (the limits are in reference CPU: a slower machine gets more time)
        tl /= getattr(self, "machine_slow", 1.0)
        return [binary("msearch", check), "--vec-size", str(self.n), "--min-sum", str(u.lo),
                "--max-sum", str(u.hi), *dfirst_args(u), "--time-limit", f"{tl:.0f}",
                "--node-limit", str(nl), "--out", out, *map(str, u.P)]

    def unit_path(self, P, lo, hi, u=None):
        # a running sequence number (globbing units/ on every launch is
        # O(files); the time stamp keeps names unique across runs)
        if not hasattr(self, "_seq"):
            self._seq = len(self.summary.files)
        seq = self._seq
        self._seq += 1
        name = f"{int(time.time())}_{seq:06d}_{p_str(P, '_')}_{lo}_{hi}{unit_tag(u)}.jsonl"
        return os.path.join(self.units, name)

    def set_machine(self, path):
        """run --machine cal.json: write a machine record at the top of each
        unit file (the summary scales its records' CPU to reference CPU),
        log the machine in the state (machines_<n>.jsonl), and scale
        --time-limit by the slower of its per-process speeds"""
        with open(path, "rb") as f:
            raw = f.read()
        cal = json.loads(raw)
        pp = cal.get("per_process") or {}
        if not pp.get("plain") or not pp.get("dfirst"):
            sys.exit(f"{path}: no per_process plain and dfirst speeds (run scripts/machine_cal.py)")
        self.machine = {"id": hashlib.sha256(raw).hexdigest()[:12],
                        "speed": {"plain": float(pp["plain"]), "dfirst": float(pp["dfirst"])}}
        self.machine_slow = min(1.0, self.machine["speed"]["plain"], self.machine["speed"]["dfirst"])
        with open(os.path.join(self.dir, f"machines_{self.n}.jsonl"), "a") as f:
            f.write(json.dumps({"t": time.time(), "id": self.machine["id"], "file": os.path.abspath(path),
                                "host": cal.get("host"), "model": (cal.get("cpu") or {}).get("model"),
                                "path": (cal.get("build") or {}).get("path"),
                                "workers": cal.get("workers"), "per_process": self.machine["speed"],
                                "plain_speed": cal.get("plain_speed"),
                                "dfirst_speed": cal.get("dfirst_speed")}) + "\n")
        if cal.get("workers") and cal["workers"] != self.args.workers:
            log(f"note: {path} recommends {cal['workers']} workers (running {self.args.workers})")
        log(f"machine {path} ({self.machine['id']}): per-process speed plain "
            f"{pp['plain']:.3f}, d-first {pp['dfirst']:.3f} reference CPU per CPU second")

    def run(self):
        args = self.args
        if getattr(args, "machine", None):
            self.set_machine(args.machine)
        plan = self.planner()
        deadline = time.time() + args.hours * 3600 if args.hours else None
        running = {}
        running_t = {}
        units_done = 0
        stopping = False
        launched = open(os.path.join(self.dir, f"launched_{self.n}.jsonl"), "a")
        nseen = len(self.summary.notable)

        def handle_sigint(sig, frame):
            nonlocal stopping
            stopping = True
            log("stopping: waiting for running units (Ctrl-C again to kill them)")
            signal.signal(signal.SIGINT, signal.default_int_handler)

        signal.signal(signal.SIGINT, handle_sigint)
        log(f"{len(self.cands)} candidate values of P, {args.workers} workers "
            f"(model analytic, pool {args.pool})")
        if self.stage1.spec:
            log(self.stage1_status())
        try:
            while True:
                for proc in list(running):
                    if proc.poll() is None:
                        continue
                    a, path = running.pop(proc)
                    running_t.pop(proc, None)
                    units_done += 1
                    self.summary.update_file(path, self.store)
                    key = self.cands.key(a)
                    for q in self.summary.notable[nseen:]:
                        if q["best_score"] >= args.announce_score:
                            tag = "MAGIC SQUARE" if q["best_score"] >= 14 else "square"
                            log(f"{tag}: P={p_str(q['P'])} S={q['S']} best={q['best_score']} "
                                f"#S={q['s_count']} #P={q['p_count']} #SP={q['sp_count']} "
                                f"{q['grid']}")
                    nseen = len(self.summary.notable)
                    plan.set_cover(a, self.summary.cover.get(key, []), self.dpart_of(key),
                                   self.calib_of(key))
                    if units_done % args.refit_every == 0:
                        self.refit()
                        log(f"refit: {describe_calib(self.calib, self.tm, self.tmd, self.lr0, self.ratio_st)}")
                        plan.build({b for b, _ in running.values()})
                    else:
                        self.scorer.set_factor(self.summary, key)
                        plan.push(a)
                if stopping or (deadline and time.time() > deadline):
                    if not running:
                        break
                    time.sleep(0.5)
                    continue
                exhausted = False
                while len(running) < args.workers:
                    u = plan.pop()
                    if u is None:
                        exhausted = True
                        break
                    path = self.unit_path(u.P, u.lo, u.hi, u)
                    fresh = self.cands.key(u.a) not in self.summary.perP
                    if getattr(self, "machine", None):
                        # (msearch appends to its --out file)
                        with open(path, "a") as f:
                            f.write(json.dumps(machine_record(u.P, u.lo, self.machine["id"],
                                                              self.machine["speed"])) + "\n")
                    proc = subprocess.Popen(self.command(u, path), stdout=subprocess.DEVNULL,
                                            env=msearch_env())
                    running[proc] = (u.a, path)
                    running_t[proc] = u.time
                    in_stage1 = self.stage1.active
                    # (stage 1: the model's predictions of the unit's records,
                    # frozen at launch for scripts/decide.py)
                    pred = unit_predictions(self.scorer, u) if self.stage1.spec else {}
                    if self.stage1.charge(u.time):
                        log(f"stage 1 done ({self.stage1.spent / 3600:.1f} reference CPU-hours "
                            f"launched): back to --dfirst {self.scorer.policy} at every N'")
                        self.scorer.stage1 = None
                        plan.build({b for b, _ in running.values()})
                    self.stage1.save()
                    launched.write(json.dumps({
                        "file": os.path.basename(path), "P": list(u.P), "lo": u.lo, "hi": u.hi,
                        "time": u.time, "squares": u.squares, "magic": u.magic, "score": u.score,
                        "fresh": fresh, "t": time.time(), "mode": u.mode, "dlo": u.dlo,
                        "dhi": u.dhi, "nd": u.nd, "calib": u.calib, "frac": u.frac,
                        "star_k": u.star_k, "star_x": u.star_x}
                        | ({"stage1": in_stage1} | pred if self.stage1.spec else {})) + "\n")
                    launched.flush()
                    log(f"start P={p_str(u.P)} S={u.lo}..{u.hi}{describe_dunit(u)} (predicted "
                        f"{u.time:.0f}s, {u.squares:.1f} squares, {u.score * 3.15e7:.3g} "
                        f"magic/CPU-year)")
                if exhausted and not running:
                    log("no more units to run (every candidate is covered or has no "
                        "predicted magic squares)")
                    break
                # poll often while the units are short (the top of the ranking
                # has units of a few seconds)
                time.sleep(0.1 if any(t < 30 for t in running_t.values()) else 0.5)
        except KeyboardInterrupt:
            for proc in running:
                proc.terminate()
            for proc, (a, path) in running.items():
                proc.wait()
                self.summary.update_file(path, self.store)
        launched.close()
        self.store.save_extra()
        self.refit()
        log("stopped")
        report_v2(self, top=10)

    def simulate(self, max_units=None, max_hours=None, detail=False):
        """greedy simulation with predicted times and squares: yields units"""
        # (stage 1: the simulated units count on a copy of the state's clock)
        st1 = Stage1.__new__(Stage1)
        st1.__dict__.update(self.stage1.__dict__)
        self.scorer.stage1 = st1.band()
        plan = self.planner()
        hours, i = 0.0, 0
        try:
            while (max_units is None or i < max_units) and (max_hours is None or hours < max_hours):
                u = plan.pop()
                if u is None:
                    break
                if detail:
                    u = plan.unit(u.a, detail=True)
                yield u
                hours += u.time / 3600
                i += 1
                plan.advance_unit(u)
                if st1.charge(u.time):
                    self.scorer.stage1 = None
                    plan.build()
                else:
                    plan.push(u.a)
        finally:
            self.scorer.stage1 = self.stage1.band()

    def stage1_status(self):
        st = self.stage1
        if not st.spec:
            return "stage 1: off"
        h, lo, hi = st.spec
        return (f"stage 1: N' {lo:.0f}-{hi:.0f} plain for the first {h:g} reference CPU-hours; "
                f"{st.spent / 3600:.2f} h launched so far ({'on' if st.active else 'done'})")


def am_s0(P):
    return _am().s0(P)


def describe_calib(calib, tm, tmd=None, lr0=None, rst=None):
    """one line: the data behind the class factors, and the current time
    laws (plain; d-first: its level; the d-first / plain ratio's level)"""
    d = calib.data
    parts = [f"{k} obs/exp {d[k][0]:.0f}/{d[k][1]:.1f}" for k in ("sq", "S", "P") if k in d]
    out = ("class-factor data: " + (", ".join(parts) or "none")
           + f"; time law (engine {tm.engine}): th=["
           + ", ".join(f"{v:.3f}" for v in tm.th) + f"] sd {tm.sd:.3f} ({tm.n} sums)")
    if tmd is not None:
        out += (f"; d-first law: {tmd.th[0]:.3f} + {tmd.th[1]:.3f} ln(N'/4000), sd {tmd.sd:.3f} "
                f"({tmd.n:.1f} sums, parts weighted)")
    if lr0 is not None:
        a1 = _am().dfirst_ratio_coefs(ENGINE)[1]
        sw = math.exp((math.log(1 + CALIB_FRAC) + lr0) / -a1) * 4000
        out += (f"; d-first/plain ratio: {lr0:+.3f} {a1:+.3f} ln(N'/4000) "
                f"({int(rst[0]) if rst else 0} pairs; auto switches at N' {sw:.0f})")
    return out


def _poisson_interval(o, e, z=1.645, phi=1.0):
    """obs/pred with an approximate 90% interval (quasi-Poisson with
    dispersion phi: the Poisson interval of o / phi, scaled back)"""
    if e <= 0:
        return float("nan"), float("nan"), float("nan")
    phi = max(float(phi), 1e-9)
    o2 = o / phi
    lo = max(0.0, o2 - z * math.sqrt(o2) + z * z / 4) if o2 > 0 else 0.0
    hi = o2 + z * math.sqrt(o2 + 1) + z * z / 2
    return o / e, phi * lo / e, phi * hi / e


def report_v2(sch, top=20, out=None):
    np = _np()
    s = sch.summary
    T = s.totals
    out = out or sys.stdout
    pr = lambda *a: print(*a, file=out)  # noqa: E731
    pr(f"\n{len(s.perP)} values of P, {T['sums']} sums, {T['cpu'] / 3600:.2f} CPU-hours, "
       f"{T['squares']} semi-magic squares ({3600 * T['squares'] / max(T['cpu'], 1):.0f}/CPU-hour)")
    if T.get("dfirst_sums") or T.get("dfirst_partial") or s.dcov:
        pr(f"d-first (in the d-first time law only): {T['dfirst_sums']} sums searched in full, "
           f"{T['dfirst_partial']} parts (dsum records of d ranges or samples), "
           f"{T['dfirst_cpu'] / 3600:.2f} CPU-hours (with the chunks of killed units), "
           f"{T['dfirst_magic']} magic squares"
           + (f", {T['dfirst_truncated']} chunks with a V_d at --node-limit"
              if T.get("dfirst_truncated") else ""))
        if T.get("dfirst_sums"):
            pr(f"  (square, SP traversal) pairs of the sums searched in full: "
               f"{T.get('dfirst_pairs', 0):.0f} found, {T.get('dfirst_est_pairs', 0.0):.1f} "
               f"estimated over every d of {T.get('dfirst_est_sums', 0)} sums"
               + (f" ({T['dfirst_star_sums']} with the star cover)"
                  if T.get("dfirst_star_sums") else ""))
        if T.get("dfirst_star_skipped"):
            pr(f"  {T['dfirst_star_skipped']} chunks in no coverage (--dfirst-star-only, or a "
               f"star cover without star_x)")
        part = [(k, int(S)) + s.dfirst_part(k, S)[:3] for k, d in s.dcov.items() for S in d]
        if part:
            pr(f"  sums searched in part ({len(part)}, to be continued from their last d): "
               + ", ".join(f"P={p_str(parse_p(k))} S={S} {100 * f:.0f}% of {nd if nd else '?'} d"
                           for k, S, nd, iv, f in sorted(part)[:top]))
    if T.get("calib_sums") or T.get("calib_dup") or T.get("calib_cpu"):
        r, lo, hi = _poisson_interval(T["calib_est"], T["calib_pred"], phi=PHI_SUM)
        pr(f"calibration streams of d-first sums (in the cells, weighted): {T['calib_sums']} sums, "
           f"{T['calib_squares']} sampled squares, {T['calib_cpu'] / 3600:.2f} CPU-hours"
           f" ({100 * T['calib_cpu'] / max(T['dfirst_cpu'], 1e-9):.0f}% of the d-first CPU, sums in "
           f"progress included); "
           f"est. squares / model x SQ12 = {T['calib_est']:.0f} / {T['calib_pred']:.1f} = "
           f"{r:.2f}" + (f"; {T['calib_dup']} repeated streams not counted" if T["calib_dup"]
                         else ""))
    if T.get("sampled_cpu"):
        pr(f"sampled records (r1-sampled, not in the fits): "
           f"{T['sampled_cpu'] / 3600:.2f} CPU-hours")
    names = {0: "0", 2: "S", 3: "P", 4: "S+S", 5: "S+P", 6: "P+P", 7: "SP",
             9: "SP+S", 10: "SP+P", 14: "MAGIC (SP+SP)"}
    pr("best pair of diagonals: " + ", ".join(
        f"{names.get(k, k)}: {v}" for k, v in sorted(s.types.items(), reverse=True)))
    pr(f"traversals: S {T['trav_S']}, P {T['trav_P']}, SP {T['trav_SP']}")
    pr(f"(sums outside the model's grid: {T['outside']}, truncated: {T['truncated']})")
    tab = s.cell_table()
    pr("\ncalibration (obs/pred with 90% quasi-Poisson intervals, dispersion from per-P "
       "factors; pred = analytic model x SQ12 for squares, 720 p for traversals, before the "
       "learned class factors):")
    groups = (("N' band", 0, NBAND_NAMES), ("k", 1, KBAND_NAMES), ("ratio", 2, RBAND_NAMES),
              ("x band", 3, XBAND_NAMES))
    parts = np.array([cell_parts(c) for c in range(NCELLS)])
    for gname, gi, labels in groups:
        pr(f"  by {gname}:")
        pr(f"    {'':10} {'sums':>8} {'squares obs/pred':>28} {'S trav obs/pred':>28} "
           f"{'P trav obs/pred':>28}")
        ee = s.cell_ee(parts[:, gi])
        for b, lab in enumerate(labels):
            sel = parts[:, gi] == b
            row = tab[sel].sum(0)
            if row[0] == 0 and row[1] == 0 and row[3] == 0:
                continue
            cols = []
            for j, (key, o, e) in enumerate((("sq", 1, 2), ("S", 3, 4), ("P", 5, 6))):
                phi = dispersion(row[e], ee[b, j], key) if b < len(ee) else 1.0
                r, lo, hi = _poisson_interval(row[o], row[e], phi=float(phi))
                cols.append(f"{row[o]:7.0f}/{row[e]:8.1f} {r:4.2f} [{lo:4.2f},{hi:4.2f}]")
            pr(f"    {lab:10} {row[0]:8.0f} " + " ".join(cols))
    tot = tab.sum(0)
    r, lo, hi = _poisson_interval(tot[7], tot[8])
    pr(f"  SP traversals: {tot[7]:.0f} observed / {tot[8]:.1f} predicted = {r:.2f} [{lo:.2f},{hi:.2f}]")
    pr("\nlearned class factors (posterior mean, se):")
    for key, name in (("sq", "squares"), ("S", "S traversals"), ("P", "P traversals")):
        b, C = sch.calib.beta[key], sch.calib.cov[key]
        pr(f"  {name:13} " + ", ".join(f"{e} {v:+.2f}({math.sqrt(C[i, i]):.2f})"
                                       for i, (e, v) in enumerate(zip(EFFECTS, b))))
    pr(f"  {describe_calib(sch.calib, sch.tm, getattr(sch, 'tmd', None), getattr(sch, 'lr0', None), getattr(sch, 'ratio_st', None))}")
    for k, tm in sorted(sch.time_models.items()):
        pr(f"    engine:mode {k}: th=[{', '.join(f'{v:.3f}' for v in tm.th)}] sd {tm.sd:.3f} "
           f"({tm.n} sums)")
    pr("\ntime obs/pred (mean ln(t / t_pred), current engine's law) by N band and label words W:")
    rows = {}
    for key, v in s.tband.items():
        eng, mode, band, W = key.split(":")
        tm = sch.time_models.get(f"{eng}:{mode}", sch.tm)
        nt = len(tm.th)
        resid = (v[nt + 1] - np.dot(v[1:nt + 1], tm.th)) / max(v[0], 1)
        rows.setdefault((int(eng), mode), []).append((NBAND_NAMES.index(band), int(W), v[0], resid))
    for (eng, mode), rs in sorted(rows.items()):
        pr(f"  engine {eng} {mode}: " + ", ".join(
            f"{NBAND_NAMES[b]} W{W}: {r:+.2f} (n {n:.0f})" for b, W, n, r in sorted(rs)))
    pr("N bias ln(nvecs_raw / N'): " + ", ".join(
        f"{b} {v[1] / v[0]:+.3f} (sd {math.sqrt(max(v[2] / v[0] - (v[1] / v[0]) ** 2, 0)):.3f}, "
        f"n {v[0]})" for b, v in sorted(s.nbias.items(), key=lambda kv: NBAND_NAMES.index(kv[0]))))
    pr("label bias ln(labels / labels_obs): " + ", ".join(
        f"{b} {v[1] / v[0]:+.3f} (n {v[0]})"
        for b, v in sorted(s.lbias.items(), key=lambda kv: NBAND_NAMES.index(kv[0]))))
    report_launched(sch, pr)
    best = sorted(s.notable, key=lambda q: (-q["best_score"], -q["sp_count"], q["S"]))[:top]
    if best:
        pr("\nbest squares:")
    for q in best:
        pr(f"  P={p_str(q['P']):18} S={q['S']:5} best={q['best_score']:2} "
           f"#S={q['s_count']} #P={q['p_count']} #SP={q['sp_count']}  {q['grid']}"
           + ("  (d-first)" if q.get("dfirst") else "")
           + ("  (calibration stream)" if q.get("calib") else ""))


def report_launched(sch, pr):
    """predicted vs observed squares of the units this scheduler launched,
    fresh P (no results before the unit) vs repeat units"""
    path = os.path.join(sch.dir, f"launched_{sch.n}.jsonl")
    if not os.path.exists(path):
        return
    rows = []
    with open(path) as f:
        for line in f:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            st = sch.summary.files.get(r["file"])
            if st and st.get("complete"):
                rows.append((r, st))
    drows = [(r, st) for r, st in rows if r.get("mode", "plain") == "dfirst"]
    rows = [(r, st) for r, st in rows if r.get("mode", "plain") != "dfirst"]
    if drows:
        o = sum(st.get("dcpu", 0.0) for r, st in drows)
        e = sum(r["time"] for r, st in drows)
        pr(f"\nlaunched d-first units (complete only): {len(drows)}, observed / predicted CPU "
           f"(d loop, calibration stream, overhead) {o:.0f} / {e:.0f} s = {o / max(e, 1e-9):.2f}")
    if not rows:
        return
    pr("\nlaunched units (complete only): observed / predicted squares (90% quasi-Poisson)")
    cut = sorted(r["score"] for r, _ in rows)[int(0.9 * (len(rows) - 1))]
    for name, sel in (("fresh P", lambda r: r["fresh"]), ("repeat", lambda r: not r["fresh"])):
        for scope, f2 in (("all", lambda r: True), ("top decile", lambda r: r["score"] >= cut)):
            o = sum(st["squares"] for r, st in rows if sel(r) and f2(r))
            e = sum(r["squares"] for r, st in rows if sel(r) and f2(r))
            nu = sum(1 for r, st in rows if sel(r) and f2(r))
            ee = sum(r["squares"] ** 2 for r, st in rows if sel(r) and f2(r))
            if nu:
                ratio, lo, hi = _poisson_interval(o, e, phi=float(dispersion(e, ee, "sq")))
                flag = "  <-- below 0.8" if name == "fresh P" and scope == "all" and ratio < 0.8 \
                    and hi < 1.0 else ""
                pr(f"  {name:8} {scope:10} {nu:6} units {o:8.0f} / {e:9.1f} = {ratio:.2f} "
                   f"[{lo:.2f},{hi:.2f}]{flag}")


# --------------------------------------------------------------------------
# v2 commands


def v2_plan(args):
    sch = SchedulerV2(args)
    if sch.stage1.spec:
        print(sch.stage1_status())
    plan = sch.planner()
    print(f"{'P':24} {'S range':>13} {'pred. time':>10} {'squares/h':>9} "
          f"{'P(magic)':>9} {'magic/CPU-year':>14}  mode (d-first: d range, calibration stride; "
          f"squares/h: those recorded, the stream's)")
    for _ in range(args.top):
        u = plan.pop()
        if u is None:
            break
        h = u.time / 3600
        pm = u.magic / max(u.squares, 1e-30) if u.mode == "plain" else float("nan")
        print(f"{p_str(u.P):24} {u.lo:6}-{u.hi:<6} {u.time:9.0f}s "
              f"{u.squares / max(h, 1e-9):9.0f} {pm:9.2e} "
              f"{u.score * 3.15e7:14.3g}  {u.mode}{describe_dunit(u)[8:]}")


def v2_emit(args):
    sch = SchedulerV2(args, quiet=True)
    units_dir = sch.units
    if os.path.abspath(units_dir).startswith(ROOT + os.sep):
        units_dir = os.path.relpath(units_dir, ROOT)  # portable plan
    for i, u in enumerate(sch.simulate(max_units=args.units)):
        out = os.path.join(units_dir,
                           f"plan{i:06d}_{p_str(u.P, '_')}_{u.lo}_{u.hi}{unit_tag(u)}.jsonl")
        print(" ".join(sch.command(u, out, check=False)[1:]))


FORECAST_TRUTHS = ("measured", "laws", "anchored")
# the "measured" truth: the measured CPU / the law's prediction per N' band
# (NBAND_NAMES index) in the calibration search in the target region
# (research/calibration-target.md 3.6, its archive/results/analyze.out
# "## Time", sum(obs) / sum(pred), as its forecast_update.py variant
# "time"): plain sums against the plain law (3-6k: 289 sums; 6-12k: 18;
# 12-24k and >= 24k: the streams' est_time of 4 and 3 d-first sums), d-first
# sums against the d-first law (7, 4 and 3 sums; at 3-6k, not measured,
# the plain factor, as the ratio law held there: 0.95); 1 below 3k (not
# measured). The plain law (engine 3) is right at 3-6k (0.97), where the
# "anchored" truth charges plain sums 1.5x too much (0.65 of it).
TIME_TRUTH_PLAIN = {2: 0.969, 3: 1.185, 4: 0.622, 5: 0.159}
TIME_TRUTH_DFIRST = {2: 0.969, 3: 0.848, 4: 0.558, 5: 0.379}


def unit_charge(sc, u, truth="laws"):
    """CPU seconds charged for a planned unit u (see unit_charge_parts)"""
    pl, df = unit_charge_parts(sc, u, truth)
    return pl + df


def unit_charge_parts(sc, u, truth="laws"):
    """(plain, d-first) CPU seconds charged for a planned unit u: the d-first
    part is the d-first search of a d-first unit, the plain part the rest
    (a plain unit, a d-first unit's calibration stream, which is a sampled
    plain search, and the unit's overhead); a machine's plain and d-first
    speeds (forecast --machine) convert each to instance time. In all: its
    predicted time (truth
    "laws"); under the "measured" truth (TIME_TRUTH_PLAIN, _DFIRST: the
    calibration search's measured CPU per N' band) a plain sum costs the
    plain law and a d-first sum the d-first law, each x the factor of its N'
    band, a calibration stream the plain cost / k; or under the "anchored"
    truth (research/scheduler-v2.md, "d-first units"): a d-first sum costs
    the d-first law, a plain sum the d-first law / r(N') at N' >= 5k (the
    measured ratio, and the d-first law back-tests at 0.94 where the plain
    law is off by 0.3-1.35x), the plain law at N' <= 3k, geometrically
    blended between; a calibration stream the plain cost / k"""
    np = _np()
    am = _am()
    if truth == "laws":
        if u.mode != "dfirst":
            return u.time, 0.0
        # (the d-first share of the predicted time: frac x the d-first law
        # of its sums, the rest is its calibration stream and overhead)
        S = (np.array([float(u.lo)]) if u.dlo is not None
             else np.arange(u.lo, u.hi + 1, dtype=float))
        td = sc.eval_modes(u.a, S)[3]
        d = min(float((u.frac * td).sum()), max(u.time - am.UNIT_OVERHEAD, 0.0))
        return u.time - d, d
    if u.dlo is not None:
        S = np.array([float(u.lo)])
    else:
        S = np.arange(u.lo, u.hi + 1, dtype=float)
    sq, m, tp, td, tc, dm, lNp, lL, cell = sc.eval_modes(u.a, S)
    if truth == "measured":
        nb = np.searchsorted(np.log(NBAND_EDGES), lNp, side="right")
        tpm = tp * np.array([TIME_TRUTH_PLAIN.get(int(b), 1.0) for b in nb])
        tdm = td * np.array([TIME_TRUTH_DFIRST.get(int(b), 1.0) for b in nb])
        Sg = sc.prof(u.a)[0]
        below = S < (Sg[0] - 1e-9) if len(Sg) else np.ones(len(S), bool)
        tpm[below] = tdm[below] = am.SUM_OVERHEAD
        if u.mode == "dfirst":
            return (float((tpm / u.calib).sum()) if u.calib else 0.0) + am.UNIT_OVERHEAD, \
                float((u.frac * tdm).sum())
        return float(tpm.sum()) + am.UNIT_OVERHEAD, 0.0
    if sc.policy == "off":
        c0, c1, c2, c3, c4, cb = sc._tcd
        k6 = int(sc.c.k[u.a])
        w3 = lL + am.LAB1 * lNp > sc._w3c - am.LAB2 * (k6 - 6)
        band = np.searchsorted(sc._tbe, lNp, side="right")
        hinge = np.maximum(0.0, lNp - LN8K)
        td = np.exp(c0 + c1 * lNp + c2 * lL + c3 * hinge + c4 * w3 + cb[band]) + am.SUM_OVERHEAD
    w = np.clip((lNp - math.log(3000)) / (math.log(5000) - math.log(3000)), 0.0, 1.0)
    tpa = np.exp((1 - w) * np.log(tp) + w * (np.log(td) - sc.log_ratio(lNp)))
    Sg = sc.prof(u.a)[0]
    below = S < (Sg[0] - 1e-9) if len(Sg) else np.ones(len(S), bool)
    tpa[below] = td[below] = am.SUM_OVERHEAD
    if u.mode == "dfirst":
        return (float((tpa / u.calib).sum()) if u.calib else 0.0) + am.UNIT_OVERHEAD, \
            float((u.frac * td).sum())
    return float(tpa.sum()) + am.UNIT_OVERHEAD, 0.0


def sp_kappa():
    """k_SP / (0.86 0.64): the model's (square, SP traversal) pairs per
    square = this x r_S r_P x 720 e^{lpSP} at f_rho = 1 (research/
    calibration-target.md 3.4, its cm.py; the SP coupling f_rho multiplies
    it, and P(magic | square) by f_rho^2)"""
    return math.sqrt(_am().KAPPA / 1.1) / (TRAV_BASE[0] * TRAV_BASE[1])


def unit_predictions(sc, u):
    """the model's predictions of the records of a planned unit u at the
    scorer's calibration and f_rho = 1: squares recorded (plain: its sums';
    d-first: its calibration stream's, 1 in k first rows), their S and P
    traversals (720 e^{lpS} x r_S, 720 e^{lpP} x r_P per square), their SP
    traversals before the class factors (720 e^{lpSP}, as the report's "SP
    traversals predicted"), their (square, SP traversal) pairs (sp_kappa()
    r_S r_P 720 e^{lpSP}, scripts/decide.py's model pairs), the pairs of a
    d-first unit's part of the d loop (found through the star cover; not
    in decide.py), and the unit's N' (of its first and last sum)"""
    np = _np()
    am = _am()
    a = u.a
    S = np.array([float(u.lo)]) if u.dlo is not None else np.arange(u.lo, u.hi + 1, dtype=float)
    sq, m, tp, td, tc, dm, lNp, lL, cell = sc.eval_modes(a, S)
    src = sc.c.src[a]
    if src >= 0:
        A = np.asarray(sc.store.arr[src], np.float64)
        v = np.asarray(sc.store.valid[src], bool)
    else:
        A, v = sc.store.get(sc.c.P(a))
        A = np.asarray(A, np.float64)
    T = NUM_TRAVERSALS[sc.n]
    if v.any():
        Sg = sc.c.S0[a] * (1 + am.GRID_U[v])
        eS, eP, eSP = (T * np.exp(np.interp(S, Sg, A[j, v])) for j in (3, 4, 5))
    else:
        eS = eP = eSP = np.zeros(len(S))
    rS, rP = np.exp(sc.ln_rS[cell]), np.exp(sc.ln_rP[cell])
    q = sp_kappa() * rS * rP * eSP
    if u.mode == "dfirst":
        rec = sq / u.calib if u.calib else np.zeros_like(sq)
        dpairs = float((u.frac * sq * q).sum())
    else:
        rec = sq
        dpairs = 0.0
    return {"pred_squares": float(rec.sum()), "pred_S": float((rec * eS * rS).sum()),
            "pred_P": float((rec * eP * rP).sum()), "pred_SP": float((rec * eSP).sum()),
            "pred_pairs": float((rec * q).sum()), "pred_pairs_dfirst": dpairs,
            "lNp": [float(lNp[0]), float(lNp[-1])]}


def forecast_truth(args):
    """the truth a forecast charges its units under: --truth, by default
    "measured" for --shipped (the shipped laws, corrected by the
    calibration search's measured CPU) and "laws" otherwise (the laws as
    learned from the state's records; their factors would correct a learned
    law twice where it has learned them. It has not learned them all: the
    online refit moves only the d-first law's level and one plain offset at
    N' >= 6k (amodel.DFIRST_TIME_LAMBDA, TIME_LAMBDA), so the measured
    per-band shape above 12k, 0.16-0.62 of the laws, is in no learned law)"""
    t = getattr(args, "truth", None)
    if t:
        return t
    return "measured" if getattr(args, "shipped", False) else "laws"


def load_machine(path):
    """(plain speed, d-first speed, cal) from scripts/machine_cal.py's
    cal.json: reference CPU-hours per instance-hour of the whole instance
    at its recommended worker count, for the plain and d-first search"""
    with open(path) as f:
        cal = json.load(f)
    ps, ds = cal.get("plain_speed"), cal.get("dfirst_speed")
    if not ps or not ds or ps <= 0 or ds <= 0:
        sys.exit(f"{path}: no plain_speed and dfirst_speed (run scripts/machine_cal.py)")
    return float(ps), float(ds), cal


def machine_desc(cal, ps, ds):
    c = cal.get("cpu", {})
    b = cal.get("build", {})
    return (f"{c.get('model') or cal.get('host', '?')}, path {b.get('path', '?')}, "
            f"{cal.get('workers', '?')} workers: plain {ps:.3g}, d-first {ds:.3g} reference "
            f"CPU-hours per instance-hour")


def v2_forecast(args):
    """simulate the planner from the current state, with posterior means
    (no optimism) and every unit taking its predicted time and finding its
    predicted squares; with --machine (scripts/machine_cal.py's cal.json)
    each unit also takes instance time, its plain part at the machine's
    plain speed and its d-first part at its d-first speed, and E and
    P(>=1) are reported at --instance-hours H"""
    np = _np()
    yr = 8766.0
    mach = getattr(args, "machine", None)
    ih = getattr(args, "instance_hours", None)
    if ih is not None and not mach:
        sys.exit("--instance-hours needs --machine cal.json")
    ps = ds = None
    if mach:
        ps, ds, cal = load_machine(mach)
    if getattr(args, "hours", None) is None:
        # (default: 1 CPU-year; with --instance-hours, through 10 CPU-years
        # and H at the faster of the two speeds)
        args.hours = yr if ih is None else max(10 * yr, ih * max(ps, ds) * 1.001)
    sch = SchedulerV2(args)
    if args.shipped:
        # the shipped calibration: prior class factors and time law, no per-P
        # factors (the state still sets the frontiers)
        sch.calib = Calibration()
        sch.tm = _am().time_prior(ENGINE)
        sch.tmd = _am().time_prior(ENGINE, mode="dfirst")
        sch.lr0 = _am().dfirst_ratio_level(None, ENGINE)
        sch.scorer.set_calibration(sch.calib)
        sch.scorer.tm = sch.tm
        sch.scorer.tmd = sch.tmd
        sch.scorer.lr0 = sch.lr0
        sch.scorer.lnfsq[:] = 0.0
        sch.scorer.lnFm[:] = 0.0
    frac = args.sample
    if frac < 1:
        # a uniform sample of the candidates with the budget scaled: the
        # same greedy, frac times smaller
        rng = np.random.default_rng(args.seed)
        drop = rng.random(len(sch.cands)) >= frac
        sch.scorer.lnfsq[drop] = -np.inf
    # (stage 1's hours scale with the budget)
    sch.stage1.scale = frac
    sch.scorer.stage1 = sch.stage1.band()
    if sch.stage1.spec:
        print(sch.stage1_status() + (f" (x{frac:g} with the sample)" if frac < 1 else ""))
    budget = args.hours * frac
    inst_budget = ih * frac if ih is not None else None
    marks = {h for h in (0.1 * yr, 0.3 * yr, yr, 3 * yr, 10 * yr, 30 * yr, 100 * yr)
             if h <= args.hours} | {args.hours}
    if args.hours < 0.1 * yr:
        marks |= {args.hours * f for f in (0.01, 0.1, 0.3)}
    marks = [m * frac for m in sorted(marks)]
    hours = squares = magic = 0.0
    touched = set()
    Ecell = np.zeros(NCELLS)
    comp = {"N' band": collections.Counter(), "k": collections.Counter(),
            "ratio": collections.Counter(), "S0": collections.Counter()}
    s0_edges = (850, 1200, 1700, 2400, 3400, 4800)
    s0_names = ("<850", "850-1.2k", "1.2-1.7k", "1.7-2.4k", "2.4-3.4k", "3.4-4.8k", ">=4.8k")
    t0 = time.time()
    last_u = None
    nunits = 0
    maxN = 0.0
    E_at = {}
    curve = [(0.0, 0.0)]
    dmagic = dhours = 0.0
    ndunits = 0
    truth = forecast_truth(args)
    inst = 0.0       # instance-hours (--machine)
    E_inst = None    # (E, CPU-hours) at --instance-hours
    if truth != "laws":
        print(f"(each unit charged under the {truth!r} truth, see --truth; the plan uses the "
              f"scheduler's laws)")
    if mach:
        print(f"machine {mach}: {machine_desc(cal, ps, ds)}")
    print(f"{'CPU-years':>10} {'squares':>12} {'magic squares':>14} {'P touched':>9} {'units':>9}"
          f" {'max N':>7}" + (f" {'inst-hours':>10} {'P(>=1)':>6}" if mach else "")
          + "  latest unit"
          + (f"   (sample {frac:g} of the candidates, scaled)" if frac < 1 else ""))

    def line(m, tail):
        E_at[m] = magic
        print(f"{m / frac / yr:10.3g} {squares / frac:12.4g} {magic / frac:14.3g} "
              f"{len(touched) / frac:9.0f} {nunits / frac:9.0f} {maxN:7.0f}"
              + (f" {inst / frac:10.4g} {1 - math.exp(-magic / frac):6.3f}" if mach else "")
              + f"  {tail}", flush=True)

    for u in sch.simulate(max_hours=budget if (truth == "laws" and inst_budget is None) else None,
                          detail=True):
        if hours >= budget and (inst_budget is None or inst >= inst_budget):
            break
        pl, df = unit_charge_parts(sch.scorer, u, truth)
        cu = pl + df
        hours += cu / 3600
        if mach:
            inst += (pl / ps + df / ds) / 3600
        squares += u.squares
        magic += u.magic
        nunits += 1
        curve.append((hours, magic))
        if u.mode == "dfirst":
            dmagic += u.magic
            dhours += cu / 3600
            ndunits += 1
        a = u.a
        touched.add(a)
        np.add.at(Ecell, u.cells, u.mvec)
        comp["k"][KBAND_NAMES[min(max(int(sch.cands.k[a]) - 5, 0), 3)]] += u.magic
        comp["ratio"][RBAND_NAMES[int(sch.cands.rbin[a])]] += u.magic
        comp["S0"][s0_names[int(np.searchsorted(s0_edges, sch.cands.S0[a]))]] += u.magic
        nb = u.cells // 36
        for b in np.unique(nb):
            comp["N' band"][NBAND_NAMES[b]] += float(u.mvec[nb == b].sum())
        if u.magic > 0:
            maxN = max(maxN, math.exp(sch.scorer.eval_sums(a, [u.hi])[3][0]))
        last_u = u
        if inst_budget is not None and E_inst is None and inst >= inst_budget:
            E_inst = (magic, hours)
        while marks and hours >= marks[0]:
            line(marks.pop(0), f"P={p_str(u.P)} S={u.lo}..{u.hi}{describe_dunit(u)}")
    marks_all = sorted(set(E_at) | set(marks))
    for m in marks:   # ran out of candidates
        line(m, "(no more units)")
    if inst_budget is not None:
        if E_inst is None:
            print(f"(ran out of units before {ih:g} instance-hours: E at {inst / frac:.4g})")
            E_inst = (magic, hours)
        Eh = E_inst[0] / frac
        print(f"on this machine: {ih:g} instance-hours = {E_inst[1] / frac / yr:.4g} reference "
              f"CPU-years (band-weighted along the plan: {E_inst[1] / max(inst_budget, 1e-12):.3g} "
              f"reference CPU-hours per instance-hour): E = {Eh:.3g}, P(>=1 magic square) = "
              f"{1 - math.exp(-Eh):.3f}")
    print("E at the marks (4 digits): " + ", ".join(
        f"{m / frac / yr:g} CPU-years {E_at[m] / frac:.4f}" for m in sorted(E_at)))
    if magic <= 0:
        return
    print(f"simulated in {time.time() - t0:.0f} s; marginal density at the end "
          f"{(last_u.score if last_u else 0) * 3.15e7:.3g} magic/CPU-year; max N' = the largest "
          f"N' of a unit with predicted magic squares")
    print("composition of E (share of the magic squares):")
    orders = {"N' band": NBAND_NAMES, "k": KBAND_NAMES, "ratio": RBAND_NAMES, "S0": s0_names}
    for name, cnt in comp.items():
        tot = sum(cnt.values())
        print(f"  {name:8} " + ", ".join(f"{k} {100 * cnt[k] / tot:.0f}%"
                                        for k in orders[name] if cnt[k]))
    sc = sch.scorer
    print(f"  mode     d-first (--dfirst {sc.policy}, N' >= {math.exp(sc.ldmin):.0f}, calibration "
          f"stream {100 * sc.calib_frac:.0f}%): {100 * dmagic / magic:.0f}% of E in "
          f"{100 * dhours / max(hours, 1e-12):.0f}% of the CPU ({ndunits / frac:.0f} units)")
    if sc.policy != "off":
        forecast_plain_only(sch, budget, frac, marks_all, E_at, curve, yr, truth)
    # uncertainty: draws of the class factors (Laplace posterior) times
    # lognormal(0, 0.2) for the pair factor and for SP+SP, and
    # lognormal(0, SELECTION_SD) for selection
    rng = np.random.default_rng(1)
    lg0, _, _, lm0 = sch.calib.rates()
    vals = []
    for lg, lm in sch.calib.draws(rng, args.draws):
        f = np.exp(lg - lg0 + lm - lm0)
        vals.append(float((Ecell * f).sum()) / frac
                    * math.exp(rng.normal(0, 0.2) + rng.normal(0, 0.2)
                               + rng.normal(0, SELECTION_SD)))
    E = magic / frac
    print(f"E({hours / frac / yr:.3g} CPU-years) = {E:.3g} at the posterior-mean log factors")
    if vals:
        q = np.percentile(vals, [5, 50, 95])
        print(f"  band from {len(vals)} draws (class-factor posterior x lognormal(0, 0.2) for PAIR "
              f"and for SP+SP x lognormal(0, {SELECTION_SD}) for selection; the plan is not "
              f"re-optimised per draw): 5% {q[0]:.3g}, median "
              f"{q[1]:.3g}, mean {np.mean(vals):.3g}, 95% {q[2]:.3g} (the draws' median is above "
              f"the point: the factors vary per cell and the sum of lognormals is skewed)")
    print("  (no selection discount: the calibration search found no winner's curse within a "
          "band; selection is in the band)")
    print(f"predicted CPU-years per magic square at this pace: {hours / magic / yr:.3g}")


def forecast_plain_only(sch, budget, frac, marks, E_at, curve, yr, truth="laws"):
    """the same simulation with the plain search only (--dfirst off): E at
    the marks with and without d-first, and the CPU d-first needs for the
    plain-only E at the end of the budget"""
    np = _np()
    sc = sch.scorer
    pol = sc.policy
    sc.policy = "off"
    try:
        hours = magic = 0.0
        todo = list(marks)
        E_off = {}
        for u in sch.simulate(max_hours=budget if truth == "laws" else None):
            if hours >= budget:
                break
            hours += unit_charge(sc, u, truth) / 3600
            magic += u.magic
            while todo and hours >= todo[0]:
                E_off[todo.pop(0)] = magic
        for m in todo:
            E_off[m] = magic
    finally:
        sc.policy = pol
    if truth == "laws":
        print("with and without d-first (the same greedy with --dfirst off), charged at the "
              "scheduler's own laws: this cannot measure the gain (the plain law under-predicts "
              "the plain search 1.35x at N' 5-8k and over-predicts it 1.5-3x above 11k); see "
              "--truth measured:")
    else:
        print(f"with and without d-first (the same greedy with --dfirst off), both charged under "
              f"the {truth!r} truth:")
    print(f"  {'CPU-years':>10} {'E d-first':>10} {'E plain':>10} {'ratio':>6}")
    for m in marks:
        e1, e0 = E_at.get(m, curve[-1][1]), E_off.get(m, magic)
        print(f"  {m / frac / yr:10.3g} {e1 / frac:10.3g} {e0 / frac:10.3g} "
              f"{e1 / max(e0, 1e-300):6.3f}")
    # CPU with d-first to reach the plain-only E of the whole budget
    h = np.array([c[0] for c in curve])
    e = np.array([c[1] for c in curve])
    target = E_off.get(marks[-1], magic) if marks else magic
    if target > 0 and e[-1] >= target:
        need = float(np.interp(target, e, h))
        print(f"  E = {target / frac:.3g} (plain only, {marks[-1] / frac / yr:.3g} CPU-years) "
              f"takes {need / frac / yr:.3g} CPU-years with d-first: "
              f"{marks[-1] / max(need, 1e-12):.2f}x less CPU")


def v2_fit(args):
    sch = SchedulerV2(args)
    sch.refit()
    print(describe_calib(sch.calib, sch.tm, sch.tmd, sch.lr0, sch.ratio_st))
    print(json.dumps(sch.calib.to_json()["beta"], indent=1))


def v2_report(args):
    sch = SchedulerV2(args, quiet=True)
    sch.save()
    report_v2(sch, top=args.top)


def v2_profile(args):
    """fill the profile store for the candidates; per-guard invalid counts"""
    np = _np()
    am = _am()
    args.profile_workers = args.workers
    t0 = time.time()
    sch = SchedulerV2(args)
    st = sch.store
    print(f"profiles for {len(sch.cands)} candidates ready ({time.time() - t0:.0f} s)")
    if st.valid is not None and len(st.valid):
        done = np.asarray(st.done, bool)
        V = np.asarray(st.valid, bool)[done]
        u3 = am.GRID_U >= 0.03
        bad = (~V[:, u3]).any(1)
        small = np.asarray(st.arr[:, 0, :], np.float32)[done] < math.log(2 * sch.n)
        ill = (~V & ~small)[:, u3].any(1)
        print(f"pool rows profiled: {done.sum()} of {len(done)}; grid points invalid: "
              f"{(~V).sum()} of {V.size}; P with an invalid point at u >= 0.03: {bad.sum()} "
              f"({100 * bad.mean():.2f}%), of which not merely N < {2 * sch.n} (ill-conditioned): "
              f"{ill.sum()} ({100 * ill.mean():.2f}%)")
        print("invalid grid points by guard (non-exclusive, counted when computed): " +
              ", ".join(f"{g} {st.guards.get(g, 0)}" for g in am.GUARDS))
        size = sum(os.path.getsize(os.path.join(st.pdir, x)) for x in os.listdir(st.pdir))
        print(f"store size {size / 1e6:.0f} MB")
    if st.extra:
        V = np.array([e["v"] for e in st.extra.values()], bool)
        print(f"extra P (outside the pool): {len(st.extra)}, grid points invalid "
              f"{(~V).sum()} of {V.size}")


def v2_pool(args):
    np = _np()
    t0 = time.time()
    pool = Pool.load(args, args.state)
    legacy = read_legacy(args.state, args.vec_size)
    idx = pool.index()
    removed = sum(1 for P, lg in legacy.items()
                  if p_str(P, "_") in idx and lg["maxS"] >= math.floor(2 * am_s0(P)))
    r = pool.ratio
    print(f"pool {pool.kind}: {len(pool)} values of P ({(r <= 1 + 1e-6).sum()} sorted, "
          f"{(r <= np.float32(1.1)).sum()} with ratio <= 1.1; {time.time() - t0:.1f} s); "
          f"{removed} removed as covered by the legacy search up to 2 S0")
    if not args.stats:
        return

    def hist(name, v, edges, fmt="{:g}"):
        b = np.searchsorted(edges, v, side="right")
        cnt = np.bincount(b, minlength=len(edges) + 1)
        labels = (["<" + fmt.format(edges[0])]
                  + [fmt.format(edges[i]) + "-" + fmt.format(edges[i + 1])
                     for i in range(len(edges) - 1)] + [">=" + fmt.format(edges[-1])])
        print(f"  {name:6} " + ", ".join(f"{l}: {c}" for l, c in zip(labels, cnt) if c))

    hist("S0", pool.S0, (600, 850, 1200, 1700, 2400, 3400, 4800))
    hist("tau", pool.tau, (2048, 4096, 8192, 16384, 32768))
    hist("k", pool.k, (5, 6, 7, 8))
    hist("ratio", pool.ratio, (1.0 + 1e-6, 1.05, 1.1, 1.15, 1.2), "{:.2f}")
    width = (pool.exps > 0).sum(0)
    print("  P using each prime: " + ", ".join(f"{p}: {w}" for p, w in zip(PRIMES, width)))


def v2_compact(args):
    """gzip finished unit files (with a done record, fully parsed)"""
    import gzip
    sch = SchedulerV2(args, quiet=True)
    n = saved = 0
    for path in glob.glob(os.path.join(sch.units, "*.jsonl")):
        st = sch.summary.files.get(os.path.basename(path))
        if not st or not st.get("done") or st["off"] != os.path.getsize(path):
            continue
        with open(path, "rb") as f, gzip.open(path + ".gz.tmp", "wb") as g:
            shutil.copyfileobj(f, g)
        os.replace(path + ".gz.tmp", path + ".gz")
        saved += os.path.getsize(path) - os.path.getsize(path + ".gz")
        os.remove(path)
        st["closed"] = 1
        n += 1
    sch.summary.save()
    print(f"compressed {n} unit files ({saved / 1e6:.1f} MB saved)")


def main():
    global UNIT_DROP, EXPLORE
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--state", default=os.path.join(ROOT, "data", "sched"))
    ap.add_argument("--vec-size", type=int, default=6)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def model_opts(p):
        p.add_argument("--model", choices=("analytic", "regression"), default=None,
                       help="analytic: the existence study's model over a wide pool (v2, "
                            "default for --vec-size 6); regression: the fitted Poisson "
                            "regression of v1 (default for other sizes)")
        p.add_argument("--pool", choices=("wide", "classic"), default=None,
                       help="candidate P: every exponent assignment within --pool-ratio of "
                            "its sorted one (wide, default with --model analytic) or "
                            "gen_candidates (classic, the only one for --model regression)")
        p.add_argument("--pool-ratio", type=float, default=1.2,
                       help="largest assignment ratio (P / P_sorted)^(1/6) in the wide pool")
        p.add_argument("--pool-s0-max", type=float, default=6000,
                       help="largest S0 = 6 P^(1/6) in the wide pool")
        p.add_argument("--pool-primes", type=int, default=10,
                       help="primes 2, 3, ... the wide pool may use")
        p.add_argument("--tau-min", type=int, default=None,
                       help="smallest number of divisors of candidate P "
                            "(default 1000 wide, 800 classic)")
        p.add_argument("--tau-max", type=int, default=12000,
                       help="largest number of divisors of candidate P (classic pool)")
        p.add_argument("--explore", type=float, default=None,
                       help="optimism of the per-P factors: posterior mean + EXPLORE x sd "
                            "(default 0 analytic, 1 regression)")
        p.add_argument("--profile-workers", type=int, default=2,
                       help="processes computing model profiles (niced)")
        p.add_argument("--no-legacy", action="store_true",
                       help="ignore imported legacy results")

    def common(p):
        model_opts(p)
        p.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1),
                       help="msearch processes at a time (default: cores - 1)")
        p.add_argument("--unit-time", type=float, default=120,
                       help="target CPU-seconds per unit of work")
        p.add_argument("--unit-drop", type=float, default=UNIT_DROP,
                       help="also end a unit where the predicted magic squares "
                            "per CPU-second fall below this fraction of its best")
        p.add_argument("--node-limit", type=int, default=20_000_000_000,
                       help="give up on a single sum after this many nodes (regression; "
                            "analytic: max(2e10, 10 x the predicted nodes of the unit's "
                            "largest sum))")
        p.add_argument("--active", type=int, default=300,
                       help="number of most promising P to consider at a time "
                            "(regression only)")
        p.add_argument("--only", help="restrict to these P, e.g. '13 6 3 2;12 7 4 2 1'")
        p.add_argument("--dfirst", choices=DFIRST_POLICIES, default="auto",
                       help="analytic: search a sum diagonal-first (msearch --diag-first, which "
                            "finds every magic square but not the other semi-magic squares) "
                            "where its predicted CPU, with its calibration stream, is below the "
                            "plain search's (auto), never (off), or always (on); at N' >= "
                            "--dfirst-min-n only")
        p.add_argument("--dfirst-min-n", type=float, default=DFIRST_MIN_NP,
                       help="analytic: smallest predicted N' of a d-first sum")
        p.add_argument("--stage1", nargs="?", const=":".join(f"{v:g}" for v in STAGE1_DEFAULT),
                       default=None, metavar="HOURS[:LO:HI]",
                       help="analytic: stage 1 (research/stage1.md): for the first HOURS of "
                            "reference CPU (the units' predicted CPU at the scheduler's laws, "
                            "counted in the state's stage1_<n>.json) search the sums with N' "
                            "in [LO, HI) plain, to learn the SP coupling f_rho from their "
                            "pairs (scripts/decide.py); default %(const)s; off = none")
        p.add_argument("--calib-frac", type=float, default=CALIB_FRAC,
                       help="analytic: CPU of a d-first sum's calibration stream (msearch "
                            "--calib-r1-stride, semi-magic squares for the models) as a share "
                            "of its predicted d-first CPU (0 = none)")

    def dispatch(v1, v2):
        return lambda a: (v2 if a.model == "analytic" else v1)(a)

    p = sub.add_parser("run")
    common(p)
    p.add_argument("--hours", type=float, default=0, help="stop after this long (0 = never)")
    p.add_argument("--machine", default=None,
                   help="analytic: this machine's cal.json (scripts/machine_cal.py): the CPU of "
                        "its records is learned as reference CPU (x its per-process speeds), "
                        "and --time-limit is scaled for a slower machine")
    p.add_argument("--refit-every", type=int, default=50, help="refit model every K units")
    p.add_argument("--announce-score", type=int, default=7,
                   help="log squares whose best diagonal pair scores at least this")
    p.set_defaults(func=dispatch(lambda a: Scheduler(a).run(), lambda a: SchedulerV2(a).run()))

    p = sub.add_parser("plan")
    common(p)
    p.add_argument("--top", type=int, default=30)
    p.set_defaults(func=dispatch(cmd_plan, v2_plan))

    p = sub.add_parser("emit")
    common(p)
    p.add_argument("--units", type=int, default=100)
    p.set_defaults(func=dispatch(cmd_emit, v2_emit))

    p = sub.add_parser("forecast")
    common(p)
    p.add_argument("--hours", type=float, default=None,
                   help="reference CPU-hours to simulate (default 8766; with --instance-hours, "
                        "through 10 CPU-years and H)")
    p.add_argument("--machine", default=None,
                   help="analytic: a machine's cal.json (scripts/machine_cal.py): also report "
                        "instance-hours (plain sums at its plain speed, d-first sums at its "
                        "d-first speed) and P(>=1)")
    p.add_argument("--instance-hours", type=float, default=None,
                   help="analytic, with --machine: E and P(>=1) after this many hours of the "
                        "instance")
    p.add_argument("--sample", type=float, default=1.0,
                   help="analytic: simulate on this uniform fraction of the candidates with "
                        "the budget scaled (the same greedy, faster)")
    p.add_argument("--seed", type=int, default=1, help="seed of --sample")
    p.add_argument("--draws", type=int, default=50,
                   help="analytic: posterior draws for the uncertainty band")
    p.add_argument("--shipped", action="store_true",
                   help="analytic: use the shipped calibration (no learned class, time or "
                        "per-P factors)")
    p.add_argument("--truth", choices=FORECAST_TRUTHS, default=None,
                   help="analytic: the CPU each planned unit is charged (the plan is always "
                        "made with the scheduler's laws): measured (the default with "
                        "--shipped) = the laws x the CPU the calibration search measured per "
                        "N' band; laws (the default otherwise) = its predicted time; anchored "
                        "= the d-first law for d-first sums and the d-first law / the measured "
                        "ratio r(N') for plain sums at N' >= 5k (the plain law at <= 3k, "
                        "blended between; it overcharges the plain sums at 3-6k 1.5x). On a "
                        "learned state (without --shipped) the laws learn only the d-first "
                        "law's level and one plain offset at N' >= 6k, so under 'laws' they "
                        "stay above the measured CPU at N' >= 12k (measured: 0.16-0.62 of the "
                        "laws there, on 7 sums); 'measured' there corrects the shipped laws, "
                        "not the learned ones")
    p.set_defaults(func=dispatch(cmd_forecast, v2_forecast))

    p = sub.add_parser("fit")
    model_opts(p)
    p.add_argument("--only", help=argparse.SUPPRESS)
    p.add_argument("--unit-time", type=float, default=120, help=argparse.SUPPRESS)
    p.add_argument("--unit-drop", type=float, default=UNIT_DROP, help=argparse.SUPPRESS)
    p.set_defaults(func=dispatch(cmd_fit, v2_fit))

    p = sub.add_parser("report")
    common(p)
    p.add_argument("--top", type=int, default=20)
    p.set_defaults(func=dispatch(cmd_report, v2_report))

    p = sub.add_parser("profile", help="analytic: compute the model profiles of the pool")
    common(p)
    p.set_defaults(func=v2_profile, workers=2)

    p = sub.add_parser("pool", help="build the candidate pool and print its size")
    model_opts(p)
    p.add_argument("--stats", action="store_true", help="sizes by S0, tau, k, ratio")
    p.set_defaults(func=v2_pool)

    p = sub.add_parser("compact", help="gzip finished unit files")
    common(p)
    p.set_defaults(func=v2_compact)

    p = sub.add_parser("ingest")
    p.add_argument("files", nargs="+")
    p.set_defaults(func=cmd_ingest)

    p = sub.add_parser("import-legacy")
    p.add_argument("stats")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--model", choices=("analytic", "regression"), default="analytic",
                   help="regression also computes vector counts over the legacy ranges "
                        "(bin/enumerate) for its fit; analytic needs none")
    p.set_defaults(func=cmd_import_legacy)

    args = ap.parse_args()
    if getattr(args, "model", "") is None:
        args.model = "analytic" if args.vec_size == 6 else "regression"
    UNIT_DROP = getattr(args, "unit_drop", UNIT_DROP)
    if hasattr(args, "pool"):
        if args.pool is None:
            args.pool = "wide" if getattr(args, "model", "analytic") == "analytic" else "classic"
        if getattr(args, "model", "analytic") == "regression" and args.pool != "classic":
            sys.exit("--model regression works on the classic pool only (--pool classic)")
        if args.tau_min is None:
            args.tau_min = 1000 if args.pool == "wide" else 800
        if getattr(args, "model", "") == "regression" and args.explore is not None:
            EXPLORE = args.explore
    args.func(args)


if __name__ == "__main__":
    main()
