#!/usr/bin/env python3
"""
Adaptive scheduler for the search over (P, S).

Instead of fixing P and sweeping S, this picks the (P, S) pairs expected to
produce the most *magic* squares per CPU-second, runs bin/msearch on them in
parallel, and keeps learning from the results.

    expected magic squares per second at (P, S)
        = (semi-magic squares per sum at (P, S)) / (CPU-seconds per sum)
        * P(a semi-magic square at (P, S) is magic)

The number of squares per sum is a Poisson regression on the number of
vectors N (piecewise linear in log N), how far S is above the smallest
possible sum S_min(P), and the divisor structure of P (number of divisors,
number of distinct primes, largest prime). It is fitted to our own runs (one
observation per sum) and to the per-P totals of the first search (one
observation per P, the sum over its range of sums). The time per sum is a
power law in N. The probability of being magic is estimated from the
traversals of the squares found: if p_S and p_P are the probabilities that a
traversal (a possible diagonal) has the magic sum / product, and rho corrects
for these not being independent, then a square has 5400 (unordered) pairs of
possible diagonals and

    P(magic) ~= 5400 * (rho * p_S * p_P)^2.

The fitted model says, roughly: p_S ~ 1/S, p_P ~ tau(P)^-0.7 (fewer primes is
better), rho ~ 2.1, CPU time per sum ~ N^5, and squares per CPU-second rise
with N up to N ~ 2000 and then fall, so the best sums are close to S_min.

All of these are multiplied by per-P factors (empirical Bayes: gamma priors
centered on the global model) learned from each P's own results, with an
optimistic bonus so new P's get explored. Each P is searched in increasing
order of S (the yield tends to decrease with S), so the scheduler keeps a
frontier S per P and repeatedly runs the next chunk of sums of the P whose
frontier has the highest score. The model is refit every few dozen units.

All state lives in a directory (default data/sched):
    units/*.jsonl    raw msearch output, one file per unit of work
    pinfo_6.json     cache of S_min(P) and the number of vectors per sum
    model_6.json     fitted model (written by `fit` and during `run`; until
                     then a built-in default model is used)
    legacy.json      per-P totals of the first search (`import-legacy`)
The state is rebuilt from units/ on every start, so runs can be interrupted
and resumed at any time; partially finished units keep their completed sums.

usage:
    scheduler.py import-legacy stats/stats_short.txt   (once; skips the sums
                                          the first search already did)
    scheduler.py run [--workers K] [--hours H] [--unit-time T]
    scheduler.py plan [--top K]           show the current ranking
    scheduler.py forecast [--hours H]     predicted squares and magic squares
                                          for the next H CPU-hours
    scheduler.py fit                      refit the model from all results
    scheduler.py report                   summary of results and best squares
    scheduler.py emit [--units K]         print msearch arguments for a static
                                          plan (e.g. SuperCloud, see
                                          submit-sc-plan.sh)
    scheduler.py ingest FILE...           add msearch outputs run elsewhere
Global options (before the command): --state DIR, --vec-size N.
Needs numpy for `fit` (and for refits during `run`).
"""
import argparse
import collections
import glob
import json
import math
import os
import shutil
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23)
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
    print(time.strftime("[%H:%M:%S] ") + msg, flush=True)


def binary(name):
    path = os.path.join(ROOT, "bin", name)
    if not os.path.exists(path):
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


class Results:
    """everything learned from msearch output files"""

    def __init__(self, n):
        self.n = n
        self.sums = {}       # (P, S) -> sum record
        self.squares = {}    # hash -> square record
        self.done_to = {}    # P -> all sums <= this have been searched
        self.legacy = {}     # P -> legacy aggregate (searched up to maxS)

    def add_file(self, path):
        try:
            f = open(path)
        except OSError:
            return
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
                if r["type"] == "sum":
                    self.sums[(P, r["S"])] = r
                    self.done_to[P] = max(self.done_to.get(P, 0), r["S"])
                elif r["type"] == "square":
                    r["P"] = P
                    self.squares[r["hash"]] = r
                elif r["type"] == "done":
                    self.done_to[P] = max(self.done_to.get(P, 0), r["last_sum"])

    def frontier(self, P, pinfo):
        """next sum to search for P"""
        f = max(self.done_to.get(P, 0), self.legacy.get(P, {}).get("maxS", 0)) + 1
        return max(f, pinfo.smin(P))

    def by_p(self):
        out = {}
        for (P, S), r in self.sums.items():
            out.setdefault(P, []).append(r)
        for v in out.values():
            v.sort(key=lambda r: r["S"])
        return out


# --------------------------------------------------------------------------
# model

# features are clamped to the range seen when fitting, so that the quadratic
# terms cannot extrapolate wildly
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
    """features for the rate of semi-magic squares at (P, S)"""
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

    # time per sum, from our own runs (it grows like N^5 in the range that
    # matters; small N are dominated by overheads and not worth modeling)
    Xt, yt = [], []
    for (P, S), r in results.sums.items():
        t = r["time"] + r["setup_time"] + r.get("enum_time", 0)
        if r["nvecs_raw"] >= min_n_time and not r["truncated"] and t > 1e-4:
            Xt.append(time_features(r["nvecs_raw"]))
            yt.append(math.log(t))
    if len(yt) < 20:
        raise SystemExit(f"not enough data to fit: {len(yt)} sums")
    lNs = [x[1] for x in Xt]
    time_range = (min(lNs), max(lNs))
    tcoef, tsd = fit_lsq(Xt, yt)
    tmodel = Model({"time": tcoef, "time_sd": tsd, "time_range": time_range})

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
    per_p = {}
    for (P, S), r in results.sums.items():
        a = per_p.setdefault(P, [0, 0.0])
        a[0] += r["squares"]
        a[1] += r["time"] + r["setup_time"]
    best_rate = max([sq / t for sq, t in per_p.values() if sq >= 20 and t > 0] or [1.0])
    # the legacy totals are dominated by large sums, far from where the
    # scheduler works, so they get a smaller weight than our own data
    wts = [1.0] * own_rows + [legacy_weight] * (len(y) - own_rows)
    sq_coef = fit_poisson(X, y, expo, groups, wts)
    model = Model({"n": n, "squares": sq_coef, "min_n": min_n, "time": tcoef,
                   "time_sd": tsd, "time_range": time_range, "rate_cap": 3 * best_rate})

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
    own_sq = len(ys)
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
                   "time_sd": tsd, "time_range": time_range, "rate_cap": 3 * best_rate,
                   "ps": ps, "pp": pp, "rho": rho,
                   "clamp": {k: list(v) for k, v in CLAMP.items()}})
    diag = {"own_sums": own_rows, "own_squares": int(own_squares),
            "own_cpu_seconds": own_time,
            "legacy_P": len(legacy_used), "legacy_squares": int(sum(y[own_rows:])),
            "traversals_S": int(sum(ys)), "traversals_P": int(sum(yp)),
            "traversals_SP": int(ysp + legacy_sp), "expected_SP_if_independent": expected_sp}
    return model, diag


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
        self.factors = {}
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

    def next_unit(self, P, unit_time, max_sums=2000):
        """the next chunk of sums of P (about unit_time CPU-seconds), with its
        predicted CPU time, semi-magic squares and magic squares; its score
        is the predicted magic squares per CPU-second"""
        lo = self.results.frontier(P, self.pinfo)
        counts = self.pinfo.counts(P, lo, lo + 200)
        total_t, total_sq, total_m, S = 0.0, 0.0, 0.0, lo
        hi = lo
        while True:
            if S not in counts:
                counts.update(self.pinfo.counts(P, S, S + 200))
            N = counts[S]
            if N >= 2 * self.n:
                t = self.model.sum_time(N)
                m, rate, _ = self.score(P, S, N)
                if total_t > 0 and total_t + t > unit_time:
                    break
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
        for path in sorted(glob.glob(os.path.join(self.dir, "units", "*.jsonl"))):
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

    def prepare(self, cands, scorer=None):
        """S_min for all candidates, then vector counts for the most promising
        ones; returns the active set of P"""
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
        active = sorted(set(active) | {P for (P, S) in self.results.sums if P in cands})
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
                            log(f"refit model on {diag['sums']} sums, {diag['squares']} squares")
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
                    proc = subprocess.Popen(self.command(u.P, u.lo, u.hi, path),
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
    for path in sorted(glob.glob(os.path.join(args.state, "units", "*.jsonl"))):
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


def simulate(sch, scorer, cands, unit_time, max_units=None, max_hours=None):
    """greedily simulate the scheduler, assuming every unit takes its
    predicted time and finds its predicted number of squares; yields units"""
    cache = {P: scorer.next_unit(P, unit_time) for P in cands}
    hours, i = 0.0, 0
    while (max_units is None or i < max_units) and (max_hours is None or hours < max_hours):
        u = max(cache.values(), key=lambda u: u.score)
        if u.score <= 0:
            break
        yield u
        hours += u.time / 3600
        i += 1
        sch.results.done_to[u.P] = u.hi
        cache[u.P] = scorer.next_unit(u.P, unit_time)


def cmd_emit(args):
    """static plan: print msearch arguments for the units the scheduler would
    run (one per line), simulating it with predicted times and squares"""
    sch = Scheduler(args)
    scorer = Scorer(sch.model, sch.results, sch.pinfo)
    cands = sch.prepare(sch.candidates(), scorer)
    units_dir = os.path.join(args.state, "units")
    if os.path.abspath(units_dir).startswith(ROOT + os.sep):
        units_dir = os.path.relpath(units_dir, ROOT)  # portable plan
    for i, u in enumerate(simulate(sch, scorer, cands, args.unit_time, max_units=args.units)):
        out = os.path.join(units_dir, f"plan{i:06d}_{p_str(u.P, '_')}_{u.lo}_{u.hi}.jsonl")
        print(" ".join(sch.command(u.P, u.lo, u.hi, out)[1:]))


def cmd_forecast(args):
    """predicted squares and magic squares for the next CPU-hours of search"""
    sch = Scheduler(args)
    cands = sch.prepare(sch.candidates(), Scorer(sch.model, sch.results, sch.pinfo))
    # the scheduler's choices use the exploration bonus, but the predicted
    # numbers of squares should not
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


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--state", default=os.path.join(ROOT, "data", "sched"))
    ap.add_argument("--vec-size", type=int, default=6)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
        p.add_argument("--unit-time", type=float, default=120,
                       help="target CPU-seconds per unit of work")
        p.add_argument("--node-limit", type=int, default=20_000_000_000,
                       help="give up on a single sum after this many nodes")
        p.add_argument("--tau-min", type=int, default=800,
                       help="smallest number of divisors of candidate P")
        p.add_argument("--tau-max", type=int, default=12000,
                       help="largest number of divisors of candidate P")
        p.add_argument("--active", type=int, default=300,
                       help="number of most promising P to consider at a time")
        p.add_argument("--only", help="restrict to these P, e.g. '13 6 3 2;12 7 4 2 1'")
        p.add_argument("--no-legacy", action="store_true",
                       help="ignore imported legacy results")

    p = sub.add_parser("run")
    common(p)
    p.add_argument("--hours", type=float, default=0, help="stop after this long (0 = never)")
    p.add_argument("--refit-every", type=int, default=50, help="refit model every K units")
    p.add_argument("--announce-score", type=int, default=7,
                   help="log squares whose best diagonal pair scores at least this")
    p.set_defaults(func=lambda a: Scheduler(a).run())

    p = sub.add_parser("plan")
    common(p)
    p.add_argument("--top", type=int, default=30)
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("emit")
    common(p)
    p.add_argument("--units", type=int, default=100)
    p.set_defaults(func=cmd_emit)

    p = sub.add_parser("forecast")
    common(p)
    p.add_argument("--hours", type=float, default=8766, help="CPU-hours to simulate")
    p.set_defaults(func=cmd_forecast)

    p = sub.add_parser("fit")
    p.add_argument("--no-legacy", action="store_true", help="ignore imported legacy results")
    p.set_defaults(func=cmd_fit)

    p = sub.add_parser("report")
    common(p)
    p.add_argument("--top", type=int, default=20)
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("ingest")
    p.add_argument("files", nargs="+")
    p.set_defaults(func=cmd_ingest)

    p = sub.add_parser("import-legacy")
    p.add_argument("stats")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.set_defaults(func=cmd_import_legacy)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
