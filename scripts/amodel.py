#!/usr/bin/env python3
"""
The analytic model of the existence study (research/existence/analytic.md,
research/existence.md sections 2-4) for scheduler v2: for one P and an array
of sums S it predicts

    N(P, S)        number of vectors (n-sets of distinct divisors of P with
                   sum S and product P), before msearch's reduction
    E_semi(P, S)   expected number of semi-magic squares with sum S
    P(magic | semi-magic square at (P, S))
    p_S, p_P, p_SP probabilities that a traversal (a possible diagonal) has
                   the magic sum / product / both
    L_raw          expected number of distinct divisors among n N draws (what
                   msearch's "labels" count, up to a fitted factor)

How. Put on the divisors d of P the tilted (max-entropy) measure
mu(d) ~ exp(-beta d / (S/n) + sum_p lambda_p e_p(d)) with E[d] = S/n and
E[e_p(d)] = a_p/n (Newton on the convex dual, vectorised over S). With H its
entropy and C the covariance of (d, e_2(d), e_3(d), ...), a local central
limit theorem with exact per-prime lattice factors (amodel_lattice.json,
from the exact exponent-matrix counts H_n, M_n) gives

    N      ~ e^{nH} (2 pi)^{-q/2} det(nC)^{-1/2} L_vec P_dist(n) / n!
    E_semi ~ e^{n^2 H} / [(2 pi)^{(2n-1)q/2} det C^{(2n-1)/2} n^{(n-1)q}]
             L_semi P_dist(n^2) / (2 (n!)^2) x K_n
    P(magic | semi) = D_n kappa p_pair,  p_pair = L_M / [(2 pi)^q det C ((n-1)^2 - c^2)^{q/2}]

where q = 1 + (number of distinct primes), P_dist(m) is the probability that
m draws from mu are distinct (saddle point), D_n = 5400 pairs of diagonals
for n = 6, and the calibration (fitted on n = 5, 6, 7 and tested on held-out
data up to N = 45k: squares obs/pred 1.04 on the first search, 0.78 at N
3.4-8.9k, 1.23 at N 15-45k) is

    log K_n = C0 + B log P_dist(n^2) + DN (n - 6),    kappa = KAPPA.

The review corrections of existence.md 4.1 (REVIEW below) and the bias of
N above ~2.5k (n_bias) are kept separate: profile() returns the model
itself, and the scheduler applies (and learns) the factors on top.

Grid. The scheduler evaluates each P once on S_j = S0 (1 + u_j), u_j in
GRID_U (24 points from 0.01 to 1, S0 = n P^{1/n} the AM-GM bound) and
interpolates in S. Points that are numerically meaningless (very close to
S_min the dual is ill-conditioned and N, E come out astronomically large)
are marked invalid: non-finite values, Newton not converged, N < 2n, or at
or before the last decrease of N along the grid.

Time. TimeModel: CPU seconds per sum as a function of the model's N (bias
corrected), L_raw and the label words, ln t = th . time_features (see
TIME_PRIOR), refit online per msearch engine (Bayesian ridge that keeps the
shape and learns an intercept, the label-word step and N' band offsets);
the d-first search has its own law (DFIRST_TIME_PRIOR, only its level
learned) and a cost profile along d (DFIRST_COST_PROFILE).
node_law(): predicted search nodes, used only for --node-limit.

Needs numpy. AMODEL_VERSION changes whenever the numbers profile() returns
or the meaning of stored statistics change (it invalidates cached
profiles and the scheduler's summary); model_hash() additionally covers
the constants and the lattice table, in case a change forgets the bump.
"""
import json
import math
import os

import numpy as np

AMODEL_VERSION = 2     # the scheduler's summary statistics (time features: 2)
PROFILE_VERSION = 1    # the numbers profile() returns (the profile store)

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53)
HERE = os.path.dirname(os.path.abspath(__file__))

# calibration of the analytic model (analytic.md 4.2-4.3)
C0, B, DN, KAPPA = -0.597, 0.1365, -2.455, 0.57
D_6 = 5400

# review corrections (existence.md 4.1): multipliers on squares (SQ12) or on
# P(magic | square) (the others); the scheduler's learning priors start here
SQ12 = 1.23     # squares per sum at N' >= 12000
NTR = 1.35      # P(magic) at N' >= 3000 (the review's value; the scheduler's prior is 1.22)
K7 = 0.85       # P(magic) at k >= 7 distinct primes
X01 = 0.92      # P(magic) at x = ln(S / smin_approx) < 0.1
PAIR = 0.87     # P(magic): the pair factor
REVIEW = dict(SQ12=SQ12, NTR=NTR, K7=K7, X01=X01, PAIR=PAIR)

# sanity cap on P(magic | square)
LOG_PM_CAP = math.log(1e-6)

# observed msearch labels / L_raw = exp(LAB0 + LAB1 ln(N'/4000) + LAB2 (k - 6))
# (sd 0.032 over 646 engine-2 sums, 10 of them at N 6.7-45k)
LAB0, LAB1, LAB2 = -0.391, -0.035, -0.016

# the model's N runs 4-6% high above 3-4k: N' = N exp(n_bias(ln N))
NBIAS_SLOPE, NBIAS_KNOT = -0.045, math.log(2500)

GRID_U = np.geomspace(0.01, 1.0, 24)
FIELDS = ("lN", "lEs", "lPm", "lpS", "lpP", "lpSP", "lLraw", "ld")

# divisors above CAP_FACTOR * max(S) / n are left out of the measure (their
# weight is negligible: |dlog E| p90 6e-5 against 80 max(S) / n on 115 P)
CAP_FACTOR = 24.0
A_MAX = 40      # the lattice table covers exponents 0..40

with open(os.path.join(HERE, "amodel_lattice.json"), "rb") as _f:
    _LAT_BYTES = _f.read()
_LAT = {int(n): v for n, v in json.loads(_LAT_BYTES).items()}


def _hash(*parts):
    import hashlib
    h = hashlib.sha1()
    for x in parts:
        h.update(x if isinstance(x, bytes) else json.dumps(x, sort_keys=True, default=float).encode())
    return h.hexdigest()[:12]


def profile_hash():
    """hash of everything profile() and the store depend on, besides the
    code: PROFILE_VERSION, the calibration constants, the grid and the
    lattice table (a cached profile is reused only when it matches)"""
    return _hash(PROFILE_VERSION, [C0, B, DN, KAPPA, D_6, CAP_FACTOR, A_MAX, LOG_PM_CAP, STORE_CLIP],
                 [float(u) for u in GRID_U], list(FIELDS), _LAT_BYTES)


def model_hash():
    """profile_hash plus the constants the scheduler's stored statistics
    depend on (review factors, label and N-bias laws, time features)"""
    return _hash(AMODEL_VERSION, profile_hash(), REVIEW, [LAB0, LAB1, LAB2, NBIAS_SLOPE, NBIAS_KNOT],
                 list(TIME_BAND_EDGES), list(TIME_FEATURES))


# --------------------------------------------------------------------------
# helpers on P (an exponent tuple over PRIMES)


def tau(P):
    t = 1
    for a in P:
        t *= a + 1
    return t


def log_p(P):
    return sum(a * math.log(p) for p, a in zip(PRIMES, P))


def s0(P, n=6):
    """AM-GM bound n P^(1/n): no vector has a smaller sum"""
    return n * math.exp(log_p(P) / n)


def smin_approx(P, n=6):
    """S_min ~ S0 (1 + 10/tau): median error 0.03%, 90% within 0.17%, max
    1.2% (5,492 P with exact S_min)"""
    return s0(P, n) * (1 + 10.0 / tau(P))


def ratio(P, n=6):
    """assignment ratio (P / P_sorted)^(1/n) >= 1, where P_sorted puts the
    same exponents in non-increasing order on 2, 3, 5, ..."""
    srt = sorted((a for a in P if a), reverse=True)
    return math.exp((log_p(P) - log_p(srt)) / n)


def num_primes(P):
    return sum(1 for a in P if a)


def ndiag_pairs(n):
    """unordered pairs of possible diagonals of a semi-magic square: 5400 for n = 6"""
    f = math.factorial(n)
    if n % 2 == 0:
        inv = 1
        for j in range(n - 1, 0, -2):
            inv *= j
    else:
        inv = n
        for j in range(n - 2, 0, -2):
            inv *= j
    return f * inv // 2


def lattice(kind, n, a):
    if a > A_MAX:
        raise ValueError(f"exponent {a} > {A_MAX}: outside the lattice table")
    return _LAT[n][kind][a]


def n_bias(lN):
    """log correction of the model's N: dN = -0.045 max(0, lN - ln 2500)"""
    return NBIAS_SLOPE * np.maximum(0.0, np.asarray(lN, float) - NBIAS_KNOT)


def labels_obs(lNp, lLraw, k):
    """predicted msearch 'labels' (distinct entries after reduction) from the
    bias-corrected lN' and L_raw"""
    return np.exp(np.asarray(lLraw, float) + LAB0 + LAB1 * (np.asarray(lNp, float) - math.log(4000))
                  + LAB2 * (k - 6))


# --------------------------------------------------------------------------
# the model


def divisors_capped(P, cap):
    """divisors d <= cap of P, their exponent matrix, and the exponents a_p
    of the primes present"""
    idx = [i for i, a in enumerate(P) if a > 0]
    ds = np.array([1.0])
    E = np.zeros((1, 0))
    for i in idx:
        p, a = PRIMES[i], P[i]
        pw = np.array([float(p) ** j for j in range(a + 1)])
        nd = ds[:, None] * pw[None, :]
        rows, cols = np.nonzero(nd <= cap)
        E = np.concatenate([E[rows], cols[:, None].astype(float)], axis=1)
        ds = nd[rows, cols]
    return ds, E, np.array([P[i] for i in idx], dtype=float)


def log_pdist(W, m, iters=60):
    """log P(m iid draws from each row of W are distinct) = log(m! e_m(w)),
    by a saddle point (existence/analytic esym2.py)"""
    R = W.shape[0]
    lo = np.full(R, -60.0)
    hi = np.full(R, 60.0)
    t = np.full(R, math.log(m))
    for _ in range(iters):
        wz = W * np.exp(t)[:, None]
        f = (wz / (1 + wz)).sum(1) - m
        fp = (wz / (1 + wz) ** 2).sum(1)
        lo = np.where(f < 0, t, lo)
        hi = np.where(f >= 0, t, hi)
        tn = t - f / np.maximum(fp, 1e-300)
        if np.abs(tn - t).max() < 1e-9:
            t = tn
            break
        bad = (tn < lo) | (tn > hi) | ~np.isfinite(tn)
        t = np.where(bad, 0.5 * (lo + hi), tn)
    wz = W * np.exp(t)[:, None]
    V = (wz / (1 + wz) ** 2).sum(1)
    return np.log1p(wz).sum(1) - m * t - 0.5 * np.log(2 * math.pi * V) + math.lgamma(m + 1)


def _maxent(X, T, S, n, k):
    """Newton (Armijo backtracking) on the dual of the max-entropy problem,
    one row per sum; returns (log partition - theta.T, weights, converged)"""
    m = len(S)
    q = 1 + k
    XX = (X[:, :, None] * X[:, None, :]).reshape(len(X), q * q)
    th = np.zeros((m, q))
    th[:, 0] = n / S

    def evalf(th):
        Z = th @ X.T
        zm = Z.max(1, keepdims=True)
        W = np.exp(Z - zm)
        Zs = W.sum(1, keepdims=True)
        W /= Zs
        return (np.log(Zs[:, 0]) + zm[:, 0]) - (th * T).sum(1), W

    val, W = evalf(th)
    scale = np.concatenate([(S / n)[:, None], np.ones((m, k))], 1)
    conv = np.zeros(m, bool)
    for _ in range(60):
        M = W @ X
        g = M - T
        conv = np.abs(g / scale).max(1) < 1e-10
        if conv.all():
            break
        Hs = (W @ XX).reshape(m, q, q) - M[:, :, None] * M[:, None, :]
        step = np.linalg.solve(Hs + 1e-14 * np.eye(q)[None], g[:, :, None])[:, :, 0]
        step[conv] = 0.0
        t = np.ones(m)
        for _ in range(30):
            nv, nW = evalf(th - t[:, None] * step)
            bad = (nv > val - 1e-4 * t * (g * step).sum(1) + 1e-12 * np.abs(val)) & ~conv
            if not bad.any():
                break
            t = np.where(bad, t * 0.5, t)
        th = th - t[:, None] * step
        val, W = nv, nW
    else:
        g = W @ X - T
        conv = np.abs(g / scale).max(1) < 1e-10
    return val, W, XX, conv


def profile(P, S, n=6, cap=None):
    """the calibrated analytic model of one P at the sums S (1-d array,
    ascending for the monotonicity guard); returns a dict of arrays over S:
    lN, lEs, lPm, lpS, lpP, lpSP, lLraw, ld (= log P_dist(n^2)), converged,
    valid, and the number of distinct primes k. No review factors applied."""
    P = tuple(int(a) for a in P)
    while P and P[-1] == 0:
        P = P[:-1]
    if any(a > A_MAX for a in P):
        raise ValueError(f"exponent > {A_MAX} in {P}")
    S = np.atleast_1d(np.asarray(S, dtype=float))
    if cap is None:
        cap = CAP_FACTOR * S.max() / n
    ds, E, a = divisors_capped(P, cap)
    k = len(a)
    q = 1 + k
    X = np.concatenate([-ds[:, None], E], axis=1)
    T = np.concatenate([-(S / n)[:, None], np.tile(a / n, (len(S), 1))], axis=1)
    with np.errstate(all="ignore"):
        H, W, XX, conv = _maxent(X, T, S, n, k)
        m = len(S)
        M = W @ X
        C = (W @ XX).reshape(m, q, q) - M[:, :, None] * M[:, None, :]
        ldC = np.linalg.slogdet(C)[1]
        ldCe = np.linalg.slogdet(C[:, 1:, 1:])[1] if k > 0 else np.zeros(m)
        L2P = math.log(2 * math.pi)
        loc = {kk: sum(lattice(kk, n, int(x)) for x in a) for kk in ("vec", "semi", "T", "M")}
        tup = n * H - 0.5 * q * L2P - 0.5 * (q * math.log(n) + ldC) + loc["vec"]
        grid = (n * n * H - 0.5 * (2 * n - 1) * q * L2P - 0.5 * (2 * n - 1) * ldC
                - (n - 1) * q * math.log(n) + loc["semi"])
        lpS = -0.5 * (math.log(2 * math.pi * (n - 1)) + np.log(C[:, 0, 0]))
        lpP = -0.5 * (k * math.log(2 * math.pi * (n - 1)) + ldCe) + loc["T"]
        lpSP = -0.5 * (q * math.log(2 * math.pi * (n - 1)) + ldC) + loc["T"]
        cdd = -1.0 if n % 2 == 0 else 0.0
        lpair = -(q * L2P + ldC + 0.5 * q * math.log((n - 1) ** 2 - cdd ** 2)) + loc["M"]
        sel = W.max(0) > 1e-16
        Wsel = W[:, sel]
        ld_n = log_pdist(Wsel, n)
        ld_n2 = log_pdist(Wsel, n * n)
        lN = tup + ld_n - math.lgamma(n + 1)
        lK = C0 + B * ld_n2 + DN * (n - 6)
        lEs = grid + ld_n2 - math.log(2) - 2 * math.lgamma(n + 1) + lK
        lPm = math.log(ndiag_pairs(n)) + math.log(KAPPA) + lpair
        Lraw = (1.0 - np.exp(-n * np.exp(lN)[:, None] * Wsel)).sum(1)
        lLraw = np.log(np.maximum(Lraw, 1.0))
    out = dict(lN=lN, lEs=lEs, lPm=np.minimum(lPm, LOG_PM_CAP), lpS=lpS, lpP=lpP, lpSP=lpSP,
               lLraw=lLraw, ld=ld_n2)
    finite = np.ones(len(S), bool)
    for key in FIELDS:
        finite &= np.isfinite(out[key])
    valid = finite & conv & (lN >= math.log(2 * n))
    # ill-conditioned points near S_min: lN falls along the grid there
    lNf = np.where(finite, lN, 1e300)
    dec = np.nonzero(np.diff(lNf) < -1e-7)[0]
    if len(dec):
        valid[:dec[-1] + 2] = False
    out.update(converged=conv, valid=valid, k=k)
    return out


def grid_sums(P, n=6, u=GRID_U):
    return s0(P, n) * (1 + np.asarray(u))


GUARDS = ("nonfinite", "not_converged", "small_N", "invalid")
STORE_CLIP = 60.0   # stored values are clipped to +-60 (e^-60 is nothing)


def profile_rows(Ps, n=6):
    """profiles of several P on the grid, as stored by the scheduler:
    float16 [m, len(FIELDS), len(GRID_U)], valid uint8 [m, len(GRID_U)],
    guard counts int32 [m, len(GUARDS)] (ProcessPoolExecutor worker)"""
    m = len(Ps)
    arr = np.zeros((m, len(FIELDS), len(GRID_U)), np.float16)
    val = np.zeros((m, len(GRID_U)), np.uint8)
    gc = np.zeros((m, len(GUARDS)), np.int32)
    for i, P in enumerate(Ps):
        try:
            pr = profile(P, grid_sums(P, n), n)
        except (ValueError, IndexError):   # exponent or prime outside the tables
            gc[i] = [len(GRID_U), 0, 0, len(GRID_U)]
            continue
        for j, key in enumerate(FIELDS):
            arr[i, j] = np.clip(np.nan_to_num(pr[key], nan=0.0, posinf=STORE_CLIP, neginf=-STORE_CLIP),
                                -STORE_CLIP, STORE_CLIP)
        val[i] = pr["valid"]
        g = guard_counts(pr)
        gc[i] = [g[key] for key in GUARDS]
    return arr, val, gc


def worker_init():
    """ProcessPoolExecutor initializer: lowest priority, one BLAS thread"""
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = "1"
    try:
        os.nice(19)
    except OSError:
        pass


def guard_counts(prof):
    """number of grid points failing each guard (for `scheduler.py profile`)"""
    finite = np.ones(len(prof["lN"]), bool)
    for key in FIELDS:
        finite &= np.isfinite(prof[key])
    return {"nonfinite": int((~finite).sum()), "not_converged": int((~prof["converged"]).sum()),
            "small_N": int((finite & (prof["lN"] < math.log(12))).sum()),
            "invalid": int((~prof["valid"]).sum())}


# --------------------------------------------------------------------------
# time and nodes per sum

UNIT_OVERHEAD = 0.05    # CPU seconds per unit (process start)
SUM_OVERHEAD = 0.002    # CPU seconds per sum

# ln t = th . time_features(lN', L_raw, k): [1, ln(N'/4000), ln(L_raw/150),
# max(0, ln(N'/8000)), [labels_obs > 128], band offsets for N' <1k, 1-3k,
# 3-6k, >=6k]. msearch engine 2, fitted (least squares on ln t, AMODEL_VERSION
# 2) with the model's N' and L_raw at each sum's S (as the scheduler uses
# them) on 3,690 full sums at N' >= 300 (pool/validate_rows, the live
# verification units, the 32 CPU-timed held-out draws) and 10 r1-sampled
# sums at N 6.7-45k (weight 30); the band offsets are 0 in the prior and
# only learned online. The step at labels > 128 is the third 64-bit word of
# the label bitsets (cost per node 0.25 -> 0.6 us); without it the law ran
# 1.2-1.3x under the sums at N' 3-6k that carry most of the forecast's E.
# Held out (fit without the draws and the live B/C units): draws CPU obs/pred
# 1.10, live C 1.03. Residual sd 0.21-0.33 by source, 0.48 on the sampled
# sums; the shipped sd 0.30 (its e^{sd^2/2} enters the mean) and the
# intercept (+0.03 over the ln-t fit) make sum(t)/sum(pred) ~1.04 at
# N' >= 3k and ~0.9 below.
TIME_BAND_EDGES = (1000.0, 3000.0, 6000.0)
TIME_FEATURES = ("intercept", "ln N'/4000", "ln L/150", "hinge N'>8k", "labels>128",
                 "band<1k", "band 1-3k", "band 3-6k", "band>=6k")
NT = len(TIME_FEATURES)
TIME_PRIOR = {"th": [3.8587, 5.4575, -3.4622, 0.9458, 0.3934, 0.0, 0.0, 0.0, 0.0],
              "sd": 0.30, "engine": 2}
# the start of a newer msearch engine's law: the previous engine's
# posterior (the shipped prior for engine 2) plus these shifts of th (see
# time_prior). Engine 3 (cx/integrated: per-r1 widths, carried bitsets up to
# 512 labels, the pretest) ran the plain search of the same sums at 0.98-1.01x
# engine 2's CPU with <= 128 labels and 0.66-0.89x (geometric mean 0.80,
# 8 sums at 137-252 labels, N 4.1-23k, paired) above, so its labels > 128
# step starts 0.227 lower (0.39 -> 0.17); above 256 labels it ran at
# 0.26x, which no feature of the law describes (3 sums at N 21-32k; refit
# the shape when engine-3 data arrives there; research/ideas.md).
# Engine 4 (round 2: the class support in the V_d searches and the star
# cover of the d loop) changes only the d-first search: the plain search's
# nodes and hashes are unchanged and bench prod's CPU is within noise, so
# its plain law starts at engine 3's posterior unshifted.
ENGINE_TIME_SHIFT = {3: {"labels>128": -0.227}}
# prior precision of th for the online refit. The shape (N, L and hinge
# slopes) is fixed: refitting it on the thousands of cheap sums at N' < 3k
# moved predictions at N' >= 3k by 1.4-2x. Learned online: a global
# intercept (sd 0.1), the W >= 3 step (sd 0.1) and one offset per N' band
# (sd 0.25), so a residual learned at small N' moves the other bands by
# about 14% of itself. TIME_SD_PSEUDO pseudo-sums at the prior sd anchor
# the residual variance.
TIME_LAMBDA = (100.0, 1e6, 1e6, 1e6, 100.0, 16.0, 16.0, 16.0, 16.0)
TIME_SD_PSEUDO = 20.0

# The d-first search (msearch --diag-first; one semi-magic search per vector
# d of the sum, on V_d = {v : |v & d| = 1}): ln t = th . time_features with
# only the intercept and the N' slope, least squares on ln t over the perf
# verifier's 11 sums at N 4.1-31.7k (process CPU per sum on the integrated
# build, engine 3; research/ideas.md, "Measurements on the integrated
# binary"), with the model's N' at each sum as the scheduler computes it:
# 3.215 (se 0.15) + 3.576 (se 0.11) ln(N'/4000), residual sd 0.21
# (leave-one-out 0.24), sum(t) / sum(pred) 1.03. The three sums with more
# than 256 labels sit on the line (residuals -0.29, +0.18, +0.07; a step
# for them fits -0.03 +- 0.19), so the law has none; adding ln(L/150) or
# the 8k hinge did not help (slopes -1.9 +- 1.0, 0.3 +- 0.5). Below N 3k
# (not in the fit) it runs above the measured sums (1.1-1.6x at N 2-3k),
# where the plain search is chosen anyway. Online only the level is learned
# (DFIRST_TIME_LAMBDA: intercept sd 0.2, the rest held; a part of a split
# sum enters with weight its share of the d loop, so a sum weighs 1).
# Back-test (engine 3, 16 distinct sums, mean prediction with
# e^{sd^2/2}): obs/pred geometric mean 0.94, sd(ln) 0.23; by N': 0.71 below
# 3k (outside the fit), 1.16 at 3-6k, 0.93 at 6-12k, 0.96 at 12-24k, 1.10
# above 24k; 1.01 at 129-256 labels, 0.96 above. Out of sample: 13 7 4 3 1
# 1 / 1950 1.24, / 2100 1.08, 12 6 3 2 1 1 / 950 0.96, six pool sums at N'
# ~6.1k 0.68-1.27 (geometric mean 0.89); the prototype's engine-2 rows
# 1.01 (sd 0.32). The scatter is mostly per P (+-25%): 13 7 4 3 1 1, 9 of
# the 16 sums, runs ~1.05-1.1x, 14 7 4 4 1 0 0 1, 9 6 4 3 1 1 1 1 and 12 6
# 3 2 1 0 1 0.6-0.75x; a label or x term on the residuals is not identified
# (slopes -0.69 +- 1.0 and -0.69 +- 0.36). If a label term is ever added,
# add the same to both laws, so that the mode choice stays label-free.
DFIRST_TIME_PRIOR = {"th": [3.215, 3.576, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                     "sd": 0.25, "engine": 3}
DFIRST_TIME_LAMBDA = (25.0, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6, 1e6)
# per-engine shifts of the d-first law after DFIRST_TIME_PRIOR's engine.
# Engine 4 (the class support in the V_d searches, gated to the V_d of >= 137
# labels, and the star cover with K = 4): ln(engine 4 / engine 3 d-first
# CPU per sum) = c0 + c1 ln(N'/4000), DFIRST_E4 = (-0.103 +- 0.023, -0.204
# +- 0.020), resid sd 0.067: the perf verifier's pooled fit of 23 paired
# sums inside the model grid (S <= S_grid max), N' 3.4-31.7k: its 13 (the
# 189b49f binary against the integrated one, identical d samples, whole-sum
# CPU = d_stride x the d loop + the fixed costs, min of 2 alternating) and
# the integrator's 10 in-grid sums (one binary, --no-class-support
# --dfirst-star 0 against the defaults). Alone, the verifier's 13 give
# -0.105 (0.032) - 0.210 (0.028); the class support alone -0.022 - 0.215;
# the integrator's first fit, -0.138 - 0.162, rested on two sums outside the
# grid (10 6 3 2 1 at S 900 and 1360, 154-194 labels at N 3.1-3.6k, where
# the class support starts at a smaller N than near S_min) and on run
# times with the whole x* choice but 1/d_stride of the loop (research/
# ideas.md, "Integration of the round-2 d-first changes"). Engine 4's law:
# 3.112 + 3.372 ln(N'/4000) (sd 0.25 kept); measured on the verifier's 13
# sums 3.181 (0.15) + 3.331 (0.13), rms 0.29 against this one. The slope
# (engine 3's, shifted) fits the ladder (mostly 13 7 4 3 1 1), but the
# calibration search's pool P at N' >= 12k came in at about 0.5x of it
# (7 sums; calibration-target.md 3.6): the online refit learns only the
# intercept (DFIRST_TIME_LAMBDA), so above 12k the planner charges pool P
# about 2x (2% of E at 1-10 CPU-years; the "measured" forecast truth
# corrects it there). The same shift goes into the ratio
# (DFIRST_RATIO_ENGINE_SHIFT), the plain search being unchanged.
DFIRST_E4 = (-0.103, -0.204)
DFIRST_E4_SD = (0.023, 0.020)
DFIRST_ENGINE_SHIFT = {4: {"intercept": DFIRST_E4[0], "ln N'/4000": DFIRST_E4[1]}}
# the CPU of the d loop is not uniform along d's index in the unreduced
# list: mean CPU per d in the deciles of u = d / N relative to the sum's
# mean, pooled over the perf verifier's 13 d logs (N 2-32k; the same shape
# below and above N 6k, se 0.01-0.04 per decile; |V_d| is flat along u).
# A range [lo, hi) of d costs W(hi / N) - W(lo / N) of the d loop, W the
# cumulative of this profile (W(0.5) = 0.524, W(0.9) = 0.961; at most 0.075
# from u). The yield per d was not measured (no pair in the samples): the
# scheduler credits a part of a sum with the same fraction of its E.
DFIRST_COST_PROFILE = (0.71, 1.08, 1.08, 1.21, 1.16, 1.21, 1.15, 1.15, 0.86, 0.39)

# The mode of a sum (plain or d-first) is chosen by the measured d-first /
# plain CPU ratio, not by the two time laws: ln r = a0 + a1 ln(N'/4000),
# pooled over the perf verifier's 13 sums with both modes (research/ideas.md,
# "Measurements on the integrated binary"; resid sd 0.17). The two laws
# were fitted separately and their errors do not cancel: the plain law
# under-predicts the plain search 1.35x at N' 5-8k (10 sums) and carries a
# label term L^-3.46 that the d-first law has not, so their quotient put the
# switch at N' ~ 9-10k (E-weighted) instead of the measured ~4.3k (with a 7%
# calibration stream: (1 + 0.07) r < 1), and chose plain on 5 of 6 pool sums
# at N' ~6.1k measured at d-first / plain 0.69-0.82. With this ratio the
# choice has regret 1.00 on all 22 engine-3 sums measured in both modes.
# Online the level a0 is learned (prior sd level_sd) from the pairs every
# d-first sum with a calibration stream gives: its d-first CPU against the
# stream's unbiased plain estimate (est_time).
DFIRST_RATIO_PRIOR = {"a": (-0.028, -0.566), "sd": 0.17, "level_sd": 0.1, "engine": 3}
# per-engine shifts (a0, a1) of the ratio after DFIRST_RATIO_PRIOR's engine
# (engine 4: the d-first law's shift, DFIRST_E4: -0.131 - 0.770
# ln(N'/4000), sd 0.17 kept, the auto switch at N' ~3.68k; measured on the
# verifier's 13 sums -0.113 (0.072) - 0.805 (0.062), rms 0.152 against it)
DFIRST_RATIO_ENGINE_SHIFT = {4: DFIRST_E4}

# ln nodes = NODE_TH . [1, ln(N'/4000), ln(labels/150), max(0, ln(N'/8000))]
NODE_TH = (17.375, 5.756, -5.765, 0.639)


def time_features(lNp, lLraw, k):
    """[..., NT] features of the time law at bias-corrected lN', L_raw and k
    distinct primes (scalar or broadcastable)"""
    lNp = np.asarray(lNp, float)
    lLraw = np.asarray(lLraw, float)
    lNp, lLraw, k = np.broadcast_arrays(lNp, lLraw, np.asarray(k, float))
    x = lNp - math.log(4000)
    w3 = (labels_obs(lNp, lLraw, k) > 128).astype(float)
    band = np.searchsorted(np.log(TIME_BAND_EDGES), lNp, side="right")
    cols = [np.ones_like(x), x, lLraw - math.log(150), np.maximum(0.0, lNp - math.log(8000)), w3]
    cols += [(band == b).astype(float) for b in range(len(TIME_BAND_EDGES) + 1)]
    return np.stack(cols, -1)


def log_nodes(lNp, lLraw, k):
    lab = labels_obs(lNp, lLraw, k)
    x = np.asarray(lNp, float) - math.log(4000)
    return (NODE_TH[0] + NODE_TH[1] * x + NODE_TH[2] * np.log(lab / 150)
            + NODE_TH[3] * np.maximum(0.0, np.asarray(lNp, float) - math.log(8000)))


def time_prior(engine, base=None, mode="plain"):
    """the prior of engine's time law for a search mode ("plain" or
    "dfirst"): base (a TimeModel of an older engine, default the shipped
    prior of the mode) with the ENGINE_TIME_SHIFT (DFIRST_ENGINE_SHIFT) of
    the engines after it, up to engine"""
    if base is None:
        base = TimeModel(mode=mode)
    shifts = DFIRST_ENGINE_SHIFT if mode == "dfirst" else ENGINE_TIME_SHIFT
    th = base.th.copy()
    for e in range(int(base.engine) + 1, int(engine) + 1):
        for f, d in shifts.get(e, {}).items():
            th[TIME_FEATURES.index(f)] += d
    return TimeModel(th, base.sd, engine, 0, mode)


def _cost_knots():
    prof = np.asarray(DFIRST_COST_PROFILE, float)
    return np.linspace(0, 1, len(prof) + 1), np.concatenate([[0.0], np.cumsum(prof) / prof.sum()])


def dfirst_cost_cum(u):
    """W(u): the share of a d-first sum's d-loop CPU in the d with index
    below u N (DFIRST_COST_PROFILE, piecewise linear)"""
    x, w = _cost_knots()
    return np.interp(np.clip(u, 0.0, 1.0), x, w)


def dfirst_cost_inv(w):
    """the inverse of dfirst_cost_cum"""
    x, wk = _cost_knots()
    return np.interp(np.clip(w, 0.0, 1.0), wk, x)


def dfirst_cost_frac(lo, hi, nd):
    """the share of a d-first sum's d-loop CPU (and of its E) in the d with
    index in [lo, hi) of nd"""
    if nd <= 0:
        return 0.0
    return float(dfirst_cost_cum(hi / nd) - dfirst_cost_cum(lo / nd))


def dfirst_ratio_coefs(engine=None):
    """the prior (a0, a1) of an msearch engine's d-first / plain ratio:
    DFIRST_RATIO_PRIOR (its engine, the default) plus the
    DFIRST_RATIO_ENGINE_SHIFT of the engines after it, up to engine"""
    a0, a1 = DFIRST_RATIO_PRIOR["a"]
    base = DFIRST_RATIO_PRIOR["engine"]
    for e in range(base + 1, int(base if engine is None else engine) + 1):
        d0, d1 = DFIRST_RATIO_ENGINE_SHIFT.get(e, (0.0, 0.0))
        a0, a1 = a0 + d0, a1 + d1
    return a0, a1


def dfirst_log_ratio(lNp, a0=None, engine=None):
    """ln(d-first CPU / plain CPU) of a sum at the model's ln N' (the
    ratio of engine, default DFIRST_RATIO_PRIOR's; level a0 if given)"""
    p0, a1 = dfirst_ratio_coefs(engine)
    return (p0 if a0 is None else a0) + a1 * (np.asarray(lNp, float) - math.log(4000))


def dfirst_ratio_level(st=None, engine=None, st_engine=None, prior=None):
    """the posterior mean of engine's ratio level a0 from pairs st = [n,
    sum x, sum y] (weighted; y = ln(d-first CPU / plain CPU), x =
    ln(N'/4000)) of engine st_engine (default engine), the slope held: a
    normal prior at prior (default st_engine's a0, dfirst_ratio_coefs) with
    sd level_sd, residual sd sd; pairs of an older engine stand in for a
    newer one's with the shift between their priors (as time_prior hands a
    time law over). To chain engines, pass as prior the level of the older
    engine's pairs shifted to st_engine (scheduler.current_ratio_level)."""
    a0, a1 = dfirst_ratio_coefs(st_engine if st_engine is not None else engine)
    shift = dfirst_ratio_coefs(engine)[0] - a0
    if prior is not None:
        a0 = float(prior)
    if not st or st[0] <= 0:
        return a0 + shift
    n, sx, sy = float(st[0]), float(st[1]), float(st[2])
    shrink = (DFIRST_RATIO_PRIOR["sd"] / DFIRST_RATIO_PRIOR["level_sd"]) ** 2
    return a0 + (sy - a1 * sx - n * a0) / (n + shrink) + shift


class TimeModel:
    """CPU seconds per sum from the model's N' and L_raw; one coefficient
    vector per msearch engine and search mode ("plain"; "dfirst": msearch
    --diag-first, whole-sum CPU). fit() is a Bayesian ridge towards the
    prior: b = (A + Lam0)^-1 (c + Lam0 b0) with A = F'F / s^2, c = F'y / s^2
    and Lam0 = diag(TIME_LAMBDA) (DFIRST_TIME_LAMBDA); s^2 is pooled with
    TIME_SD_PSEUDO pseudo-sums at the prior sd."""

    def __init__(self, th=None, sd=None, engine=None, n=0, mode="plain"):
        prior = DFIRST_TIME_PRIOR if mode == "dfirst" else TIME_PRIOR
        self.mode = mode
        self.th = np.array(prior["th"] if th is None else th, float)
        if self.th.shape != (NT,):
            raise ValueError(f"time law with {self.th.size} coefficients, expected {NT}")
        self.sd = float(prior["sd"] if sd is None else sd)
        self.engine = prior["engine"] if engine is None else engine
        self.n = n  # rows in the last fit

    def log_time(self, lNp, lLraw, k):
        return time_features(lNp, lLraw, k) @ self.th

    def time(self, lNp, lLraw, k):
        """expected CPU seconds per sum"""
        return np.exp(np.minimum(self.log_time(lNp, lLraw, k), 60.0) + 0.5 * self.sd ** 2) + SUM_OVERHEAD

    @staticmethod
    def stats(F, y):
        """sufficient statistics of rows (F, y): n, F'F, F'y, y'y"""
        F = np.asarray(F, float).reshape(-1, NT)
        y = np.asarray(y, float)
        return {"n": len(y), "FF": F.T @ F, "Fy": F.T @ y, "yy": float(y @ y)}

    def fit(self, st, prior=None, iters=4):
        """posterior from sufficient statistics st (see stats), starting from
        prior (a TimeModel; default the shipped prior)"""
        prior = prior or TimeModel(mode=self.mode)
        b0, s0_ = prior.th, prior.sd
        lam = np.diag(DFIRST_TIME_LAMBDA if self.mode == "dfirst" else TIME_LAMBDA)
        s2 = s0_ ** 2
        th = b0.copy()
        for _ in range(iters):
            A = st["FF"] / s2 + lam
            c = st["Fy"] / s2 + lam @ b0
            th = np.linalg.solve(A, c)
            rss = st["yy"] - 2 * th @ st["Fy"] + th @ st["FF"] @ th
            s2 = (max(rss, 0.0) + TIME_SD_PSEUDO * s0_ ** 2) / (st["n"] + TIME_SD_PSEUDO)
        return TimeModel(th, math.sqrt(s2), self.engine, st["n"], self.mode)

    def to_json(self):
        return {"th": [float(v) for v in self.th], "sd": self.sd, "engine": self.engine, "n": self.n,
                "mode": self.mode}

    @staticmethod
    def from_json(d):
        return TimeModel(d["th"], d["sd"], d.get("engine"), d.get("n", 0), d.get("mode", "plain"))
