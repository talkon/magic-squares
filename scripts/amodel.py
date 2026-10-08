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
shape and learns an intercept, the label-word step and N' band offsets).
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
# prior precision of th for the online refit. The shape (N, L and hinge
# slopes) is fixed: refitting it on the thousands of cheap sums at N' < 3k
# moved predictions at N' >= 3k by 1.4-2x. Learned online: a global
# intercept (sd 0.1), the W >= 3 step (sd 0.1) and one offset per N' band
# (sd 0.25), so a residual learned at small N' moves the other bands by
# about 14% of itself. TIME_SD_PSEUDO pseudo-sums at the prior sd anchor
# the residual variance.
TIME_LAMBDA = (100.0, 1e6, 1e6, 1e6, 100.0, 16.0, 16.0, 16.0, 16.0)
TIME_SD_PSEUDO = 20.0

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


class TimeModel:
    """CPU seconds per sum from the model's N' and L_raw; one coefficient
    vector per msearch engine (and search mode). fit() is a Bayesian ridge
    towards the prior: b = (A + Lam0)^-1 (c + Lam0 b0) with A = F'F / s^2,
    c = F'y / s^2 and Lam0 = diag(TIME_LAMBDA); s^2 is pooled with
    TIME_SD_PSEUDO pseudo-sums at the prior sd."""

    def __init__(self, th=None, sd=None, engine=None, n=0):
        self.th = np.array(TIME_PRIOR["th"] if th is None else th, float)
        if self.th.shape != (NT,):
            raise ValueError(f"time law with {self.th.size} coefficients, expected {NT}")
        self.sd = float(TIME_PRIOR["sd"] if sd is None else sd)
        self.engine = TIME_PRIOR["engine"] if engine is None else engine
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
        prior = prior or TimeModel()
        b0, s0_ = prior.th, prior.sd
        lam = np.diag(TIME_LAMBDA)
        s2 = s0_ ** 2
        th = b0.copy()
        for _ in range(iters):
            A = st["FF"] / s2 + lam
            c = st["Fy"] / s2 + lam @ b0
            th = np.linalg.solve(A, c)
            rss = st["yy"] - 2 * th @ st["Fy"] + th @ st["FF"] @ th
            s2 = (max(rss, 0.0) + TIME_SD_PSEUDO * s0_ ** 2) / (st["n"] + TIME_SD_PSEUDO)
        return TimeModel(th, math.sqrt(s2), self.engine, st["n"])

    def to_json(self):
        return {"th": [float(v) for v in self.th], "sd": self.sd, "engine": self.engine, "n": self.n}

    @staticmethod
    def from_json(d):
        return TimeModel(d["th"], d["sd"], d.get("engine"), d.get("n", 0))
