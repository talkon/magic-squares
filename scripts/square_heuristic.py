#!/usr/bin/env python3
"""
First-principles per-square predictions of the rarer events of a 6x6
semi-magic square: how many traversals (possible diagonals) have the magic
sum S, the magic product P or both, how many of the 5400 partner pairs of
traversals are of each type (S+S ... SP+SP, i.e. the pair scores of
square_diag_stats in src/c/square.c), and the probability of each
best_score. Nothing is fitted: the prediction uses only the square's 36
entries (moments, cumulants and integer lattices), never which traversals
actually hit S or P (except in mode="conditional", see below).

    predict(grid, S, P, mode="moment") -> {event: expected count}
    observe(grid, S, P)                -> {event: exact count}  (same keys)
    components(grid, S, P) / combine(comp, mode=..., ew=..., pair_lattice=...)
        (the expensive part once, then any variant; ~0.1 s per square)

grid: 6x6 ints (an msearch "square" record's "grid"); S: the magic sum; P:
the record's exponent list (e.g. [12, 6, 3, 2, 1, 0, 1]) or the integer.

Keys: S, P, SP (traversals with sum S incl. SP, product P incl. SP, both),
S_only, P_only, none; the unordered partner-pair types 0+0 S+0 P+0 S+S S+P
P+P SP+0 SP+S SP+P SP+SP (exclusive classes, scores 0,2,3,4,5,6,7,9,10,14);
best=k and best>=k (probabilities, k in those scores); magic = SP+SP.

The model (the derivation and its validation on all squares found in October
2026 are in research/calibration.md):

* A traversal sigma gives Z(sigma) = sum_i z[i, sigma(i)] with
  z[i, j] = (a_ij, v_p1(a_ij), ..., v_pk(a_ij)); sigma has sum S iff X = Z_0 =
  S, product P iff Y = Z_1.. = (v_p(P)). For a semi-magic square E[Z] =
  t = (S, v(P)) exactly over a uniformly random sigma, and Hoeffding's formula
  gives Cov Z = 1/(n-1) sum_ij d_ij d_ij^T, d_ij = z_ij - t/n.
* Z lives on Z(id) + L, L the lattice of the 2x2 interaction contrasts. Point
  probabilities at the mean come from the lattice local CLT,
      P(Z_I = t_I) ~ [t_I in Z(id)_I + L_I] covol(L_I) / ((2 pi)^(r/2) sqrt(det Cov Z_I)),
  times first-order Edgeworth factors (ew="pair": every coordinate's
  univariate factor and every pair of coordinates' bivariate factor, in a
  Kirkwood product, each clamped to [0.25, 4]).
* Pairs: W = (Z(sigma), Z(sigma o tau)) for each of the 15 fixed-point-free
  involutions tau, with the exact covariance [[Sigma, C_tau], [C_tau, Sigma]],
  C_tau = 1/(n-1) sum_i sum_c d_ic d_tau(i)c^T, and the exact lattice of the
  pair. That lattice has >= k + 1 - 4 independent congruences mod 3 (k = number
  of primes of P), so a pair type tau can be magic only if (S, v(P)) mod 3 lies
  in a subspace depending on the square, and is then 3^d times likelier
  (pair_lattice=False ignores this). Class indicators (none, S only, P only,
  SP) are inclusion-exclusion sums of point events.
* best_score: Poisson clumping (approximate; calibrate.py also evaluates the
  best pair exactly on the partner graph under i.i.d. traversal classes).

mode="moment" (default): everything from moments and lattices.
mode="conditional": one rung at a time. Takes the square's realized
inclusive S and P counts; SP = n_S n_P / 720 * rho_square (the moment model's
S-P correlation of this square); pair rungs = the own-count null on those
classes (with the predicted SP) times the moment model's per-tau pair
correlation. Never uses the realized SP or pair counts.

Approximations, explicitly: uniform random traversal with exact moments and
lattices; local CLT at the mean with only 6 terms; first-order Edgeworth kept
to pairwise order and multiplied (an ad-hoc truncation chosen because it fits
the P+P rung and the 12-cell product event, and keeps the P rung's ~12%
deficit constant in the number of primes); the 3^d congruence boost assumes
smoothness on the coarse lattice; Poisson clumping for best_score.

usage: square_heuristic.py FILE.jsonl ...   (prints observed vs predicted
       for each square record)
"""
import itertools
import json
import math
import sys

import numpy as np

PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
N = 6
PERMS = np.array(list(itertools.permutations(range(N))), dtype=np.int64)  # lexicographic
PERM_INDEX = {tuple(p): i for i, p in enumerate(PERMS.tolist())}
ROWS = np.arange(N)


def _partners(n=N):
    """fixed-point-free involutions of {0..n-1}, in the order of
    gen_partners() in src/c/square.c"""
    out = []

    def rec(cur):
        a = next((i for i in range(n) if cur[i] < 0), None)
        if a is None:
            out.append(tuple(cur))
            return
        for b in range(a + 1, n):
            if cur[b] < 0:
                cur[a], cur[b] = b, a
                rec(cur)
                cur[a] = cur[b] = -1

    rec([-1] * n)
    return out


PARTNERS = np.array(_partners(), dtype=np.int64)       # 15 x 6
# PIDX[s, t] = index of the permutation p2 with p2[i] = perm_s[partners[t][i]]
PIDX = np.array([[PERM_INDEX[tuple(PERMS[s][PARTNERS[t]])] for t in range(len(PARTNERS))]
                 for s in range(len(PERMS))], dtype=np.int64)
NT = len(PERMS)            # 720 traversals
NTAU = len(PARTNERS)       # 15 partner involutions
NPAIRS = NT * NTAU // 2    # 5400 unordered partner pairs

# traversal classes 0 none, 1 S only, 2 P only, 3 SP; scores as square.c
SCORE = np.array([0, 2, 3, 7])
PAIR_NAMES = {(0, 0): "0+0", (0, 1): "S+0", (0, 2): "P+0", (1, 1): "S+S", (1, 2): "S+P",
              (2, 2): "P+P", (0, 3): "SP+0", (1, 3): "SP+S", (2, 3): "SP+P", (3, 3): "SP+SP"}
PAIR_KEYS = list(PAIR_NAMES.values())
PAIR_SCORE = {"0+0": 0, "S+0": 2, "P+0": 3, "S+S": 4, "S+P": 5, "P+P": 6,
              "SP+0": 7, "SP+S": 9, "SP+P": 10, "SP+SP": 14}
SCORES = [0, 2, 3, 4, 5, 6, 7, 9, 10, 14]


# ---------------------------------------------------------------------------
# exact counts


def p_exps(P):
    """record P (exponent list) or int -> (primes, exponents) of the primes
    dividing P"""
    if isinstance(P, (list, tuple)):
        return [p for p, e in zip(PRIMES, P) if e], [e for e in P if e]
    pr, ex = [], []
    for p in PRIMES:
        e = 0
        while P % p == 0:
            P //= p
            e += 1
        if e:
            pr.append(p)
            ex.append(e)
    assert P == 1, "P has a prime factor > 31"
    return pr, ex


def cell_vectors(grid, P):
    """Z[i, j] = (a_ij, v_p1(a_ij), ..., v_pk(a_ij)) (6x6x(k+1) int64) and
    a = (0, v_p1(P), ..., v_pk(P))"""
    pr, ex = p_exps(P)
    g = np.array(grid, dtype=np.int64)
    Z = np.zeros((N, N, 1 + len(pr)), dtype=np.int64)
    Z[:, :, 0] = g
    for t, p in enumerate(pr):
        x = g.copy()
        e = np.zeros_like(g)
        while True:
            m = (x % p == 0)
            if not m.any():
                break
            e += m
            x = np.where(m, x // p, x)
        Z[:, :, 1 + t] = e
    return Z, np.array([0] + ex, dtype=np.int64)


def traversal_values(Z):
    """W[s] = sum_i Z[i, perm_s(i)], shape 720 x (k+1)"""
    return Z[ROWS[None, :], PERMS].sum(axis=1)


def classes(Z, S, a):
    W = traversal_values(Z)
    isS = W[:, 0] == S
    isP = np.all(W[:, 1:] == a[1:], axis=1)
    return isS.astype(np.int64) + 2 * isP.astype(np.int64)   # 0 none, 1 S, 2 P, 3 SP


def observe(grid, S, P):
    """exact counts with the keys of predict() (best=k / best>=k are 0/1);
    as square_diag_stats: 720 traversals, pairs {sigma, sigma o tau}"""
    Z, a = cell_vectors(grid, P)
    c = classes(Z, S, a)
    out = {"S": int(np.sum((c == 1) | (c == 3))), "P": int(np.sum(c >= 2)),
           "SP": int(np.sum(c == 3)), "S_only": int(np.sum(c == 1)),
           "P_only": int(np.sum(c == 2)), "none": int(np.sum(c == 0))}
    c1 = np.repeat(c[:, None], NTAU, axis=1)
    c2 = c[PIDX]
    code = np.minimum(c1, c2) * 4 + np.maximum(c1, c2)
    cnt = np.bincount(code.ravel(), minlength=16)
    best = -1
    for (x, y), nm in PAIR_NAMES.items():
        v = int(cnt[x * 4 + y] // 2)   # each unordered pair is seen twice
        out[nm] = v
        if v:
            best = max(best, PAIR_SCORE[nm])
    for k in SCORES:
        out[f"best={k}"] = int(best == k)
        out[f"best>={k}"] = int(best >= k)
    out["best"] = best
    out["magic"] = out["SP+SP"]
    return out


# ---------------------------------------------------------------------------
# integer lattices (echelon / Hermite form with Python ints)


def pivots(basis):
    return [next(j for j, x in enumerate(r) if x) for r in basis]


def echelon(rows):
    """row-echelon (Hermite) basis, positive pivots, of the lattice generated
    by the integer rows"""
    A = [list(map(int, r)) for r in rows if any(r)]
    if not A:
        return []
    d = len(A[0])
    basis = []
    col = 0
    while A and col < d:
        nz = [r for r in A if r[col] != 0]
        rest = [r for r in A if r[col] == 0]
        if not nz:
            col += 1
            continue
        while len(nz) > 1:
            nz.sort(key=lambda r: abs(r[col]))
            p = nz[0]
            new = [p]
            for r in nz[1:]:
                q = r[col] // p[col]
                r2 = [x - q * y for x, y in zip(r, p)]
                if r2[col] != 0:
                    new.append(r2)
                elif any(r2):
                    rest.append(r2)
            nz = new
        p = nz[0]
        if p[col] < 0:
            p = [-x for x in p]
        basis.append(p)
        A = rest
        col += 1
    piv = pivots(basis)
    for j in range(len(basis)):
        c, h = piv[j], basis[j][piv[j]]
        for i in range(j):
            q = basis[i][c] // h
            if q:
                basis[i] = [x - q * y for x, y in zip(basis[i], basis[j])]
    return basis


def reduce_vecs(V, basis):
    """reduce the rows of the integer array V by an echelon basis (a copy;
    Python ints when the entries are big)"""
    big = max((abs(x) for r in basis for x in r), default=0) > 2 ** 31 or \
        (V.dtype != object and np.abs(V).max(initial=0) > 2 ** 31)
    V = V.astype(object) if big else V.copy()
    for r, c in zip(basis, pivots(basis)):
        row = np.array(r, dtype=object if big else np.int64)
        q = V[:, c] // row[c]
        V = V - q[:, None] * row[None, :]
    return V


def lattice_of(V):
    """echelon basis of the lattice generated by the rows of the int array V
    (incremental: only rows not yet in the lattice are added)"""
    V = np.unique(V[np.any(V != 0, axis=1)], axis=0)
    if not len(V):
        return []
    rng = np.random.default_rng(0)
    d = V.shape[1]
    first = rng.choice(len(V), min(len(V), 2 * d), replace=False)
    basis = echelon(V[first].tolist())
    while True:
        R = reduce_vecs(V, basis)
        bad = np.any(R != 0, axis=1)
        if not bad.any():
            break
        basis = echelon(basis + R[bad][:4].tolist())
        V = V[bad]
    return basis


def in_lattice(v, basis):
    R = reduce_vecs(np.array([v], dtype=np.int64), basis) if basis else np.array([v])
    return not np.any(R != 0)


def project(basis, idx):
    return echelon([[r[i] for i in idx] for r in basis])


def sat_index(B):
    """index of the lattice spanned by the rows of B in its saturation
    (span_R(B) intersected with Z^m) = gcd of the maximal minors"""
    r, m = len(B), len(B[0])
    if r == m:
        return abs(int(np.prod([row[c] for row, c in zip(B, pivots(B))])))
    Bm = np.array(B, dtype=float)
    g = 0
    for cols in itertools.combinations(range(m), r):
        g = math.gcd(g, int(round(abs(np.linalg.det(Bm[:, cols])))))
        if g == 1:
            break
    return g


def lclt_point(basis_full, idx, cov_full, offset, saturate=False):
    """lattice local-CLT probability that the coordinates idx of the integer
    vector are at their mean, given the lattice of its differences, its
    covariance and offset = mean - value at the identity (an integer vector).
    0 if the mean is not in the coset. saturate=True ignores the congruences
    (the lattice is replaced by the integer points of its real span)."""
    if not idx:
        return 1.0
    B = project(basis_full, idx)
    off = [int(round(offset[i])) for i in idx]
    if not B:
        return 1.0 if not any(off) else 0.0
    if not saturate and not in_lattice(off, B):
        return 0.0
    Bm = np.array(B, dtype=float)                 # r x m
    r = Bm.shape[0]
    M = Bm.T @ np.linalg.inv(Bm @ Bm.T)           # m x r: lattice coordinates c = x M
    Cc = M.T @ cov_full[np.ix_(idx, idx)] @ M
    sign, ld = np.linalg.slogdet(Cc)
    if sign <= 0:
        return float("nan")
    p = math.exp(-0.5 * (r * math.log(2 * math.pi) + ld))
    return p / sat_index(B) if saturate else p


def lattice_factor(basis_full, idx, offset):
    """[mean in coset] * [saturation : lattice] for the coordinates idx"""
    if not idx:
        return 1.0
    B = project(basis_full, idx)
    off = [int(round(offset[i])) for i in idx]
    if not B:
        return 1.0 if not any(off) else 0.0
    return sat_index(B) if in_lattice(off, B) else 0.0


# ---------------------------------------------------------------------------
# moments and Edgeworth factors


def moments(grid, S, P):
    Z, a = cell_vectors(grid, P)
    k1 = Z.shape[2]
    t = a.astype(float).copy()
    t[0] = S
    D = Z.astype(float) - t / N
    Dr = D.reshape(N * N, k1)
    Sigma = Dr.T @ Dr / (N - 1)
    # the adjacent 2x2 interaction contrasts generate the lattice
    G = (Z[:-1, :-1] + Z[1:, 1:] - Z[:-1, 1:] - Z[1:, :-1]).reshape(-1, k1)
    L = echelon(G.tolist())
    Zid = Z[ROWS, ROWS].sum(axis=0)
    off = t - Zid
    Ctau = [np.einsum("icx,icy->xy", D, D[tau]) / (N - 1) for tau in PARTNERS]
    return dict(Z=Z, t=t, a=a, D=D, Sigma=Sigma, L=L, off=off, Ctau=Ctau, k1=k1)


# class indicator as a signed sum over constraint sets: classes 0 none, 1 S,
# 2 P, 3 SP; constraint sets 0 = {}, 1 = {X = S}, 2 = {Y = a}, 3 = both
COEF = {0: {0: 1, 1: -1, 2: -1, 3: 1}, 1: {1: 1, 3: -1}, 2: {2: 1, 3: -1}, 3: {3: 1}}
EW_CLAMP = (0.25, 4.0)


def coords(cset, k1, shift=0):
    out = []
    if cset & 1:
        out.append(shift)
    if cset & 2:
        out += list(range(shift + 1, shift + k1))
    return out


def edgeworth_pairs(V, mean):
    """F[i, i] = univariate and F[i, j] = bivariate first-order Edgeworth
    factors of the density at the mean, for the uniform distribution on the
    rows of V:
      1 + (1/8) k_abcd G_ab G_cd - (1/8) k_abc k_def G_ab G_cd G_ef
        - (1/12) k_abc k_def G_ad G_be G_cf    (G = inverse covariance)"""
    X = np.asarray(V, float) - mean
    n, m = X.shape
    C = X.T @ X / n
    X2, X3 = X * X, X * X * X
    M21 = X2.T @ X / n          # E[x_i^2 x_j]
    M31 = X3.T @ X / n          # E[x_i^3 x_j]
    M22 = X2.T @ X2 / n         # E[x_i^2 x_j^2]
    d = np.diag(C)
    m3 = np.diag(M21)
    m4 = np.diag(M31)
    F = np.ones((m, m))
    with np.errstate(divide="ignore", invalid="ignore"):
        k3 = m3 / d ** 1.5
        k4 = m4 / d ** 2 - 3
        f1 = 1 + k4 / 8 - 5 * k3 ** 2 / 24
    f1 = np.where(d > 0, f1, 1.0)
    iu, ju = np.triu_indices(m, 1)
    s11, s22, s12 = d[iu], d[ju], C[iu, ju]
    det = s11 * s22 - s12 ** 2
    good = det > 1e-9 * np.maximum(s11 * s22, 1e-300)
    npair = len(iu)
    K3 = np.zeros((npair, 2, 2, 2))
    t30, t21, t12, t03 = m3[iu], M21[iu, ju], M21[ju, iu], m3[ju]
    for (a, b, c), v in {(0, 0, 0): t30, (0, 0, 1): t21, (0, 1, 1): t12, (1, 1, 1): t03}.items():
        for perm in set(itertools.permutations((a, b, c))):
            K3[(slice(None),) + perm] = v
    S2 = np.zeros((npair, 2, 2))
    S2[:, 0, 0], S2[:, 1, 1], S2[:, 0, 1], S2[:, 1, 0] = s11, s22, s12, s12
    raw = {0: m4[iu], 1: M31[iu, ju], 2: M22[iu, ju], 3: M31[ju, iu], 4: m4[ju]}  # by # of 2nd index
    K4 = np.zeros((npair, 2, 2, 2, 2))
    for idx in itertools.product((0, 1), repeat=4):
        a, b, c, e = idx
        K4[(slice(None),) + idx] = raw[sum(idx)] - (S2[:, a, b] * S2[:, c, e] + S2[:, a, c] * S2[:, b, e] +
                                                     S2[:, a, e] * S2[:, b, c])
    G = np.zeros((npair, 2, 2))
    dd = np.where(good, det, 1.0)
    G[:, 0, 0], G[:, 1, 1] = s22 / dd, s11 / dd
    G[:, 0, 1] = G[:, 1, 0] = -s12 / dd
    A = np.einsum("pabcd,pab,pcd->p", K4, G, G)
    B = np.einsum("pabc,pdef,pab,pcd,pef->p", K3, K3, G, G, G)
    Cc = np.einsum("pabc,pdef,pad,pbe,pcf->p", K3, K3, G, G, G)
    f2 = 1 + A / 8 - B / 8 - Cc / 12
    f2 = np.where(good, f2, f1[iu] * f1[ju])
    lo, hi = EW_CLAMP
    F[iu, ju] = F[ju, iu] = np.clip(f2, lo, hi)
    F[np.arange(m), np.arange(m)] = np.clip(f1, lo, hi)
    return F


def ew_factors(F, J, k1=None):
    """Edgeworth factors of the coordinate set J (0..k1-1 = first diagonal
    (X, Y_p...), k1..2k1-1 = second diagonal), as (marg, pair, mixed):
      marg  = prod_i f_i
      pair  = marg * prod_{i<j} f_ij / (f_i f_j)          (all pairs)
      mixed = marg * prod over the pairs that are not two exponent
              coordinates of the same diagonal of f_ij / (f_i f_j)"""
    J = sorted(J)
    f = np.diag(F)
    marg = float(np.prod(f[J])) if J else 1.0
    pw = mx = marg
    m = k1 if k1 else F.shape[0]
    for a in range(len(J)):
        for b in range(a + 1, len(J)):
            i, j = J[a], J[b]
            r = F[i, j] / (f[i] * f[j])
            pw *= r
            if not ((i // m) == (j // m) and (i % m) != 0 and (j % m) != 0):
                mx *= r
    return marg, float(pw), float(mx)


def sub_sets(k1):
    """coordinate sets I of the sub-events (0 = the sum X, 1..k = the
    exponents of the k primes of P): the k cyclic windows of d consecutive
    exponent coordinates for d < k, all k for d = k, each alone and together
    with X. Deterministic and balanced (every prime is in d of the k windows
    of size d). The full sets are rungs of the ladder: (0,) = S, (1..k) = P,
    (0..k) = SP (one diagonal), and on both diagonals of a partner pair
    "both sums S", "both products P" and magic."""
    k = k1 - 1
    wins = []
    for d in range(1, k + 1):
        for j in (range(k) if d < k else [0]):
            w = tuple(sorted(1 + (j + i) % k for i in range(d)))
            if w not in wins:
                wins.append(w)
    return wins + [(0,)] + [(0,) + w for w in wins]


SUB_VARIANTS = ("gauss", "marg", "pair", "mixed")   # order of e1 / e2 in the sub-event rows
CONG_CLASSES = ("forbidden", "no extra congruence", "extra congruences")


def pair_congruence_class(L2, idx2, off, lam1):
    """congruence class of a partner involution for the coordinates idx2 of
    both diagonals: 0 if the target is not in the coset of the pair's lattice
    (no pair of this type can hit it), 2 if it is and the pair lattice is
    finer than the product of the two single-diagonal lattices (index
    sat_index / lam1^2 > 1: the heuristic boosts these pairs by that index),
    else 1"""
    B = project(L2, idx2)
    if not B:
        return 1
    if not in_lattice([int(round(off[i])) for i in idx2], B):
        return 0
    return 2 if lam1 and sat_index(B) > lam1 * lam1 else 1


def components(grid, S, P, sub_events=False):
    """the expensive, variant-independent part: Gaussian local-CLT point
    probabilities (one traversal; and per tau for pairs, with and without the
    pair congruences) and the Edgeworth factor triples.

    sub_events=True also returns comp["sub"]: for every coordinate set I of
    sub_sets(), the exact number of traversals with Z_I = t_I (n1) and of
    unordered partner pairs with both diagonals so (n2), and the predictions
    e1, e2 (one value per SUB_VARIANTS: plain Gaussian local CLT, univariate,
    pairwise (the default) and mixed Edgeworth factors) and e2n (pairwise
    Edgeworth, without the congruences of the pair), and split by the
    congruence class of each tau (pair_congruence_class: forbidden / no extra
    congruence / extra congruences): the number of taus ntc, the observed
    pairs n2c, e2c and e2nc. These events pin fewer coordinates than the
    rungs, so they are 10-10^5x commoner: they test the local-CLT machinery,
    and the two parts of the pair factor (congruences, correlation), below
    the magic rung."""
    m = moments(grid, S, P)
    k1, t = m["k1"], m["t"]
    Wall = traversal_values(m["Z"])
    F1 = edgeworth_pairs(Wall, t)
    Qg = {c: lclt_point(m["L"], coords(c, k1), m["Sigma"], m["off"]) for c in (1, 2, 3)}
    Qg[0] = 1.0
    ew1 = {c: ew_factors(F1, coords(c, k1), k1) for c in range(4)}
    lam = {c: lattice_factor(m["L"], coords(c, k1), m["off"]) for c in (1, 2, 3)}
    lam[0] = 1.0
    sub = None
    if sub_events:
        ti_int = np.round(t).astype(np.int64)
        sub, sub_lam, sub_hit = [], [], []
        for I in sub_sets(k1):
            I = list(I)
            hit = np.all(Wall[:, I] == ti_int[I], axis=1)
            g1 = lclt_point(m["L"], I, m["Sigma"], m["off"])
            f1 = ew_factors(F1, I, k1)
            sub.append({"I": I, "n1": int(hit.sum()), "n2": int((hit[:, None] & hit[PIDX]).sum() // 2),
                        "e1": [NT * g1] + [NT * g1 * f for f in f1], "e2": [0.0] * 4, "e2n": 0.0,
                        "ntc": [0, 0, 0], "n2c": [0, 0, 0], "e2c": [0.0] * 3, "e2nc": [0.0] * 3})
            sub_lam.append(lattice_factor(m["L"], I, m["off"]))
            sub_hit.append(hit)
    pairs = []
    tt = np.concatenate([t, t])
    for ti in range(NTAU):
        W = np.concatenate([Wall, Wall[PIDX[:, ti]]], axis=1)
        L2 = lattice_of(W - W[0])
        C = m["Ctau"][ti]
        cov = np.block([[m["Sigma"], C], [C, m["Sigma"]]])
        off = tt - W[0]
        Q2g, Q2n = {(0, 0): 1.0}, {(0, 0): 1.0}
        for c1 in range(4):
            for c2 in range(4):
                if (c1, c2) in Q2g:
                    continue
                if (c2, c1) in Q2g:        # symmetric under swapping the diagonals
                    Q2g[(c1, c2)] = Q2g[(c2, c1)]
                    Q2n[(c1, c2)] = Q2n[(c2, c1)]
                    continue
                idx = coords(c1, k1) + coords(c2, k1, k1)
                Q2g[(c1, c2)] = lclt_point(L2, idx, cov, off)
                # without the congruences of the pair (those of each diagonal kept)
                Q2n[(c1, c2)] = lclt_point(L2, idx, cov, off, saturate=True) * lam[c1] * lam[c2]
        F2 = edgeworth_pairs(W, tt)
        ew2 = {(c1, c2): ew_factors(F2, coords(c1, k1) + coords(c2, k1, k1), k1)
               for c1 in range(4) for c2 in range(4)}
        pairs.append({"Q2g": Q2g, "Q2n": Q2n, "ew2": ew2, "index": sat_index(L2),
                      "ok": bool(in_lattice([int(round(x)) for x in off], L2))})
        if sub is not None:
            for row, lm, h in zip(sub, sub_lam, sub_hit):
                idx2 = row["I"] + [k1 + i for i in row["I"]]
                g = lclt_point(L2, idx2, cov, off)
                f = ew_factors(F2, idx2, k1)
                row["e2"][0] += NT / 2 * g
                for j in range(3):
                    row["e2"][j + 1] += NT / 2 * g * f[j]
                gn = lclt_point(L2, idx2, cov, off, saturate=True) * lm * lm if lm else 0.0
                row["e2n"] += NT / 2 * gn * f[1]
                # this tau's congruence class (exact): 0 = the target is not in
                # the pair's coset (no pair can hit), 1 = no extra congruence,
                # 2 = extra congruences (the heuristic boosts these pairs)
                c = pair_congruence_class(L2, idx2, off, lm)
                row["ntc"][c] += 1
                row["n2c"][c] += int((h & h[PIDX[:, ti]]).sum()) // 2
                row["e2c"][c] += NT / 2 * g * f[1]
                row["e2nc"][c] += NT / 2 * gn * f[1]
    cls = classes(m["Z"], S, m["a"])
    out = {"k1": k1, "f": [float(x) for x in np.diag(F1)], "Qg": Qg, "ew1": ew1, "pairs": pairs,
           "n": [int(np.sum(cls == u)) for u in range(4)], "varX": float(m["Sigma"][0, 0])}
    if sub is not None:
        out["sub"] = sub
    return out


def class_probs(Q):
    return {u: sum(cf * Q[c] for c, cf in COEF[u].items()) for u in range(4)}


def pair_class_probs(Q2):
    return {(u, v): sum(cu * cv * Q2[(c1, c2)] for c1, cu in COEF[u].items()
                        for c2, cv in COEF[v].items()) for u in range(4) for v in range(4)}


# ---------------------------------------------------------------------------
# best_score: Poisson clumping


def best_probs(trav, pairs, sp_partner):
    """P(best_score = k) from expected counts (Poisson clumping).

    trav: expected traversal counts by class; pairs: expected pair counts by
    name; sp_partner: v -> P(a given partner of an SP traversal is of class v).
    SP levels (7, 9, 10, 14): Poisson number of SP traversals, each with 15
    partners of class v independently with probability sp_partner[v]. Lower
    levels (2..6), given no SP traversal: Poisson number of clumps (a
    traversal of the higher class with >= 1 qualifying partner: mu pairs over
    n carriers -> n (1 - exp(-mu / n)) clumps)."""
    lam_sp = trav[3]
    q_s, q_p, q_sp = (sp_partner.get(v, 0.0) for v in (1, 2, 3))
    ge = {7: 1 - math.exp(-lam_sp),
          9: 1 - math.exp(-lam_sp * (1 - max(0.0, 1 - q_sp - q_p - q_s) ** NTAU)),
          10: 1 - math.exp(-lam_sp * (1 - max(0.0, 1 - q_sp - q_p) ** NTAU)),
          14: 1 - math.exp(-lam_sp * (1 - max(0.0, 1 - q_sp) ** NTAU) / 2)}

    def clumps(names, n):
        mu = sum(pairs[nm] for nm in names)
        if mu <= 0:
            return 0.0
        if n <= 1e-12:
            return mu
        return n * (1 - math.exp(-mu / n))

    lam6 = clumps(["P+P"], trav[2])
    lam5 = lam6 + clumps(["S+P"], trav[2])
    lam4 = lam5 + clumps(["S+S"], trav[1])
    lam3 = lam4 + clumps(["P+0"], trav[2])
    lam2 = lam3 + clumps(["S+0"], trav[1])
    for k, lam in ((6, lam6), (5, lam5), (4, lam4), (3, lam3), (2, lam2)):
        ge[k] = ge[7] + (1 - ge[7]) * (1 - math.exp(-lam))
    ge[0] = 1.0
    for i in range(len(SCORES) - 2, -1, -1):   # monotone
        ge[SCORES[i]] = max(ge[SCORES[i]], ge[SCORES[i + 1]])
    out = {}
    for i, k in enumerate(SCORES):
        nxt = ge[SCORES[i + 1]] if i + 1 < len(SCORES) else 0.0
        out[f"best={k}"] = max(0.0, ge[k] - nxt)
        out[f"best>={k}"] = ge[k]
    return out


# ---------------------------------------------------------------------------
# predictions


def cond_base(x, y, trav, nS, nP, rho):
    """conditional mode: P(a random ordered partner pair has classes (x, y))
    if the square's classes were placed at random, estimated without bias
    from the realized inclusive counts nS, nP (products of rates by falling
    factorials; SP = rho (nS/720) (nP/720), the realized SP is never used).
    The classes are exclusive: an SP diagonal with an S_only partner is (SP,
    S inclusive) minus (SP, SP), and likewise for P_only."""
    N1, N2 = NT, NT * (NT - 1)
    if x != 3 and y != 3:
        if x == y:
            return trav[x] * max(trav[x] - 1, 0.0) / N2
        return trav[x] * trav[y] / N2
    spsp = rho ** 2 * nS * max(nS - 1, 0) * nP * max(nP - 1, 0) / N2 ** 2
    if x == 3 and y == 3:
        return spsp
    o = y if x == 3 else x
    if o == 1:
        return max(rho * nS * max(nS - 1, 0) * nP / (N2 * N1) - spsp, 0.0)
    if o == 2:
        return max(rho * nS * nP * max(nP - 1, 0) / (N1 * N2) - spsp, 0.0)
    return trav[3] * trav[0] / N2


def combine(comp, mode="moment", ew="pair", pair_lattice=True, details=False):
    """predictions from components(). mode: "moment" or "conditional"; ew:
    "pair" (default), "marg", "mixed" or None (plain Gaussian); pair_lattice:
    use the exact lattice of the pair of diagonals (with its congruences mod
    3) or only those of each diagonal"""
    k1 = comp["k1"]
    if ew is True:
        ew = "pair"
    w = {"marg": 0, "pair": 1, "mixed": 2}.get(ew)

    def fac1(c):
        return 1.0 if w is None else comp["ew1"][c][w]

    Q = {c: comp["Qg"][c] * fac1(c) for c in range(4)}
    P2 = []
    for pr in comp["pairs"]:
        Q2 = {cc: q * (1.0 if w is None else pr["ew2"][cc][w])
              for cc, q in pr["Q2g" if pair_lattice else "Q2n"].items()}
        P2.append(pair_class_probs(Q2))
    pu = class_probs(Q)
    rho = Q[3] / (Q[1] * Q[2]) if Q[1] > 0 and Q[2] > 0 else 0.0
    if mode == "moment":
        trav = {u: NT * pu[u] for u in range(4)}
        pc = {}
        for (u, v), nm in PAIR_NAMES.items():
            s = sum(p[(u, v)] for p in P2)
            if u != v:
                s += sum(p[(v, u)] for p in P2)
            pc[nm] = s * NT / 2
    elif mode == "conditional":
        n = comp["n"]
        nS, nP = n[1] + n[3], n[2] + n[3]
        sp = nS * nP / NT * rho
        trav = {1: max(nS - sp, 0.0), 2: max(nP - sp, 0.0), 3: sp}
        trav[0] = NT - trav[1] - trav[2] - trav[3]
        pc = {}
        for (u, v), nm in PAIR_NAMES.items():
            tot = 0.0
            for p in P2:
                for (x, y) in ({(u, v), (v, u)} if u != v else {(u, u)}):
                    den = pu[x] * pu[y]
                    kap = p[(x, y)] / den if den > 0 else 1.0
                    tot += kap * cond_base(x, y, trav, nS, nP, rho)
            pc[nm] = tot * NT / 2
    else:
        raise ValueError(mode)
    sp_partner = {}
    if pu[3] > 0:
        for v in (1, 2, 3):
            sp_partner[v] = float(np.mean([p[(3, v)] for p in P2]) / pu[3])
            if mode == "conditional" and pu[v] > 0:
                sp_partner[v] *= (trav[v] / NT) / pu[v]
    out = {"S": trav[1] + trav[3], "P": trav[2] + trav[3], "SP": trav[3],
           "S_only": trav[1], "P_only": trav[2], "none": trav[0]}
    out.update(pc)
    out.update(best_probs(trav, pc, sp_partner))
    out["magic"] = out["SP+SP"]
    if details:
        out["_rho"] = rho
        out["_f"] = list(comp["f"])
        out["_kappa_SPSP"] = float(np.mean([p[(3, 3)] for p in P2]) / pu[3] ** 2) if pu[3] > 0 else float("nan")
        out["_pair_ok"] = sum(pr["ok"] for pr in comp["pairs"])
    return out


def predict(grid, S, P, mode="moment", ew="pair", pair_lattice=True, details=False):
    """expected count of every event for one semi-magic square"""
    return combine(components(grid, S, P), mode=mode, ew=ew, pair_lattice=pair_lattice,
                   details=details)


def main():
    keys = ["S", "P", "SP", "S+S", "S+P", "P+P", "SP+S", "SP+P", "SP+SP"]
    print("hash      P                     S  " + " ".join(f"{k:>14}" for k in keys) + "   (observed / predicted)")
    for path in sys.argv[1:]:
        with open(path) as f:
            for line in f:
                if '"type":"square"' not in line:
                    continue
                r = json.loads(line)
                if r.get("n") != 6 or "grid" not in r:
                    continue
                o = observe(r["grid"], r["S"], r["P"])
                e = predict(r["grid"], r["S"], r["P"])
                print(f"{r['hash'][:8]}  {' '.join(map(str, r['P'])):20} {r['S']:4d}  " +
                      " ".join(f"{o[k]:>5} / {e[k]:<6.3g}" for k in keys))


if __name__ == "__main__":
    main()
