# How many additive-multiplicative magic squares should exist? An analytic heuristic, calibrated

Lens: analytic number-theory heuristic (exponent matrices + local limit theorem), calibrated on the
searches' data. Orders n = 5, 6, 7, 8. Everything here was computed in
`analytic/` of the study (scripts in `research/existence/scripts.tar.xz`, listed in section 11; the JSON data are not in the repo). Enumeration/search CPU used by this lens
(both runs): ~25 min, one process at a time; model evaluation (Python) another ~1.5 CPU-h.

## 0. Summary

Definitions: N_semi(X) and N_magic(X) are the expected numbers of *primitive* n x n
additive-multiplicative semi-magic squares (rows and columns) and magic squares (rows, columns
and both diagonals) with magic sum S <= X, counted up to row/column permutations and transposition.
"Expected" means under the heuristic of sections 2-5, calibrated in section 4.

| | n = 5 | n = 6 | n = 7 | n = 8 |
|---|---|---|---|---|
| N_magic(inf) (total) | **5e-6** (x/÷ ~10) | **~2.0** (68% range 0.8-4) | infinite in practice (1e14 by S = 2e4, still growing) | infinite in practice (6e30 by S = 1e4) |
| P(at least one exists) | ~1e-5 | **~0.86** (0.55-0.98) | 1 | 1 |
| X where N_magic = 1 | never | **S ~ 5,900** (quad calibration 6,650; ~2,900 to never within the uncertainty) | **S ~ 430** (390-470) | **S ~ 377-410** (360-425) |
| median of the smallest magic S | - | 4,400 (if one exists: ~3,700) | ~422 | ~375 |
| local growth g = dlnN_magic/dlnX | 2.5 (X=400) -> 0 by 10^4 | 20 (300), 8.5 (500), 5.7 (1000), 2.6 (2000), 1.2 (4000), 0.4 (10^4), 0 (>3x10^4) | 29 (400), 23 (465), 12 (1000), 3.4 (10^4) | 69 (392), 45 (500), 21 (1000), 8 (10^4) |
| known smallest | none known | none known | Miquel 2016, S = 465 (model: E[# with S <= 465] ~ 7) | Boyer 2005, S = 600; Horner 1955, S = 840 (model: ~10^9 / 10^14 with S below these) |
| X where N_semi = 1 | S ~ 380 (found 498 by exhaustive search, see 8.1) | S ~ 294 (known: Morgenstern S+P type, 289) | S ~ 300 (known: Shirakawa S+S type, 310) | S ~ 331 |
| largest entry / S (data) | 0.39-0.50 | 0.38 (0.31-0.55) | 0.47-0.54 | 0.31-0.375 (Horner, Boyer) |

Main conclusions:

1. **The structure.** For a given product P = prod p^{a_p}, multiplicative (semi-)magicness factorises
   over primes into (semi-)magic exponent matrices, counted exactly by Ehrhart (quasi-)polynomials
   H_n(a), M_n(a) (computed exactly, section 2). The additive conditions are handled by a joint
   maximum-entropy / local-limit model for the vector (size, exponents) of the entries (section 3).
   It gives, with no free parameter, the number of vectors per (P, S) to ~10% and, with one smooth
   calibration factor (2 parameters + a per-n offset), the number of semi-magic squares: the first
   search's 758,949 squares over 1,136 P are predicted to 4%, the near-S_min data within a factor 0.64-1.21 in every
   bin of N up to N = 8,000 (0.92 for the 3 sums with N = 6-8k), and **an exhaustive out-of-sample search of all 3,966 sorted-exponent P for 5x5
   squares with S <= 700, run for this lens, found 13 semi-magic squares where the model predicted 12.8**.
2. **P(magic | semi-magic square) = D_n kappa p_pair ~ C / (S^2 prod_p (a_p/n)(1 + a_p/n))**, with
   D_6 = 5,400 diagonal pairs, the size-exponent correlation rho = (1 - R^2)^{-1/2} ~ 2.2 and per-prime
   lattice factors. This is the empirical rule 5400 (rho p_S p_P)^2, derived; it matches the observed
   S, P and SP traversal counts to 0.64-0.97 and gives 1 in 2.8x10^7 (arithmetic mean) near S_min.
3. **Convergence.** For every n >= 3 the expected total number of primitive AM-magic squares is finite
   (Dirichlet-series argument, section 6), but the tail is (log X)^{M_n(1)-1} X^{-n} with M_6(1) = 96,
   so finiteness is invisible in practice except for n = 5 and 6. For n = 6 the calibrated model
   saturates: 50% of N_magic(inf) lies below S ~ 5,800, 90% below S ~ 20,000.
4. **n = 6: ~2 expected magic squares in total, most of them out of reach.** Only 0.33 of the ~2.0 come
   from sums with N <= 8,000 vectors (where the model is verified on squares); 0.37 from N 8-16k,
   0.39 from 16-32k, 0.30 from 32-64k, 0.4 beyond. With the current search time t(N) ~ N^4.7
   (N^2.7 at 8k by the few large-N data), an ideal scheduler over the whole universe collects
   E ~ 0.07-0.11 by 1 CPU-year, 0.16-0.22 by 10, 0.29-0.46 by 100 and 0.48-0.77 by 1,000 CPU-years
   (section 9). Better N-scaling helps most: t ~ N^2.5 beyond N = 2,000 gives 0.13 / 0.33 / 0.63
   at 1 / 10 / 100 CPU-years.
5. **The empirical forecast is probably too pessimistic** (section 9.1): on its own candidate pool
   the analytic model gives E = 0.16-0.20 for S <= 1.5-1.8 S_min against the forecast's E_max = 0.022.
   Two causes, both testable: (a) its squares-per-sum fit under-predicts the observed squares at
   N = 5,000-8,000 by 5.5x (88 predicted, 484 observed; analytic 558), and (b) its diagonal model
   under-predicts the observed SP traversals of our squares by 2.5x (analytic: 0.97x).
6. **n = 5: no.** N_magic(inf) ~ 5x10^-6 (sampling se 5%, calibration uncertainty x/÷ 10): 5x5
   AM-magic squares almost surely do not exist, even though ~9,000 primitive 5x5 semi-magic squares
   are expected (smallest near S ~ 380; we found one with S = 498).
7. **n = 7, 8: consistent with the known examples.** Smallest 7x7 magic predicted at S ~ 430 (Miquel's
   465 is the smallest known: the model expects ~7 primitive 7x7 magic squares with S <= 465, so a
   smaller one probably exists). Smallest 8x8 predicted at S ~ 375-410 (Boyer's 600 is far above: the
   model expects ~10^9 with S <= 600). Because N_magic grows like X^25-X^70 there, these crossing
   points are insensitive to the calibration (a factor 10 moves them by 3-9%).

## 1. Objects and counting conventions

* An n x n **additive-multiplicative semi-magic square** (AM-semi) is a grid of n^2 *distinct*
  positive integers whose n rows and n columns all have the same sum S and the same product P.
  It is **magic** if, in addition, both main diagonals have sum S and product P.
* Squares are counted **up to the 2 (n!)^2 row permutations, column permutations and
  transposition** (this is what the search's canonical forms count; "a semi-magic square" below is
  such a class). A class contains n! "traversals" (permutation patterns) that can be moved onto the
  main diagonal, and D_n distinct unordered pairs (diagonal, anti-diagonal) that can be produced by
  row/column permutations: D_n = n! (n-1)!!/2 for even n, n! n (n-2)!!/2 for odd n, i.e.
  D_5 = 900, **D_6 = 5,400**, D_7 = 264,600, D_8 = 2,116,800. A class is magic iff one of these
  D_n pairs has both diagonals with sum S and product P; the expected number of magic squares
  (counted as magic classes; a class with several magic pairs is rare and counts once to first order) is
  E_magic = E_semi x D_n x p_pair.
* Scaling: if M is AM-magic then so is cM (sum cS, product c^n P). The heuristic below treats every
  (P, S) as an independent random trial, which is right for primitive squares (gcd of entries 1);
  non-primitive copies of a primitive square are deterministic consequences and are discussed in
  section 6.
* P = prod p^{a_p}, written by its exponents (a_2, a_3, a_5, ...); k = number of primes with
  a_p > 0, q = 1 + k, tau(P) = prod (a_p + 1). S_min(P) = smallest possible line sum (smallest sum of n
  distinct divisors with product P); numerically S_min = (1.000-1.05) x n P^{1/n} (median 1.005, checked
  on 12,400 P for n = 6) -- AM-GM gives n P^{1/n} as the lower bound.

## 2. Multiplicative structure: exponent matrices

For a fixed P, a grid of divisors of P has all row and column products equal to P iff, for every
prime p, the matrix of p-exponents is a nonnegative integer matrix with all row and column sums
a_p (a "semi-magic exponent matrix"); it is multiplicatively magic iff these matrices also have both
diagonal sums a_p. Since the primes are independent, the number of ordered multiplicatively
semi-magic grids (repetitions allowed) is

  f_n(P) = prod_p H_n(a_p),      and multiplicatively magic grids: f^M_n(P) = prod_p M_n(a_p),

where H_n(a) is the number of n x n nonnegative integer matrices with line sums a (Ehrhart polynomial
of the Birkhoff polytope B_n, degree (n-1)^2) and M_n(a) the number that also have both diagonal sums a
(an Ehrhart quasi-polynomial of degree (n-1)^2 - 2 = n^2 - 2n - 1). By Birkhoff-von Neumann
every semi-magic exponent matrix is a sum of a_p permutation matrices, so H_n(1) = n! and a
multiplicatively semi-magic square is an entrywise product over a multiset of (prime, permutation) pairs.

**Exact counts** (cell-by-cell transfer-matrix DP, `c/ehr.c`, 52 s total; long double, exact below
2^64). Selected values (T = one diagonal, M = both diagonals):

| n | a=1 H / M | a=2 H / M | a=3 H / M | a=4 H / M | a=5 H / M |
|---|---|---|---|---|---|
| 5 | 120 / 20 | 6,210 / 449 | 153,040 / 6,792 | 2,224,955 / 67,063 | 22,069,251 / 484,419 |
| 6 | 720 / 96 | 202,410 / 14,763 | 20,933,840 / 957,936 | 1.0476e9 / 3.3177e7 | 3.0768e10 / 7.1851e8 |
| 7 | 5,040 / 656 | 9,135,630 / 657,370 | 4.6629e9 / 2.1212e8 | 9.3667e11 / 2.9782e10 | 9.4162e13 / - |
| 8 | 40,320 / 5,568 | 5.4501e8 / 3.9928e7 | 1.5791e12 / 7.3671e10 | 1.4559e15 / - | 5.6930e17 / - |

(H_6(a) for a <= 10, H_5 for a <= 12, H_7 for a <= 7, H_8 for a <= 5 are in `c/ehr.txt`.) With
reciprocity H_n(-a-n) = (-1)^{(n-1)^2} H_n(a) and H_n(-1..-(n-1)) = 0 these values determine the
full polynomials for n = 5, 6 (`ehrpoly.py`); leading coefficients (relative volumes of B_n)
2.255e-7 (n = 5) and 9.455e-13 (n = 6).

* **Canfield-McKay**: H_n(a) ~ C(a+n-1, n-1)^{2n} / C(na+n^2-1, n^2-1) x e^{1/2} is within
  8-20% of the exact counts for every n = 5..8 and every a >= 1 (ratio 1.08-1.20; e.g. n = 6:
  1.11 at a = 1, 1.17 at a = 10).
* **Gaussian (max-entropy) approximation** of H_n(a) (cells iid truncated-geometric with mean a/n,
  2n-1 independent line constraints, local CLT) overestimates by 3.5-5.8x at a = 1, 2.3-2.7x at
  a = 2 and ~1.35-1.5x at a = 10. These "lattice factors" L_semi(a) = H_n(a)/Gauss_n(a) (and
  L_T, L_M for one/two diagonals) are applied per prime below; for a beyond the exact range they are
  extrapolated with log L = L_inf - A/(a + B).
* **Diagonal fraction**: M_n(a)/H_n(a) = 0.13-0.14 at a = 1 for all n = 5..8 and decays like
  c/a^2 (a^2 M/H = 0.13, 0.29, 0.41, 0.51, 0.58 for n = 6, a = 1..5): for one prime the chance that a
  random semi-magic exponent matrix has both diagonal sums right is ~ (2 pi Var)^{-1} ~ a^{-2}.
* Dirichlet series: sum_P f_n(P) P^{-s} = prod_p (sum_a H_n(a) p^{-as}) = zeta(s)^{n!} G(s) with
  G holomorphic for Re s > 1/2, so sum_{P <= Y} f_n(P) ~ c_n Y (log Y)^{n!-1}/(n!-1)!; for magic
  exponent matrices the exponent is M_n(1) - 1 (95 for n = 6). These enormous log powers are what
  makes the problem "pre-asymptotic" in every practical range (section 6).

## 3. Joint maximum-entropy model for one (P, S)

**Single-cell measure.** For fixed (P, S) put on the divisors d of P the tilted (max-entropy) measure

  mu(d) proportional to exp( -beta d/(S/n) + sum_p lambda_p e_p(d) ),

with (beta, lambda) chosen so that E_mu[d] = S/n and E_mu[e_p(d)] = a_p/n. (`maxent.py`, Newton on
the convex dual.) Let H = entropy of mu, C = the q x q covariance of v(d) = (d, e_2(d), e_3(d), ...)
under mu.

**Vectors** (n-sets of distinct divisors with sum S and product P), local CLT for the q-vector sum
of n iid cells, distinctness and ordering:

  N(P,S) ~ e^{nH} (2 pi)^{-q/2} det(nC)^{-1/2} x prod_p L_vec(a_p) x P_dist(n) / n!

**Ordered grids** with all 2n line sums S and products P: the 2n line constraints on the q-vector are
2n-1 independent (row totals = column totals), and the local CLT gives

  G(P,S) ~ e^{n^2 H} / [ (2 pi)^{(2n-1)q/2} det(C)^{(2n-1)/2} n^{(n-1)q} ] x prod_p L_semi(a_p)

(the n^{(n-1)q} is the determinant of the line-sum lattice). For the exponent coordinates alone this
is exactly the Gaussian approximation of prod_p H_n(a_p) corrected by the exact lattice factors, so
the multiplicative part is exact up to the coupling with the size constraint.

**Distinct entries and symmetry:**

  E_semi(P,S) = G(P,S) x P_dist(n^2) / (2 (n!)^2),  P_dist(m) = m! e_m(mu)

(probability that m iid mu-draws are distinct; computed by a saddle point, `esym2.py`).

**Diagonals.** Conditional on all rows and columns, a given traversal has q-vector deviation with
covariance (n-1) C, so

  p_S  = [2 pi (n-1) Var_mu(d)]^{-1/2}
  p_P  = prod_p L_T(a_p) x [(2 pi (n-1))^k det C_ee]^{-1/2}
  p_SP = prod_p L_T(a_p) x [(2 pi (n-1))^q det C]^{-1/2} = rho p_S p_P,  rho = (det C_dd det C_ee / det C)^{1/2} = (1 - R^2)^{-1/2}

where R^2 is the squared multiple correlation of the size d with the exponent vector (sizes and
exponents are strongly correlated because log d = sum e_p log p; the model gives rho = 1.8-3.6).
For the two diagonals of one arrangement the deviations have covariance C (x) [[n-1, c],[c, n-1]] with
c = -1 for even n (disjoint diagonals) and 0 for odd n, so

  p_pair = prod_p L_M(a_p) / [ (2 pi)^q det C ((n-1)^2 - c^2)^{q/2} ]

and

  **P(magic | semi-magic square) = D_n kappa p_pair  ~  D_n kappa (rho p_S p_P)^2 x [(n-1)^2/((n-1)^2-c^2)]^{q/2}**

(kappa is the empirical calibration of section 4). This is the empirical rule
P(magic | square) = 5400 (rho p_S p_P)^2 of the forecast, derived. Its scaling: Var_mu(d) =
c_S^2 (S/n)^2 with c_S ~ 0.3-0.6 and Var_mu(e_p) ~ (a_p/n)(1 + a_p/n) (truncated geometric), so

  P(magic | semi) ~ C_n kappa prod_p l_M(a_p) / [ S^2 (1-R^2) prod_p (a_p/n)(1+a_p/n) ]

-- the "C / (S^2 prod a_p^2)" form for large a_p (and ~1/(S^2 prod a_p) for small a_p, where the lattice
factors l_M take over). Equivalently p_S ~ n/(c_S sqrt(2 pi (n-1)) S) (the report's "1/S" rule) and
p_P ~ prod_p 1/sqrt(2 pi (n-1) Var e_p) (close to the report's "1/tau(P)" rule only when all a_p are
small). Within one P the effective decline is much steeper than S^-2 near S_min, because the entry
distribution broadens quickly as S leaves S_min (Var_mu(d)/(S/n)^2 grows): the model gives
d log P(magic|semi)/d log S ~ -8 between 1.3 and 2 S_min, matching the forecast's S^-7..-9.

**Worked numbers** (calibrated model of section 4):

| case | n | P | S (S_min) | vectors N | E_semi per sum | P(magic \| semi) |
|---|---|---|---|---|---|---|
| Morgenstern 0-type semi-magic | 6 | 10 4 3 2 | 327 (171) | 446 (obs 452) | 2.3e-4 | 1 in 5.0e6 |
| Morgenstern S+P semi-magic | 6 | 11 4 3 1 1 | 289 (206) | 450 | 9.2e-4 | 1 in 3.2e6 |
| best P region | 6 | 13 7 4 3 2 1 | 2820 (2562) | 7,540 | 414 | 1 in 6.1e8 |
| Miquel magic | 7 | 10 7 3 2 1 | 465 (279) | 8,445 | 3.9e3 | 1 in 1.35e6 |
| Horner magic | 8 | 7 8 3 1 0 1 1 1 1 1 | 840 (658) | 50,447 | 8.6e5 | 1 in 6.0e7 |
## 4. Calibration against the data

Data used (all exhaustive per searched sum): n = 6 near S_min (19,834 sums, 6,498 squares: the
forecast's units plus the verifiers' runs), n = 6 low-S survey (`runs/n6low.jsonl`: 30 P, every sum from
S_min to 400, 4,241 sums, 124 squares, smallest S = 304), n = 5 survey (`runs/n5_survey.jsonl`: 200 P with S_min 1,830-2,700, every
sum up to S = 3,000 within 20 s per P: 109,193 sums, 20 squares; a weighted subsample of 3,020 sums enters the fit), n = 7 (53 sums of 2 P near S_min, 15 squares), and the first search's per-P
totals (`stats/stats_short.txt`: 1,136 P, 758,949 squares, S up to max_S, 519 SP traversals).

### 4.1 Vectors (one line): the max-entropy model is accurate

log(N_obs / N_pred) (sums with N >= 20):

| data | mean | sd | by N (n = 6 near S_min) |
|---|---|---|---|
| n = 6 near S_min | +0.06 | 0.095 | +0.27 (N 20-100), +0.14 (100-300), +0.05 (300-1000), +0.01 (1000-3000), -0.03 (>3000) |
| n = 6 low S | +0.07 | 0.08 | |
| n = 5 survey | -0.09 | 0.06 | |
| n = 7 | +0.13 | 0.02 | |

So the single-line count (local CLT in dimension q = 1 + k with exact per-prime lattice factors and
the iid distinctness factor) is right to ~10% (1% for N > 1000), e.g. (10 4 3 2, 327): 446 predicted,
452 actual; (11 4 3 1 1, 289): 450 vs 484.

### 4.2 Semi-magic squares: the raw model needs a factor, and the factor is smooth

The raw square model (section 3, no free parameters) versus observed squares, log(obs/pred):

| data | squares | pred.-weighted ld = log P_dist(n^2) | raw (Gauss + lattice) | raw + 1st-order Edgeworth |
|---|---|---|---|---|
| n = 5 survey | 20 | -2.2 | +1.46 | +2.34 |
| n = 6 low S | 124 | -10.2 | -1.78 | -0.98 |
| n = 6 near S_min | 6,498 | -8.7 | -1.78 | -0.86 |
| n = 7 | 15 | -22.6 | -6.28 | -5.69 |

Why it is off: the local CLT is used far from its asymptotic regime. Each line has only n = 5..8
terms, the cell law is strongly skewed (sizes truncated near 2-3 S/n, exponents with mean a_p/n < 1)
and 2n-1 constraints are imposed at once. The size of the corrections shows it: the per-prime
lattice factors contribute -3.8 to -7.0 to log E, and the first-order Edgeworth correction alone is
-4.9 to -7.6 (so higher orders are not negligible). Second, the iid distinctness factor misses the
correlations that the line constraints induce: the residual falls with ld = log P_dist(n^2).

A Poisson regression of the observed counts on log E_raw (offset) with covariates gives
(sampling-weighted; deviance in brackets):

| correction log K | coefficients (se) | obs/fit: n5, n6 low, n6 near, n7 |
|---|---|---|
| c + b ld | b = 0.209 | 6.6, 1.38, 1.00, 0.21 (fails across n) |
| c + b ld + e ld^2 | | 14.2, 1.29, 0.99, 0.99 (fails at n = 5) |
| **c + b ld + d (n - 6)** | **c = -0.597 (0.078), b = 0.1365 (0.009), d = -2.455 (0.18)** | **0.89, 1.23, 1.00, 0.86** |
| c_n + b ld (free c_5, c_7) | c_5 - c_6 = +2.35 (0.23), c_7 - c_6 = -2.62 (0.29), b = 0.135 | 1.00, 1.22, 1.00, 1.00 |

The free per-n offsets are linear in n (+2.35, 0, -2.62), so the adopted correction is

  **log K_n = -0.597 + 0.1365 ld - 2.455 (n - 6)**,   E_semi(P,S) = K_n x E_raw(P,S).

Equivalently, squares behave as if P_dist entered with power 1.14 (stronger than iid), with an extra
factor e^{-2.46} per unit of n. Within n = 6 a per-P fixed-effects fit gives b = 0.107 (0.012), and an
earlier within-n = 6 quadratic form log K = -1.087 - 0.0091 ld^2 (+ ad hoc offsets
+2.6 / -0.58 / -1.2 for n = 5 / 7 / 8) fits equally well; it is carried along as the alternative
"quad" calibration below (it differs mainly for n = 8).

**Out-of-sample check on the first search** (1,136 P, 758,949 squares, mostly at S well above the
calibration data): predicted/observed totals = 1/1.04 (adopted), 1/1.03 (quad), but 1/0.23 for the raw
model; by S_min bins (0-250, 250-500, ..., 1250-1500) the calibrated ratio stays within 0.96-1.10.
Per-P scatter: sd of log(obs/pred) over P with >= 10 squares is 0.28, of which 0.17 is Poisson, i.e.
~25% genuine per-P scatter. The n = 6 low-S survey (not used in the within-n = 6 quad fit)
is reproduced to 1.23x.

### 4.3 Diagonals (traversals)

Observed / model numbers of traversals with the magic sum (S), product (P), both (SP):

| data | squares | S | P | SP |
|---|---|---|---|---|
| n = 6 near S_min (per square) | 6,598 | 8,285 / 9,631 = 0.86 | 2,037 / 3,170 = 0.64 | 14 / 14.4 = 0.97 |
| n = 6 first search (integrated) | 758,949 | 0.89 | 0.67 | 519 events: 0.68 |
| n = 7 near S_min | 12 | 216 / 230 = 0.94 | 14 / 19.3 = 0.73 | 0 / 0.13 |
| n = 5 | 20 | 2 / 2.0 | 1 / 0.3 | 0 / 0.003 |

The model's own S-P correlation (rho = (1 - R^2)^{-1/2} ~ 2.2 for n = 6) reproduces the empirical
rho ~ 2.15-2.6, so the SP rate is right to within 0.68-0.97. For two diagonals the forecast found
S+S pairs 1.1-1.2x more often than independence. Adopted:

  **kappa = kappa_SP^2 x 1.1 = 0.72^2 x 1.1 = 0.57**  (plausible range 0.46-1.2),

so P(magic | semi) = D_n x 0.57 x p_pair. Near S_min for n = 6 this gives, averaged over the observed
squares, 1 in 2.8e7 (arithmetic mean; geometric mean 1 in 4.1e7), inside the forecast's "1 in 20-90
million".
## 5. Summing over P and S

N_semi(X) = sum_P sum_{S <= X} E_semi(P,S),  N_magic(X) = sum_P sum_{S <= X} E_magic(P,S).

**Which (P, S) matter.** For fixed P the semi-magic yield per sum rises steeply from S_min (where
P_dist(n^2) is tiny: few divisors are available), peaks at 1.3-2 S_min and falls; the magic yield per
sum peaks at 1.1-1.3 S_min because P(magic|semi) falls ~S^-8 (section 3); both are negligible beyond
4 S_min (the model is integrated over S_min < S <= min(X, 4 S_min) on a 22-26 point log grid; the
log-linear quadrature underestimates relative to a 48-point grid by 1-4% for n = 5, 6 and by up to
28% (n = 7) / 2x (n = 8) at the smallest X; corrected in section 7). Across P,
yields need tau(P) large compared to n^2 (for n = 6 the yield is < 1e-9 per P for tau < 1,000 and peaks
for tau ~ 4,000-30,000) and many divisors near P^{1/n}.

**Universe and sampling** (`gs3.py`). Candidates: all exponent vectors over the primes 2..31 that are
nonincreasing over the primes used (one zero "gap" allowed after the 4th prime), 3 <= k <= 10 primes,
at most 6 exponents equal to 1, tau >= tau_min (1,000 / 300 / 2,000 / 4,000 for n = 6 / 5 / 7 / 8) and
scale n P^{1/n} <= X_hi (50,000 / 200,000 / 20,000 / 10,000): 447k / 389k / 533k / 614k P. They are
stratified by (half-octave of scale, floor log2 tau); in each scale bin a proportional first phase
(40-50 P) is followed by a Neyman-allocated second phase (80-110 P, allocation by the first-phase RMS of
the magic yield). Each sampled P is evaluated on its own S grid (exact S_min from `bin/enumerate` when
P < 2^62, else 1.005 n P^{1/n}). Totals are stratified (Horvitz-Thompson) estimates; standard errors
include the finite-population correction.

**Exponent permutations.** The sorted universe misses P whose exponents are not nonincreasing in p
(e.g. 13 7 4 3 1 2, or Horner's 7 8 3 1 0 1 1 1 1 1). Direct evaluation of all permutations of three
P shows they add 0.4-1.7x of the sorted P's yield, almost all from permutations whose scale is within
1.2x of the sorted one (13 7 4 3 2 1: x2.66; 10 4 2 2 1 1: x1.41; 12 6 4 2 2 1: x1.89 for the magic
yield). A 10-12% subsample of the sampled P is evaluated with all its permutations of scale <= 1.3x,
and the totals are multiplied by the resulting ratio estimator (per X).
## 6. Growth with X, convergence of the total

### 6.1 Why the total converges (heuristically) for every n >= 3

Drop distinctness and the requirement that entries be balanced (both only lower the count). A grid
with product P chosen uniformly among the f_n(P) = prod_p H_n(a_p) multiplicatively semi-magic grids
has entries of geometric mean P^{1/n}; its line sums are spread over a range >~ n P^{1/n}, so the 2n
line sums all coincide (for some S) with probability <~ c (n P^{1/n})^{-(2n-2)}. Hence

  E[# primitive AM-semi-magic squares of all sizes] <~ c sum_P f_n(P) P^{-(2-2/n)} = c D_n(2 - 2/n),

and for magic squares (2n+1 additive equalities, multiplicative diagonals: f^M_n = prod_p M_n(a_p))

  E[# primitive AM-magic squares] <~ c sum_P f^M_n(P) P^{-2} = c D^M_n(2).

Both Dirichlet series have abscissa of convergence 1 (Euler factors 1 + n! p^{-s} + ..., resp.
1 + M_n(1) p^{-s} + ...; D_n(s) = zeta(s)^{n!} G(s)), and 2 - 2/n > 1 for n >= 3. **So the expected
number of primitive AM-(semi-)magic squares of each order n >= 3 is finite; the cumulative counts
N(X) converge as X -> infinity.** (Counting non-primitive squares as well, any primitive magic square M
generates the multiples cM, so the full count would grow linearly, N_all(X) ~ sum_M X/S_M; these are not
independent events and are excluded throughout.)

Where the bound's mass sits (Euler products evaluated with Canfield-McKay H_n and
M_n/H_n = 1.18/(a+2.1)^2, which fits the exact ratios for all n = 5..8 to ~5%): the tail of the
convergent sum is negligible only beyond

| n | semi: log10 D, mass at log10 S | magic: log10 D, mass at log10 S |
|---|---|---|
| 5 | 21.5, ~7 | 6.5, ~2 |
| 6 | 51.1, ~14 | 20.0, ~4.1 |
| 7 | 117, ~31 | 49.4, ~8.4 |
| 8 | 270, ~75 | 112, ~19 |

(the constants c are not determined; this only locates where the bound peaks). For semi-magic squares
the bound peaks at astronomically large S, so in any searchable range N_semi(X) keeps growing
(pre-asymptotic regime, governed by the huge log power (log Y)^{n!-1} of the Dirichlet series). For
**magic** squares the extra two additive and k multiplicative diagonal constraints (P^{-2} instead of
P^{-5/3} for n = 6, and M_n(a) ~ H_n(a)/a^2) move the peak down to S ~ 10^4 for n = 6 -- the same place
where the calibrated model below finds the magic yield per octave peaking. For n = 7, 8 the magic peak
of the crude bound lies far above the searched range, i.e. the magic counts are still in their
growing phase at S ~ 10^3; for n = 5 the crude bound is uninformative (it ignores that 25 distinct
balanced divisors need tau(P) in the hundreds or more, i.e. S >~ 300-500; the smallest 5x5 square we
found has S = 498 and tau = 1,680, section 8.1).

### 6.2 Asymptotic form of the convergence (formal)

Partial sums of a Dirichlet series with a pole of order m at s = 1, evaluated at s = alpha > 1:
sum_{P > Y} f(P) P^{-alpha} ~ c Y^{1-alpha} (log Y)^{m-1}/((alpha - 1)(m-1)!). With Y ~ (X/n)^n this gives
the formal tails N(inf) - N(X) ~ X^{-(n-2)} (log X)^{n!-1} (semi) and ~ X^{-n} (log X)^{M_n(1)-1}
(magic), valid only once log X is large compared to the log powers (n!, M_n(1) = 96 for n = 6), i.e.
never in practice. In the practical range the relevant description is the effective local exponent
g(X) = d log N / d log X computed from the calibrated model (section 7).

## 7. Results: N_semi(X), N_magic(X) for n = 5, 6, 7, 8

Stratified estimates of section 5 (`gs3_n*.json`, `final_tables.py`). They include the
exponent-permutation multiplier, and the quadrature-grid correction measured by re-evaluating 20 P per n
on a 48-point S grid (`gridbias.py`). That correction is x1.01-1.04 for n = 5, 6. For n = 7 it is
x1.28 at X ~ 390, falling to 1.02 at the top of the grid; for n = 8 it is x2.0 at X ~ 325, falling to 1.06.
"se" is the sampling standard error only (section 10 has the full uncertainty). "quad" is the alternative
calibration of section 4.2. g = d ln N / d ln X over X/1.15..1.15X.

### n = 6 (universe: 447k sorted-exponent P with scale 6 P^{1/6} <= 50,000 and tau >= 1,000; 2,809 evaluated; permutation multiplier 1.0-1.77)

| X | N_semi(S<=X) | g_semi | N_magic(S<=X) (se) | g_magic | N_magic quad | P(>=1 magic with S<=X) |
|---|---|---|---|---|---|---|
| 289 | 0.69 | 23 | 1.1e-7 | 21 | 7.1e-8 | 1e-7 |
| 327 | 9.1 | 19 | 1.2e-6 | 16 | 8.4e-7 | 1e-6 |
| 400 | 211 | 14 | 1.5e-5 (2.7e-6) | 11 | 1.2e-5 | 1.5e-5 |
| 500 | 3.4e3 | 11.5 | 1.25e-4 (2.2e-5) | 8.5 | 1.1e-4 | 1.3e-4 |
| 700 | 1.1e5 | 9.6 | 0.0014 (0.0002) | 6.3 | 0.0014 | 0.0014 |
| 1,000 | 3.1e6 | 9.2 | 0.012 (0.002) | 5.7 | 0.012 | 0.012 |
| 1,500 | 6.7e7 | 6.2 | 0.071 (0.009) | 3.5 | 0.070 | 0.069 |
| 2,000 | 3.7e8 | 5.6 | 0.17 (0.02) | 2.6 | 0.16 | 0.15 |
| 3,000 | 3.3e9 | 5.1 | 0.42 (0.04) | 1.7 | 0.40 | 0.34 |
| 4,000 | 1.3e10 | 4.2 | 0.63 (0.06) | 1.2 | 0.58 | 0.47 |
| 6,000 | 5.6e10 | 3.5 | 1.02 (0.08) | 1.0 | 0.92 | 0.64 |
| 10,000 | 3.0e11 | 3.1 | 1.44 (0.11) | 0.4 | 1.27 | 0.76 |
| 20,000 | 2.1e12 | 2.7 | 1.76 (0.12) | 0.2 | 1.53 | 0.83 |
| 50,000 | 1.7e13 | 2.3 | 1.96 (0.12) | 0.0 | 1.68 | 0.86 |

Crossings (lin / quad): N_magic = 0.1 at S ~ 1,660 / 1,670; ln 2 at 4,370 / 4,660; 1 at 5,880 / 6,650;
never 3. Semi-magic: N_semi = 0.1 / ln 2 / 1 / 3 at S = 267 / 289 / 294 / 309. The tail beyond
the universe adds ~0.03-0.05 to N_magic: the per-half-octave yield falls ~0.6x per half-octave of
scale beyond 16,000, and P with tau < 1,000 contribute 0 (checked on a 30% sample of the 3,330 such P with scale <= 6,000, `lowtau.py`). Up to the
calibration, the universe restrictions (k <= 10 primes up to 31, at most 6 exponents equal to 1) are
immaterial: P with 9-10 primes carry 1.2% of the total.

Where the n = 6 magic yield sits (sorted-exponent P, before the permutation multiplier):
- **scale** 6 P^{1/6}: per half-octave 0.08 (1,000), 0.11 (1,450), 0.13 (2,050), 0.10 (2,900),
  0.22 (4,100), 0.13 (5,800), 0.10 (8,200), 0.07 (11,600), 0.06 (16,400), 0.04 (23,000), 0.02 (32,800);
- **number of primes** k = 5: 15%, 6: 44%, 7: 28%, 8: 11%;
- **tau(P)**: 2^12-2^13: 13%, 2^13-2^14: 25%, 2^14-2^15: 30%, 2^15-2^16: 19%, 2^16-2^17: 9%;
- **S/S_min**: 1.05-1.4 (P(magic | semi) falls ~S^-8 within a P);
- **vectors per sum N** (per-sum accounting, `frontier_an.py`): N < 2k: 0.023, 2-4k: 0.087, 4-8k:
  0.22, 8-16k: 0.37, 16-32k: 0.39, 32-64k: 0.30, 64-128k: 0.22, 128k-1M: 0.20, > 1M: 0.01
  (total 1.83 for S <= 2 S_min). Mean P(magic | semi) falls from 1 in 2x10^7 (N < 1k) to 1 in 1.4x10^9
  (4-8k) and 1 in 3x10^10 (16-32k).

Top single P (sorted, magic yield over all S): 15 6 4 3 2 1 (S_min 2,687), 11 8 4 2 2 1 (1,766),
13 7 4 3 1 (1,122), 16 8 4 3 2 (2,838), 15 6 4 3 2 (1,754); each ~3x10^-4. **No single P holds more than
0.02% of the total.** The expectation is spread over ~10^5 P, which is why it cannot be
collected cheaply.

### n = 5 (389k P, scale <= 200,000, tau >= 300; 2,753 evaluated)

| X | 300 | 400 | 500 | 700 | 1,000 | 2,000 | 5,000 | 10^4 | 2x10^4 | 2x10^5 |
|---|---|---|---|---|---|---|---|---|---|---|
| N_semi | 0.26 | 1.3 | 3.5 | 13.6 | 56 | 330 | 1,480 | 3,660 | 5,820 | 9,190 |
| N_magic | 9.7e-8 | 2.3e-7 | 3.6e-7 | 8.5e-7 | 1.5e-6 | 2.7e-6 | 4.0e-6 | 4.8e-6 | 5.0e-6 | 5.1e-6 |

N_semi = 1 at S ~ 380 (median smallest 355; 90% range 266-480). N_magic converges to 5.1x10^-6
(se 0.24x10^-6; quad 5.6x10^-6). The magic yield peaks at scale 500-4,000 (k = 5-7, tau 2^10-2^13). P with
tau < 300 contribute nothing (checked).

### n = 7 (533k P, scale <= 20,000, tau >= 2,000; 1,715 evaluated) and P with tau 300-2,000

| X | 300 | 310 | 350 | 400 | 435 | 465 | 500 | 600 | 1,000 | 2,000 | 10^4 | 2x10^4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N_semi | 1.3 | 8.0 | 2.2e3 | 3.1e5 | 4.0e6 | 2.2e7 | 1.4e8 | 1.0e10 | 3.5e13 | 8.5e17 | 2.9e24 | 3.3e26 |
| N_magic (tau >= 2000) | 2e-6 | 1.1e-5 | 0.0017 | 0.15 | 1.24 | 6.0 | 32 | 1.2e3 | 8.3e5 | 1.4e9 | 1.7e13 | 1.1e14 |
| + P with tau 300-2,000 | | | | +0.05 | +0.17 | +0.43 | +1.0 | | | | | |
| N_magic quad | 5e-7 | 3.7e-6 | 0.0015 | 0.23 | 2.7 | 14.5 | 88 | 4.3e3 | 4.4e6 | 9.3e9 | 1.0e14 | 6.0e14 |

(The low-tau supplement was evaluated with a 40-point grid on all 5,022 such P with scale <= 500; tau < 300
gives 0.) Combined: **N_magic = 1 at S ~ 429 (quad ~418), median smallest S ~ 422, 90% range 391-449**
(N = 0.1 and 3); N_semi = 1 at S ~ 300. Growth: g_magic = 29 at 400, 23 at 465, 12 at 1,000, 7 at 2,000,
3.4 at 10^4, 2.4 at 2x10^4. N_magic never stops growing in the searchable range. The crude bound of
section 6.1 puts the turnover near S ~ 10^8.

### n = 8 (614k P, scale <= 10,000, tau >= 4,000; 1,453 evaluated) and P with tau 500-4,000

| X | 330 | 360 | 375 | 392 | 420 | 450 | 500 | 600 | 840 | 1,000 | 2,000 | 10^4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N_semi | 0.62 | 5.9e3 | 1.1e5 | 1.9e6 | 1.2e8 | 6.5e9 | 1.8e12 | 1.0e16 | 1.5e21 | 2.6e23 | 5.4e30 | 1.6e42 |
| N_magic (tau >= 4000) | 1.4e-6 | 0.0065 | 0.10 | 1.7 | 114 | 6.4e3 | 8.4e5 | 1.6e9 | 1.9e14 | 9.8e15 | 1.3e22 | 6.2e30 |
| + P with tau 500-4,000 | | +0.06 | +0.62 | +6.0 | +129 | +1.6e3 | | | | | | |
| N_magic quad (tau >= 4000) | 9e-11 | 2.5e-5 | 8.6e-4 | 0.027 | 3.9 | 433 | 3.1e5 | 3.7e9 | 1.4e15 | 1.4e17 | 4.6e23 | 2.4e32 |

The tau_min = 4,000 universe misses most of the yield at the smallest X: 9,937 P with tau 500-4,000 add
3.5x the main total at X = 392. P with tau < 500 add 0. Combined (lin): **N_magic = 1 at S ~ 377,
median smallest ~375, 90% range 362-385**. With the quad calibration the crossing is ~410. N_semi = 1 at S ~ 331.
Growth: g_magic = 69 at 392, 45 at 500, 21 at 1,000, 16 at 2,000, 8 at 10^4.

### 7.1 Largest entry instead of S

All entries of a line are distinct and positive with sum S, so the largest entry is at least S/n.
The observed ratio r = (largest entry)/S (`maxent_ratio.py`):

| data | r: 5% / median / 95% |
|---|---|
| n = 6, 6,498 squares near S_min | 0.31 / 0.39 / 0.55; median 0.30 at S/S_min 1.0-1.1, 0.36 at 1.1-1.2, 0.45 at 1.3-1.4, 0.55 at 1.6-1.7 |
| n = 6 low-S survey (124 squares) | 0.38 / 0.47 / 0.58 |
| n = 5 survey (20 squares, S > 2,200) | 0.32 / 0.39 / 0.46 |
| n = 5 low S (13 squares of 8.1, S <= 700) | 0.37 / 0.50 / 0.61 |
| n = 7 (12 squares; Miquel's magic square 252/465 = 0.54) | 0.42 / 0.47 / 0.52 |
| n = 8 magic (Horner 261/840 = 0.31, Boyer 225/600 = 0.375) | |

The magic yield sits at S/S_min ~ 1.05-1.4, where r ~ 0.33-0.45. So N_magic(max entry <= Y) ~
N_magic(S <= Y / r) with r ~ 0.38 for n = 6. This is a change of scale only: the growth exponents g
are unchanged. Expected smallest largest entry: n = 6 ~ 2,200 (where N = 1; 1,700 at the median), n = 7 ~ 200
(Miquel's: 252), n = 8 ~ 125-150 (Boyer's: 225).

## 8. Checks against data that were not used in the calibration

### 8.1 New: exhaustive 5x5 search at low S (out of sample)

The n = 5 calibration used only the survey at S = 1,834-3,000 (20 squares). The model predicts 5x5 AM
semi-magic squares already at S ~ 380-500, far below any searched range. To test this I ran
`bin/msearch --vec-size 5 --max-sum 700` on **all 3,966 sorted-exponent P with scale 5 P^{1/5} <= 700**
(the gen2 universe with tau >= 100: k >= 3 primes up to 31, at most one zero gap). That is 637,438 sums
with N <= 451, 5 s of search CPU (`runs5/`). Every 5x5 square with S <= 700 whose product lies in that class (nonincreasing exponents up to one gap,
tau >= 100, at most 6 exponents equal to 1) is therefore found.

| S <= | 400 | 500 | 600 | 700 |
|---|---|---|---|---|
| predicted (same P set, lin calibration) | 0.99 | 3.07 | 6.91 | 12.8 |
| observed | 0 | 1 | 9 | **13** |

Poisson-consistent at every threshold (P(<= 1 | 3.07) = 0.19). This tests the n = 5 offset of the
calibration and the extrapolation in S by a factor of 3-4, out of sample. The smallest square found has
S = 498 and P = 2^6 3^4 5^2 7 11 13 19 (largest entry 252):

    133 195  40 108  22
    135  20  78 209  56
     88  57 252  26  75
     12 154  95 120 117
    130  72  33  35 228

Among the 13 squares (13 x 120 traversals) there are 17 S-type, 4 P-type and 1 SP-type traversals.
The SP-type one is at S = 570, P = 2^5 3^5 5^3 7^2 11 13, best diagonal score 7. No magic square,
as expected (model: 7x10^-7).

### 8.2 Smallest known examples

| n | known | model |
|---|---|---|
| 6 semi | Morgenstern 2007, S+P type, S = 289; 0-type S = 327; our low-S survey of 30 P: smallest 304 | N_semi(289) = 0.69, N_semi(327) = 9; smallest expected at 267-309 (90%), median 289 |
| 7 semi | Shirakawa 2010, S+S type, S = 310; our P 10 4 2 2 1 1 run: smallest 310 | N_semi = 1 at 300, 90% range 289-304; N_semi(310) = 8 |
| 7 magic | Miquel 2016, S = 465 (P = 2^10 3^7 5^3 7^2 11, tau = 2,112) | smallest expected ~422 (391-449); E[# with S <= 465] ~ 6.4 (quad 15) |
| 8 magic | Horner 1955 S = 840; Boyer 2005 S = 600 (largest entry 225) | smallest expected ~375 (lin) - 410 (quad); E[# with S <= 600] ~ 10^9 |

Every known record lies at or above the predicted smallest value, never below it. The near-records for
semi-magic squares (289, 310) are where the model puts the smallest ones. A record found far below the
predicted crossing would have refuted the model; none exists. The magic records for n = 7, 8 were found
by non-exhaustive (and, for n = 8, constructive) methods, so they are expected to sit above the true
minimum.

### 8.3 Large N (where the n = 6 total comes from)

Observed squares in all exhaustively searched n = 6 sums, by N (`nbins.py`, `schedobs.py`):

| N | sums | observed | analytic (calibrated) | scheduler fit (model_6.json) |
|---|---|---|---|---|
| < 1,000 | 11,571 | 106 | 166 | 240 |
| 1,000-2,000 | 7,994 | 4,182 | 4,328 | 6,355 |
| 2,000-3,000 | 260 | 1,625 | 1,343 | 1,831 |
| 3,000-5,000 | 5 | 101 | 122 | 206 |
| 5,000-8,000 | 4 | 484 | 558 | 88 |

Individual large sums: (12 6 3 2 1 1, S = 988, N = 6,671): 245 observed, 298 analytic, 12 scheduler fit;
(12 8 4 2 1, 1,302, 5,995): 87 / 126 / 11; (14 7 5 3, 2,050, 6,206): 24 / 22 / 2;
(14 7 4 4 1 0 0 1, 3,231, 7,901): 128 / 112 / 64. The analytic model has no free parameter in N, and it
holds up to the largest N searched. The first search's P, binned by the N of their predicted squares, give
obs/pred 0.84 / 1.06 / 1.03 / 1.04 for N < 1k / 1-2k / 2-3k / 3-5k. No sum with N > 8,000 has been
searched exhaustively, so 82% of the n = 6 total rests on extrapolating the model in N (section 10).

## 9. Comparison with the empirical forecast; what speed and N-complexity buy

### 9.1 Same candidate pool

The forecast (`forecast/forecast.md`) uses the scheduler's fitted model on its pool of ~10^4 P
(scheduler.gen_candidates: 2^8-18 3^3-9 5^2-6 7^1-4 nonincreasing, at most two of 11-19, tau
800-12,000). It finds E_max = 0.022 for all unsearched sums up to 1.5-1.8 S_min (central 0.015;
optimistic scenarios up to 0.2). The analytic model on the same pool (`pool_cmp.py`, 1,209 of 10,075 P, stratified):

| sums | E_magic | E_semi | mean P(magic \| semi) |
|---|---|---|---|
| S <= 1.3 S_min | 0.083 (se 0.005) | 1.3e8 | 1 in 1.6e9 |
| S <= 1.5 S_min | 0.155 (0.009) | 6.2e8 | 1 in 4.0e9 |
| S <= 1.8 S_min | 0.197 (0.010) | 1.5e9 | 1 in 7.7e9 |
| S <= 4 S_min | 0.209 (0.011) | 2.5e9 | 1 in 1.2e10 |

Split by N (S <= 1.8 S_min, approximate): N < 3k: 0.014; 3-8k: 0.052; 8-20k: 0.055; 20-60k: 0.025. So
**even the part with N <= 8,000 (0.066) is 3x the forecast's total.** Point by point on 200 sums of the
pool (`schedcmp.py`), analytic/scheduler is as follows:
- squares per sum: x1.1 (median) for N < 3k, x1.6 for 3-6k, x27 for 6-10k, x500 for 10-20k;
- P(magic | semi): x3.7 for N < 3k, x1.9 for 3-6k, x1.1-1.35 above.

Which is right:
- **Squares at large N**: table 8.3. On the only data beyond N = 5,000 (4 sums, 484 squares) the
  scheduler fit is 5.5x low and the analytic model 1.15x high. The scheduler's square count is a
  polynomial in log N fitted mostly on N <= 2,500. Beyond that it bends over, and the data do not.
- **Diagonals**: on our 5,821 squares the forecast's own check found SP traversals 2.47x above its model
  (14 vs 5.67), with a 90% CI of 1.5-3.9x. The analytic diagonal model gives 14.4 (0.97x). The analytic
  P(magic | semi) is close to the forecast's scenario (b) ("SP rate of our 5821 squares"), which raised
  E_max to 0.137.
- **Time model**: the forecast's time model (quadratic in log N, local exponent 5.2 at N = 2k and 7.3 at 8k)
  predicts 838 s and 2,111 s for the sums with N = 6,671 and 7,901. With the current build they took 291 s
  and 348 s. A power law fitted to 20,039 current-build sums (N = 800-7,901) gives t = e^{-34.21} N^{4.72}
  (sd 0.47 in log t). It is driven by N < 2,500. The 6 sums with N > 4,000 lie 4-10x below it. A
  quadratic fit bends to a local exponent of 2.7 at N = 8k.

### 9.2 E(CPU) under the analytic model (whole universe; ideal greedy scheduler)

`frontier_an.py` / `frontier_eval.py` work as follows. Every sum S in (S_min, 2 S_min] of every P of the n = 6
universe (874 P subsampled from the stratified sample; x1.77 for exponent permutations) gets E_magic and
a search time t(N). Sums are taken in decreasing E/t order. N here is the model's vector count / 1.05.
The table is cumulative E_magic after C CPU:

| time model t(N) | 100 h | 1000 h | 1 CPU-yr | 10 CPU-yr | 100 CPU-yr | 1000 CPU-yr | C for E = 0.5 | C for E = 1 |
|---|---|---|---|---|---|---|---|---|
| power law N^4.72 (current build fit) | 0.010 | 0.030 | 0.072 | 0.16 | 0.29 | 0.48 | 1,300 CPU-yr | 2x10^5 CPU-yr |
| data-curved (N^2.7 at 8k, then N^2.7) | 0.010 | 0.033 | 0.090 | 0.22 | 0.46 | 0.77 | 140 CPU-yr | 5,000 CPU-yr |
| forecast's time model | 0.020 | 0.052 | 0.11 | 0.20 | 0.33 | 0.51 | 930 CPU-yr | 3x10^5 CPU-yr |
| 10x faster at all N (N^4.72) | 0.030 | 0.076 | 0.16 | 0.29 | 0.48 | 0.70 | 130 CPU-yr | 2x10^4 CPU-yr |
| N^3.5 beyond N = 2,000 | 0.010 | 0.034 | 0.094 | 0.22 | 0.43 | 0.69 | 190 CPU-yr | 1.2x10^4 CPU-yr |
| N^2.5 beyond N = 2,000 | 0.010 | 0.041 | 0.13 | 0.33 | 0.63 | 0.98 | 40 CPU-yr | 1,200 CPU-yr |

At 1 CPU-year the marginal sum has N ~ 3,500-5,000, and the marginal cost is 16-37 CPU-years per magic
square. At 100 CPU-years the marginal sum has N ~ 6,500-20,000 and the cost is 700-1,700 CPU-years per
magic square. To exhaust every sum of an N band (power law / data-curved time model) costs:

| N band | E_magic | CPU-years to exhaust the band | CPU-years per magic square |
|---|---|---|---|
| < 2k | 0.023 | 12 | 520 |
| 2-4k | 0.087 | 410 / 260 | 4,800 / 2,900 |
| 4-8k | 0.22 | 1.3x10^4 / 3.2x10^3 | 5.8x10^4 / 1.4x10^4 |
| 8-16k | 0.37 | 6.8x10^5 / 4.0x10^4 | 1.8x10^6 / 1.1x10^5 |
| 16-32k | 0.39 | 2.4x10^7 / 3.7x10^5 | 6x10^7 / 9x10^5 |

What follows for speed versus complexity in N:
- Under the analytic model E(C) does not saturate as in the forecast. It grows by ~0.1-0.15 per decade
  of CPU, because each octave of N from 8k to 64k carries a similar ~0.3-0.4 of expectation. The cost
  per expected magic square in a band grows roughly like N^alpha: like N^3.2-5.0 for the power law
  and N^2.3-3.1 for the data-curved model.
- A constant-factor speedup s shifts the curve left by s. 10x lifts E(1 CPU-yr) from 0.07 to 0.16.
- Lowering the exponent alpha matters more at large budgets. alpha 4.7 -> 2.5 beyond N = 2,000 doubles
  E(100 CPU-yr) to 0.63 and cuts the CPU for E = 0.5 from ~1,300 to ~40 CPU-years. The few large-N
  timings suggest that the current build is already closer to the "data-curved" row than to the
  power law. Measuring t(N) at N = 10-30k (a handful of sums) is the cheapest way to pin this down.
- These E(C) values assume an ideal targeting of the best sums under the analytic model, with sums
  searched independently (`msearch --sums`). Under this model the sums already searched hold
  E_magic = 0.006 (first search: 729,748 predicted squares at ~1 in 1.2x10^8 each; the report's 1 in
  8.6x10^8 omits the S-P correlation rho^2 ~ 5). Our near-S_min runs hold 0.0002.

## 10. Uncertainty budget and caveats (n = 6 total)

| source | factor on N_magic(inf) | basis |
|---|---|---|
| sampling (stratified, 2,809 P) | x/÷ 1.06 | se |
| permutation multiplier (365-P ratio estimator) | x/÷ 1.1 | spread between X bins |
| quadrature grid, universe tail, low tau | x1.02-1.05 | measured |
| semi-magic calibration K (lin vs quad, coefficient se, n6-low residual 1.23) | x0.85-1.25 | section 4.2 |
| diagonal calibration kappa = 0.57 (0.46-1.2) | x0.8-2.1 | traversal counts, section 4.3 |
| squares at N > 8k (82% of the total; validated only to N ~ 8k) | x0.3-1.5 on that part | obs/pred 0.87 at 5-8k; no free parameter in N |
| P(magic \| semi) at large N, diagonal correlations | x0.5-1.5 | SP checked to N ~ 5k (0.68-0.97) |
| **combined** | **N_magic(inf) ~ 2.0, 68% range 0.8-4, 90% range 0.4-6** | |

Hence P(a 6x6 AM magic square exists) = 1 - e^{-N}: ~0.86 central, 0.55-0.98 (68%), 0.33-1 (90%).

Caveats and approximations, explicitly:
1. **Independence.** Each (P, S) is treated as an independent random trial, and squares in different
   (P, S) as independent. Deterministic families are counted separately or not at all: multiples cM of a
   primitive square are excluded, and structured constructions (products of smaller squares, which make
   n = 8 easy) are not modelled. They can only add squares, so for n = 8 the model is a lower bound in
   spirit. For n = 6 no such construction is known.
2. **Local CLT far from its regime.** Lines have only n = 5-8 terms and the cell laws are skewed. The
   raw model is therefore off by e^{-0.6} to e^{-10} (section 4.2) and is corrected by a smooth empirical
   factor K_n(ld). The factor is fitted on n = 5, 6, 7 near S_min with N <= 8k, so n = 8 relies on its
   linear-in-n extrapolation. The crossing point is insensitive to it because g ~ 70: lin and quad give
   377-410. The n = 6 total relies on extrapolating in N.
3. **Diagonals.** p_pair uses the conditional Gaussian for two diagonals with exact lattice factors
   L_M(a) per prime and one empirical factor kappa. The diagonal-pair correlation (S+S pairs 1.1-1.2x
   more frequent than independent) is included only as a constant.
4. **Distinctness** enters through P_dist(n^2) of the iid tilted measure plus the calibrated power
   (1.14). The real constraint is correlated with the line conditions. This is the largest source of the
   raw-model error, and it is absorbed by K.
5. **Universe.** Sorted exponents (with a measured permutation multiplier), primes <= 31, k <= 10, at most
   6 exponents equal to 1. All are checked to be immaterial for n = 5, 6. For n = 7, 8, low-tau P were
   added near the crossing only.
6. **Formal asymptotics vs practice.** The finite totals of section 6 are rigorous consequences of the
   heuristic, but the formal tail ~X^{-n} (log X)^{M_n(1)-1} only takes over at astronomically large X.
   The practical description is the local exponent g(X) of section 7.

## 11. Files (`analytic/` in `research/existence/scripts.tar.xz`; data not in the repo)

- Model: `maxent.py` (tilted measure), `batch2.py` (vectorised per-P model over S), `model.py`,
  `corr.py` (lattice factors), `esym2.py` (P_dist saddle point), `gmodel.py` (calibrated model; C0, B, DN,
  KAPPA), `c/ehr.c` + `c/ehr.txt` (exact H_n, M_n), `ehrpoly.py`.
- Calibration: `calib3.json` (19,834 n = 6 sums), `xn_comp.json` + `xn_fit.py` (cross-n Poisson fit),
  `legacy_pred.json` / `legacy_an.py` / `legacy_nbins.py` (first search), `travcal.py` (traversals),
  `nbins.py`, `schedobs.py` (by N, vs scheduler fit).
- Global sums: `gs3.py`, `gs3_est.py`, `gs3_report.py`, `gs3_bins.py`, `gs3_n{5,6,7,8}.json`,
  `gridbias.py`, `lowtau.py`, `final_tables.py` -> `final_tables.json`.
- New checks: `n5low.py` -> `n5low_{400,500,600,700}.json`; `runs5/` (exhaustive 5x5 search,
  `squares5_le700.jsonl`), `n5cmp.py`; `maxent_ratio.py`.
- Forecast comparison: `pool_cmp.py` -> `pool_cmp.json`, `pool_split.py`, `schedcmp.py` -> `schedcmp.json`,
  `frontier_an.py` -> `frontier_segs.json`, `frontier_eval.py` -> `frontier_eval.json`, `timefit.json`.
