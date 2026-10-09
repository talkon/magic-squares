# Do 6x6 additive-multiplicative magic squares exist, and how large is the smallest one?

A reconciled heuristic estimate, written 2026-10-08 and revised the same day after two reviews (a
skeptical mathematician's and a skeptical data scientist's; their scripts are in `skeptic/` and
`datasci/`). §10 lists what changed. The estimate combines three investigations, each with its own
write-up in `research/existence/`:

* `analytic.md`: an analytic heuristic. Exact exponent-matrix counts are combined with a
  max-entropy / local-limit model and calibrated on the search data.
* `model.md`: an empirical model ("model lens"). Squares are counted per sum with a structural
  "random hypergraph" model plus a large-N correction measured by sampled searches at N = 15-32k.
* `literature.md`: the known results for each order, and checks on normal squares.

The synthesis's own checks are in `synth/`, and the post-review recomputation is in `synth/rev/`
(§4.1). The study's scripts (Python, plus the exact H_n counter `analytic/c/ehr.c` and the r1-sampling patch `model/sample/sampling.patch`) are archived in `research/existence/scripts.tar.xz`, with paths relative to the study's directory. The data they read (msearch output, caches) are not in the repo. Every number is an expectation under a heuristic, not a search result, unless it is marked
"observed". "Fit" marks data a parameter was fitted on; "test" marks data no parameter saw.

## 1. Answer

The definitions:

* N_magic(X) is the expected number of primitive 6x6 AM-magic squares with S <= X, over all P.
  The squares have 36 distinct positive integers, and rows, columns and both diagonals have equal
  sum S and equal product P. Primitive means the gcd of the entries is 1.
* Squares are counted up to row and column permutations and transposition.
* **Analytic** is the central curve: the analytic model, corrected for exponent-assignment coverage
  (§3.2) and for the review findings (§4.1).
* **Model lens** is the empirical model's central curve. It is an alternative model form, not an
  upper bound: it is lower than the analytic curve below X ~ 2,300.

| X | 1,000 | 1,500 | 2,000 | 3,000 | 4,000 | 6,000 | 10^4 | 2x10^4 | 5x10^4 | inf |
|---|---|---|---|---|---|---|---|---|---|---|
| **N_magic(X), analytic** | 0.014 | 0.10 | 0.27 | **0.74** | 1.2 | 2.0 | 3.2 | 3.9 | 4.4 | **4.6** |
| 68% band (model error, §4.1) | 0.007-0.03 | 0.05-0.2 | 0.13-0.5 | 0.37-1.5 | 0.6-2.5 | 1.0-4.4 | 1.5-7.5 | 1.9-11 | 2.1-12 | 2.2-13 |
| local growth g = dlnN/dlnX | 5.8 | 4.1 | 2.9 | 2.1 | 1.5 | 1.1 | 0.5 | 0.2 | ~0 | |
| N_magic(X), model lens | 0.0065 | 0.058 | 0.22 | 1.2 | 3.2 | 11 | 62 | 400 | (diverges) | |
| P(smallest S <= X), analytic with model error | 0.02 | 0.12 | 0.27 | 0.53 | 0.68 | 0.82 | 0.90 | 0.93 | 0.94 | 0.95 |
| previous version (analytic x coverage) | 0.014 | 0.095 | 0.23 | 0.60 | 0.96 | 1.55 | 2.3 | 2.9 | 3.2 | ~3.2 |

1. **Growth with the maximum magic sum.**
   * N_magic grows steeply and then ever more slowly: like X^6 at X ~ 1,000, X^3 at 2,000 and
     X^2 at 3,000.
   * For X <= 3,000 the two models differ by 0.45-1.6x. Within the sums with N <= 32k, the range
     where square counts are measured, they differ by 1.3x at X = 3,000 and 1.6x at 4,000.
   * Beyond ~6,000 they part ways:
     * The analytic curve saturates. Its local exponent is ~1 at 6,000, 0.5 at 10^4 and ~0
       beyond 3x10^4. Its total is **N_magic(inf) ~ 4.6**, with a 68% range of 2.2-13 and a 90%
       range of 1.3-24, skewed upward.
     * The model lens keeps growing like X^3.
   * The largest entry is r S with r ~ 0.36-0.45 at the relevant sums. In terms of the largest
     entry Y, the curves are the same with X = Y/r.
2. **Infinitely many?**
   * For primitive squares with distinct entries, the independence heuristic says no. A
     Dirichlet-series bound makes the expected total finite for every n >= 3 (§4.2), and the
     calibrated n = 6 model saturates at a few squares by S ~ 2-5x10^4.
   * **This is a property of the heuristic, not something the data show.** 35% of the analytic total
     lies at sums with N > 45k, beyond anything measured. The large-N trend of the square counts is
     +0.12 ± 0.26 per ln N. At +1 sigma the total grows ~1.25x, and at +2 sigma ~2x, but it stays
     finite.
   * The model lens diverges, but it is disfavoured at ~2 sigma by the two sums measured at
     N = 39-45k (§3.3).
   * The heuristic excludes structured families:
     * The multiples cM of any magic square are trivially infinite in number.
     * Latin-square constructions at n = 8, 9 and 16 may give infinitely many primitive squares.
       None is known at n = 6.
3. **The smallest one.**

   | | median S | 10-90% range | where N_magic(X) = 1 |
   |---|---|---|---|
   | analytic, with model error | ~2,870 | 1,450-~9,900 | ~3,600 |
   | model lens, with model error | ~2,640 | 1,630-4,000 | ~2,900 |
   | **summary** | **~2,800** (2,650-2,950 across models) | **~1,500 to 6,000-10,000** | 2,900-3,600 |

   * Measured by the largest entry, the median smallest is ~1,000-1,300.
   * The smallest square is expected at:
     * P with tau(P) ~ 5-15k divisors and 5-7 distinct primes;
     * S_min ~ 1,500-2,500 and S/S_min ~ 1.1-1.6;
     * sums with N ~ 3k-32k vectors (91% of N_magic(3,000)).
   * **P(at least one primitive 6x6 AM-magic square exists) ~ 0.95** (90% range 0.75-1). This is
     conditional on no structural obstruction. Such obstructions exist for n = 3 and 4, and no
     counting heuristic can see them.
4. **Calibration.**
   * **Two-sided tests that no parameter saw** (§5) are all 6x6 semi-magic or 5x5 data, plus the
     first search's traversal and SP-pair rates:
     * the first search's 758,949 squares: obs/pred 1.04;
     * the full sums at N = 3.4-8.9k: 0.78;
     * the sampled sums at N = 15-45k: 1.23;
     * the exhaustive 5x5 search with S <= 700: 13 squares against 12.8 predicted.
   * **P(magic | semi-magic)** is calibrated on the S, P and SP traversals (kappa; 14 SP
     traversals) and checked rung by rung by the ladder of `research/calibration.md`. No SP+SP event
     has been seen.
   * **Records of other orders are consistency checks with almost no power** (7x7 magic at 465; 8x8
     at 600). The predicted counts grow like X^23-70 there, so even a 30-100x error in the level
     would leave the records unsurprising.
5. **Search.**
   * Expected magic squares in everything searched so far: 0.004-0.006. A model-free lower bound from
     the first search's SP counts gives 0.0037.
   * Assume the current build, an ideal scheduler, and that the model is right. Then the expected
     number found is **0.16 in 1 CPU-year, 0.34 in 10, 0.64 in 100 and 1.1 in 1,000** (68% band
     x/÷2).
   * **The CPU needed for E = 1 is undetermined.** It is ~650 CPU-years at the central model, but
     40 to 6x10^4 CPU-years at ±1 sigma, and anywhere from ~8 to ~10^8 at 90%. E rises by
     only 0.2-0.5 per decade of CPU, so any constant factor on E moves the CPU for E = 1 by orders
     of magnitude.
   * Speedups are best judged by ratios at a fixed budget, which no constant factor on E affects:
     * a 10x faster search multiplies E by 2.2 / 1.9 / 1.7 / 1.5 at 1 / 10 / 100 / 1,000 CPU-years;
     * t ~ N^2 beyond N = 4,000 multiplies it by 1.11 / 1.35 / 1.60 / 1.68;
     * t ~ N^3 multiplies it by 1.04 / 1.11 / 1.21 / 1.25 (§6).

## 2. The heuristic

For a product P = prod_p p^{a_p} (k primes, tau(P) = prod (a_p + 1)) and a line sum S, every
(P, S) is treated as an independent random trial:

    E_magic(P, S) = E_semi(P, S) x P(magic | semi-magic; P, S)
    N_magic(X)    = sum_P sum_{S <= X} E_magic(P, S),     P(smallest S <= X) = 1 - exp(-N_magic(X)).

**Multiplicative part, exact.** A grid of divisors of P has all row and column products P iff, for
every prime p, its p-exponent matrix has all line sums a_p. The number of such grids is
prod_p H_6(a_p) (Birkhoff polytope Ehrhart polynomial; H_6(1) = 720, H_6(2) = 202,410). With both
diagonals it is prod_p M_6(a_p) (M_6(1) = 96, M_6(a)/H_6(a) ~ 1.18/(a + 2.1)^2). The review
recomputed H_6(1..4) and M_6(1..5) independently. The asymptotic regime is never reached
(H_6(10) is 665x its leading term), so the exact values and lattice factors are needed.

**Additive part, local CLT.** Take the max-entropy measure mu on the divisors of P with
E[d] = S/6 and E[e_p(d)] = a_p/6. Let H be its entropy and C the covariance of
v(d) = (d, e_2(d), e_3(d), ...), which has q = 1 + k dimensions. Then:

    N(P,S)    ~ e^{6H} (2 pi)^{-q/2} det(6C)^{-1/2} prod_p L_vec(a_p) P_dist(6) / 6!         (vectors; 10% accuracy)
    E_semi    ~ K(ld) e^{36H} / [(2 pi)^{11q/2} det(C)^{11/2} 6^{5q}] prod_p L_semi(a_p) P_dist(36) / (2 (6!)^2)
    log K     = -0.597 + 0.1365 ld,   ld = log P_dist(36)                                    (fitted on n = 5, 6, 7)
    P(magic | semi) = D_6 kappa p_pair = 5400 x 0.57 x prod_p L_M(a_p) / [(2 pi)^q det C 24^{q/2}]
                    ~ 5400 (rho p_S p_P)^2 x (pair factor),   rho = (1 - R^2)^{-1/2} ~ 2.2-2.6

In these formulas:

* L_* are exact per-prime lattice factors.
* P_dist(m) is the probability that m draws from mu are distinct.
* 5,400 is the number of diagonal/anti-diagonal pairs reachable by row and column permutations
  (720 x 15/2). It matches 2(6!)^2/192, where 192 is the order of the group preserving magic
  squares.
* R is the multiple correlation of an entry's size with its exponent vector.
* The pair factor is 1.1 x (25/24)^{q/2} = 1.23-1.27: the S+S excess, times the conditional
  correlation -1/5 of two disjoint diagonals.

Rough scalings: N ~ e^{6H}/(S prod a_p)^{1/2}, E_semi ~ N^6 x (density)^5 x P_dist(36)^{1.14},
and P(magic | semi) ~ 1/(S^2 prod_p (a_p/6)(1 + a_p/6)).

**Within one P, P(magic | semi) falls steeply with x = ln(S/S_min)**, much more steeply near S_min
than the forecast's S^-7..-9. For 13 7 4 3 1 1, 1/P(magic | semi) is:

| S/S_min | 1.05 | 1.1 | 1.2 | 1.3 |
|---|---|---|---|---|
| 1/P(magic \| semi) | 7.3x10^7 | 6.0x10^8 | 4.8x10^9 | 1.6x10^10 |

* The local slope is -20 to -25 for x < 0.25, about -12 for x = 0.2-0.4, and -8 only at
  1.3-2 S_min. The agreement with the forecast holds only for x > 0.26.
* The data say the model is somewhat too steep. The P-traversal ratio rises from 0.52 (x = 0.05-0.10)
  to 0.68 (x = 0.3-0.5), about 1.7x per diagonal pair across that range.
* The effect on N_magic is under ±20%, because kappa is averaged over x and the yield peaks at
  x = 0.1-0.25. The x < 0.1 correction is applied in §4.1.

The model lens uses the same P(magic | semi) form, with p_S and p_P taken from the scheduler's
log-linear fit. For squares per sum it uses a regression on the random-hypergraph count:

    mu = e^{-11.02} H(N,L)^{0.46} e^{-0.58 x} N^{2.56} x h(N),   H = (L)_36 (N/C(L,6))^12 / (2 x 720^2),
    h(N) = 0.64 (N/10^4)^{-0.27} for N >= 3,000

Here L is the number of distinct labels (entries in use).

## 3. Why the two quantitative estimates differ, and which to trust

The two families agree on most ingredients:

* **Vectors N(P, S).** Each lens matches the actual counts to within 10%. At the sums checked, their
  two predictions differ by 3-15%.
* **P(magic | semi) at the sampled large-N sums.** They agree within 1.5x, for example 1 in
  3.46x10^9 vs 1 in 3.47x10^9 at (13 7 4 3 1 1, S = 2,200).
* **Squares per sum on the measured sums.** On the 10 full sums at N = 3.4-8.9k the observed/predicted
  ratio is 0.78 for the analytic model and 1.03 for the model lens. On the 7 sampled sums at
  N = 15-32k it is 1.28 and 0.73. The data sit between the two models.

They disagree in three places:

| source | effect on N_magic | resolution |
|---|---|---|
| (a) **Exponent assignments not covered by the analytic universe.** The analytic lens took sorted exponents with at most one zero gap, plus permutations within the same primes. It left out spreads onto larger primes, e.g. 13 7 4 3 0 0 1 1, and primes >= 53. | analytic too low by x1.14 (X = 1,000), 1.46 (3,000), 1.79 (inf) | computed (§3.2), applied |
| (b) **Model form for squares per sum within the measured range (N <= 32k).** | the models differ by 1.3x at X = 3,000 and 1.6x at 4,000, using only N <= 32k | the data sit between the two models (above). The analytic curve is scaled by its measured 1.23 at N >= 12k (§4.1). The X <= 4,000 gap is carried in the error |
| (c) **Squares per sum beyond N ~ 32k.** The analytic model saturates at ~5x10^4 squares per sum, while the model lens grows like N^2.5-3. | only 4% of the analytic N_magic(3,000) and 10% of N_magic(4,000) lie at N > 32k; at X >= 6,000 it is 21% or more | two sums at N = 39-45k favour the analytic form at ~2 sigma (§3.3). The model lens is disfavoured beyond X ~ 6,000, not beyond 3,000 |

### 3.2 Coverage of exponent assignments

I enumerated every distinct assignment of an exponent multiset to the first 15 primes with scale at
most 1.35x the sorted one, and summed the analytic E(S <= X) over them. The estimator uses **24
multisets drawn with probability proportional to their magic yield** from the analytic lens's
stratified sample (`synth/assign_an.py`, `ratio.py`; 4 hand-picked multisets were run but not used).

| X | 1,000 | 1,500 | 2,000 | 3,000 | 4,000 | 6,000 | >= 10^4 |
|---|---|---|---|---|---|---|---|
| all assignments / analytic coverage, R(X), to scale 1.35 | 1.14 | 1.33 | 1.39 | 1.44 | 1.53 | 1.52 | 1.62 |
| all assignments / analytic universe members | 1.33 | 1.82 | 2.07 | 2.33 | 2.63 | 2.69 | 2.92 |

* **Scales beyond 1.35.** At scale 1.30-1.35 the yield density is still 1.4% per 0.05 of scale at
  X = 3,000 and 3.5% at infinite X, decaying ~1.5x per step. Extrapolating, larger scales add ~1.5%
  at X = 3,000 and ~5% at infinite X, so **R(inf) ~ 1.70**. This is an extrapolation; nothing was run
  beyond scale 1.35.
* **Primes >= 53.** No universe member or R(X) assignment covers them. On a 40-P sample drawn by
  yield, adding one prime p >= 53 to P adds 5.3% ± 0.3% of the yield (`skeptic/addprime.py`). The
  share is 3% for base P with S_min 1-2k and 8-10% at 1-3x10^4. The yield of P·p sits at
  ~2x the sums of P, so it adds only 0.5% at X = 3,000 but 5.4% at infinity.
* Per multiset, R(inf) is 1.06-2.16 (median 1.58). It is largest for multisets with several
  exponents equal to 1 or 2.
* The model lens's own all/monotone factor was 4.6 at X = 3,000. The analytic model gives 2.1-5.1
  per multiset; the review's spot check got 3.99 and 5.07 against the 3.65 and 4.69 recorded here.
  So the two lenses agree on this factor.

### 3.3 Squares per sum at large N

The observed/predicted ratio of squares per sum is fitted over all 19 sums with N >= 3.4k by a
Poisson regression with overdispersion (`datasci/bigN_trend.py`):

    obs/analytic = exp(-0.19 + 0.12 ln(N/10^4)),   slope se 0.26,   overdispersion 7.4

The slope is consistent with zero, so the analytic model's N-shape is right within the errors. Extrapolated, the ratio
is 1.10 (1-sigma 0.60-2.0) at N = 10^5 and 1.25 (0.52-3.0) at 3x10^5.

| sums | observed | analytic (no free parameter in N) | model lens H_corr (fitted to the 15-32k sums) |
|---|---|---|---|
| 10 full sums, N 3.4-8.9k (test) | 860 | 1,104 (obs/pred 0.78) | 837 (1.03) |
| 7 sampled sums, N 15-32k (squares in the samples) | 60 | 47.0 (1.28) | 82.4 (0.73) |
| 2 new sampled sums, N 39-45k (test) | 8 | 8.2 | 43.5 |
| all 9 sampled sums, N >= 12k | 68 | 55.2 (**1.23**, 1-sigma 0.99-1.53) | 126 (0.54) |

The two new sums lie beyond N ~ 32k, at larger P. They were sampled with `msearch_sample`, which
searches every k-th first row; its estimate is unbiased for both counts and time.

| P (tau, k) | S (x) | N, labels | first rows sampled | CPU (s) | squares in sample | analytic | model lens | model lens with observed labels |
|---|---|---|---|---|---|---|---|---|
| 14 8 5 3 2 1 1 (38,880, 7) | 8,000 (0.10) | 44,732, 330 | 23 of 44,732 | 750 | **1** | 3.85 | 19.6 | 1.4 |
| 13 7 5 3 2 1 1 (32,256, 7) | 6,015 (0.11) | 39,309, 310 | 26 of 39,309 | 540 | **7** | 4.37 | 23.9 | 1.7 |

* **Estimated full-sum counts.** The first sum has 1,945 ± 1,945 squares (predicted: 7,480
  analytic, 38,190 model lens). The second has 10,583 ± 4,111 (predicted: 6,611 and 36,203).
* **Estimated full-sum CPU.** 1.46x10^6 s (17 days) and 8.1x10^5 s (9 days).
* **Significance.** Against pure Poisson noise, the model lens's 19.6 → 1 and 23.9 → 7 would be
  overwhelming (P = 6x10^-8 and 5x10^-5). Per-sum counts scatter lognormally with sd ~1
  (overdispersion 8-17), and with that scatter **the two sums favour the analytic form at about 2-2.5
  sigma**. The 7 sums at 15-32k do not discriminate between the models.
* **Why the model lens grows (a hypothesis, not shown).**
  * Its count goes as H(N, L)^0.46, which is proportional to L^-16.6 at fixed N.
  * Its population takes L from a label model fitted on small P. That label model underestimates
    the labels at large P: by 5-12% at the 7 sampled sums and by 15% at the two new sums (279 vs 330
    and 262 vs 310).
  * With the observed labels, the model lens predicts 3.1 squares against the 8 observed, so the
    H form is fragile in both directions.
  * Rerunning its population with labels refitted on the large-N sums would test this, and was not
    done.

**Decision.** The central curve is the analytic model with the corrections of §3.2 and §4.1:

* It has a structural form with no free parameter in N.
* It is tested out of sample up to N = 45k, tau = 39k and k = 7.
* Its per-P S-profile is confirmed: exhaustive P (12 6 3 2 to 20 S_min; 11 4 3 2 1 to 10 S_min)
  give 37 and 67 predicted against 42 and 74 observed. The model lens predicts 25 and 26.

The model lens is kept as an alternative form for X <= 4,000, where the data do not separate the
two, and enters the smallest-S summary. Its divergence beyond ~6,000 is disfavoured at ~2 sigma.

### 3.4 Diagonals by N, k and x

The analytic per-traversal rates were compared with the observed S-, P- and SP-traversals of the
6,957 distinct squares from our runs. These are the 6,949 of the near-S_min, large-N and verifier
runs plus the 8 in the two new samples (`synth/trav_an.py`, `datasci/trav_bins.py`). The squares
are counted differently elsewhere: the analytic lens's per-square calibration used the 6,598
near-S_min squares, and `research/calibration.md` used 7,021, deduplicated over all msearch output.

| subset | squares | S obs/model | P obs/model | SP obs/model | rS·rP vs global |
|---|---|---|---|---|---|
| all | 6,957 | 0.87 | 0.64 | 14 / 14.5 | 1 |
| N < 1k | 113 | 0.82 | 0.55 | 1 / 0.55 | 0.81 |
| N 1-3k | 5,916 | 0.82-0.88 | 0.62-0.66 | 13 / 13.6 | 0.95-1.04 |
| N 3-6k | 188 | 0.99 | 0.75 | 0 / 0.12 | 1.35 |
| N 6-9k | 672 | 0.93 | 0.63 | 0 / 0.25 | 1.06 |
| N 12-50k | 68 | 1.17 | 0.94 | 0 / 0.01 | 1.99 (x/÷1.5) |
| **N >= 3k, pooled** | 928 | 0.95 | 0.68 | 0 / 0.4 | **1.17 (x/÷1.10)** |
| x < 0.10 | 335 | 0.66-0.83 | 0.52-0.57 | 0 / 0.6 | |
| x 0.10-0.20 | 3,244 | 0.84-0.85 | 0.57-0.65 | 6 / 7.3 | |
| x 0.20-0.50 | 2,956 | 0.86-0.91 | 0.64-0.68 | 7 / 5.9 | |
| k = 4 / 5 / 6 | 6,920 | | 0.66 / 0.64 / 0.61 | | |
| k >= 7 | 37 | | 0 / 2.4 (P = 0.09) | | |

* **The diagonal calibration kappa = 0.57 is valid for k <= 6, N <= 9k (full sums) and x >= 0.1.**
  The yield lies largely outside that range, so §4.1 applies four explicit factors:
  * The S and P ratios rise with N, by +0.036 ± 0.025 and +0.056 ± 0.053 per ln N. If the SP rate
    scales like rS·rP, the pair rate at N >= 3k is 1.17^2 = 1.37x the global value. **Factor x1.35
    at N >= 3k (range 1.0-1.7).** The N >= 12k bin alone suggests ~4x but is only 1.5 sigma. There
    is no SP or pair calibration at N > 12k, nor at k >= 7.
  * The P ratio drifts down with k, and k >= 7 has 0 P-traversals against 2.4 predicted.
    **Factor x0.85 for k >= 7 (range 0.7-1.0).** This class holds 15% of N_magic(3,000) and 37%
    of the total.
  * **Factor x0.92 for x < 0.1 (range 0.86-0.97).** This class holds 3% of N_magic(3,000) and 18%
    of the total.
  * **Pair factor: x0.87 (range 0.75-1.0).** The model's 1.23-1.27 is high against the ladder of
    `research/calibration.md`:
    * on sub-events (65 to 12 million pairs) the pair factor given the counts stays at 1.00, and
      in pair types without extra congruences it falls to 0.88 at 4 pinned exponents;
    * S+S pairs are 1.19 [1.04, 1.37];
    * the first search's SP+S and SP+P pairs are 10 observed against 12.4 expected at global rates
      (9 vs 9.3 and 1 vs 3.0), a ratio of 0.81 [0.44, 1.37]. The report's 4.58 and 1.46 omit a
      factor 2 for unordered mixed pairs.
* **Net effect.** The four factors give x1.12 on N_magic(3,000) and x1.08 on the total, so they
  mostly cancel. Their range is skewed upward in the yield region.
* **Cross-check.** The analytic model's E[magic] over the 6,598 near-S_min squares is 2.4x10^-4.
  The ladder's heuristic, with its pair factor corrected by the sub-events, gives 1.8x10^-4
  (1.6-2.3x10^-4) over 7,021 squares. Its SP rung is 14 observed against 8.8 predicted, an excess
  `research/calibration.md` attributes to seed P inside the first search's range.
* **Literature lens.** For normal squares of orders 5 and 6, P(magic | semi-magic)/p_S^2 is 0.97 and
  0.92 with no fitted parameter, which supports the diagonal-independence assumption.

## 4. Reconciled growth, finiteness, smallest size

### 4.1 Revised curve and its error

`synth/rev/segs.py` re-evaluates all 2,809 P of the analytic stratified sample. It uses the full
S-profile of each P, without the S <= 5x10^4 cap, and keeps N, x and k for every S segment. On the
uncorrected model it reproduces the published sorted totals within 1% (0.266 vs 0.264 at 3,000;
1.10 vs 1.10 at 5x10^4). `rev.py` then applies the corrections segment by segment and scales the
published N_magic(X) x R(X) by the ratio:

| correction | factor at X = 3,000 | at inf | basis |
|---|---|---|---|
| assignment scales beyond 1.35 (§3.2) | x1.015 | x1.05 | extrapolated decay |
| primes >= 53, at ~2x the S of the base P (§3.2) | x1.005 | x1.054 | 40-P sample by yield |
| P beyond scale 5x10^4 | x1 | x1.039 | the per-half-octave yield decays 0.70x from scale 8k to 46k (fit); a geometric tail |
| squares per sum x1.23 at N >= 12k (§3.3) | x1.08 | x1.17 | 9 sampled sums, test |
| diagonals: x1.35 at N >= 3k, x0.85 at k >= 7, x0.92 at x < 0.1, pair x0.87 (§3.4) | x1.12 | x1.08 | traversals, ladder |
| **total** | **x1.23 (0.60 → 0.74)** | **x1.45 (3.2 → 4.6)** | |
| all review factors at their low / high ends | N = 0.44 / 1.2 | N = 2.3 / 8.7 | `variants.out` |

The skeptic's corrections alone (scales, primes, tail) give N_magic(inf) = 3.7. Adding the
squares-per-sum factor gives 4.3, and adding the diagonal factors gives 4.6.

**Where the revised analytic N_magic(X) sits, by N of the sum** (`rev.out`):

| X | N < 3k | 3-12k | 12-32k | 32-45k | > 45k |
|---|---|---|---|---|---|
| 3,000 | 6% | 53% | 38% | 3% | 1% |
| 4,000 | 5% | 44% | 41% | 6% | 4% |
| 6,000 | 3% | 34% | 42% | 10% | 11% |
| 10^4 | 2% | 26% | 37% | 11% | 22% |
| inf | 2% | 21% | 32% | 11% | 35% |

The median smallest S therefore depends barely at all on extrapolating beyond the measured N. The
total and P(exists) do: 35% of the total lies beyond the largest N measured (45k), and 46% beyond
32k.

**Model error.** It is a lognormal factor lambda on N_magic, with a sigma that depends on X.

* **At X <= 4,000, sigma = 0.7.** The components are:
  * the K calibration (x0.85-1.25);
  * the base kappa (0.46-1.2, from 14 near-S_min SP traversals against the first search's 0.68);
  * the review factors above (90% corner range x0.6-1.6);
  * the model form (the model lens is 1.3-1.6x higher on N <= 32k);
  * sampling and coverage (x/÷1.06, x/÷1.1).

  They combine to ~0.55 in quadrature, rounded up for what is not modelled.
* **At X >= 2x10^4, sigma = 0.75 downward and 1.0 upward**, interpolated in log X from 4,000. The
  extra upward error is the extrapolation beyond N = 45k:
  * the large-N slope of +0.12 ± 0.26 changes the decay per half-octave from 0.70 to 0.70·1.8^b;
  * the diagonal N-trend is untested beyond 12k.

  The 68% band in §1 uses these sigmas.

* **Update ([calibration-target.md](calibration-target.md) §5).** A pre-registered search of 321 sums
  at N' 3-45k measures the squares level (to ±4%) and the S/P rates (to ±10% on (f_S f_P)^2) at
  3-12k. It finds 7 SP pairs against 5.46, which leaves kappa at ±0.27 (κ-only prior). Model form,
  the pair factor and sampling are unchanged. On the same components the quadrature sum goes from
  0.52 to 0.45, so sigma at X <= 4,000 becomes ~0.57 (x/÷1.8) instead of 0.7 with the same
  round-up. The widening at X >= 2x10^4 is not reduced: that search has 7 sums at 12-45k.

### 4.2 Finite or infinite

**The coincidence bound needs distinct (generic) entries.** Without distinctness the expected count
diverges because of degenerate grids:

* For P = p prime, all 96 multiplicatively magic grids are additively magic with S = p + 5.
* Any grid built from pairwise-disjoint permutation patterns has equal line sums identically.

For grids with distinct entries, a uniformly random multiplicatively magic grid with product P has
all 2n + 2 line sums equal with probability <~ c (n P^{1/n})^{-2n}. So

    E[# primitive AM-magic with distinct entries] <~ c sum_P prod_p M_n(a_p) P^{-2} = c zeta(2)^{M_n(1)} G(2) < infinity

for every n >= 3. Distinctness and balance only lower this. The semi-magic analogue at
s = 2 - 2/n > 1 is finite too.

**The bound shows convergence, not the size of the total.** With exact Euler factors the
review gets the following:

* Magic bound: log10 D = 19.7, median scale 1.0x10^4, 90% below 6.8x10^4. The yield per
  half-octave decays from 0.97x down to 0.7x past the peak.
* Semi-magic bound: log10 D = 50.3, median scale 4x10^13.
* The formal tail X^{-6} (log X)^95 would need P with ~96 log log P prime factors, so it is
  irrelevant at any searchable size.

**The calibrated n = 6 yield.**

* The magic yield per half-octave of P-scale peaks at scale ~2-6k.
* It decays ~0.70x per half-octave from 8k to 46k (untruncated). Beyond ~2x10^4 this decay is an
  extrapolation of the model, tested at 9 sums.
* So the n = 6 total saturates at ~5.
* For n = 7 and 8 the totals are still growing like X^3-8 at X = 10^4, and are astronomically large
  (10^14 and 10^30 by X = 2x10^4 and 10^4).

**Assumption that matters:** (P, S) cells are independent trials over generic grids.

* Deterministic families are excluded: multiples cM, and Latin-square constructions (Horner and
  others at n = 8, 9, 16, which may give infinitely many primitive squares).
* None is known at n = 6, and order 6 has no Graeco-Latin square to build one from.

### 4.3 The smallest square

Poisson first occurrence: P(min S <= X) = 1 - E_lambda[exp(-lambda N_magic(X))], with lambda
lognormal as in §4.1 and the same draw at every X (`synth/rev/quant.py`):

| curve | X(N = 1) | 10% | 25% | median | 75% | 90% | P(exists) |
|---|---|---|---|---|---|---|---|
| analytic, no model error | 3,590 | 1,510 | 2,060 | 2,920 | 4,440 | 6,810 | 0.99 |
| **analytic, with model error** | 3,590 | 1,450 | 1,950 | **2,870** | 4,840 | 9,900 | **0.95** |
| model lens, sigma 0.7 | 2,900 | 1,630 | 2,080 | 2,640 | 3,240 | 3,970 | ~1 |
| 50/50 mixture of the two | | 1,540 | 2,020 | 2,710 | 3,710 | 5,790 | |
| previous version (analytic, no error) | 4,150 | 1,560 | 2,200 | 3,290 | 5,450 | 9,840 | 0.96 |

* The lower quantiles agree well across models. The upper tail depends on the large-N behaviour.
* Summary: **median ~2,800; 25% below ~2,000; 10% below ~1,500.**
* **P(no primitive 6x6 magic square exists) ~ 5%.** This is the error-averaged analytic value. It
  is 11% at the error distribution's 16th percentile and 26% at its 5th, so the 90% range is
  0-26%. All of this is conditional on no structural obstruction.
* The largest entry: r = max/S is 0.31-0.55 for our squares, with median 0.36-0.45 at
  S/S_min = 1.1-1.4. The median smallest largest entry is therefore ~1,000-1,300.

## 5. Calibration

Two-sided tests of 6x6 semi-magic counts and traversals:

| test | fit / test | heuristic | observed | note |
|---|---|---|---|---|
| near-S_min sums, by N: <1k / 1-2k / 2-3k / 3-5k / 5-8k | **fit** (K fitted on all 19,834 sums) | 166 / 4,328 / 1,343 / 122 / 558 | 106 / 4,182 / 1,625 / 101 / 484 | the bin spread (0.64-1.21) is level scatter between P |
| within-P residual slope in N, 81 P | residual of the fit | 0 | −0.02 ± 0.04 per ln N | the strongest evidence that the N-shape is right |
| first search, 1,136 P, 758,949 squares, N <= 5.5k | test | 729,748 (obs/pred 1.04; 0.96-1.10 by S_min bin, 0.84-1.06 by N bin) | 758,949 | overlaps the near-S_min data on 34 seed P (4,741 squares, 0.6% of the total) |
| 10 full sums, N 3.4-8.9k | test | 1,104 | 860 (0.78) | |
| 9 sampled sums, N 15-45k | test (2 new) | 55.2 (model lens 126) | 68 (1.23, 1-sigma 0.99-1.53) | now applied as x1.23 at N >= 12k |
| traversals S / P / SP, 6,957 squares | **fit** (kappa) | calibrated 0.87 / 0.64 / 14.5 | 0.87 / 0.64 / 14 | drifts with N, k and x (§3.4) |
| first search's traversals S / P / SP | test | | 0.89 / 0.67 / 519 events (0.68) | flat across S_min bins |
| first search's SP+S, SP+P pairs | test | 9.3 / 3.0 | 9 / 1 (combined 0.81 [0.44, 1.37]) | no pair excess at SP |
| 5x5 semi-magic, all 3,966 sorted P, S <= 400 / 500 / 600 / 700 | test, exhaustive | 0.99 / 3.07 / 6.91 / 12.8 | 0 / 1 / 9 / 13 | tests the S-extrapolation by 3-4x |

Consistency checks with little or no power, mostly one-sided:

| check | heuristic | observed | what it can and cannot exclude |
|---|---|---|---|
| 6x6 smallest semi-magic S | N_semi(289) = 0.69 | 289 (Morgenstern 2007) | the semi-magic level within ~10x (g_semi ~ 20) |
| 7x7 smallest semi-magic | N_semi(310) = 8 | 310 (Shirakawa 2010) | within ~10-100x |
| 5x5 smallest semi-magic | N_semi = 1 at ~380 (all P) | 476 (P not in the sorted class) | |
| 7x7 magic | N_magic(465) ~ 6; N = 1 at ~430 | 465 (Miquel 2016, ~600 PC-hours, not exhaustive) | a record at 465 is a 6-18% event even if the model is 30-100x too high, so it excludes only a model >~100x too high; it says nothing if the model is too low |
| 8x8 magic | 1.6x10^9 with S <= 600 | 600 (Boyer), 840 (Horner), constructive | nothing below a 10^9-fold error |
| 5x5 magic | 5x10^-6 in total with the n = 6 kappa; 5x10^-5 to 10^-2 if the n = 5 traversal rates are used | none known (Boyer prize open) | not a calibration success: on the 33 5x5 squares the diagonal rates are under-predicted ~3x per diagonal (P 6 vs 2.0, SP 1 vs 0.012, S 19 vs 6). "None expected" stands |
| 3x3, 4x4 | | impossible (algebraic proofs) | structure the heuristic cannot see |

**None of the cross-order records tests P(magic | semi) at the 2-5x precision that matters for
n = 6.** That rests on the traversal rungs (§3.4) and the ladder of `research/calibration.md`.
Squares with a best pair of SP+S or better: 2 observed against 0.19-0.53 expected under the ladder's
predictors. The first search's SP+S and SP+P pairs: 10 observed against 12.4. No SP+SP has been
seen.

## 6. Implications for the search

* **What has been searched:** E = 0.004-0.006 in both models, nearly all of it from the first
  search. A model-free lower bound from the first search's per-P SP counts,
  5400 Σ n_SP(n_SP − 1)/(720^2 n_P), gives 0.0037 (0.0041 with a 1.1 pair factor). So the absence
  of a magic square so far says nothing.
* **Where the first magic square is expected:**
  * P with tau ~ 5-15k, 5-7 primes and S_min ~ 1,500-2,500;
  * 55-70% of the yield from exponent assignments outside the analytic universe (orders that are
    not non-increasing, spreads onto larger primes), plus ~5% (at large S) from primes >= 53;
  * S ~ 1.1-1.6 S_min, at sums with N ~ 3k-32k;
  * ~10^3 semi-magic squares per sum, each magic with probability ~10^-9-10^-10.5.

  The yield is spread over ~10^5 P, and no single P holds more than ~0.1% of it.
* **Cost with the current build** (`synth/rev/frontier_rev.py`). The frontier is built as follows:
  * It takes an ideal greedy order over all sums with S <= 2 S_min, using the revised analytic
    E(P, S).
  * It assumes the extra exponent assignments have the same E/t distribution (m = 2.7). The P·p
    sums are costed per prime.
  * Time per sum is the scheduler's model below N = 4,000. Above it is 51 s (N/4000)^3.92,
    steepening to slope 5.2 above 2x10^4. This reproduces the measured sums at N = 31.7k, 39.3k and
    44.7k within 1.3x.

  | time model | 1 CPU-yr | 10 | 100 | 1,000 | 10^4 | CPU-yr for E = 1 at F = 2 / 1 / 0.5 |
  |---|---|---|---|---|---|---|
  | **current build: E** (68% band x/÷2) | **0.16** (0.08-0.31) | **0.34** (0.17-0.68) | **0.64** (0.32-1.3) | **1.1** (0.55-2.2) | 1.6 | 39 / 650 / 5.9x10^4 |
  | t ~ N^3 beyond N = 4,000 | 0.16 | 0.38 | 0.78 | 1.37 | 2.0 | 24 / 250 / 8.7x10^3 |
  | t ~ N^2 beyond N = 4,000 | 0.17 | 0.46 | 1.03 | 1.83 | 2.6 | 13 / 91 / 1.6x10^3 |
  | 2x faster at all N | 0.20 | 0.41 | 0.76 | 1.25 | 1.8 | 19 / 325 / 2.9x10^4 |
  | 10x faster at all N | 0.34 | 0.64 | 1.09 | 1.61 | 2.1 | 3.9 / 65 / 5.9x10^3 |
  | 100x faster at all N | 0.64 | 1.09 | 1.61 | 2.11 | 2.6 | 0.4 / 6.5 / 590 |

  * F is a constant factor on E, i.e. the model error. F = 2 or 0.5 is ±1 sigma. At F = 3.2 or
    1/3.2 (90%), E = 1 needs 8 or ~10^8 CPU-years. At F = 0.25 it is never reached within
    S <= 2 S_min.
  * **So the CPU for E = 1 is undetermined between ~8 and ~10^8 CPU-years (90%), with 40 to
    6x10^4 at ±1 sigma.** E(C) rises only 0.18 / 0.30 / 0.45 / 0.52 per decade of CPU from 1 to
    1,000 CPU-years, which is why a constant factor on E moves the CPU for E = 1 so much.
  * The machine has 4 cores, so 1 CPU-year is ~3 months of the whole machine.
  * The model lens's frontier (monotone P only) gives E = 0.03 / 0.08 / 0.21 / 0.48 at
    1 / 10 / 100 / 1,000 CPU-years. With the other assignments (x4.5) its E = 1 falls at ~500-3,000
    CPU-years, inside the band above.
* **Speedups, judged by ratios at a fixed budget.** These do not depend on F:

  | change | 1 CPU-yr | 10 | 100 | 1,000 |
  |---|---|---|---|---|
  | 2x faster at all N | x1.28 | x1.22 | x1.19 | x1.14 |
  | 10x faster at all N | x2.17 | x1.90 | x1.71 | x1.47 |
  | 100x faster at all N | x4.1 | x3.2 | x2.5 | x1.9 |
  | t ~ N^3 beyond N = 4,000 | x1.04 | x1.11 | x1.21 | x1.25 |
  | t ~ N^2 beyond N = 4,000 | x1.11 | x1.35 | x1.60 | x1.68 |

  A constant speedup helps most on small budgets, and a better exponent in N helps most on large
  ones.
* **Cost per expected magic square by N band** (current build, base assignments):

  | N band | < 2k | 2-4k | 4-8k | 8-16k | 16-32k |
  |---|---|---|---|---|---|
  | E available | 0.03 | 0.14 | 0.39 | 0.70 | 0.85 |
  | CPU-years per expected magic square | 160 | 1.4x10^3 | 1.4x10^4 | 1.9x10^5 | 4.5x10^6 |

  The cost rises x9-24 per octave of N, i.e. a cost per magic square of ~N^3.2-4.6, while each
  octave from 4k to 64k holds a similar 0.4-0.9 expected squares. **So the N-exponent of the search
  time decides how much of the existing expectation is reachable.**
* **For per-version retrospectives, use these ratios, not absolute E.** The time model is also
  less firm than its fit suggests:
  * It came from `msearch_sample`, whose sources were copied at 23:10 on Oct 7, before the round-2
    commit 0f6fa1f (8-9% faster natively and 12-15% on cascadelake at N <= 3k, README).
  * It measures wall-clock time on a shared 4-core machine, where node rates within one sum vary
    0.34-1.6 M nodes/s.

  Compare builds by node counts or by controlled benchmarks, then convert to E with the ratio table.
* **The scheduler as fitted would under-collect:**
  * Its square model has no growth beyond N = 5,540. It predicted 0.08 squares where 60 were found
    at N = 15-32k, and 88 where 484 were found at N = 5-8k.
  * Its diagonal model under-predicts SP traversals ~2.3x (`research/calibration.md`).
  * Its candidate pool omits most exponent orders.

  The forecast's E_max = 0.022 reflects these limits, not the expected number of squares.
  Replacing the scheduler's E(P, S) by the analytic model and widening the pool would make the
  frontier above attainable. The analytic model is per-(P, S), cheap (~0.05-1 s per P) and tested
  to N = 45k.

## 7. What would test the heuristic further

1. **SP traversal rates at k = 7-8 and N = 10-40k.** The diagonal calibration has no data where 37%
   of the yield lies. The sum + 3-4 exponent sub-events of `research/calibration.md` give 10-50x
   more events than SP, so a few thousand squares from such sums would pin the N-trend and the
   k >= 7 factor.
2. **Sampled sums at N >= 10^5** (tau 3-5x10^4, k = 7-8; ~10 sums at ~10 CPU-min each). This
   regime decides how large the saturated total is (35% of it lies beyond 45k) and the large-N
   slope.
3. **The model lens's population with labels refitted on the 9 large-N sums.** This would test
   whether its divergence comes from the label model.
4. **5x5 exhaustive search at S <= 1,000 over all exponent assignments.** It predicts ~56
   semi-magic squares for a few CPU-minutes and tests the growth exponent, not just the level.
5. **The timing law at N = 10-50k for each build, in node counts.** The exponent of t(N) sets the
   slope of E(C).
6. **7x7 at S <= 440** for the best P near the predicted crossing. Under the model a 7x7 magic
   square smaller than Miquel's probably exists (E ~ 6 with S <= 465). Finding one would be a
   positive test of the whole machinery, though a weak test of its level (§5).

## 8. Assumptions and approximations

1. Each (P, S) is an independent Poisson trial over generic grids with distinct entries.
   * Multiples cM and structured constructions are excluded.
   * A structural obstruction, as for n = 3 and 4, is invisible to the heuristic.
2. The local CLT is used far from its regime (6-term lines, skewed cells). One smooth factor K(ld)
   absorbs this.
   * It was fitted on n = 5-7 with N <= 8k, and is tested out of sample to N = 45k (n = 6) and to
     S = 700 (n = 5).
   * At N >= 12k the analytic model is scaled by the measured 1.23.
   * Beyond N = 45k (35% of the total) the square counts are an extrapolation.
3. P(magic | semi) = D_6 kappa p_pair treats the two diagonals as conditionally Gaussian with
   per-prime lattice factors.
   * kappa = 0.57 is calibrated for k <= 6, N <= 9k and x >= 0.1. The explicit factors of §3.4
     extend it, with a range skewed upward at large N.
   * No SP+SP event has been observed, so the final step is untested.
4. **Coverage.** The universe is sorted exponents with at most one gap, times R(X):
   * R(X) comes from 24 multisets drawn by yield (x/÷1.1) and is extrapolated beyond scale 1.35;
   * primes to 47 are covered, and primes >= 53 are added from a 40-P sample;
   * only P with tau >= 1,000 are counted; smaller tau contributes nothing, which was checked.
5. **Corrections applied to the perm-multiplied curve.** The corrections were computed on the
   sorted-universe sample and applied as ratios to the published curve, which includes the
   permutation multiplier. This assumes the permuted assignments have the same N, k and x mix.
6. **The frontier.** It assumes ideal ordering, independent sums and E_all(C) = m E_U(C/m) for the
   extra assignments. The time per sum has a scatter of ~x/÷1.4.

## 9. Files (in `research/existence/scripts.tar.xz`)

| file | content |
|---|---|
| `synth/perP.py`, `synth/persum.py` | both models evaluated on the same P and S |
| `synth/bigN_an.py` -> `bigN_an.json` | analytic vs observed for the 17 sums with N >= 3k |
| `synth/assign_an.py`, `pick_ms.py`, `ratio.py`, `assign_an_*.json` | the exponent-assignment coverage R(X) |
| `synth/trav_an.py`, `datasci/trav_bins.py` | traversal rates by x, N and k |
| `synth/cand.py`; `synth/runs/p8000.*`, `m8000.*`, `q6015.*` | the two new sampled sums |
| `synth/rev/segs.py` -> `segs.json` | full S-profiles of the 2,809 sampled P, with N, x and k |
| `synth/rev/rev.py` -> `rev.out`, `variants.out`, `curve_*.json` | the revised curve, the factor table, the N-band shares and the low/high corners |
| `synth/rev/n53.py` | where the P·p yield sits (same N as for P, at ~2x the S) |
| `synth/rev/quant.py` -> `quant.out` | the smallest-S quantiles, P(exists) and the bands |
| `synth/rev/frontier_rev.py` -> `frontier_rev.out` | the cost frontier, speedup ratios and N-band costs |
| `skeptic/*`, `datasci/*` | the reviews' recomputations (Ehrhart counts, Dirichlet bound, added primes, trends, N shares) |

CPU used: the first synthesis used ~22 CPU-min of sampled search, ~0.5 min of enumeration and ~20 min
of model evaluation. This revision used ~3 CPU-min of model evaluation and no search, one heavy
process at a time.

## 10. Changes after review

* **The central numbers moved.**
  * N_magic(3,000) went from 0.60 to 0.74, and N_magic(inf) from ~3.2 to ~4.6 (68% range 2.2-13).
  * The median smallest S went from ~3,000-3,300 to ~2,800 (2,650-2,950 across models).
  * X(N = 1) went from 2,900-4,150 to 2,900-3,600.
  * P(exists) is now one value, ~0.95 (90% range 0.75-1), conditional on no obstruction.
  * The sources are listed in §4.1:
    * coverage beyond scale 1.35;
    * primes >= 53;
    * the untruncated universe tail (decay 0.70, not 0.6, per half-octave);
    * x1.23 squares at N >= 12k;
    * the explicit diagonal factors.
* **Cost.** The point estimate "E = 1 needs 2x10^3-10^4 CPU-years" was replaced by E at fixed
  budgets with bands, plus a statement that the CPU for E = 1 is undetermined (~8 to ~10^8
  CPU-years at 90%). Speedups are now given as ratios that no constant factor on E affects, and the
  provenance of the time model is noted for per-version retrospectives.
* **The models' disagreement was re-attributed.** At X <= 4,000 the gap is model form within the
  measured N range, not squares per sum beyond 32k. The model lens is disfavoured beyond X ~ 6,000,
  not beyond 3,000, at ~2 sigma rather than ~3. Its label explanation is now called a hypothesis.
* **Finiteness.** The bound is restated for distinct-entry (generic) grids; without distinctness
  degenerate grids make the count diverge. Saturation is labelled a property of the heuristic.
  Latin-square families at n = 8, 9 and 16 are noted. The exact-Euler-factor Dirichlet numbers are
  added.
* **Calibration.**
  * Results are split into fit and test, and two-sided tests are separated from consistency checks.
    The n = 7 and 8 magic records are relabelled as nearly powerless.
  * The 5x5 magic "5x10^-6" is no longer claimed as a success.
  * The within-P N slope and the first search's SP+S and SP+P pair test are added; the report's
    4.58 and 1.46 lack a factor 2.
  * `research/calibration.md` is cited.
* **Diagonals.** kappa's validity range is stated. The N-, k- and x-drifts and the pair factor are
  carried explicitly. The within-P slope of P(magic | semi) is given by x range: -20 to -25 near
  S_min, not -8.
* **Consistency fixes.**
  * R(X) uses 24 multisets, not 28.
  * Square counts are given per dataset (6,598 / 6,949 / 6,957 / 7,021).
  * The curves are renamed "analytic" and "model lens".
  * The model error depends on X (sigma 0.7 at X <= 4,000, rising to 0.75 / 1.0 down / up).
  * P(exists) is stated once.
