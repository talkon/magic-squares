# How many 6x6 additive-multiplicative magic squares exist with magic sum S <= X? (model-based lens)

Written 2026-10-07. This is the second version. The first version (kept as `model_v1_before_largeN.md`)
had to extrapolate the number of semi-magic squares per sum from N <= 8,891 vectors to the
N = 2x10^4-10^5 where most of the expected magic squares are. This version measures that
range directly (section 2): 7 sums with N = 15,199-31,743, using a sampled version of the
search. The central model was recalibrated on those sums.

Everything below is an *expected* count from the models described here. None of it is a
search result, apart from the calibration runs in sections 2 and 4.
Scripts: `model/` in `research/existence/scripts.tar.xz` (section 12); the data are not in the repo.

## 1. Bottom line

**N_magic(X)** is the expected number of 6x6 magic squares with magic sum S <= X, over all P.
Central model "H_corr mid" (sections 3-4):

| X | 1,000 | 1,500 | 2,000 | 2,500 | 3,000 | 4,000 | 5,000 | 7,000 | 10^4 | 2x10^4 | 3x10^4 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **N_magic(X), central (all P)** | 0.0065 | 0.058 | 0.22 | 0.54 | **1.2** | 3.2 | 6.2 | 18 | 62 | 400 | 990 |
| low / high calibration at large N (all P) | 0.007 / 0.006 | 0.06 / 0.06 | 0.19 / 0.25 | 0.43 / 0.71 | 0.84 / 1.7 | 2.0 / 5.5 | 3.4 / 12 | 8 / 42 | 24 / 180 | 110 / 1,600 | 230 / 4,800 |
| central, counting only sums with N <= 32,000 (the measured range) | 0.0065 | 0.058 | 0.22 | 0.51 | 0.94 | 1.8 | 2.4 | 3.4 | 5.4 | 7.9 | 9.4 |
| central, P with non-increasing exponents only | 0.0025 | 0.018 | 0.058 | 0.13 | 0.25 | 0.67 | 1.4 | 3.7 | 9.9 | 49 | 106 |
| cross-check: scheduler model + anchored N^2.25 growth, x1.18 at large N (all P) | 0.006 | 0.05 | 0.16 | 0.43 | 1.05 | 4.0 | 12 | 48 | 170 | 1,500 | 4,700 |
| scheduler model as fitted (clamped at N = 5,540, no growth), monotone P | 0.0016 | 0.004 | 0.006 | 0.007 | 0.009 | 0.013 | 0.016 | 0.022 | 0.030 | 0.044 | 0.052 |
| first version's central (H x 0.8, all P) | 0.0056 | 0.072 | 0.26 | 0.69 | 1.6 | 5.1 | 13 | 48 | 160 | 1,400 | 4,100 |

* **Growth with X.** The local exponent d ln N_magic / d ln X is 5.5 at X = 1,200, 4.9 at
  1,500, 4.2 at 2,000 and 4.0 at 3,000. It falls to 3.0-3.2 at 5,000-10^4 and to ~2.2 at
  2x10^4. That last value is too low, because the population is truncated at
  S_min(P) <= 3x10^4.
  * Doubling X near X_1 multiplies the count by ~16.
  * Near X_1, N_magic(X) ~ 1.2 (X / 3,000)^4.
* **N_magic = 1 at X_1 ~ 2,900** (central). Read as the first point of a Poisson process, the
  smallest magic square has its magic sum at a **median ~2,650**, with a Poisson 10-90% range of
  1,700-3,600.
  * The large-N calibration band (low/high rows) moves X_1 to 2,700-3,150 and the median to
    2,500-2,850.
  * All model factors combined (section 10), including that calibration band, are a factor
    ~2.5 (1 sigma) on E. This moves X_1 to 2,300-3,700 (1 sigma) and to ~2,000-4,400 (90%).
  * Mixing that model error into the Poisson spread gives a median smallest S of **~2,600**, a
    10-90% range of **1,600-4,200** and a 5-95% range of 1,350-4,900.
  * Given the model, the probability that a magic square exists is ~24% for S <= 2,000, ~66%
    for S <= 3,000, ~96% for S <= 5,000 and ~100% for S <= 10^4.
* **What changed since the first version.** The first version's central extrapolation (H x 0.8)
  over-predicted the measured large-N sums by 2.1x.
  * Recalibrating lowers N_magic(3,000) from 1.6 to 1.2. X_1 moves from 2,700 to 2,900.
  * X_1 now depends only weakly on extrapolation. Near X_1, **84% of the expected squares are
    in sums with N <= 32,000**, and squares per sum have now been measured up to that N.
  * Counting only those sums, N_magic(3,000) = 0.94 and N_magic(5,000) = 2.4.
* **The scheduler's fitted model (no growth beyond N = 5,540) is ruled out.** On the 7 sampled
  sums with N = 15k-32k it predicts 0.08 squares in the samples, against **60 found**.
  * The "same-x N^2.25" variant predicts 0.9.
  * The H model predicts 154; it over-predicts by 2.6x and grows a little too fast with N.
  * The anchored N^2.25 growth model (region verifier) predicts 51, but with a per-sum scatter
    of 0.16-3.3x.
* **Divergence.** In every variant that fits the data, N_magic(X) still grows like X^2-3 at
  X = 10^4-3x10^4. The total over all P diverges, because the per-P totals E(P) keep growing
  up to P ~ 10^26 (S_min ~ 10^5) and then fall only slowly, while the number of P grows faster.
  * The divergence lives entirely at N >> 3x10^4, which nothing has measured.
  * Counting only sums with N <= 32k, N_magic levels off near ~10.
* **Where they are (central, X ~ X_1):**
  * P with tau = 3,000-15,000 (93%), 5-7 distinct primes (k = 6 most often) and
    S_min = 1,000-2,500 (92%). About 3/4 of the yield comes from P whose exponents are not
    in decreasing order on 2, 3, 5, ... (section 7).
  * Sums at S = 1.2-1.65 S_min (x = ln(S/S_min) = 0.2-0.5: 74%), with N = 3k-50k
    (N 10k-20k: 35%, 20k-32k: 23%, 32k-50k: 13%, 3k-10k: 24%).
  * Each such sum holds ~10^3 semi-magic squares (250-7,200 measured), each magic with
    probability ~10^-10 (geometric mean 1.9x10^-10).
  * There are ~2x10^10 semi-magic squares with S <= 2,900 in all. The two searches have found
    ~7.6x10^5 semi-magic squares in total.
* **Cost with the current search.** Measured here, a sum costs t(N) ~ 51 s (N/4000)^3.9,
  steepening to ~N^5.2 above 2x10^4. That is 7,000-50,000 CPU-s at N = 15k-23k and
  ~375,000 CPU-s at N = 31,743. Each semi-magic square costs 4-80 CPU-s at N = 15k-32k.
  * The frontier assumes an ideal order of sums and P with non-increasing exponents only. It
    reaches E = 0.03 after 1 CPU-year, 0.08 after 10, 0.21 after 100 and 0.47-0.49 after 1,000.
  * **E = 1 takes ~7,000-14,000 CPU-years**: 2,500 with the high calibration and 34,000 with the
    low one. The marginal cost there is 1 magic square per ~20,000-40,000 CPU-years.
  * Adding the other exponent assignments could at best cut this ~15x, to ~500 CPU-years, if
    their cells were as cheap per unit E as the monotone ones. They are not, so the estimate
    is **~10^3-10^4 CPU-years (central ~3,000) for one expected magic square**. That is about
    10^9-10^10 semi-magic squares searched, at ~25 CPU-s each.
  * The best single sums measured have a density of one magic square per ~500-900 CPU-years:
    13 7 4 3 1 1 at S = 2,200 and 11 6 4 3 2 1 at S = 2,250. But each such sum holds only
    ~10^-6 expected magic squares.

## 2. New measurement: semi-magic squares per sum at N = 15k-32k (sampled search)

**Method.**
* The arrangement search splits into independent subproblems, one per first row r1 (the root loop
  in `arrange_core.h`). A square is found exactly once, under its lowest-index vector. The
  vectors are sorted by their largest label, so the other rows and columns come from later
  vectors.
* A copy of the C sources in `sample/src/` (the repo is unchanged) runs only
  r1 = off, off + k, off + 2k, ..., with a random offset, and logs squares, nodes and time per
  r1. k times the sample count is an unbiased estimate of the sum's square count and CPU time.
  The standard errors use the per-r1 counts (simple random sampling formula).
* One sum (S = 2650) was sampled in two strata: r1 < 5,164 at 1/600 and the rest at 1/3,000.
  The 1/600 run was stopped after 8 subproblems, because it would have exceeded the CPU budget.
  The cost per r1 turned out to fall ~1,000x from the first r1 to the last.
* **Validation.**
  * With k = 1 the output is identical to `bin/msearch`: 81 sums of 12 6 3 2, the same 6
    squares and 70,030,459 nodes.
  * On a fully searched sum (12 6 3 2 1 1, S = 988, N = 6,671: 245 squares in 290 s), two
    independent 1/10 samples estimated 230 and 270 squares and 277 s.

**Results** (predictions are for the whole sum; "est." is the sample count / sampling fraction):

| P | S | x | N | labels | fraction | squares in sample | est. squares | H raw | H x0.8 (v1 central) | anchored N^2.25 x0.84 | scheduler as fitted (x0.35) | est. CPU-s for the sum |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 7 4 3 1 1 | 2200 | 0.25 | 15,199 | 211 | 1/100 | 28 | 2,800 +- 670 | 5,500 | 4,400 | 940 | 1.8 | 12,400 |
| 12 9 6 2 1 1 | 3500 | 0.25 | 17,142 | 256 | 1/250 | 1 | 250 +- 250 | 900 | 720 | 1,570 | 1.8 | 7,200 |
| 14 7 4 4 1 0 0 1 | 3600 | 0.24 | 19,346 | 247 | 1/300 | 3 | 890 +- 660 | 4,330 | 3,460 | 2,080 | 1.8 | 14,900 |
| 11 6 4 3 2 1 | 2250 | 0.28 | 20,825 | 234 | 1/400 | 18 | 7,200 +- 2,350 | 12,400 | 9,900 | 2,170 | 0.75 | 41,800 |
| 9 6 4 3 1 1 1 1 | 2700 | 0.13 | 20,896 | 259 | 1/400 | 4 | 1,600 +- 780 | 4,600 | 3,680 | 6,220 | 19 | 49,900 |
| 13 7 4 3 1 1 | 2400 | 0.33 | 22,992 | 245 | 1/600 | 3 | 1,800 +- 1,020 | 11,900 | 9,530 | 2,390 | 0.23 | 46,700 |
| 13 7 4 3 1 1 | 2650 | 0.43 | 31,743 | 279 | 1/600 + 1/3000 | 3 | 4,600 +- 3,400 | 29,100 | 23,300 | 4,950 | 0.02 | 375,000 |
| **pooled (expected count in the samples)** | | | | | | **60** | | 154 | 123 | 51 | 0.08 | |

**What this shows.**
* **Squares per sum keep growing past N = 10^4.** Over all 17 sums with N >= 3,000 (10 full,
  7 sampled), squares go as N^(3.0 +- 0.4) across P. There is no x trend (-0.05 +- 0.95) and
  large P-to-P scatter: overdispersion 21, i.e. a factor ~2-3 per sum.
  * The scheduler model's clamp at N = 5,540 is off by ~700x at these N.
  * Its "same-x N^2.25" form is off by ~70x.
* **H (the first version's central model) over-predicts at large N**, increasingly with N.
  * In-sample the ratio obs/H was 0.16-1.4 for N < 3,000 (rising with N) and 0.62-0.67 at
    N = 3k-9k.
  * Over all 17 sums with N >= 3k, a Poisson fit with overdispersion gives
    **obs/H_raw = 0.64 (N/10^4)^(-0.27)**, with se 0.17 on ln of the level and 0.35 on the
    slope. An extra x term is not significant (+0.45 +- 0.92).
  * This is 0.77 at N = 5k, 0.53 at 2x10^4, 0.41 at 5x10^4 and 0.34 at 10^5 (x/ 1.2, 1.4,
    2.0, 2.5).
  * The per-sum ratios at N >= 12k are 0.15-0.58. Within 13 7 4 3 1 1, squares grow only from
    ~2,800 to ~4,600 between x = 0.25 and 0.43, while N doubles.
* **The traversal model holds at large N with new P** (tau up to 22,400, k up to 8, with the
  unclamped p_S, p_P used here).
  * On the 60 sampled squares: S-traversals 33 vs 27.3 predicted, P-traversals 7 vs 5.0.
  * On all squares at N >= 3k: S 649 vs 568, P 126 vs 121.
  * The clamped variant predicts S 23.7 and P 7.1 on the 60 squares, so these data do not
    separate the two.
* **Time per sum.** Over the 17 sums: ln t = 3.34 + 4.33 ln(N/4000), with residual sd 0.35.
  * With the scheduler's t(4000) = 51.3 s fixed, the exponent is a = 3.92.
  * Above N = 12k the local exponent is 5.2. Sums with > 256 distinct numbers, which use the
    N x N matrix path, are slower: N = 31,743 took 375,000 CPU-s, while 3.92 alone predicts
    172,000.
  * Seconds per semi-magic square: 1-6 at N < 9k, **4-80 at N = 15k-32k**.
* CPU used in this version: ~900 s of sampled search and validation, plus a few seconds of
  enumeration. Together with the first version's ~1,300 s of msearch and ~2 min of
  enumerate, the lens used ~39 CPU-min of search and enumeration.

## 3. Pipeline

N_magic(X) = sum over P of sum over S <= X of E(P, S), with E(P, S) = squares(P, S) x 5400 (rho p_S p_P)^2.

1. **Population of P.**
   * All exponent tuples over the primes up to 43 with tau >= 500 and 6 P^(1/6) <= 30,000,
     with exponents non-increasing in the prime ("monotone"): 105,345 tuples, of which 102,514
     have N >= 300 at some sum.
   * Then all other assignments of the same exponent multisets to the first 16 primes, with
     P up to 64x larger (S_min up to 2x), by importance sampling (section 7).
   * S_min(P) = 6 P^(1/6) (1 + 0.2-1%) analytically. This matches bin/enumerate within +-1 on
     the P checked.
2. **Vector counts N(S): analytic** (`nmodel.py`).
   * Take a random ordered 6-tuple of divisors with product P. Its log-sizes are approximately
     Gaussian, with variance s2 = sum_i (ln p_i)^2 5 e_i (e_i + 6)/252 and correlation -1/5.
     Then N(S) = c (V/720) dF/dS, with V = prod_i C(e_i + 5, 5) and F the tabulated
     distribution of ln mean e^v.
   * The correction c (distinctness etc.) is fitted on 10,703 windows of 3,000 P. In practice
     it is 0.67-0.81, with **residual sd 4.5%**. It is checked against full counts of 7 P
     (`counts/`) and the counts of the sampled sums: analytic vs actual N is 15,090 vs 15,199,
     19,349 vs 20,825, 20,332 vs 20,896 and 33,070 vs 31,743.
   * Over 158 P with counts to 2.5 S_min, E(P) from the analytic N vs the real counts has a
     median ratio of 1.10 and sd(log) 0.18 (Jensen).
3. **Labels L(S)** (`lmodel.py`): distinct divisors near S/6, calibrated to the searched
   sums, with **residual sd 3.8%**.
4. **Semi-magic squares per sum.**
   * **H model** (`hmodel.py`): ln mu = -11.02 + 0.46 ln H(N, L) - 0.58 x + 2.56 ln N, where
     H = (L)_36 / (2 x 720^2) (N / C(L,6))^12 is the random-hypergraph count. It was fitted by
     Poisson ML on 21,474 sums with N <= 8k.
   * **Central "H_corr mid"** = H_raw x h(N), with h = 0.8 for N < 3,000 and
     **h = 0.64 (N/10^4)^(-0.27) for N >= 3,000** (section 2).
     * "low": h = 0.54 (N/10^4)^(-0.62). "high": h = 0.76 (N/10^4)^(+0.08). These are
       roughly the 1-sigma corners of the fit.
     * In the region that matters (N = 10^4-5x10^4), h is 0.41-0.64 (low 0.2-0.54, high
       0.76-0.86). The first version used 0.8.
   * **Cross-check "base_cal"**: the scheduler's model (model_6.json, tau/k/p clamped) x 0.35
     for fresh P. Beyond N = 4,000 it grows as (N/4000)^2.25, anchored where N first reaches
     4,000, x 0.84. Here it is multiplied by a further 1.18 (its pooled ratio at N >= 12k).
5. **P(magic | square) = 5400 (rho p_S p_P)^2**, with rho = 2.6.
   * p_S and p_P come from the scheduler's log-linear fit, **not clamped in tau, k, p** (x
     clamped at 1.1).
   * Evidence:
     * k = 8 P (25 squares): 0 P-traversals vs 0.6 unclamped and 4.4 clamped.
     * tau > 6,500: 120 S-traversals vs 117 unclamped and 94 clamped.
     * The per-prime argument: each extra prime with exponent 1 multiplies p_P by ~264/720;
       the fit gives 0.36.
     * The new large-N data are neutral.
   * Clamping (as the scheduler does) gives 2-4x more E at X ~ 2-3k (section 10).
6. **Large S.**
   * For S > 3 S_min the fit gives p_S ~ S^-1.13 (~c/S) and p_P ~ S^0.18 (flat). The legacy
     exhaustive P confirm the traversal totals: S 1.06, P 1.04.
   * This never matters. The exhaustive runs (12 6 3 2 to 20x S_min, 11 4 3 2 1 to 10x)
     found no squares beyond x ~ 1.0, although N is near its peak there. The models put < 1%
     of E beyond x = 0.8.

## 4. What the models were checked against

| check | data | result |
|---|---|---|
| exhaustive runs (msearch, all S up to 20x / 10x S_min) | P = 12 6 3 2 (71 CPU-s) and 11 4 3 2 1 (192 CPU-s) | 42 and 74 squares, identical to the legacy totals. All at S/S_min = 1.35-2.7, none beyond x = 1.0 |
| legacy per-P totals (1,136 P, 758,949 squares; analytic N over each P's searched range) | out of sample for H | obs/pred: scheduler 1.30 (trained on them), H 1.44; per-P sd(log) 0.38 / 0.40 |
| 40 fresh scheduler units (near S_min) | 901 squares | scheduler 0.35, H 0.88 |
| 5 fresh P with tau 10.9-22.4k (x < 0.06) | 92 squares | scheduler 0.25, H 0.50 |
| 10 full sums at N 3,418-8,891 | 860 squares | scheduler clamped 2.9 (obs/pred), same-x N^2.25 2.7, anchored N^2.25 0.84, H_raw 0.73 |
| **7 sampled sums at N 15,199-31,743 (section 2)** | 60 squares in the samples | **scheduler clamped 740, same-x N^2.25 69, anchored N^2.25 x0.84 1.18 (per sum 0.16-3.3), H_raw 0.39 (per sum 0.15-0.58), H_corr mid 0.74** |
| traversal rates (all squares at N >= 3k, unclamped model) | 649 S-, 126 P-traversals | S 1.14, P 1.04; the sampled squares alone 1.21, 1.39 |
| legacy exhaustive P (293 P, 53,501 squares) | traversal totals | S 1.06, P 1.04; SP 37 vs 53 predicted at rho 2.6 (rho_eff ~2.2) |
| expected magic squares in what has been searched | legacy + ours | 0.004-0.006 (legacy) + 0.0001 (ours) + 1.3x10^-8 (the 60 sampled squares here). Finding none had probability ~99.5% |

H_corr is H with the large-N factor. The factor was fitted on the 17 N >= 3k sums, so the
0.74 (0.73 exactly) is in-sample. Over all 17 sums the fit matches the total (920 vs 919).
At N >= 12k it is below 1 because the power law in N also has to fit the high-count full
sums at N = 3k-9k.

## 5. Per-P expected counts (sample across sizes, central H_corr mid)

E(P) is summed over all S. x50 is where half of E(P) is reached. "share > 1.6" is the
fraction of E(P) above 1.6 S_min.

| P | tau | k | S_min | N_peak | **E(P)** | S50/S_min | N at x50 | share > 1.6 | semi-magic squares (all S) | P(magic/square), E-weighted |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 4 3 2 (Morgenstern) | 660 | 4 | 171 | 484 | 5.9e-9 | 1.66 | 369 | 0.65 | 0.04 (one known) | 1e-7 |
| 12 6 3 2 | 1,092 | 4 | 310 | 1,560 | 1.1e-6 | 1.56 | 893 | 0.41 | 25 (obs 42) | 4e-8 |
| 11 4 3 2 1 | 1,440 | 5 | 285 | 2,170 | 3.7e-7 | 1.50 | 1,120 | 0.32 | 26 (obs 74) | 1e-8 |
| 13 7 3 2 1 | 2,688 | 5 | 620 | 9,190 | 2.8e-5 | 1.54 | 3,990 | 0.37 | 2.0e4 | 1e-9 |
| 14 7 5 3 | 2,880 | 4 | 1,104 | 10,100 | 9.1e-5 | 1.59 | 4,310 | 0.49 | 9.4e4 | 1e-9 |
| 12 6 3 2 1 1 | 4,368 | 6 | 706 | 20,700 | 6.8e-5 | 1.47 | 8,160 | 0.29 | 2.3e5 | 3e-10 |
| 13 7 4 2 0 1 | 3,360 | 5 | 834 | 13,700 | 5.9e-5 | 1.56 | 6,130 | 0.43 | 1.0e5 | 6e-10 |
| 12 8 4 2 1 | 3,510 | 5 | 868 | 15,200 | 8.9e-5 | 1.52 | 6,070 | 0.36 | 1.5e5 | 6e-10 |
| 13 7 4 3 1 1 | 8,960 | 6 | 1,718 | 83,100 | 1.3e-3 | 1.53 | 32,300 | 0.38 | 6.1e7 | 2e-11 |
| 11 6 4 3 2 1 | 10,080 | 6 | 1,693 | 87,800 | 1.1e-3 | 1.52 | 35,900 | 0.36 | 5.0e7 | 2e-11 |
| 12 9 6 2 1 1 | 10,920 | 6 | 2,731 | 106,000 | 2.1e-4 | 1.58 | 44,900 | 0.47 | 4.2e7 | 5e-12 |
| 14 7 4 4 1 0 0 1 | 12,000 | 6 | 2,841 | 124,000 | 7.6e-4 | 1.58 | 52,000 | 0.45 | 1.5e8 | 5e-12 |
| 9 6 4 3 1 1 1 1 | 22,400 | 8 | 2,360 | 264,000 | 2.3e-4 | 1.47 | 94,600 | 0.27 | 3.3e8 | 7e-13 |
| 10 5 3 2 1 1 1 1 1 | 25,344 | 9 | 2,057 | 261,000 | 2.8e-5 | 1.44 | 100,000 | 0.21 | 1.1e8 | 3e-13 |
| 14 8 5 3 2 1 1 | 38,880 | 7 | 7,239 | 1.1e6 | 0.030 | 1.52 | 4.2e5 | 0.36 | 2.1e11 | 1e-13 |
| 16 9 5 3 2 1 1 1 | 97,920 | 8 | 17,890 | 5.3e6 | 0.042 | 1.49 | 1.6e6 | 0.35 | 1.3e13 | 3e-15 |
| 18 10 6 4 3 2 1 1 | 351,120 | 8 | ~112,000 | 4.6e7 | 0.20 | 1.53 | 1.5e7 | 0.40 | 7.4e15 | 3e-17 |
| 20 12 7 5 3 2 2 1 1 | 1.9e6 | 9 | ~995,000 | 6.6e8 | 0.11 | 1.55 | 2.2e8 | 0.42 | 4.5e18 | 3e-20 |
| 24 14 8 6 4 3 2 2 1 1 | 1.7e7 | 10 | ~2.7e7 | 1.9e10 | 0.011 | 1.62 | 6.1e9 | 0.55 | 5.7e21 | 2e-24 |

The S = 2,650 sum of 13 7 4 3 1 1 sampled in section 2 is the median-yield sum of that P:
x50 = 0.42. Its measured ~4,600 squares, each magic with probability 6.4x10^-11, carry
3x10^-7 expected magic squares for ~375,000 CPU-s.

**How E(P) scales with the size of P** (monotone population, 102,514 P with S_min <= 3x10^4):
* OLS: **ln E = -46.7 + 16.4 ln tau - 9.85 ln S_min - 6.5 k** (residual sd 3.2). Equivalently
  E ~ tau^16.4 P^-1.64 e^(-6.5 k).
  * At fixed size, more divisors help enormously.
  * At fixed tau and size, each extra distinct prime costs a factor ~650, through p_P^2 and
    through L.
* The largest E(P) per S_min class rises steeply, then flattens:
  * 6x10^-7 (S_min 150-300), 1.4x10^-5 (300-500), 2.2x10^-4 (500-1,000), 9.5x10^-4
    (1,000-1,500), 1.7x10^-3 (1,500-2,000), 6.4x10^-3 (2,000-3,000), 0.023 (3-5k), 0.067
    (5-10k), 0.15 (10-20k), 0.21 (20-30k).
  * Along the family in the sample table it reaches ~0.2 at tau ~ 3.5x10^5 (S_min ~ 10^5,
    P ~ 10^26). It then falls slowly: 0.11 at S_min 10^6, 0.011 at 2.7x10^7. This is
    extrapolated far beyond any measured N.
* The median E(P) by tau is 3x10^-10 (tau 500-1,000), 1x10^-9 (1-4k), 3x10^-8 (4-8k),
  1x10^-6 (8-16k), 2x10^-5 (16-32k) and 2x10^-4 (> 32k). The top 1% of P in each class
  hold 18-26% of the class's E.
* The probability that a semi-magic square is magic falls from 10^-7 (S_min ~ 200) to 10^-24
  (S_min ~ 3x10^7). The number of semi-magic squares per P grows even faster: from 0.04 to
  6x10^21.
* **Number of P** (monotone, N >= 300) per S_min class:
  * 718 (150-300), 1,246 (300-500), 1,272 (500-700), 1,860 (700-1,000), 2,953 (1-1.5k),
    2,787 (1.5-2k), 5,190 (2-3k), 9,272 (3-5k), 8,305 (5-7k), 11,250 (7-10k), 16,600
    (10-15k), 14,813 (15-20k), 26,207 (20-30k).
  * That is roughly X^1.1 cumulatively. All assignments of the exponents to primes add
    400-3,000 variants per multiset, which carry 2.6-9.4x the monotone E (section 7).
* Class totals E (all S) per unit ln S_min are still rising at the top of the population:
  * 1.9 (S_min 2-3k), 7.2 (3-5k), 32 (5-10k), 117 (10-20k), 247 (20-30k), monotone.
  * Hence the divergence of the model total. It comes from P with N >> 3x10^4 at their good
    sums, where h(N) is extrapolated.

## 6. N_magic(X) by variant

The full table is in section 1. Monotone-P curves for the first version's variants are in
`model_v1_before_largeN.md`, section 5. Of those, "scheduler clamped (nogrowth)",
"samex" and "H, only N <= 10^4" are now excluded by the large-N data.

| variant (all P unless noted) | X at N_magic = 0.105 / ln 2 / 1 / 2.303 |
|---|---|
| H_corr mid (central) | 1,700 / 2,650 / **2,900** / 3,600 |
| H_corr low | 1,700 / 2,850 / 3,150 / 4,250 |
| H_corr high | 1,650 / 2,500 / 2,700 / 3,200 |
| H_corr mid, only N <= 32k | 1,700 / 2,750 / 3,050 / 4,800 |
| base_cal x1.18 | 1,800 / 2,750 / 2,950 / 3,550 |
| first version (H x0.8) | 1,600 / 2,500 / 2,700 / 3,300 |
| H_corr mid, monotone P only | 2,350 / 4,050 / 4,550 / 5,950 |

The two independent model families (structural H with the measured large-N factor; the
scheduler model with the measured anchored growth) agree to within 10% on X_1. They differ
by 3-5x on N_magic at X >= 10^4, where the extrapolation dominates.

## 7. Non-monotone P (all exponent assignments)

The monotone list misses P such as 13 7 4 2 0 1 or 14 7 4 4 1 0 0 1. These have the same
exponent multiset on other primes, the same V and tau, and an S_min a few % larger.

**Method** (`assign3.py`, run on H_corr mid):
* In each S_min class, sample 25 monotone P with probability proportional to E(P).
* Enumerate every distinct assignment of each one's exponents to the first 16 primes, with P
  up to 64x larger.
* Accumulate E(S <= X) for each assignment, and form the unbiased importance-sampling
  estimate of the class total. As a check, the same estimator reproduces the monotone totals
  within 0.95-1.10.

**Result.**
* All-assignment total / monotone total is 3.4 (S_min 150-300), 3.7, 5.0, 5.4 (700-1,000),
  5.8, 8.6 (1,500-2,000), 7.5 (2-3k), 10 (3-5k), 12, 14 (7-10k). The larger classes are
  truncated by the P limit.
* In terms of X: x2.6 (X = 1,000), 3.7 (2,000), **4.6 (3,000)**, 4.5 (5,000), 6.3 (10^4),
  9.4 (3x10^4).
* The E-weighted S_min shift is x1.09-1.22. Sampling error is ~+-20-30% per class.
* The low/high/capped variants use the mid ratio, and base_cal uses the first version's H
  ratio. The two ratios agree within ~20% at X <= 5,000.

## 8. Where the expected squares are

E-weighted shares of N_magic(X) (monotone P; the all-P picture is similar but shifted ~10-20%
to larger S_min):

| | X = 2,000 | **X = 2,900 (~X_1)** | X = 5,000 | X = 10^4 |
|---|---|---|---|---|
| N of the sum | 3-10k 51%, 10-20k 34%, 20-32k 7% | 3-10k 24%, **10-20k 35%, 20-32k 23%**, 32-50k 13%, 50-100k 3% | 10-20k 14%, 20-32k 18%, 32-50k 21%, 50-100k 31%, 0.1-0.3M 9% | 50-100k 23%, 0.1-0.3M 49%, > 0.3M 10% |
| x = ln(S/S_min) | 0.2-0.5: 77% | 0.1-0.2: 9%, 0.2-0.3: 21%, 0.3-0.4: 32%, 0.4-0.5: 21%, 0.5-0.6: 10%, > 0.6: 6% | 0.2-0.5: 71% | 0.2-0.5: 68% |
| tau(P) | 3k-6.5k 57%, 6.5k-10k 30% | 3k-6.5k 33%, 6.5k-10k 36%, 10k-15k 24% | 10k-15k 32%, 15k-30k 37% | 15k-30k 50%, > 30k 32% |
| S_min(P) | 1,000-1,500 65% | 1,000-1,500 32%, 1,500-2,000 46%, 2,000-2,500 15% | 2-3k 49%, 3-5k 33% | 3-5k 37%, > 5k 53% |
| k (distinct primes) | 5: 40%, 6: 41% | 5: 33%, **6: 43%**, 7: 17% | 6: 45%, 7: 24% | 6: 42%, 7: 32% |
| semi-magic squares with S <= X (monotone; all P ~4-6x) | 2.9e8 | 4.5e9 | 1.9e11 | 1.3e13 |
| geometric-mean P(magic/square) | 6.6e-10 | 1.9e-10 | 3.2e-11 | 3.7e-12 |

* Top monotone P near X_1: 13 7 4 3 2, 13 7 4 3 1 1, 12 7 3 3 2 1, 11 6 4 3 2 1, 12 7 4 3 2,
  12 6 4 3 2 1, 12 7 4 2 2 1. Each has ~0.8-0.9x10^-3 expected magic squares with
  S <= 2,900, ~0.4% of the monotone total.
* The yield is very spread out: hundreds of P each carry 10^-4-10^-3.

## 9. Cost with the current search

**Frontier** (`sample/frontier2.py`).
* It takes every (P, S) cell of the monotone population, sorted by expected magic squares
  per CPU-second, in any order of sums. That is ideal; the scheduler goes in increasing S
  per P.
* Time per sum is the scheduler's model up to N = 4,000, then the measured
  51.3 s (N/4000)^3.92 ("a = 3.92"). The "steep" variant uses slope 5.2 above N = 2x10^4,
  which matches the 375,000 s measured at N = 31.7k. There is a floor of 1 s per semi-magic
  square.

| model | E after 1 CPU-yr | 10 | 100 | 1,000 | 10^4 | CPU-yr for E = 0.1 / 0.5 / 1 | marginal cost at E = 1 | cells used for E = 1 (E-weighted) |
|---|---:|---:|---:|---:|---:|---|---|---|
| H_corr mid, a = 3.92 | 0.029 | 0.079 | 0.21 | 0.49 | 1.15 | 22 / 1,260 / **7,300** | 1 per 18,000 CPU-yr | N ~ 22k, x 0.30, tau ~ 11.7k, S_min ~ 3,300, S ~ 4,900 |
| H_corr mid, steep | 0.029 | 0.079 | 0.21 | 0.47 | 0.94 | 22 / 1,540 / **13,700** | 1 per 40,000 CPU-yr | N ~ 21k, S ~ 5,100 |
| H_corr low, a = 3.92 | 0.031 | 0.079 | 0.17 | 0.36 | 0.75 | 23 / 3,150 / 34,000 | 1 per 10^5 CPU-yr | N ~ 25k |
| H_corr high, a = 3.92 | 0.029 | 0.089 | 0.26 | 0.72 | 1.9 | 18 / 510 / 2,500 | 1 per 5,000 CPU-yr | N ~ 20k |
| first version, H x0.8, a = 3.5 | 0.034 | 0.11 | 0.35 | 1.04 | 2.9 | 11 / 280 / 1,200 | | |

* **Including non-monotone P.** With m ~ 4.5x as much E at similar cost densities,
  E_all(C) = m E_mono(C/m). E = 1 would then need ~520 CPU-yr (mid, a = 3.92), and ~550
  with the steep time model, which bites only above N = 2x10^4.
  * This is optimistic. The extra assignments are 400-3,000 per multiset, most with tiny E,
    so their cost per unit E is higher.
  * Realistic: **~10^3-10^4 CPU-years (central ~3,000) for one expected magic square**. The
    low-calibration corner needs 34,000 CPU-yr for monotone P.
* **The first CPU-years are cheap and stay near S_min.** E = 0.1 needs ~20 CPU-yr, in cells with
  N ~ 7k, x ~ 0.28, S ~ 2,600.
  * This frontier (ideal ordering) gives E(1 CPU-yr) ~ 0.03. The forecast's scheduler
    simulation gives 0.015-0.018.
* **Search time per semi-magic square at the relevant N**: 4-80 CPU-s at N = 15k-32k
  (median ~25 s), not ~1 s. The finds per CPU-hour fall from ~600-10,000 near S_min to
  ~40-800 at N = 15k-32k.
* The magic-square density of the best measured single sums: 13 7 4 3 1 1 at S = 2,200
  (2,800 squares x 2.9x10^-10, 12,400 s) and 11 6 4 3 2 1 at S = 2,250 (7,200 x 2.0x10^-10,
  41,800 s). Each gives one magic square per ~500-900 CPU-years. That is the marginal density
  near E ~ 0.1 on the frontier.

## 10. Uncertainty on N_magic(X) at fixed X (central H_corr mid)

| source | factor on E near X_1 | basis |
|---|---|---|
| large-N calibration h(N) (level and slope) | x0.72-1.47 | low/high rows, section 2 fit; P-to-P scatter included through overdispersion |
| rho (2.2-3.7 around 2.6) | x0.7-2.0 | legacy exhaustive SP (rho_eff 2.2) vs our 14 SP (~3.7) |
| p_S p_P at x ~ 0.2-0.5 (squared) | x0.55-1.75 | obs/pred: 1.14 and 1.04 at N >= 3k; 1.21 and 1.39 on the sampled squares; 1.06 and 1.04 legacy |
| clamping p_S, p_P in tau, k (scheduler's choice) | x1-3.8 (one-sided) | the first version's travclamp variant; data mildly favour unclamped |
| multiplicity sampling | x0.75-1.3 | importance-sampling spread |
| analytic N (Jensen) | x1.0-1.15 | 158-P comparison |
| h at N < 3,000 | negligible near X_1 (2.5% of E) | |
| model form: diagonal-pair correlation 1.1-1.44 (not included) | x1.1-1.4 (one-sided) | forecast |

Combined, this is a factor **~2.5 (1 sigma)**, i.e. sigma(ln E) ~ 0.9, and a factor ~4.5
at 90%. It is skewed upward by the one-sided items.
* N_magic ~ X^4 near X_1, so a factor 2.5 moves X_1 by x1.25: **X_1 ~ 2,900 (central),
  2,300-3,700 (1 sigma), ~2,000-4,400 (90%)**.
* The distribution of the smallest magic sum, with that error mixed in: median ~2,600,
  10-90% 1,600-4,200, 5-95% 1,350-4,900.
* P(a magic square with S <= 2,000) ~ 0.24, P(S <= 3,000) ~ 0.66, P(S <= 5,000) ~ 0.96.
* Counting only sums inside the measured N range (<= 32k): P(S <= 3,000) ~ 0.6,
  P(S <= 5,000) ~ 0.83-0.91, P(S <= 10^4) ~ 0.95.
* **Not covered by any of this**: a structural obstruction. For n = 3 and 4, add-mult magic
  squares are impossible for algebraic reasons, which no counting heuristic sees. The model
  also assumes that the SP+SP event is as predicted by the product of single-traversal rates.
  It has been checked only up to SP and SP+S / SP+P frequencies, which are few.

## 11. Assumptions and approximations (explicit)

* **Squares per sum.**
  * The validation range is now N <= 31,743, x <= 1, tau <= 22,400, k <= 8, primes <= 23.
  * Near X_1, 84% of E is inside it. At X = 10^4 only ~10% is, so N_magic(10^4) and above
    rest on the extrapolated h(N) = 0.64 (N/10^4)^-0.27.
  * The fit uses 17 sums from 9 P, and the sampled sums have 1-28 squares each. The per-sum
    scatter is large (factor 2-3). The population average is assumed to follow the fitted
    mean.
* **Sampling.** Systematic r1 sampling with a random offset is unbiased for the count and the
  time. Squares cluster weakly by r1: at most 3 per r1, and the SRS se is 1.0-1.4x the
  Poisson se.
  * The S = 2,650 sum used two strata. The cost per r1 falls ~1,000x from early to late r1,
    so a time-limited sample must not be truncated; that run was stopped and restratified.
* P(magic | square) = 5400 (rho p_S p_P)^2 treats the two diagonals as independent given the
  square, apart from rho. The pair correlation (1.1-1.44) is not included.
* E counts magic arrangements of semi-magic squares (SP-SP diagonal pairs), up to row,
  column and transpose symmetry. That is about the number of essentially different magic
  squares while it is small.
* Sums are treated as a continuum on an x grid of step 0.004. Only sums with N >= 300 can hold
  squares; the smallest known square is at N = 452.
* The population uses primes up to 43 (monotone) and 53 (assignments). Primes beyond 37 carry
  < 0.1%. tau >= 500, and S_min(P) <= 3x10^4, so curves above X ~ 1.5x10^4 are lower bounds.
* **Time model:** 51.3 s (N/4000)^3.92 beyond N = 4,000, with or without slope 5.2 beyond
  2x10^4. It is fitted to 17 sums with residual sd 0.35 (factor 1.4) per sum. The cost
  frontier assumes ideal ordering.

## 12. Files (`model/` in `research/existence/scripts.tar.xz`; data not in the repo)

* First version (unchanged):
  * `mlib.py`: scheduler-model evaluation.
  * `nmodel.py`, `Ftab.npz`, `ncal.*`: analytic N.
  * `lmodel.py`, `lfit.py`, `lcal.json`: labels.
  * `hmodel.py`, `hcoef.json`: the H model.
  * `epop.py`, `variants.py`, `popgen.py`, `runpop.py`, `pop_*_30000.json`: populations.
  * `assign3.py`: non-monotone multiplicity.
  * `frontier.py`: cost frontier.
  * `legcheck.py`, `calib.py`, `sampletab.py`: checks.
  * `runs/`: exhaustive runs and the N = 8,891 sum. `counts/`: full vector counts.
  * `model_v1_before_largeN.md`: the first write-up.
* This version, in `sample/`:
  * `src/` (patched copy of `src/c` at 23a34f7; only `arrange_core.h` changed, see
    `sampling.patch`: r1 stride/offset/min and a per-r1 log, via env SAMPLE_STRIDE,
    SAMPLE_OFFSET, SAMPLE_MIN, SAMPLE_LOG) and `msearch_sample`.
  * `t*.jsonl` / `t*.log`, `s2650_137_A.*`, `u2650_137_B.*`: the sampled runs.
  * `v988_*`: the validation samples. `ref_12632.jsonl`, `s1_12632.jsonl`: the stride-1
    identity check.
  * `est.py`: estimates from per-r1 logs.
  * `predict.py`: all model predictions for one sum.
  * `bigN.py`, `bigN_3000.txt`, `bigN_rows.json`: all N >= 3k sums vs models, the h(N) fit
    and the traversal check.
  * `corrfit.py`: correction forms.
  * `corr.py`: runs any model script with the H_corr variant.
  * `pop_Hc_{mid,low,high,cap32k}_30000.json`, `logs_Hc_*.txt`, `assign3_Hc_mid.txt`,
    `assign2_H_corr.json`: populations and multiplicity.
  * `final2.py`, `final2_curves.json`: the section 1 table.
  * `pertab2.py`, `pertab_mid.txt`: per-P table.
  * `percoef.py`: E(P) scaling.
  * `where2.py`, `where_mid.txt`: where the yield is.
  * `frontier2.py`, `front_*.txt`: costs with the measured time model.
