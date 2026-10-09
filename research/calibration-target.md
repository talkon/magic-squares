# Calibration search in the target region (October 2026)

**Question.** Have we run a quick search in the range of N where the first magic square is expected,
to calibrate the estimates?

**Short answer.** Not before this search; now we have. 321 pre-registered sums were searched in the
region where scheduler v2 puts its expected yield (N' 3-45k, 321 different P, 312 of them never
searched before). The search took 4.93 CPU-hours.

* **Squares per sum and the S and P traversal rates are now calibrated where most of the forecast's E
  sits.** At N' 3-6k the error is a few percent; at 6-12k it is about 10%.
  * Squares are 1.04 of prediction at 3-6k and 1.28 at 6-12k.
  * S traversals per square are 0.94-0.97, P traversals 1.02-1.08.
  * The shipped x0.8 correction of squares at 6-12k mostly goes away. The refit puts it at x0.93.
* **The SP rate is still not calibrated, as the pre-registration said it would not be.** There were 7
  (square, SP diagonal) pairs against 5.46 expected. The coupling beyond the S and P rates is f_rho =
  1.06 (68% 0.84-1.34). A factor-2 deficit is excluded (p = 0.02); a factor-2 excess is not
  (p = 0.15).
* **Forecast (main update):**
  * 1 CPU-year: E = **0.126**, 90% band 0.052-0.31.
  * 10 CPU-years: E = **0.28**, 90% band 0.12-0.73.
  * The shipped forecast was 0.12 and 0.25, quoted as about 0.1 and 0.2 after a x0.8 selection
    discount. The data show no selection effect within a band, so the discount is dropped and the
    figures to quote are about **0.13 and 0.28**.
  * Charged at the measured CPU instead of the anchored time truth, E is about 5% higher: 0.13 and
    0.30.
* **The band is narrower but still wide.**
  * The class-factor part shrank from ln sd 0.56-0.59 to 0.10-0.12.
  * What is left is the SP coupling: f_rho^2 has ln sd 0.47.
  * The shipped band was x8 between its 5% and 95% points. It carried only ln sd 0.2 for the SP term;
    with an honest SP term it would have been about x16. The updated band is x6.
* **The existence study's level factor F** goes from x/÷2.0 to x/÷1.8 (68%), compared component by
  component. This search does not reduce the extra error at large X (sums beyond N = 12k).

The first write-up of this analysis (not committed) quoted 0.159 / 0.361 with a x3.9 band and F
x/÷1.5. Two reviews found that it had counted the SP evidence twice and that several intervals were
too narrow. This document is the corrected version; "Changes after review" lists the fixes.

## 1. Design (pre-registered)

The predictions were frozen before the first run.

* `PREREGISTERED.txt`, registered 2026-10-09T04:16:57Z, while runs/ was still empty. The first run
  started at 04:17:52Z.
* `predictions.json` sha256 `c1e673b595d7de9452a5c2c01113967c4c6579220c885a3d25e15ed270b06c08`.
* `plan.json` sha256 `4be7288d…`, `run.sh` sha256 `60565a8c…`, `runner.py` sha256 `c8145b30…`.
* `bin/msearch` sha256 `4bf407c2…`, engine 3, built at e0ace54.

`analyze.py` and `prepare.sh` check all of these hashes before using the files.

**Model.** Scheduler v2 at its shipped calibration on a fresh state: `Calibration()` priors, the
engine-3 time laws, ratio level -0.028. The prediction formulas are in `archive/work/cm.py` and in
`predictions.json` under "model".

**Where E sits.** A greedy simulation of v2's plan (`forecast --shipped --truth anchored`, 10% of
candidates) put the E of 1 / 10 CPU-years at:

| N' band | 1 CPU-yr | 10 CPU-yr |
|---|---:|---:|
| < 3k | 34% | 18% |
| 3-6k | 61% | 54% |
| 6-12k | 5% | 26% |
| 12-24k | 0% | 2% |

The existence study's N_magic puts 13 / 21 / 39 / 28% of its N' 3-45k part in 3-6k / 6-12k / 12-24k /
24-45k.

**Strata.** Four N' bands × 3 assignment classes (sorted, ratio ≤ 1.1, ratio > 1.1) × 3 k bins (k ≤ 5,
6, ≥ 7). At 3-12k, sums were drawn from the 10 CPU-year plan in proportion to predicted magic; at
12-45k, from all pool sums with S ≤ 2 S0, also in proportion to predicted magic. There is one sum per
P, and (P, S) pairs searched before were excluded.

| band | sums | plain | d-first | d stride k_d | stream stride k_cal | predicted CPU-h (anchored) |
|---|---:|---:|---:|---|---|---:|
| 3-6k | 289 | 289 | 0 | – | – | 3.58 |
| 6-12k | 25 | 18 | 7 (whole d loop) | 1 | 5-20 | 2.44 |
| 12-24k | 4 | 0 | 4 | 5, 11, 16, 24 | 5-37 | 0.80 |
| 24-45k | 3 | 0 | 3 | 354, 592, 876 | 236-672 | 0.60 |

(An earlier text summary of the design said "17 plain, 8 d-first (k_cal 4-20)" and listed other
strides. The runs were made as `plan.json` and `work/design.out` give them; the table above follows
those files.)

**Pre-registered expectations (`plan.json` "verdict").**

* About 5.5 SP pairs were expected. Even at their expectation, that gives an interval of [0.36,
  1.93] on obs/pred. The verdict was "sp_count_useful: false": about 45 pairs, roughly 39 CPU-hours at
  3-6k, would be needed for ±25%.
* The SP+S and SP+P rungs were expected to have about 0.07 and 0.02 events, i.e. no information.
* The analysis was to lean on the squares, the S and P traversals, the sub-events and the time laws.
* The plan also said the sub-events would test the SP rate indirectly "with 10-50x the events". That
  premise turned out wrong at j ≥ 4 (section 4).

**Limits.** At most 2 msearch processes, all at `nice -n 19`. Each process had RLIMIT_CPU of 3 × its
prediction + 90 s, and the run had an 8 CPU-hour hard cap.

## 2. What ran

* All 321 runs finished with status "done", rc 0, and none was truncated.
* The ledger never shows more than 2 runs at once.
* **Search CPU** was 4.93 h (process CPU), against 7.41 h predicted under the anchored truth.

  | band | CPU-h | obs/pred |
  |---|---:|---:|
  | 3-6k | 2.34 | 0.65 |
  | 6-12k | 1.97 | 0.81 |
  | 12-24k | 0.41 | 0.51 |
  | 24-45k | 0.21 | 0.36 |
  | all | 4.93 | 0.665 |

* **The analysis CPU** was outside the search budget, about 1.8 CPU-h in all, all at nice 19:
  * the first pass, mostly calibrate.py's heuristic on 11,121 squares (2 processes, 40 min);
  * this revision, about 0.2 h: two forecast simulations run one after the other, and calibrate
    reruns from the heuristic cache.

  In the first pass, a forecast simulation ran alongside the 2-process heuristic, which made 3 busy
  processes on the shared machine for about 10 minutes. This revision ran at most 2.
* **Records.** There were 10,677 plain squares and 444 calibration-stream csquares (weight =
  stride; total weight 18,470), with no duplicate hashes. All grids re-verify: the recomputed counts
  match the records exactly. There were 0 dsquare records: the 7 full d loops at 6-12k found no SP
  pair (0.71 expected).

## 3. Observed vs predicted

**Conventions.**

* O/E pools over sums. **HH** is the design-weighted (Hansen-Hurwitz) mean of per-sum ratios. It is the
  estimator that matches the draw design; it agrees with the pooled ratio within 0-4%.
* Totals that include stream estimates (est_squares, k × traversals) are weighted. Their **Poisson**
  interval uses the effective count O/s against E/s, with s = Σ w_i E_i / Σ E_i, so a stream's events
  carry the uncertainty of one event each.
* **boot** resamples sums within bands.
* **od** adds the measured between-sum sd 0.41 over the E-weighted effective number of sums.
* All intervals are 90%.

### 3.1 Squares per sum

| band | sums | O / E | HH | effective count | Poisson | boot (sums) | od (between 0.41) |
|---|---:|---|---:|---:|---|---|---|
| 3-6k | 289 | 8,739 / 8,443 = **1.035** | 1.03 | 8,739 | [1.02, 1.05] | [0.98, 1.09] | [0.98, 1.10] |
| 6-12k | 25 | 5,229 / 4,100 = **1.28** | 1.29 | 679 | [1.20, 1.36] | [1.16, 1.39] | [1.04, 1.56] |
| 6-12k plain | 18 | 1,938 / 1,501 = 1.29 | 1.31 | 1,938 | [1.24, 1.34] | [1.10, 1.47] | [1.08, 1.54] |
| 6-12k d-first streams | 7 | 3,291 / 2,599 = 1.27 | 1.26 | 284 | [1.14, 1.39] | – | [0.94, 1.71] |
| 12-24k | 4 | 2,558 / 2,087 = 1.23 | 1.45 | 77 | [1.00, 1.47] | – | [0.78, 1.92] |
| 24-45k | 3 | 1,944 / 3,220 = 0.60 | 0.41 | 4.2 | [0.20, 1.33] | – | [0.24, 1.51] |
| 12-45k | 7 | 4,502 / 5,307 = 0.85 | 0.96 | 15 | [0.52, 1.29] | – | [0.50, 1.43] |

* **6-12k.** Plain sums and d-first streams agree (1.29 and 1.27). These ratios are against
  predictions that include the shipped x0.80 at 6-12k; against the uncorrected model this band is at
  about 1.03.
* **Comparison with the earlier held-out live units.** Those units gave 0.78 against the uncorrected
  model at N' ≥ 6k, about 2 sd from this search's 1.03. The GLM refit with the measured dispersion
  (A_SQ 6, section 5) combines the prior from those units with these sums: the 6-12k effect is
  -0.08 ± 0.09 (x0.93) in place of the shipped -0.22 ± 0.15 (x0.80). With A_SQ 15 it is
  -0.05 ± 0.06, and with A_SQ 4.7 -0.09 ± 0.09. The x0.8 shrinks to about x0.93; it is not refuted.
* **12-45k.** These sums are effectively 15 squares, 4.2 of them at 24-45k: R004 sampled 1 csquare,
  R100 4 and R204 0 (stride 236-672). The SQ12 = 1.23 in the predictions would become about 1.04,
  with a 90% range of roughly [0.62, 1.76]. That is uninformative, so it is not pooled with the 9
  earlier sums.
* **Between-sum spread** (plain sums). The sd of ln(O/E) is 0.46 (0.40 after stratum means), and
  Pearson φ on raw counts is 5.4 at 3-6k. The scheduler assumes 0.26 (A_SQ 15). With a per-sum 90%
  interval built from 0.26 and φ 2.5, 15-18% of sums fall outside it.

### 3.2 S and P traversals per square

| band | S / square | Poisson (eff.) | boot | P / square | Poisson (eff.) | boot |
|---|---|---|---|---|---|---|
| 3-6k | 5,408 / 5,759 = **0.939** | [0.918, 0.960] | [0.912, 0.968] | 1,544 / 1,517 = **1.02** | [0.975, 1.06] | [0.974, 1.06] |
| 6-12k | 0.973 (eff. 342) | [0.89, 1.07] | [0.82, 1.08] | **1.08** (eff. 101) | [0.91, 1.28] | [0.89, 1.31] |
| 3-12k | **0.950** (eff. 2,443) | [0.919, 0.982] | [0.90, 1.00] | **1.04** (eff. 792) | [0.98, 1.10] | [0.97, 1.11] |
| 12-24k | 0.64 (eff. 18) | [0.41, 0.95] | – | 1.21 (eff. 6) | [0.50, 2.25] | – |
| 24-45k | 2 raw traversals | – | – | 0 raw | – | – |

* **The S deficit is real.** At 3-6k it is 0.939, with Poisson p(≤) = 1.5e-6. The traversal counts are
  overdispersed (φ 1.41), and adjusted for that p = 5e-5. Over 3-12k, p = 0.016.
* **The effect on P(magic) without SP**, (f_S f_P)², is 0.91 at 3-6k, 1.10 at 6-12k and 0.98 over
  3-12k.

### 3.3 Strata and subgroups

Per stratum (band × class × k), the squares O/E ranges from 0.78 to 1.76. Strata with more than
about 10 sums sit within 0.93-1.14 at 3-6k (`archive/results/analyze.out`).

**Subgroups at 3-12k, against the overall 3-12k rate.** There are 48 contrasts: squares, S and P,
over the splits class, k, S/S_min, labels and actual N. Each is a bootstrap ratio over sums with a
two-sided p.

* 3 have p < 0.05; 2.4 are expected by chance. None survives Holm.
* The 3 are all in the actual-N split, restating the 6-12k excess:

  | actual N | squares | relative to 3-12k | p |
  |---|---:|---:|---:|
  | < 4k | 0.94 | 0.85 | 0.003 |
  | 4-6k | 1.02 | 0.92 | – |
  | 6-9k | 1.36 | 1.22 | 0.0015 |

* Two subgroups have low S rates that are **not** significant against the overall 0.95:
  * S/S_min < 1.1: S 0.886, relative 0.93, p 0.10;
  * k ≥ 7: S 0.89, relative 0.94, p 0.23.
* Squares by class: sorted 1.11, ratio ≤ 1.1 1.13, ratio > 1.1 0.98 (relative 0.88, p 0.31).

### 3.4 SP pairs (the key number)

| | O | E | O/E | Poisson | boot (sums) | pre-registered overdispersed |
|---|---:|---:|---:|---|---|---|
| pairs found, per sum | 7 | 5.46 | **1.28** | [0.60, 2.41] | [0.54, 2.15] | [0.63, 2.62] |
| per square (SP rate), 3-12k | 7 | 5.92 | **1.18** | [0.56, 2.22] | [0.50, 1.99] | – |
| per square, plain only | 7 | 5.02 | 1.40 | [0.65, 2.62] | [0.59, 2.31] | – |

* **Where the pairs came from.**
  * All 7 are sp_count 1, best_score 7 squares in plain runs: R017, R154, R224, R227 and R270 at
    3-6k, R111 and R317 at 6-12k.
  * By band, 3-6k had 5 against 4.23 and 6-12k had 2 against 1.21.
  * The d-first loops found 0 against 0.73, and the streams 0 against 0.09.
* **What the count can exclude.**
  * p(≥ 7 | 5.46) = 0.31.
  * A factor-2 deficit is excluded: p(≥ 7 | 2.73) = 0.022.
  * A factor-2 excess is not: p(≤ 7 | 10.9) = 0.15.
* **SP coupling.** f_rho = f_SP / (f_S f_P) is the SP rate beyond what the refit S and P rates imply.
  Here E = 5.92 × 0.950 × 1.037 = 5.86, so the MLE is 1.19.

  | prior on ln f_rho | f_rho median (68%) | f_rho² median / mean | ln sd of f_rho² |
  |---|---|---|---:|
  | **N(0, 0.29²): κ calibration 0.145 plus extrapolation 0.25 (main, count only)** | **1.06 (0.84-1.34)** | **1.12 / 1.24** | **0.47** |
  | N(0, 0.145²): κ only | 1.02 (0.89-1.17) | 1.04 / 1.08 | 0.27 |
  | flat | 1.14 (0.76-1.63) | 1.30 / 1.63 | 0.78 |

  The extrapolation sd of 0.25 (carrying κ from near S_min to N' 3-12k) is an assumption. The κ-only
  and flat rows are sensitivity analyses.

### 3.5 Rungs above SP

The plan counted these from the grids: the partners of each SP diagonal, as unordered pairs.

| rung | observed | expected |
|---|---:|---:|
| SP+S | 0 | 0.073 |
| SP+P | 0 | 0.019 |
| magic | 0 | 2.9e-5 |

The best_score distribution over all squares is 0: 5,239; 2: 4,118; 3: 1,689; 4: 42; 5: 21; 6: 5;
7: 7. As pre-registered, these rungs carry no information.

### 3.6 Time per sum (process CPU, Σ obs / Σ pred)

| sums | plain law | d-first law | ratio law (d-first / plain vs r) | anchored truth |
|---|---|---|---|---|
| plain 3-6k (289) | **0.97** (gm 0.87, sd 0.29) | – | – | 0.65 |
| plain by N' | 3-4k 0.83, 4-5k 0.95, 5-6k 1.01, 6-8k 1.19 | – | – | 6-8k 0.75 |
| d-first 6-12k (7, full d loops) | streams' est. 0.99 | **0.85** (sd 0.18) | **0.94** (sd 0.06) | – |
| d-first 12-24k (4) | 0.62 | 0.56 | 0.92 | – |
| d-first 24-45k (3) | 0.16 | 0.38 | 1.03 | – |
| all 14 d-first sums | – | – | gm **0.95**, sd 0.11 | – |

* The plain law holds at 3-6k, while the anchored truth charges these sums about 1.5x too much.
* The ratio law holds: a0 moves from -0.028 to about -0.08, so the d-first switch moves from N' 4.3k
  to about 3.9k.
* Both laws over-predict above 12k, but that rests on 7 sums.
* N / N' has a median of 1.00 at 3-6k and 1.02-1.05 above.

## 4. Sub-events and the per-square heuristic

calibrate.py was run on the 10,677 plain squares at 3-12k, unweighted
(`archive/results/calibrate_unw/`). The weighted run on all 11,121 squares is in `calibrate/`. Its rungs
are dominated by a few heavy 12-45k csquares (the 1% heaviest squares carry 27% of the weight), so it
is not used for the figures below.

**Rungs** (O/E, bootstrap 90%):

| predictor | S | P | SP | S+S | S+P | P+P |
|---|---|---|---|---|---|---|
| v2 (pre-registered rates) | 0.937 [0.92, 0.96] | 1.026 [0.99, 1.06] | 1.40 [0.60, 2.39] | 0.96 | 0.88 | 1.74 (6 / 3.45) |
| heuristic | 0.97 [0.95, 0.99] | 1.18 [1.14, 1.23] | 1.74 [0.74, 2.99] | 1.00 | 1.05 | 1.62 |
| heuristic, conditional | – | – | 1.56 | 1.05 | 0.96 | 1.25 |

**Sum-product coupling.** This is "X + j exponents pinned" over "j exponents alone". The SP
traversals satisfy every window of the square's exponents, so at j ≥ 4 the sub-events are largely the
SP count itself:

* at j = 5, 32 of the 35 events are the 7 SP traversals;
* at j = 4, 40 of the 109 are.

They are not independent evidence on the SP rate, contrary to the "10-50x the events" premise of the
plan. `subevents_nonsp.py` removes them from both sides:

| j | coupling | without the SP traversals |
|---:|---:|---:|
| 1 | 0.98 | 0.98 |
| 2 | 0.99 | 0.98 |
| 3 | 1.06 | 1.02 |
| 4 | 1.37 (109 events) | **1.25** (69 events; about [0.92, 1.70]) |
| 5 | 2.20 (35) | not estimable (3 events; the SP part is about 95% of E) |

The j = 4 excess without SP points the same way as the count, but it is not significant on its own.
Going from j = 4 to j = k would need a model of how the coupling grows with j, so it is not combined
with the count.

**Pair residual and E[magic].**

* The pair residual is r(d) = 0.99, 0.96, 0.93, 0.88, 0.68 for d = 1-5, as in calibration.md, with
  β = 0.14.
* E[magic] over these squares:
  * v2: 2.7e-5;
  * heuristic: 3.7e-5;
  * heuristic × r(k): 2.4e-5;
  * heuristic with congruences only: 2.6e-5.

## 5. Forecast update

**Method.** The class-factor GLMs (squares, S, P) are refit on the 321 sums with `Calibration.fit`
(`forecast_update.py`). The fit uses the shipped priors and the scheduler's stream weights, with A_SQ
= 6 from the measured between-sum sd (0.40 after stratum means). The scheduler ships A_SQ = 15
(0.26); with 15, E changes by less than 1%.

The forecast is then re-simulated as `forecast --shipped --truth anchored` (10% of candidates, seed 1,
budget scaled). The point is E(updated) × the f_rho² median (1.12). The band is built from 4,000
draws of:

* the class-factor posterior;
* × lognormal(0, 0.2) for the pair factor;
* × the f_rho² posterior (ln sd 0.47);
* × lognormal(0, 0.12) for selection.

The plan is not re-optimised per draw.

**Selection.** The x0.8 discount is replaced by the lognormal(0, 0.12) term. Within 3-6k, the
quartiles by predicted magic per CPU-second are measured against the stratum-adjusted model, on
squares × (S P)². The top quartile holds 54% of the band's E.

| quartile | ratio |
|---|---|
| Q1 (top) | 0.96 [0.82, 1.11] |
| Q2 | 1.05 |
| Q3 | 0.96 |
| Q4 | 1.09 |

The slope of ln(O/E) on ln density is +0.020 ± 0.027 (244 sums). There is no detectable winner's
curse within the band. The 0.7-0.9 seen on earlier held-out sets is consistent with the 6-12k
squares correction, which this refit now carries.

| E | 1 CPU-yr | 10 CPU-yr |
|---|---|---|
| shipped (pre-registered) | 0.120 [0.047-0.37] (x7.8) | 0.254 [0.096-0.84] (x8.7) |
| shipped, with an honest SP term (f_rho² prior, ln sd 0.58) | 0.120 [0.034-0.53] (x15.5) | 0.254 [0.072-1.15] (x16) |
| updated class factors only | 0.113 [0.071-0.19] (x2.7) | 0.254 [0.16-0.42] (x2.7) |
| **main: + f_rho² (count only) + selection term** | **0.126 [0.052-0.31]** (x6.0) | **0.284 [0.12-0.73]** (x6.0) |
| main, charged at the measured CPU per N' band | 0.132 | 0.30 |
| f_rho² variants: κ-only prior / flat prior | 0.117 / 0.146 | 0.265 / 0.33 |
| (rejected, double-counts SP: count + "heuristic route") | (0.159) | (0.359) |

* The bands are 90% (5-95%).
* The previous quote was about 0.1 and 0.2: the shipped point × 0.8. The quote is now **about 0.13 at
  1 CPU-year and 0.28 at 10**.
* **Class-factor part of the band:** ln sd 0.56 → 0.115 at 1 CPU-year and 0.59 → 0.097 at 10.

**Where E moved** (updated / shipped, before f_rho):

| N' | 1 CPU-yr | 10 CPU-yr | basis |
|---|---|---|---|
| 1-3k | x0.85 | x0.86 | main effects learnt at 3-12k; not measured here |
| 3-6k | x0.91 | x0.91 | |
| 6-12k | x1.8 | x1.28 | |
| 12-24k | – | x1.2 | 2% of E |

**Level factor F of the existence estimate** (68% ln sd; existence.md §4.1). The comparison is
component by component on existence.md's basis, with no extrapolation term for κ:

| component | before | after | basis of "after" |
|---|---:|---:|---|
| squares level (K calibration, x0.85-1.25) | 0.12 | 0.04 | 3-6k od interval |
| S/P rates (diagonal review factors) | 0.25 | 0.10 | 3-12k bootstraps of S and P, squared |
| SP coupling (base κ 0.46-1.2) | 0.29 | 0.27 | 7 pairs, κ-only prior |
| pair factor | 0.11 | 0.11 | no magic-rung events |
| model form | 0.30 | 0.30 | the magic-rung form is not tested by squares |
| sampling and coverage | 0.11 | 0.11 | unchanged |
| in quadrature | 0.52 | 0.45 | |
| σ (existence.md rounds 0.55 up to 0.7; the same factor applied) | 0.70 (x/÷2.0) | **0.57 (x/÷1.8)** | |

* With the SP extrapolation term on both sides, the two totals are 0.73 and 0.59.
* existence.md's wider σ at X ≥ 2×10^4 (0.75 down, 1.0 up) comes from extrapolation beyond the measured
  N. This search does not reduce it: it has 7 sums at 12-45k with an effective 15 squares.

## 6. What remains uncalibrated

1. **SP coupling at the target.** It is x/÷1.6 on P(magic) (f_rho² ln sd 0.47) and is now the largest
   part of the band. About 45 SP pairs would give ±25% on f_rho, roughly 39 CPU-hours of plain search
   at 3-6k. The cheaper j = 4 sub-event route needs a model of the coupling's growth with j first.
2. **The magic rung itself.** The pair factor is v2's 1.07, against the heuristic's pair residual. No
   SP+S, SP+P or SP+SP event has been seen.
3. **N' ≥ 12k.** There are 7 sums with an effective 15 squares. Squares are 0.85 [0.50, 1.43] there, the
   traversals are a handful, and SQ12 stays at 1.23 (1.04 here, uninformative). This band holds 2% of
   E at 10 CPU-years but 28-46% of the existence total.
4. **N' < 3k.** 34% of E at 1 CPU-year sits here. It moved only through the GLM's main effects
   (k, ratio, x) learnt at 3-12k.
5. **Time laws at ≥ 12k.** They measured 0.16-0.6 of prediction on 7 sums.
6. **Between-sum dispersion of squares.** It is 0.41-0.46, against 0.26 in the shipped A_SQ. The
   forecast refit uses A_SQ 6, but the scheduler still ships 15.
7. **k ≥ 7, ratio > 1.1 and x < 0.1.** These are calibrated only at main-effect level (26, 42 and few
   sums).

## 7. Reproduce

Everything is in `research/calibration-target/`. `archive/` (2 MB) holds:

* `runs.tar.xz`: all 321 run outputs and the runner ledger;
* the frozen `predictions.json`, `plan.json`, `run.sh` and `runner.log` (xz);
* `PREREGISTERED.txt` and `runner.py`;
* `work/`: the design scripts and the design printout;
* `results/`: every analysis output (results.json, per_run.tsv, analyze.out, bands.json,
  glm_update.json, the forecast logs, subevents_nonsp, and both calibrate ladders with their
  results.json).

The msearch binary is not archived; it is built from e0ace54.

```sh
cd research/calibration-target
./prepare.sh                         # unpacks archive/ into work/, checks the pre-registered hashes
export CALIB_TARGET_DIR=$PWD/work
nice -n 19 python3 analyze.py > work/analysis/analyze.out          # 10 s; identical to archive/results
# the ladder and sub-events (heuristic ~45 CPU-min, cached afterwards; S_min from bin/enumerate
# when work/state has no pinfo cache):
mkdir -p work/plain312 && for r in $(awk -F'\t' '$5=="plain"{print $1}' archive/results/per_run.tsv); do
  ln -sf ../runs/$r.jsonl work/plain312/; done
nice -n 19 python3 ../../scripts/calibrate.py --state work/state --heuristic-cond --unweighted \
  --plugin v2=$PWD/v2plugin.py --jobs 2 --out work/calibrate_unw plain312=work/plain312
nice -n 19 python3 subevents_nonsp.py
# forecasts: need a filled scheduler state (profile store) in work/state (~1 CPU-h to profile, or copy one)
for v in shipped updated time; do nice -n 19 python3 forecast_update.py $v; done   # ~4 CPU-min each
nice -n 19 python3 forecast_update.py glm
nice -n 19 python3 bands.py
```

The scripts are:

* `analyze.py`: observed vs predicted per run and stratum;
* `forecast_update.py`: GLM refit and re-simulated forecast;
* `bands.py`: f_rho posterior, bands, F budget;
* `subevents_nonsp.py`;
* `v2plugin.py`: the pre-registered rates as a calibrate.py predictor;
* `cm.py`: the prediction formulas;
* `ctpaths.py`: paths.

Run them one at a time.

## 8. Changes after review

The first analysis was reviewed twice before this write-up. The fixes:

* **SP evidence counted twice (high).** The first posterior multiplied the 7-pair count by a "heuristic
  route" likelihood, a hard-coded N(ln 1.35, 0.25) taken from the j = 4-5 sub-events, which contain
  the same SP traversals. Its effect was f_rho² 1.41 and E 0.159 / 0.361 with a x3.9 band. The main
  result is now the count alone (pre-specified), and the sub-event route is recomputed without the SP
  traversals as a cross-check only.
* **Effective counts (medium).** Intervals on totals mixing plain and stream data used Σ(est²/se²), and
  the S/P intervals used raw counts. They now use the event scale s = Σ w E / Σ E, and φ is computed on
  raw counts. The 6-12k, 12-45k and S/P-at-6-12k intervals widened 2-4x. The "all" pooled ratio is no
  longer quoted, since it mixes differently defined populations.
* **F (medium).** It is now compared component by component on existence.md's basis, with model form
  kept at 0.30. The result is x/÷2.0 → x/÷1.8, not x/÷1.5.
* **d-first records (medium; no effect here, 0 dsquares).** calibrate.py and analyze.py now count
  magic squares by unique hash. The estimate with a d stride is Σ k_d / 2. SP+S, SP+P and SP+SP are
  judged per record from its own SP diagonal (`calibrate.d_rooted`). The d_stride falls back to the
  dchunk records when a run was killed before its dsum.
* **calibrate.py weighting.** The pair-correction totals are now weighted. Bootstraps resample whole
  sums when csquares carry weights, and the ladder states how far a few heavy csquares dominate it.
* **A_SQ.** The GLM refit uses A_SQ 6 (measured), with 15 and 4.7 as sensitivity.
* **Selection discount.** It is measured within band against the stratum-adjusted model, and carried
  as a band term instead of a point x0.8 (or x0.9).
* **Multiple comparisons.** Subgroups are reported as contrasts against the overall rate, with the
  number of tests and a Holm count. The S deficit carries a dispersion-adjusted p.
* **Rungs above SP.** They are counted from the grids as pre-registered.
* **Wording.** "Fresh P" is now "fresh (P, S); 312 of 321 P never searched", and the design summary
  is corrected (section 1).
* **Reproducibility.** The data are archived in the repository, and the paths are set by
  `CALIB_TARGET_DIR`.
