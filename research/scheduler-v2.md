# Scheduler v2: the analytic model as a scheduler (October 2026)

`scripts/scheduler.py` (default `--model analytic`) chooses the next
(P, S) to search by expected magic squares per CPU-second, using the
existence study's analytic model ([existence.md](existence.md) sections 2-6,
[existence/analytic.md](existence/analytic.md)) instead of v1's fitted
Poisson regression. This note records what it does, what it was checked
against, and how much to trust its forecast. The previous scheduler is still
available as `--model regression --pool classic` (and is the default for
`--vec-size` other than 6). Since this revision v2 also searches sums
diagonal-first where that is cheaper ("d-first units" below).

## Short answer

* **d-first units** (this revision; "d-first units" below). v2 picks plain
  or `msearch --diag-first` per sum from two time laws, splits d-first
  sums of hours into units of d, merges and resumes them, and feeds a
  calibration stream's squares to its models. Under its own laws d-first
  adds x1.00 / 1.02 / 1.07 / 1.14 to E at 1 / 10 / 100 / 1,000 CPU-years,
  and the plain-only E of 1,000 CPU-years needs 1.8x less CPU. The gain
  grows with the budget, as in retrospective.md section 7, but is smaller
  there. The numbers in the bullet below are plain-only, from before this
  revision.
* Forecast for the current build (`forecast --shipped`, a fresh state, model
  at its shipped calibration):

  | CPU time | 100 h | 1 CPU-yr | 10 CPU-yr | 100 CPU-yr |
  |---|---:|---:|---:|---:|
  | E, point | 0.016 | 0.119 | 0.25 | 0.39 |
  | 90% band from 50 draws | 0.005-0.053 | 0.038-0.32 | 0.075-0.62 | 0.10-1.0 |
  | E x 0.8 (selection discount) | 0.013 | 0.095 | 0.20 | 0.31 |
  | v1's choices, scored by the same model | 0.0087 | 0.059 | 0.096 | (0.18 *) |

  10 and 100 CPU-years are simulated on 10% and 1% of the candidates with
  the budget scaled (on 1-3% samples the earlier version came out 5-7% lower
  than on 10%). (*) v1 at 100 CPU-years is the earlier time law's number.
  On all existing data the learned calibration gives 0.122 at 1 CPU-year
  (band 0.09-0.22), which matches the shipped prior.
* **Update (calibration search in the target region,
  [calibration-target.md](calibration-target.md)).** The GLMs are refit on
  321 pre-registered sums at N' 3-45k (A_SQ 6) and the SP coupling is taken
  from 7 SP pairs. The forecast becomes **0.126 at 1 CPU-year (90% band
  0.052-0.31) and 0.28 at 10 (0.12-0.73)**. No selection discount applies:
  there is no winner's curse within band, so it is a band term of ln sd
  0.12. **Quote about 0.13 and 0.28.** The class-factor part of the band
  shrank from ln sd 0.56-0.59 to 0.10-0.12. Most of what is left is the
  SP coupling (f_rho^2 ln sd 0.47), and the 90% band is x6 wide. The
  shipped band's x8 had an SP term of only 0.2; with an honest one it would
  have been x16. The text below is from before this update.
* The quoted figure was the **discounted** one, about **0.1 in 1 CPU-year and
  0.2 in 10**. Units that v2 itself ranks highest came in at 0.7-0.9 of their
  predicted magic density on held-out data. This is a winner's curse, mostly
  through squares and P traversals (below).
* Gain over v1 under the same model: x1.9 at 100 h, x2.0 at 1 CPU-year,
  x2.6 at 10. v1 forecast 0.015-0.018 for itself (research/forecast.md).
  That "10x less" came mostly from v1's own pessimistic model (squares per
  sum flat beyond N ~ 5.5k), not from its choices.
* The existence study's "ideal" frontier (0.16 / 0.34 / 0.64) is **not**
  comparable to these numbers. It uses another time law, about 2.4x slower
  than v2's on the sums v2 picks, and assumes a x2.7 multiplier for
  non-sorted exponent assignments. The enumerated pool supports about x1.4
  of that x1.8 implied gain at a fixed budget. Re-planned under the
  frontier's time law, v2 reached 62-66% of the ideal (verification of the
  first version, before the time-law fix). Under its own law v2 is within
  1% of an exhaustive greedy over the same pool.

## The model per (P, S)

    magic per CPU-s = E_sq(P,S) g_sq(c) f_sq(P) x min(P_m(P,S) m(c) F_m(P), 1e-6) / t(P,S)

* `E_sq`, `P_m` (= 5400 kappa p_pair), the traversal probabilities, N and
  the raw label count come from `scripts/amodel.py`, a port of the
  calibrated model. It is evaluated once per P on 24 sums
  S0 (1 + u), u = 0.01..1, kept in a float16 store (`profile`, about 0.9
  CPU-hours for the whole pool), and interpolated in S. Grid points the
  model cannot evaluate (not converged, N < 12, at or below the last fall
  of N) are invalid, and no squares are predicted there.
* `c` is the calibration cell: 6 N' bands x 4 classes of the number of
  primes k x 3 assignment-ratio bins x 3 bins of x = ln(S / S_min), 216
  cells. `g_sq(c)` and `m(c)` = PAIR x (r_S r_P / base)^2 come from three
  main-effects GLMs (squares, S traversals, P traversals). Their priors
  carry the review corrections of existence.md 4.1, revised here on
  held-out data:
  * squares: x1.23 at N' >= 12k (SQ12, in E_sq); x0.80 at N' 6-12k;
    x0.78 below 1k; x0.90 at x < 0.1;
  * P(magic): x1.22 at N' >= 3k (the review's x1.35 lowered: S traversals
    at N' >= 3k ran 0.91 of the prior on the live units, 0.94 on all
    data); x0.85 at k >= 7; x0.92 at x < 0.1; the pair factor x0.87.
* `f_sq`, `F_m` are per-P gamma factors (empirical Bayes). For squares:
  (15 + o / 2.5) / (15 + e / 2.5). The shape 15 is the measured between-P
  sd of 0.26 in ln(obs/pred) per unit, and 2.5 is the overdispersion of
  squares between the sums of one P (per-sum Pearson chi2/df). For
  traversals the shape is 30. The observed counts are taken over the same
  sums as the expectations: inside the model's grid, not truncated.
* The GLMs are quasi-Poisson. Each cell's dispersion is
  2.5 + (1/15) sum_P e_P^2 / sum_P e_P for squares (traversals: 1 +
  (1/30) ...). Counts of one P move together, so one or two P can no
  longer pin a class factor. The `report` intervals use the same
  dispersion. The verification measured a between-unit chi2/df of 9.4;
  the formula gives about 9 for P with ~100 expected squares.

## The pool

`scripts/pool.py` enumerates every exponent assignment over 2..29 with
S0 = 6 P^(1/6) <= 6000, tau >= 1000, 4-10 primes, exponents <= 31 and
(P / P_sorted)^(1/6) <= 1.2: 377,908 P (25,498 sorted). The generator is an
exact branch-and-bound; it is checked against brute force in the tests and
cached in `pool_6.npz`, keyed by its parameters and POOL_VERSION. It is for
n = 6 only, and `--pool wide` refuses other sizes. Widening the pool further
adds little at 1-100 CPU-years: ratio 1.2-1.35 adds +0.7% / 1.7% / 3%, and
S0 6000-9000 adds +2.2% / 3.9% / 6%.

## The time model

    ln t = th . [1, ln(N'/4000), ln(L_raw/150), max(0, ln(N'/8000)),
                 [labels_obs > 128], band offsets (N' <1k, 1-3k, 3-6k, >=6k)]

* **The label-word step is new.** msearch's search is compiled per number
  of 64-bit words of the label bitsets. Above 128 labels the third word
  raises the cost per node from about 0.25 to 0.6 us. Without the step, the
  first law under-predicted the N' 3-6k sums, which carry 61% of the
  1-CPU-year E, by 1.2-1.3x (verification: 1.20 [1.05, 1.35] on 32 new
  CPU-timed sums, 1.31 on the older W = 3 sums).
* **Prior** (AMODEL_VERSION 2). th = (3.859, 5.458, -3.462, 0.946, 0.393,
  0, 0, 0, 0), sd 0.30. It is fitted by least squares on ln t with the
  model's N' and L_raw at each sum, as the scheduler uses them. The data
  are 3,690 full engine-2 sums at N' >= 300 (validation rows, the live
  verification units A/B/C, the 32 held-out draws) and 10 r1-sampled sums at
  N 6.7-45k (weight 30). Held out, the same law fitted without the draws and
  the live B/C units gives sum(t)/sum(pred) 1.10 on the draws, 1.03 on live
  C and 1.04 on live B. The previous law gave 1.20 and 1.10 there. As
  shipped, obs/pred is 1.04 at N' 3-6k and 1.06 at 6-12k, but 0.87-0.90
  below N' 3k: small sums are slightly over-predicted, which costs them a
  little priority. The residual sd is 0.21-0.33 by source and 0.48 on the
  sampled sums.
* **Engine 3** (cx/integrated: per-r1 widths, carried bitsets up to 512
  labels, the pretest; msearch writes `"engine":3`). Its plain search costs
  0.98-1.01x engine 2's up to 128 labels, 0.66-0.89x at 129-256 labels
  (geometric mean 0.80) and 0.26x above 256 labels (paired, 13 sums at
  N 2-32k; ideas.md, "Measurements on the integrated binary"). Its law
  starts from engine 2's posterior with the label-word step 0.227 lower
  (`amodel.ENGINE_TIME_SHIFT`, applied by `amodel.time_prior`), and learns
  its level online as above. The prior's fixed shape (W3 step, N slope,
  hinge) no longer matches this build above 128 labels or 11k vectors:
  against the engine-2 prior the integrated plain search came out at
  obs/pred 0.24-0.72 at N >= 11k. Refit the shape for engine 3 when large-N
  data arrives; the 13 measured sums can serve as anchors.
* **The online refit learns the level, not the shape.** The first version
  refit all coefficients by ridge regression on every sum. Thousands of
  cheap sums at N' < 3k then moved the label slope from -2.2 to -3.1..-3.5,
  and the N' 3-12k sums came out 1.4-1.9x cheaper than they are (three
  independent checks). Now the N, L and hinge slopes are held at the prior
  (precision 1e6). Learned are: a global intercept (sd 0.1), the label-word
  step (sd 0.1) and one offset per N' band (sd 0.25). A residual learned in
  one band moves the others by about 14% of itself, and a test asserts that
  3,000 sums at N' < 3k with a different shape move predictions at
  N' >= 3k by less than 10%. Refit on all existing data, the band offsets
  come out at -0.09 / -0.07 / +0.10 / +0.17 and the residuals are within
  +-0.08 in every band and word class.
* **CPU, not wall time.** msearch now writes `"cpu"` (process CPU time of
  the reduction, the setup, the search and the sum's share of the
  enumeration) next to the wall times, and the scheduler learns from it
  when present. The smoke run, made during other benchmarks, had read 1.28x
  slow at W = 2 on wall time. Sampled research outputs (records with a
  `sample`/`stride` field) are skipped. The r1-sampled files of the
  existence study carry no such field and must not be put in `units/`.

## Back-test: what was fitted on what

Most of the data that exists was used to fit something. The numbers the
first version reported as back-tests on these data (seed36 1.085, sched40
0.92, the time law at 1-3k) are **fit residuals**:

| constant | fitted on |
|---|---|
| C0, B, DN (squares) | n = 6 sums near S_min: seed36 + sched40 + verify-forecast (19,834 sums, 6,498 squares) + existence n6low |
| KAPPA, TRAV_BASE 0.86 / 0.64 | the squares of the same set |
| SQ12 | 9 r1-sampled sums at N >= 12k |
| NTR, K7, X01 | existence full and sampled squares |
| time prior, LAB0-2 | full engine-2 sums (validation rows, live A/B/C, the 32 draws) + 10 sampled sums |
| GLM prior revisions (N' < 1k, 6-12k, x < 0.1, S trav N' >= 3k) | all data + the live verification units |
| A_SQ, PHI_SUM | per-unit residuals of the live units |

Held-out evidence (not used for any constant when it was measured):

* **The first search's per-P totals** (legacy; 326 P searched up to 2 S0,
  574,001 squares): squares 1.05 [1.04, 1.06] (P bootstrap); 1.05 for the
  209 P whose median square is at N' 3-6k. S traversals 1.02, P traversals
  1.05, SP 0.73 (consistent with kappa_SP).
* **v2's own top units** (A11 smoke run, 46 units, 5,924 sums): squares
  0.92 [0.86, 0.98] (v1 on the same: 1.43). S traversals 0.99, P traversals
  0.84 [0.74, 0.94], (r_S r_P)^2 0.69 [0.53, 0.90].
* **Draws from v2's 1-CPU-year plan** (32 sums, chosen in proportion to
  predicted magic, run in full): squares 651/693 = 0.94 [0.75, 1.18],
  E-weighted 0.85 [0.67, 1.03] (v1 on the same: 0.56). Sorted assignments
  1.01; ratio <= 1.1 0.74 [0.57, 0.93] (12 P). P traversals 1.06.
* **Live units, predictions saved before each run**: stage A (v2's top 53
  units on a fresh state) squares 0.97 [0.92, 1.03]. Stages B + C (15
  E-weighted units, 9 with non-monotone exponents): squares 0.87
  [0.84, 0.91] (0.94 below N' 3k, 1.00 at 3-6k, 0.78 at >= 6k), S
  traversals 0.92, P traversals 1.03. The implied magic-density factor is
  0.72 pooled (unit bootstrap 0.53-1.01), or 0.85 design-weighted.
* **Time**: see above (held-out 1.03-1.10 with the new law).
* **Calibration search in the target region** (321 pre-registered sums,
  [calibration-target.md](calibration-target.md)). The sums were drawn in
  proportion to predicted magic from the 10 CPU-year plan at N' 3-12k and
  from the pool at 12-45k, and predictions.json was frozen before any run.
  Results:
  * **Squares.** 3-6k: 1.035 (boot over sums [0.98, 1.09], design-weighted
    1.03). 6-12k: 1.28 [1.16, 1.39], i.e. ~1.03 against the model without
    its x0.8 at 6-12k. 12-45k: 0.85 [0.50, 1.43], an effective 15 squares,
    uninformative on SQ12.
  * **Traversals per square.** S 0.939 at 3-6k (dispersion-adjusted p 5e-5)
    and 0.95 at 3-12k. P 1.02 at 3-6k and 1.04 at 3-12k.
  * **SP pairs.** 7 / 5.46 (1.28 [0.60, 2.41]). SP+S, SP+P and magic were
    all 0, against 0.07, 0.02 and 3e-5.
  * **Between-sum sd of squares obs/pred.** 0.41-0.46, against the 0.26
    behind A_SQ 15.
  * **Selection.** Within 3-6k the top density quartile (54% of E) is 0.96
    [0.82, 1.11] on squares x (S P)^2 against the stratum-adjusted model,
    and the slope is +0.02 +- 0.03.
  * **Time.** The plain law is 0.97 at 3-6k. The anchored truth over-charges
    these sums 1.5x. The ratio law is 0.95 (sd 0.11) over 14 d-first sums.
  * **Refit.** The GLM refit with A_SQ 6 puts the 6-12k squares effect at
    -0.08 +- 0.09 (x0.93) instead of -0.22. The earlier live units (0.78 at
    >= 6k) and these sums (1.03) differ by ~2 sd, so the correction is
    reduced, not dropped.

So the squares model holds on held-out data within about +-10% on average.
What v2 selects came in 5-15% low on squares and somewhat low on P(magic)
(0.7-1.05). The **selection discount of 0.8** that `forecast` prints is the
middle of that range. The calibration search above finds no such effect
within band once the 6-12k squares correction is refit. Its update
therefore drops the discount and carries selection as a band term instead.
`forecast` still prints the x0.8 (SELECTION_DISCOUNT in scheduler.py is
unchanged).

## Caveats and open points

* **Overdispersion.** Unit-level obs/pred has sd 0.26, so any statement
  from a handful of units is weak. Some earlier intervals (e.g. "fresh-P
  0.92 [0.87, 0.97]") were Poisson-only and are about 3x too narrow at the
  unit level.
* **Non-sorted assignments.** Held-out evidence is mixed: the draws give
  0.74 at ratio <= 1.1, the live units and the legacy data about 1.05, and
  the learned all-data effect is +0.05 (0.09). The prior stays at 0, and
  the GLM learns it under the quasi-Poisson likelihood. 50% of E at
  1 CPU-year is in sorted assignments, 40% at ratio <= 1.1 and 9% above.
  If the 0.74 held for all non-sorted P, E(1 CPU-year) would be about 15%
  lower.
* **Winner's curse.** The discount is applied to the quoted forecast but
  not to the ranking. Shrinking each (P, S)'s predicted P-traversal rate
  toward its class mean before ranking was considered and not done. The
  evidence for a per-cell effect is split (smoke 0.84, top decile of data
  v2 did not choose 0.80, draws 1.06, live B + C 1.03), and the planner's
  choices are robust to a uniform factor.
* **The forecast's point lies below the median of its draws.** The point
  uses the posterior-mean log factors per cell, and the draws vary
  independent cell effects (a sum of lognormals). `forecast` prints both.
  The band does not re-plan per draw.
* **Large N.** E at 1 CPU-year stays below N' 8k (max N' of a unit with
  predicted magic squares). Beyond 12k the time law rests on 10 sampled sums
  (sd 0.48) and the squares on SQ12 from 9 sums. These matter from about 10
  CPU-years on (8% of E at 100).
* **Not done.** Sharded `units/` and a compact summary store. The summary
  is JSON and is rewritten at every refit, which will be slow beyond ~10^5
  P (about 620k units per CPU-year). The per-P prior strengths have not
  been re-estimated beyond A_SQ. (The d-first mode, listed here before, is
  now in: "d-first units" below.)

## d-first units (this revision)

`msearch --diag-first` searches a sum as one semi-magic search per vector
d of the sum, on V_d = {v : |v & d| = 1}. It finds every magic square
(twice) but not the other semi-magic squares, and above N ~ 4-5k it is
cheaper than the plain search: 0.57x at N = 7.6k, 0.42-0.47x at 15-21k,
0.29-0.43x at 23-32k (ideas.md, "Measurements on the integrated binary").
v2 now launches it.

**Mode per sum** (`--dfirst auto`, the default; `off`; `on`). A sum is
searched d-first where the measured d-first / plain CPU ratio with the
calibration stream (below) is below 1, (1 + 0.07) r(N') < 1, and its N' is
at least `--dfirst-min-n` (default 2000). `on` takes d-first at every N' >=
`--dfirst-min-n` (for tests). A d-first unit runs `msearch --diag-first
--diag-first-min-n 0`, so msearch never overrules the scheduler's choice. A
unit holds one mode only, and the launched record (`launched_6.jsonl`) has
`mode`, `dlo`, `dhi`, `nd`, `calib` and `frac`. `emit` prints the same
arguments.

    ln r = a0 - 0.566 ln(N'/4000),   a0 = -0.028 (prior, sd 0.1)

(`amodel.DFIRST_RATIO_PRIOR`, pooled over the perf verifier's 13 sums with
both modes, resid sd 0.17; ideas.md "Measurements on the integrated
binary"). With the prior level `auto` switches at N' = 4.29k for every P.
The level a0 is learned online from pairs (below), the slope held.

* The first version compared the two time laws instead, t_d + t_cal < t_p.
  That put the switch at N' ~ 9-10k weighted by E, where the
  measurements put it at ~4.3k. The two laws were fitted separately and
  their errors do not cancel:
  * the plain law under-predicts the plain search at N 5-8k generally
    (obs/pred geometric mean 1.35 over 10 engine-3 sums, range 1.02-1.73,
    9 of 10 at >= 1.10), while the d-first law is at 0.94 there;
  * the quotient t_d / t_p inherits the plain law's label term L^-3.46,
    which the d-first law has not: a sum with labels 10% above typical
    got a predicted ratio ~1.4x higher. The pool's E-heavy sums at N'
    4-12k have model labels 7-14% above the measured sums'.
* Six pool sums at N' ~6.1k (13 8 4 2 2 1 / 2469, 13 6 4 3 2 1 / 2368,
  13 9 4 3 1 1 / 2750, 14 7 3 4 2 1 / 3315, 13 7 4 2 1 1 / 1492, 12 7 4 3 2
  / 1854; 154-168 labels), measured by the review as d-first at d-stride
  16 plus a stream at stride 12: d-first / plain 0.69-0.82, flat in labels,
  on the pooled law (r(6.1k) = 0.77). The laws' quotient predicted
  0.93-1.49 and chose plain on 5 of the 6 (regret 1.21-1.33x), and on 13
  7 4 3 1 1 / 1950 (measured 0.78, predicted 1.16).
* On all 22 engine-3 sums measured in both modes (the 13, 13 7 4 3 1 1 /
  1850 and 1950, 12 6 3 2 1 1 / 950 and the six) the ratio law's choice
  has regret 1.00. Between 4k and 12k the laws' quotient agreed with it on
  only 32% of the plan's E.
* The time laws now only price the sums (the planner's density, the
  units' sizes). If a label term is ever added to one law, the same goes
  into the other, so that pricing stays consistent.
* `--dfirst on --dfirst-min-n 5000` reproduces msearch's own threshold
  (and retrospective.md 7's assumption) instead.

**The d-first time law** (`amodel.DFIRST_TIME_PRIOR`, its own TimeModel
per engine, mode "dfirst"):

    ln t = 3.215 (se 0.15) + 3.576 (se 0.11) ln(N'/4000),   sd 0.25

* Fitted by least squares on ln t over the 11 measured sums at N 4.1-31.7k
  (process CPU per sum of the integrated build, d-sampled), with the
  model's N' at each sum as the scheduler computes it (the profile
  interpolated in S). The residual sd is 0.21 (leave-one-out 0.24), and
  sum(t) / sum(pred) is 1.03.
* The three sums with more than 256 labels sit on the line (residuals
  -0.29, +0.18, +0.07; a step for them fits -0.03 +- 0.19). A label term
  (slope -1.9 +- 1.0) or the plain law's 8k hinge (0.3 +- 0.5) did not
  help, so the law has neither.
* Below N 3k it is not fitted and runs 1.1-1.6x above the measured sums
  (N 2-3k), where the plain search is chosen anyway.
* Online only its level is learned (prior sd 0.2; `DFIRST_TIME_LAMBDA`
  holds the rest), as for the plain law. Its data are the "dsum" records
  with d_stride 1, no truncated V_d, at least 32 d, and at least
  `DFIRST_TIME_MIN_N` = 2000 vectors (below that the overheads outside the
  d loop dominate: the end-to-end test's sums of 300-400 vectors ran 5x
  the law). A part [d_lo, d_lo + nd) of a sum is scaled to the whole sum:
  its overhead (cpu - time + index_time) plus its d loop (time -
  index_time) divided by the part's share of the d loop. The row's
  weight is that share, so a sum split into m parts weighs 1, not m. (With
  one row per part, a 25k sum split into 150 units of 120 s moved the
  level from +0.19 to -0.12 against 50 whole sums at 6k, and shrank the
  sd: the parts share the sum's N' error, small parts are noisy, sd
  0.26-0.43 in ln for parts of 100-1,000 d, and the tail parts biased.)
* The d-first records stay out of the plain law. The calibration streams
  enter it (below).
* Back-test (the review; engine 3, 16 distinct sums, mean prediction
  with e^{sd^2/2}): obs/pred geometric mean 0.94, sd(ln) 0.23. By N':
  0.71 below 3k (outside the fit, where plain is chosen), 1.16 at 3-6k,
  0.93 at 6-12k, 0.96 at 12-24k, 1.10 above 24k; 1.01 at 129-256 labels
  and 0.96 above. Out of sample: 13 7 4 3 1 1 / 1950 1.24, / 2100 1.08,
  12 6 3 2 1 1 / 950 0.96, the six pool sums at ~6.1k 0.68-1.27
  (geometric mean 0.89), the prototype's engine-2 rows 1.01 (sd 0.32).
  Most of the scatter is per P (about +-25%): 13 7 4 3 1 1, 9 of the 16
  sums, runs at ~1.05-1.1, 14 7 4 4 1 0 0 1, 9 6 4 3 1 1 1 1 and 12 6 3 2
  1 0 1 at 0.6-0.75. A label or x term on the residuals is not identified.
  With the switch at ~4.3k this law sets the density of every sum above
  it.

**The cost along d is not uniform.** In the perf verifier's 13 d logs
(N 2-32k), the mean CPU per d by decile of d's index u = d / N, relative
to the sum's mean, is

    0.71  1.08  1.08  1.21  1.16  1.21  1.15  1.15  0.86  0.39

with the same shape below and above N 6k (se 0.01-0.04 per decile) and a
flat |V_d| along u. Single d vary with CV 0.5-0.67. A range [lo, hi) costs
W(hi/N) - W(lo/N) of the d loop, W the cumulative of the profile
(`amodel.DFIRST_COST_PROFILE`, `dfirst_cost_frac`). W is at most 0.075
from u (W(0.5) = 0.524, W(0.9) = 0.961). The yield per d was not measured
(no (square, d) pair in the samples). The scheduler credits a part with
the same share of the sum's E as of its CPU, so every part of a sum has
the sum's density. Under uniform yield instead, the E credited to a sum
searched in part would differ by at most 7.5% of the sum's E, and a
complete sum's E is the same either way.

**Units of d.** A d-first sum predicted to take more than 1.5 x
`--unit-time` gets units of d, one sum each, with `--d-range lo:hi`:

* each unit is about `--unit-time` of d loop, by the cost profile;
* the unit that carries the sum's calibration stream (normally the first)
  gets `--unit-time` minus the stream of d loop, at least a quarter of a
  unit (the stream runs after the unit's d loop, in one piece). Before,
  its d loop alone was a unit, so with the stream it ran 2.6x / 7.8x /
  23.5x `--unit-time` 120 at N' 14.7k / 22.1k / 30.9k;
* a unit runs from the first d not searched to at most the next part
  already searched;
* a tail below a quarter of a unit is merged into the last unit, so a
  unit takes at most ~1.5 units (tested at N' 30k);
* while the sum's number of d is only predicted (N' predicts nvecs_raw),
  the last unit runs to the end (`lo:`). After the first unit has run,
  the dchunk records give it exactly, and the later units are priced at
  the sum's own N (the d-first law's slope applied to nvecs_raw / N'):
  the review's live run had parts at obs/pred 1.21 / 1.14 / 1.05 on a sum
  with nvecs_raw / N' = 1.03 (x1.105 in time) and 0.63-0.88 on one at
  0.98 (x0.92);
* a d-range unit passes `--d-chunk` = 1/8 of its d (at least 8), so it
  writes 8 checkpoints and `--time-limit` can stop it between them (with
  the default 256, a unit of < 256 d, every unit above N' ~15k at 120 s,
  was one chunk).

Smaller d-first sums share units, all d-first, like plain sums. Sums below
the model's grid go with the first sum inside it. If that sum gets units of
d, they are skipped, as `Planner._walk` already skipped such gaps.

**The parts merge into coverage** (`Summary._dcover`, `Summary.dcov`).

* msearch's "dchunk" records now carry `nvecs_raw` and `truncated`. Every
  chunk with d_stride 1 is merged as a half-open interval of d into
  `dcov[P][S]`. This holds also for killed units (a unit keeps every chunk
  it finished; a partial last line is ignored), duplicates (a unit run
  twice) and overlaps.
* Chunks of a `--d-stride` sample are ignored.
* A sum is covered once its intervals cover [0, nvecs_raw), or by one
  complete "dsum". It is then counted once (`dfirst_sums`), also when
  searched again.
* A chunk with a V_d at `--node-limit` counts as searched, as a truncated
  plain sum does (`dfirst_truncated` in the report).
* The "done" record of a d-range unit still leaves the sum out (its
  incomplete dsum is a hole).
* Records of a sum already searched in full (a duplicate unit, e.g. one
  still running when the scheduler restarted, or the dsum of the unit
  whose chunks completed the sum, read in a later pass) are skipped, so
  they no longer recreate a `dcov` entry that `report` showed as
  "searched in part" forever; incremental and from-scratch summaries
  agree.
* A plain search of a sum searched in part replaces its parts; with
  `--dfirst off` such a sum is searched plain in a normal unit (before,
  the parts mask stopped the unit before it, and a P whose sums below
  the grid went with it was dropped with a zero-score unit).
* `dfirst_cpu` counts each chunk's d loop as it completes (also those of
  killed units, whose d count as searched) and the rest of a dsum
  (index, reduction, enumeration share) with it.
* The planner keeps `dpart` per P (from the summary in `run`, simulated in
  `forecast` and `emit`). It continues the frontier sum from its first gap
  and credits each unit with its share (`frac`) of the sum's E and CPU.
* While a unit of a P runs, the P is busy, so the units of one sum run one
  after the other, each planned from the records of the previous ones.

**The star cover** (msearch `--dfirst-star K`, research/ideas.md "The star
cover of the d loop"; integrated with the units in round 2). A d-first sum
may skip the d through one number x* except every K-th (by rank over the
whole sum) and still find every magic square, because the two diagonals of
an even square share no number. That holds for one x* only: two parts
searched with different x* could each skip one diagonal of the same magic
square. So:

* Every d-first unit gets `--dfirst-star K` explicitly (K = 4 =
  `DFIRST_STAR_K` for even n, 0 for odd n), so that its cost and records
  do not depend on msearch's default.
* msearch's "dchunk" records carry `star_x` and `star_k` (and `star_only`,
  `nd_star_skipped`, `nd_other_skipped`, `pairs_star`).
* `Summary._dcover` keeps one coverage per (x*, K) of a sum ("-" for
  chunks without a star cover, e.g. engine 3's). Chunks of different
  groups never merge; the sum is covered once one group covers every d.
  `--dfirst-star-only` chunks (the star d only) and star chunks without
  `star_x` (the c2/star prototype) count towards no coverage
  (`dfirst_star_skipped` in the report).
* The planner continues the group with the most of the sum's d loop (by
  the cost profile; `Summary.dfirst_part`) and passes that group's K and,
  with `--dfirst-star-x`, its x*, so that the new chunks join it even if
  msearch would now choose another x*. Fresh sums get msearch's x*
  (deterministic over the whole unreduced list, so every unit of a sum
  agrees). The launched record has `star_k` and `star_x`.
* The pairs of a covered sum and their estimate over every d come from its
  group's chunks: a chunk's estimate is its pairs on d without x* plus K x
  those of the star d it searched (none with K = -1: no estimate), each
  chunk weighted by the share of its range not yet covered in the group
  (a duplicate adds nothing). Disjoint chunks add up exactly to msearch's
  own `est_pairs` of the whole sum (tested on msearch output). A sum
  covered by one complete "dsum" takes its `est_pairs` (`dsum_est_pairs`).
* The d-first time law learns "as run": from engine 4 on only from "dsum"
  records with v2's K (engine 3 records have no star cover), never from
  star-only ones. A part's share of the sum is that of [d_lo, end) with
  end = d_lo + nd + the d the star filter skipped (`dsum_span`; nd alone
  counts only the d searched, and d_hi of a stopped unit is beyond).

**The calibration stream** (`--calib-r1-stride k`, `--calib-frac`, default
0.07). msearch runs it after the d loop of every d-first sum of the run. v2
passes it to one unit per sum: the first planned while the summary has no
stream of that sum, normally the unit from d 0.

* k = round(1 / (0.07 r(N'))) by the ratio law, i.e. the plain CPU taken
  as t_dfirst / r: ~15 at the switch, ~18 at 6k, ~41 at 25k. The stream
  records ~1/k of the sum's semi-magic squares. (Before, k came from the
  plain law, which under-predicts at 5-8k, so the stream cost ~9.5%
  instead of 7% there, and 10.8% at 13 7 4 3 1 1 / 2000 in the live run;
  above 11k, where the plain law over-predicts up to 3x, it cost < 3% and
  recorded ~3x fewer squares.)
* k is raised so that the stream takes at most one `--unit-time`
  (`DFIRST_STREAM_MAX`): at 120 s units that binds above N' ~15k, where
  the stream then costs less than 7%.
* msearch skips the stream of a sum whose d loop `--time-limit` (or
  `--total-nodes`) stopped (it ran 7.2 s past a 1-s limit before); the
  summary then has no stream of the sum, and the unit that continues it
  carries it.
* The planner's density charges it to the sum as 0.07 t_d. In a split
  sum the unit that carries it has it in its time, and the other units'
  scores leave it out.
* Its squares enter the cells, i.e. the squares and traversal GLMs and
  the per-P factors (`Summary._ingest`, csum mode "calib"):
  * squares as w x the sampled count against w x E_sq / k;
  * traversals as w x the csquares' counts against w x their expectations
    per square, as for a plain sum's squares;
  * the weight w = min(1, PHI_SUM / (phi_s + PHI_SUM / k)) discounts the
    first-row clustering. phi_s = se_squares^2 / (k est_squares) is the
    stream's own dispersion where it sampled >= 5 squares, and 2.5
    otherwise.
* So at N' >= 5k, where no plain sums will be run any more, the class and
  per-P factors keep learning from ~20-60 squares per sum.
* No double counting:
  * a second stream of a sum, or a stream of a sum with a plain "sum"
    record, is not counted (`calib_dup`);
  * a plain search after a stream subtracts the stream's contribution
    (kept per sum in `perP[P]["calS"]`);
  * the csquares wait for their csum across incremental reads, and those
    of r1-sampled research runs stay out.
* Its est_time (plus the reduction and the enumeration share, as a "sum"
  record's cpu) is an unbiased estimate of the sum's plain CPU (relative
  se ~7%). It is a row of the plain law, weighted by its precision
  (1 / (1 + (se/est)^2 / 0.3^2)), so the plain law's band offsets learn at
  and above the switch, where no plain sums run (in the review's live run
  the streams measured plain at 1.2-1.65x the plain law). The sum's N and
  label bias are recorded with it, once per sum.
* With the sum's d-first CPU (its parts' whole-sum estimates, weighted,
  once they cover at least half the d loop and the sum is searched in
  full) it is one pair y = ln(t_dfirst / t_plain) of the ratio law
  (`Summary.ratio`, per engine). The level a0 is their posterior mean
  with prior sd 0.1 and residual sd 0.17 (`amodel.dfirst_ratio_level`):
  `refit` and `report` show it and the switch N' it implies.
* Its CPU is reported as the stream's (`calib_cpu`).
* `report` prints the streams' est_squares against the model (x SQ12,
  before class factors); a plain search that replaces a stream now
  removes it from these totals too.
* A notable csquare (best >= 7) joins the notable squares once (flagged
  calib, and `report` says so). So does a "dsquare" with best >= 7 (SP,
  SP+S, SP+P: in a d-first sum the only record of these), once per (P,
  S, hash) and flagged d-first; before, only magic ones were kept.

**A check on real output** (msearch, 12 6 3 2 1 1 / S = 950, N = 5,704,
not among the sums of the cost profile; 100 CPU-s):

* The sum was run as three `--d-range` units: d 0-1,901 with
  `--calib-r1-stride 20`, then d 1,901-3,802 killed after two chunks of
  128, then d 3,802 to the end.
* The summary held [0, 2,157) of 5,704 after the first two units. The
  resumed unit ran 2,157-3,802, and the sum was then covered and counted
  once.
* The stream cost 6.6 CPU-s, 7.1% of the 93.5 CPU-s d loop, and sampled
  4 squares (est. 80).
* Each part's whole-sum estimate for the time law:

  | part | share of d | share of CPU, cost profile | measured share of CPU | whole-sum CPU, by the profile | if CPU were uniform in d |
  |---|---:|---:|---:|---:|---:|
  | d 0-1,901 | 0.333 | 0.327 | 0.337 | 96 | 95 |
  | d 2,157-3,802 | 0.288 | 0.340 | 0.369 | 102 | 120 |
  | d 3,802-5,704 | 0.333 | 0.278 | 0.250 | 84 | 70 |

  The true whole-sum CPU is 93.5 s. The profile's estimates (84-102) are
  closer than the uniform ones (70-120).

**Forecast with and without d-first.** `forecast` simulates with the mode
choice and then repeats the same greedy with `--dfirst off`. It prints E
at the marks for both and the CPU d-first needs to reach the plain-only E
of the budget. Both arms are planned with the scheduler's laws, and by
default each unit is also charged its planned time (`--truth laws`). That
cannot measure the gain: the plain law under-predicts the plain search 1.35x
at N' 5-8k (so plain looks cheaper than d-first exactly where the ratio
says d-first wins) and over-predicts it up to 3x above 11k. `--truth
anchored` charges each planned unit instead under the review's
measurement-anchored truth (T2): a d-first sum costs the d-first law (it
back-tests at 0.94), a plain sum the d-first law / r(N') at N' >= 5k, the
plain law at <= 3k, blended geometrically between; a stream costs the
plain cost / k.

Fresh state, shipped calibration and laws at their priors, 10% of the
candidates with the budget scaled (seed 1), after the fixes below:

| CPU-years | 0.1 | 1 | 10 |
|---|---:|---:|---:|
| E with d-first, `--truth anchored` | 0.0467 | 0.117 | 0.246 |
| E plain only, `--truth anchored` | 0.0466 | 0.115 | 0.231 |
| ratio, anchored | 1.004 | 1.021 | 1.066 |
| ratio, the scheduler's own laws as truth | 0.990 | 0.959 | 0.957 |

* Under the anchored truth the plain-only E of 10 CPU-years takes 8.06
  CPU-years with d-first: 1.24x less CPU. d-first units hold 26% of E in
  38% of the CPU at 1 CPU-year, 54% of E in 75% at 10.
* The review's truth simulation (truthsim.py; T2, 10% sample seed 1 at 1
  CPU-year, seed 2 at 10, 1% at 100) gives for this label-free choice
  x1.017 / 1.054 / 1.138 at 1 / 10 / 100 CPU-years and 1.25x / 1.78x less
  CPU at 10 / 100; for the first version (the laws' quotient) x1.004 /
  1.019 / 1.092 and 1.09x / 1.50x. Retrospective.md section 7's "d-first
  on top" of the plain search is x1.03 / 1.07 / 1.13 (its x1.11 / 1.14 /
  1.19 are against fdb77fc and include the plain search's own speed-up,
  which v2's engine-3 plain prior carries in both arms).
* So the d-first gain v2 realises is about what the retrospective
  predicted once the mode is chosen by the measured ratio. The first
  version's smaller figures (x1.002 / 1.019 / 1.072 / 1.140 at 1 / 10 /
  100 / 1,000 CPU-years under its own laws) came mostly from its mode
  choice, not from v2's plan; and the own-law ratio at 100 CPU-years
  (x1.08) lay between the review's two truths (1.07, 1.09), so it was
  not an upper bound either.
* Under its own laws the forecast now shows d-first as a loss (0.96):
  those laws price d-first at 4.3-8k above plain. That is the reason the
  choice no longer uses them, and why `--truth laws` (the default, for
  E itself) should not be read as the gain.
* Charged under the anchored truth, E itself is 2-3% lower at 1-10
  CPU-years than under the planning laws: the plain law under-predicts
  the plain sums at 3-5k.

**Fixes after the review of the d-first scheduler.** Two verifiers (a
back-test of the laws and the mode choice with 6 new pool sums, an
adversarial review, and a live run with kills and restarts) found:

* the mode choice by the two laws' quotient switched at ~9-10k instead of
  ~4.3k (above, "Mode per sum"): now the ratio law, label-free, its level
  learned from the streams;
* the d-first time law took one row per part, so a sum split into m units
  counted m times: now weighted by the part's share (above);
* the first unit of a split sum ignored its stream (up to 23.5x
  `--unit-time` at N' 31k), units of < 256 d had no checkpoint, and the
  stream ran past `--time-limit`: now sized with the stream (capped at one
  unit), `--d-chunk` 1/8 of the unit, and msearch skips the stream after a
  stop;
* the plain law never saw the streams' est_time and the stride came from
  the plain law: now the streams are plain-law rows and ratio pairs, and k
  comes from the ratio;
* minors: records after a sum was covered recreated a dcov entry;
  `--dfirst off` dropped a P whose first grid sum was searched in part;
  a replaced stream stayed in the report's totals; SP-type squares found
  d-first were dropped; killed units' CPU was lost; later parts used the
  model's N' with the sum's true N known; the report's stream rows
  showed "sums 0" and no "(calib)" flag, and the N bias ignored d-first
  sums.
* Not taken: the review's suggestion to keep one row per sum in the
  d-first law at the sum's cover (weighting by share is simpler and
  teaches the law before a large sum completes).

**Tests** (`scripts/test_scheduler.py`, ~45 s in all):

* `test_dfirst_merge`: duplicates, overlaps, a killed unit, a sampled
  part, chunks without nvecs_raw, incremental reads and save/load, the
  weighted time-law rows of parts, the CPU of a killed unit's chunks, and
  a plain search after parts.
* `test_dfirst_plan`:
  * the mode choice by the ratio law (one switch N', ~4.29k, unchanged
    when the plain law moves) and against 16 measured sums (the 13,
    1850 / 1950 along 13 7 4 3 1 1, the six pool sums at ~6.1k);
  * grid_density against eval_sums;
  * the d-first prior and an online level with the slope held;
  * a sum split into units of ~600 s whose fractions add to 1 and whose
    E adds to the sum's, the first unit's d loop plus its stream one
    unit, the stride from the ratio law;
  * the stream on the first unit only;
  * resume from records with gaps, priced at the sum's known N, with
    `--d-chunk`;
  * every unit of a sum at N' 30k (120-s units) within 1.5 units;
  * `--dfirst off` on a sum searched in part, also as the first sum in
    the grid with sums below it;
  * whole-sum d-first units;
  * upper bounds of the lazy heap with d-first.
* `test_calib_cells`: the streams' cell contributions, weights,
  duplicates, replacement by a plain search (cells and report totals),
  r1-sampled records, the plain-law row of a stream, and the GLM learning
  from streams.
* `test_dfirst_learning`: one sum in 1 or 40 parts gives the same d-first
  law; a stream's plain-law row and, once the sum is covered, its ratio
  pair and the scheduler's learned level; records after a sum is covered
  (no dcov entry; incremental equals a full parse); a d loop stopped by
  `--time-limit` (nd < d_hi - d_lo); a truncated chunk; SP-type and magic
  dsquares in the notable list; a multi-sum unit stopped inside a sum
  resumed from its last chunk with the stream.
* `test_dfirst_star`: synthetic star records through both schedulers
  (a complete K = 4 sum's pairs by its estimate, a duplicate K = 0 run, K =
  -1, a star-only part, magic squares found once and twice, the d-first law
  from v2's K only) and msearch's units of one sum (one x*, chunks adding
  up to the whole sum's est_pairs, `--dfirst-star-x`).
* `test_dfirst_star_cover`: the coverage per (x*, K): two x* never cover a
  sum together, star-only and star_x-less chunks count towards none, the
  planner continues the largest group with its x* and K, duplicate and
  overlapping chunks count once in the pairs and their estimate, K = -1
  covers without an estimate, engine-3 chunks are a group of their own.
* `test_dfirst_e2e`: `run` with msearch, `--dfirst on --dfirst-min-n 0` and
  tiny units, so that sums are split into 2 d-range units. The second unit
  starts where the first one's records stop and knows the sum's number of
  d, the stream comes with the first unit, and the sums are merged into
  coverage; the chunks of every covered sum are disjoint and cover [0,
  nvecs_raw), with one (x*, K): K 4 on every unit, msearch's x* on the
  first unit and the same x* passed to the units that continue the sum.
  Then `report`, `emit` and `forecast` with d-first.

**Not done.**

* The plain law's shape for engine 3 (obs/pred 1.35 at 5-8k, 0.65-0.69 at
  8-24k, 0.30-0.34 above 256 labels): the streams now teach its band
  offsets, but one offset for N' >= 6k cannot fix a shape. With the
  switch at ~4.3k it prices only `--dfirst off` and the plain sums below
  the switch.
* The ratio's slope is held; only its level learns. The P-level scatter
  of the d-first law (+-25%) is not modelled.
* The yield per d is assumed proportional to its cost.
* A sum's units run one after the other (one unit per P at a time). A
  sum of 11 CPU-hours at the default `--unit-time` 120 is ~330 units with
  ~0.3 s of enumeration and index each (0.3%).
* A stream's plain-law row stays when a later plain search of the sum
  replaces its squares in the cells (the sum then has two plain-time rows).
* v1 (`--model regression`) still never launches d-first units, and counts
  a sum searched in parts as unsearched.

## Fixes after verification

* Time law: label-word step, prior refit with sd 0.30, shape fixed online,
  band offsets, CPU time from msearch.
* Quasi-Poisson GLMs and report intervals; per-P squares factor with
  A_SQ = 15 and tempering by 2.5. Per-P observed counts are taken from the
  cells, the same sums as the expectations.
* GLM priors: squares N' < 1k -0.25, 6-12k -0.22, x < 0.1 -0.10; S
  traversals N' >= 3k +0.04 (was +0.09).
* Summary: a square counts only once its sum record arrives. A unit killed
  mid-sum used to leave orphan squares, which its rerun then counted
  again. Sampled records are skipped.
* `run` exits when no unit is left, polls faster while units are short,
  and defaults to cores - 1 workers. `emit` no longer needs a local msearch.
  `--model` defaults to regression for n != 6. Unit file names use a
  counter instead of globbing `units/` on every launch. A coverage gap of
  negligible predicted density (< 1e-3 of the P's best) no longer holds
  back a P's later sums.
* Caches: the profile store is keyed by PROFILE_VERSION and a hash of the
  model constants and lattice table; the summary by AMODEL_VERSION, a hash
  of the constants it depends on and the cell edges; the pool by
  POOL_VERSION.
* `forecast` prints the draws' median and mean, labels the point, and
  prints the selection-discounted E.

## Files

Scratch (not in the repo): `sched-v2/fix/` (tdata.py builds the time
dataset; tfit.py and tfit2.py fit the candidate laws; tcheck.py runs the
held-out and refit checks; fc_*.out holds the forecasts; report_all.out is
the all-data report). The verification scripts are in `sched-v2/backtest/`,
`sched-v2/vfc/` and `sched-v2/live/work/`. For the d-first units:
`sched-dfirst/` (tfit.py fits the d-first law, dprofile.py and dhomog.py
the cost profile along d from the perf verifier's d logs in
`integ/perf/runs/`, xover.py the switch points, realsplit.py the check on
real output, fc_*.out the forecasts).
