# Scheduler v2: the analytic model as a scheduler (October 2026)

`scripts/scheduler.py` (default `--model analytic`) chooses the next
(P, S) to search by expected magic squares per CPU-second, using the
existence study's analytic model ([existence.md](existence.md) sections 2-6,
[existence/analytic.md](existence/analytic.md)) instead of v1's fitted
Poisson regression. This note records what it does, what it was checked
against, and how much to trust its forecast. The previous scheduler is still
available as `--model regression --pool classic` (and is the default for
`--vec-size` other than 6).

## Short answer

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
* The quoted figure is the **discounted** one, about **0.1 in 1 CPU-year and
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

So the squares model holds on held-out data within about +-10% on average.
What v2 selects comes in 5-15% low on squares and somewhat low on P(magic)
(0.7-1.05). The **selection discount of 0.8** that `forecast` prints is the
middle of that range.

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
  P (about 620k units per CPU-year). There is also no d-first mode in the
  scheduler (it launches plain msearch units; the Summary reads the records
  of `msearch --diag-first` runs, covering their complete sums and keeping
  them out of the fits, see ideas.md, "Integration of cx/wide, cx/pretest
  and cx/dfirst"), and the per-P prior strengths have not been re-estimated
  beyond A_SQ.

## Fixes after verification (this revision)

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
`sched-v2/vfc/` and `sched-v2/live/work/`.
