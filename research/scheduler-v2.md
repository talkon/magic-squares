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
searched d-first where the d-first law's CPU plus its calibration stream
(below) is below the plain law's, and its N' is at least `--dfirst-min-n`
(default 2000). `on` takes d-first at every N' >= `--dfirst-min-n` (for
tests). A d-first unit runs `msearch --diag-first --diag-first-min-n 0`, so
msearch never overrules the scheduler's choice. A unit holds one mode only,
and the launched record (`launched_6.jsonl`) has `mode`, `dlo`, `dhi`,
`nd`, `calib` and `frac`. `emit` prints the same arguments.

* With the shipped laws (engine 3), `auto` picks the mode the
  measurements pick on all 13 measured sums: plain at N = 2.0k, 3.0k and
  4.1k, d-first from 6.7k on.
* Along 10 pool P it switches to d-first at these N' (from where every
  later sum is d-first; 13 6 3 2 never gets there, N' <= 1.6k up to
  2 S0):

  | P | N' |
  |---|---:|
  | 12 6 3 2 1 1 | 4.8k |
  | 10 6 4 2 1 1, 12 6 3 2 1 0 1 | 5.3-5.4k |
  | 13 7 4 3 1 1, 12 7 4 2 2 1, 11 6 4 3 2 1 | 7.1-7.9k |
  | 16 8 4 3 1 1, 14 7 4 4 1 0 0 1 | 9.0-9.3k |
  | 9 6 4 3 1 1 1 1 | 10.7k |

* This is later than the measured crossover (3.8-4.9k) for some P.
  The plain law under-predicts the integrated build's plain search
  between ~5k and 8k for them: obs/pred 1.68 at 13 7 4 3 1 1 / 2000 (N
  7.6k), against 1.01 at 12 6 3 2 1 1 / 988 (6.7k). Its shape for engine
  3 is the open point of "The time model". In such a band d-first would
  save up to ~1.7x per sum.
* `--dfirst on --dfirst-min-n 5000` reproduces msearch's own threshold
  (and retrospective.md 7's assumption) instead.
* At N >= 11k the plain law over-predicts (0.31-0.94). That makes d-first
  look better than it is there, but it is the cheaper mode there anyway.

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
  index_time) divided by the part's share of the d loop.
* d-first records stay out of the plain law, as before.

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
* a unit runs from the first d not searched to at most the next part
  already searched;
* a tail below a quarter of a unit is merged into the last unit;
* while the sum's number of d is only predicted (N' predicts nvecs_raw),
  the last unit runs to the end (`lo:`). After the first unit has run,
  the dchunk records give it exactly.

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
* A plain search of a sum searched in part replaces its parts.
* The planner keeps `dpart` per P (from the summary in `run`, simulated in
  `forecast` and `emit`). It continues the frontier sum from its first gap
  and credits each unit with its share (`frac`) of the sum's E and CPU.
* While a unit of a P runs, the P is busy, so the units of one sum run one
  after the other, each planned from the records of the previous ones.

**The calibration stream** (`--calib-r1-stride k`, `--calib-frac`, default
0.07). msearch runs it after the d loop of every d-first sum of the run. v2
passes it to one unit per sum: the first planned while the summary has no
stream of that sum, normally the unit from d 0.

* k = round(t_plain / (0.07 t_dfirst)), with both times from the laws:
  15-19 where `auto` switches. The stream records ~1/k of the sum's
  semi-magic squares.
* Where the plain law over-predicts (N' > 11k), k comes out larger (~45
  at a law ratio of 0.3). The stream then costs less than 7% and records
  fewer squares.
* The mode choice and the planner's density charge it to the sum: (t_d +
  t_cal) against t_p, with t_cal = t_p / k. In a split sum the first unit
  carries it (that unit's time is `--unit-time` plus the stream), and the
  later units' scores leave it out.
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
* Their CPU is reported as the stream's (`calib_cpu`), not in the plain
  time law.
* `report` prints the streams' est_squares against the model (x SQ12,
  before class factors).
* A notable csquare (best >= 7) joins the notable squares once (flagged
  calib).

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
of the budget. Fresh state, shipped calibration and both laws at their
priors, as in "Short answer" (10, 100 and 1,000 CPU-years on 10%, 1% and
0.1% of the candidates with the budget scaled; one run each):

| CPU-years | 1 | 10 | 100 | 1,000 |
|---|---:|---:|---:|---:|
| sample of the candidates | 10% | 10% | 1% | 0.1% |
| E with d-first (`--dfirst auto`) | 0.126 | 0.266 | 0.436 | 0.603 |
| E plain only (`--dfirst off`) | 0.126 | 0.261 | 0.406 | 0.529 |
| ratio | 1.002 | 1.019 | 1.072 | 1.140 |
| CPU for the plain-only E, plain / with d-first | 1.01 | 1.07 | 1.37 | 1.81 |
| d-first units: share of E / of CPU | 4% / 7% | 13% / 23% | 34% / 56% | 52% / 85% |

* The 1% sample gives x1.008 at 1 and x1.032 at 10 CPU-years.
* The 0.1% sample is noisy: its E at 1 CPU-year is 0.085 against 0.126.
* The calibration streams cost 6.5% of the d-first CPU. That is 0.5% of
  all CPU at 1 CPU-year and 4-6% at 100-1,000.

* The gain grows with the budget, as in retrospective.md section 7, but
  starts later and stays smaller than there ("d-first on top" in T10:
  x1.03 / 1.07 / 1.13 / 1.16 at 1 / 10 / 100 / 1,000 CPU-years). There
  are two reasons:
  * v2's plan is not the ideal frontier. At 1 CPU-year 94% of its E is
    at N' < 6k, where d-first wins little or nothing.
  * The plain law under-predicts the plain search at N' 6-8k (obs/pred
    1.0-1.7), where the retrospective used the measured ratios.
* At 100-1,000 CPU-years the comparison leans on the plain law at
  N' >= 11k, which over-predicts the integrated build's plain search
  (obs/pred 0.31-0.94). That inflates the plain-only CPU there, so these
  ratios are upper bounds. Refitting the plain law's shape for engine 3
  above 11k (open, see "The time model") would tighten them.

**Tests** (`scripts/test_scheduler.py`, 43 s in all):

* `test_dfirst_merge`: duplicates, overlaps, a killed unit, a sampled
  part, chunks without nvecs_raw, incremental reads and save/load, the
  time-law rows of parts, and a plain search after parts.
* `test_dfirst_plan`:
  * the mode choice against both laws and against the measured sums;
  * grid_density against eval_sums;
  * the d-first prior and an online level with the slope held;
  * a ~11k CPU-s sum split into 20 units of ~600 s whose fractions add to
    1 and whose E adds to the sum's;
  * the stream on the first unit only;
  * resume from records with gaps;
  * `--dfirst off` on a sum searched in part;
  * whole-sum d-first units;
  * upper bounds of the lazy heap with d-first.
* `test_calib_cells`: the streams' cell contributions, weights,
  duplicates, replacement by a plain search, r1-sampled records, and the
  GLM learning from streams.
* `test_dfirst_e2e`: `run` with msearch, `--dfirst on --dfirst-min-n 0` and
  tiny units, so that sums are split into 2 d-range units. The second unit
  starts where the first one's records stop and knows the sum's number of
  d, the stream comes with the first unit, and the sums are merged into
  coverage. Then `report`, `emit` and `forecast` with d-first.

**Not done.**

* The calibration streams' `est_time` could teach the plain law above the
  crossover, where no plain sums run. One band offset for N' >= 6k cannot
  fix its shape there (obs/pred 1.7 at 7.6k, 0.3 at 32k), so they are only
  reported.
* The yield per d is assumed proportional to its cost.
* A sum's units run one after the other (one unit per P at a time). A
  sum of 11 CPU-hours at the default `--unit-time` 120 is ~330 units with
  ~0.3 s of enumeration and index each (0.3%).
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
