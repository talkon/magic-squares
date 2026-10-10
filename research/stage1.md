# Stage 1: learn the SP coupling before deciding how much to spend (October 2026)

**Question.** No measured speed lever moves E at a fixed budget by more than
about 10%, so the user's real decision is whether to spend 1, 3 or 10
CPU-years at all. The forecast's 90% band is about x6 wide, and most of it
is the SP coupling f_rho (f_rho^2, ln sd 0.47; calibration-target.md 3.4).
Can a short, cheap first stage narrow that band before the money goes?

**Short answer.**

* **The pairs are cheap; stage 1 itself adds almost nothing.** Scheduler v2's
  own first 60 reference CPU-hours already search 134 units at N' 3-6k
  plain (the d-first switch is at N' ~3.7k and the best sums sit near
  S_min). The model predicts **60.5 (square, SP traversal) pairs in that
  band in 60 h, and 45 after 51.4 h** (55.4 h without `--stage1`). Without
  stage 1 only 9 units (0.28 CPU-h) of the first 60 h are d-first, so
  turning the band plain changes 0.5% of the CPU and E not at all
  (forecast --shipped: 0.1202 / 0.1843 / 0.2827 at 1 / 3 / 10 CPU-years,
  with and without stage 1, to 4 digits).
* The premise "45 pairs = 39 CPU-hours of plain search at 3-6k" was 10x off
  for the scheduler's order: the top 3-6k sums give **~14 pairs per CPU-hour**
  (model), against 1.8 for the calibration search's sums (drawn across the
  whole 10-year plan). 45 band pairs need only ~3.2 CPU-h of 3-6k units;
  the scheduler interleaves them with ~48 h of N' < 3k units, which bring
  ~900 more pairs.
* **What the 45 pairs buy depends on how clustered the pairs are.** With
  Poisson-like counts (dispersion phi 1-2.1) decide.py's 90% band of E at 1
  CPU-year narrows from **x5.8 to x3.3-3.8** (simulation, shipped prior).
  But the pairs found so far are strongly clustered: outside the band, on
  the 166 pairs of all earlier records, phi = 6.8 (whole P with 9-20 pairs
  against 0.06-0.8 expected). At phi 6.8 the band at 45 pairs is **x4.6**
  (x4.4 at 60), and x4 needs ~120 band pairs (~100 reference hours).
* **First measurement (203 units, 0.32 CPU-h):** squares 0.979 [0.957,
  1.000] of the frozen prediction, S 0.98, P 0.93; pairs 8 against 12.2
  (0.66 [0.33, 1.19]), all at N' < 3k (no 3-6k unit yet). 24.9 pairs per
  measured CPU-hour against the model's 41 per reference hour; CPU 1.08x
  the reference prediction. Projected: 45 band pairs after ~56 CPU-h at f_rho
  = 1, ~70 CPU-h if the band's pairs come at the 0.66 seen at N' < 3k.
* **E gain: x1.00, as expected.** The value is the decision. Example
  (decide.py's predictive, shipped prior): P(>=1) at 10 CPU-years is 28%
  before stage 1; it would read ~43% if the band came in at f_rho = 1.4
  and ~13% at f_rho = 0.7 (f_rho known exactly; 25% at 1).

## 1. What was built

* `scheduler.py run / plan / emit / forecast --stage1 [HOURS[:LO:HI]]`
  (default 60:3000:6000). For the first HOURS of reference CPU (the units'
  predicted CPU at the scheduler's laws, i.e. the fast x86 build's clock,
  kept in the state's `stage1_6.json` and counted from the start of stage
  1) the sums with N' in [LO, HI) are searched plain (`AnalyticScorer._choose`);
  the planner's heap is rebuilt when the clock runs out. `plan`, `emit` and
  `forecast` count their simulated units on a copy of the clock (`forecast
  --sample f` scales HOURS by f). Under stage 1, `run` writes the model's
  predictions of each unit's records at launch into `launched_6.jsonl`
  (`unit_predictions`: squares, S and P traversals, SP traversals before
  the class factors, pairs = k_SP r_S r_P 720 e^{lpSP} per square at f_rho
  = 1). `forecast` now also prints E at 0.3, 3 and 30 CPU-years and E at
  every mark to 4 digits.
* `scripts/laptop/make_plan.py --stage1` writes the same plan statically
  (for rented machines without the scheduler state), with those predictions
  per unit, and a PREREGISTERED file: the plan's and the decision table's
  sha256, the code commit, and the totals by N' band and mode. The plan of
  this note is `stage1/plan-stage1.jsonl.xz` (21,298 units, 60.0 h,
  E[magic] 0.0137, 999.7 pairs of which 60.5 at N' 3-6k), registered in
  `stage1/PREREGISTERED.txt` before any of its units ran. Its state is the
  laptop plan's, so its first units are research/laptop's.
* `scripts/decide.py STATE` (or `--plan PLAN --units-dir DIR` for a static
  plan run by scripts/laptop/run.py): reports, never decides.
  * the records: by default the units launched under stage 1;
  * checks: observed / predicted squares, S and P traversals per N' band
    (learned and shipped class factors), and against the frozen per-unit
    predictions;
  * the f_rho posterior: the f1/p1 fit as written (scratchpad
    `p1-learn-the-sp-coup/frho_block.py`; p1 is not committed, so its fit
    is copied, not its Summary changes): per P the band's plain pairs
    Y_P against M_P, plus the d-first sums searched in full in one record
    (their est_pairs; their calibration streams left out), quasi-Poisson
    with a gamma prior, phi = the Pearson dispersion of the per-P totals
    floored at 1 (`--phi-min` raises the floor). Priors: shipped (the
    calibration search's posterior, median 1.06, f_rho^2 ln sd 0.47) for new
    records; pre (its prior) to replay its own records; flat. Pairs outside
    the band are a check only (pooling them needs a model of f_rho's N
    dependence, which is the extrapolation this stage measures), and their
    dispersion is shown as a sensitivity;
  * the predictive E and P(>=1) at 1 / 3 / 10 CPU-years: E_ship(C)
    (`stage1/ecurve.json`, forecast --shipped) x f_rho^2 draws x
    lognormal(class factors 0.115-0.097, pair 0.2, selection 0.12), and
    the CPU for P(>=1) = 25 / 50%;
  * the decision table (`stage1/decision-table.json`, frozen with the
    plan): each row's P(>=1) against the user's threshold X (null until the
    user sets it);
  * `--simulate PLAN`: the coverage gate below.

## 2. Gates

| gate | result |
| --- | --- |
| `--stage1 off` (or none) reproduces today's plan | emit 3,000 units and a 120-hour make_plan (38,646 units, 72 of them d-first): sha256 identical to 80d15cc's |
| forecast --shipped 0.120 / 0.185 / 0.283 | 0.1202 / 0.1843 / 0.2827 (the 3-year value is 0.184 to 3 digits; the target's 0.185 was rounded) |
| `--stage1` E at 1 CPU-year drops <= 1% | 0.1202 -> 0.1202 (0.0%); also 3 and 10 CPU-years unchanged to 4 digits |
| replay: the calibration search's records, 'pre' prior -> median 1.06, f_rho^2 ln sd 0.47 within 0.02 | median 1.051, f_rho^2 ln sd 0.460 (7 pairs / 5.74: 4.85 from plain squares, 0.89 from its 7 whole d loops, which found 0); without the d loops the median is 1.099 (fails) |
| simulation (1,000 replicates per f_rho 0.5 / 1 / 2, per-P gamma clustering): 90% interval covers 85-95%, posterior median unbiased within 5% | flat prior, phi 1.0 / 1.41 / 2.1: covers 86.2-91.3%, bias -3.6% to +0.2%: **passes**. phi 2.5: 85.5-87.2%. phi 6.8 (measured on the older pairs): 78.8-85.6%, bias -1.9 to -7.7%: **fails** (phi_hat 5-6 underestimates it) |
| (same, shipped prior) | f_rho = 1: covers 93.8-96.3%, bias -0.3 to +1.5%; f_rho drawn from the prior: 86-89%; at f_rho 0.5 / 2 the prior pulls the median +32-60% / -14-34% (45 pairs against the prior's weight of ~19): by design, but a far-off truth needs more pairs to show |
| band <= x4 at 45 pairs (shipped prior, f_rho = 1) | x3.32 / 3.50 / 3.80 at phi 1 / 1.41 / 2.1: passes; x4.60 at phi 6.8: fails (x4.40 at 60 pairs) |
| ctest -R fast_ | 51/51 passed; test_scheduler.py all ok (new test_stage1) |
| measurement: squares within Poisson + 20% | 5,695 / 5,820 = 0.979 [0.957, 1.000]: passes |

## 3. Measurement (one process, nice 10, 20 minutes, 0.32 CPU-hours)

`scheduler.py run --stage1 --workers 1 --hours 0.333` on a copy of the
plan's state: 203 units, all at N' < 3k (the first 3-6k unit of the plan
comes later), 0.297 reference hours predicted, 0.321 CPU-h measured.

| | observed | predicted (frozen at launch, f_rho = 1) | ratio [90%] | per CPU-hour (measured) | model per reference hour |
| --- | ---: | ---: | --- | ---: | ---: |
| squares | 5,695 | 5,820 | 0.979 [0.957, 1.000] | 17,700 | 19,600 |
| S traversals | 7,313 | 7,460 | 0.980 [0.962, 0.999] | | |
| P traversals | 1,732 | 1,855 | 0.934 [0.897, 0.972] | | |
| pairs (= SP traversals) | 8 | 12.16 (17.6 before the class factors) | 0.66 [0.33, 1.19] | 24.9 | 41 |

One of the 8 is an SP+P square (P = 2^12 3^7 5^4 7^2 11^2, S = 1220,
best_score 10). The 203 units are the plan's densest: the plan averages 16.7
pairs per reference hour over 60 h and 14.2 in the 3-6k band.

Projected CPU for 45 band pairs: 51.4 reference hours at f_rho = 1, x1.08
measured / reference = 55.6 CPU-h; at the 0.66 seen at N' < 3k the 45
pairs need 68 predicted ones, which the scheduler's order reaches at ~65
reference hours (~70 CPU-h). Both are within the 80-hour criterion (<= 1%
of a 1 CPU-year budget).

## 4. The clustering of pairs, and what it means

All records of the plan's state (1,527 files, 38k squares) give at N' < 3k
166 pairs against 38 predicted (f_rho 4.3, dispersion 6.8), but the first
203 stage-1 units give 8 against 11.6 (0.66). The difference is a few P:
2^16 3^5 5^4 7^2 has 11 pairs against 0.06 expected, 2^11 3^6 5^2 7 11 13 9
against 0.07, 2^10 3^6 5^3 7^2 11 20 against 2.9. Pairs come in clusters
(a P, and an SP vector d that sits in many squares of a sum), so the count
behaves like ~Y / 7 independent events there. In the 3-6k band the
dispersion is not measured yet (12 pairs in all records: phi 2.1).
decide.py estimates phi from the band's own pairs; with ~45 pairs over ~75 P
that estimate runs low when the truth is ~7 (simulation: 5-6), and the
interval under-covers. decide.py therefore prints the fit at the
dispersion measured outside the band as a sensitivity line; `--phi-min`
makes it the floor.

What 45-60 band pairs are worth, as the band of E at 1 CPU-year (x5.8
before stage 1):

| phi of the band's pairs | 45 pairs (51 h) | 60 pairs (60 h) | pairs for x4 |
| ---: | ---: | ---: | ---: |
| 1.0 | x3.3 | | 18 |
| 1.41 | x3.5 | | 25 |
| 2.1 | x3.8 | x3.6 | 37 |
| 6.8 | x4.6 | x4.4 | ~120 (~100 reference hours) |

## 5. Using it

```
python3 scripts/scheduler.py run --stage1          # 60 h, then on as usual (resumable)
python3 scripts/decide.py data/sched               # any time: the report
# or, on rented machines:
python3 scripts/laptop/make_plan.py --state data/sched --stage1 --out plan.jsonl
git add plan.jsonl.PREREGISTERED.txt && git commit   # before the first unit runs
python3 scripts/laptop/run.py plan.jsonl --out runs/
python3 scripts/decide.py data/sched --plan plan.jsonl --units-dir runs/
```

Set the thresholds X in `stage1/decision-table.json` before stage 1 runs
(the plan registers its sha256). Read the report's sensitivity line: if the
band's pairs are as clustered as the older ones, extend stage 1 (`--stage1
120`) rather than read 45 pairs as x3.5.

Files: `scripts/decide.py`, `scripts/scheduler.py` (Stage1, unit_predictions),
`scripts/laptop/make_plan.py`, `research/stage1/` (plan, PREREGISTERED.txt,
decision-table.json, ecurve.json). Scratch: `scratchpad/f1/stage1/`
(fc_*.txt forecasts, run.log and state_run the measurement, decide_*.txt).
