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
  S_min). The model predicts **60.5 (square, SP traversal) pairs in those
  134 units in 60 h, and 45 of them after 51.4 h** (55.4 h without
  `--stage1`). Counted by sum, the band holds **91.6 pairs**: 243 units that
  start below N' 3k reach into the band with 31.1 more, and 45 band pairs
  come after **40.3 h** (make_plan's recount; the 51.4 h above counts only
  the 134 units). Without stage 1 only 9 units (0.28 CPU-h) of the first
  60 h are d-first, so turning the band plain changes 0.5% of the CPU and E
  not at all (forecast --shipped: 0.1202 / 0.1843 / 0.2827 at 1 / 3 / 10
  CPU-years, with and without stage 1, to 4 digits).
* **The band is measured at its low end.** The 134 band units' sums lie at
  N' 3.0-4.2k (first sums 3.0-4.0k; 100 of them end below 3.6k, with 77% of
  their predicted pairs). The 4.2-6k part of the band is not reached in 60 h.
* The premise "45 pairs = 39 CPU-hours of plain search at 3-6k" was 10x off
  for the scheduler's order: the top 3-6k sums give **~14 pairs per CPU-hour**
  (model), against 1.8 for the calibration search's sums (drawn across the
  whole 10-year plan). 45 band pairs need only ~3.2 CPU-h of 3-6k units;
  the scheduler interleaves them with ~48 h of N' < 3k units, which bring
  ~900 more pairs.
* **What the 45 pairs buy.** With Poisson-like counts (dispersion phi 1-2.1)
  decide.py's 90% band of E at 1 CPU-year narrows from **x5.8 to x3.3-3.8**
  (simulation, shipped prior). The pairs found so far are not clustered
  once each sum is counted once (section 4): outside the band phi 1.05, in
  it 1.0. The first version of this note read phi 6.8 and f_rho ~4 from
  the state's raw records, where benchmark and regression sums searched
  8-16 times each count every pair 8-16 times; that risk (x4.6 at 45 pairs,
  ~120 pairs needed) is withdrawn.
* **Measurements.** (a) 203 units of the scheduler's own order (0.32 CPU-h,
  all at N' < 3k): squares 0.979 [0.957, 1.000] of the frozen prediction,
  S 0.98, P 0.93; pairs 8 against 12.2 (0.66 [0.33, 1.19]). (b) 10 random
  band units of the committed plan (verification, 1,323 CPU-s): squares
  1.09, S 1.04, P 1.08, pairs 6 against 4.58 (1.31 [0.57, 2.59]; 3 of them
  in one P), CPU 1.15x the reference; decide.py gives f_rho 1.08 median
  (shipped prior) on them.
* **E gain: x1.00, as expected.** The value is the decision. Example
  (decide.py's predictive, shipped prior): P(>=1) at 10 CPU-years is 28.7%
  before stage 1; it would read ~43% if the band came in at f_rho = 1.4
  and ~13% at f_rho = 0.7 (f_rho known exactly; 25% at 1).

## 1. What was built

* `scheduler.py run / plan / emit / forecast --stage1 [HOURS[:LO:HI]]`
  (default 60:3000:6000). For the first HOURS of reference CPU (the units'
  predicted CPU at the scheduler's laws, i.e. the fast x86 build's clock,
  kept in the state's `stage1_6.json` and counted from the start of stage
  1) the sums with N' in [LO, HI) are searched plain (`AnalyticScorer._choose`);
  the planner's heap is rebuilt when the clock runs out. HOURS is the
  total: the same band keeps its clock, so `--stage1 120` after a finished
  `--stage1 60` runs 60 more hours (120 in all); another band starts a new
  clock; units launched after stage 1 ended are not counted. `plan`, `emit`
  and `forecast` count their simulated units on a copy of the clock
  (`forecast --sample f` scales HOURS by f). Under stage 1, `run` writes
  the model's predictions of each unit's records at launch into
  `launched_6.jsonl` (`unit_predictions`: squares, S and P traversals, SP
  traversals before the class factors, pairs = k_SP r_S r_P 720 e^{lpSP}
  per square at f_rho = 1). `forecast` now also prints E at 0.3, 3 and 30
  CPU-years and E at every mark to 4 digits. `run --stage1` works with
  `--machine cal.json` (the clock counts reference CPU; the unit files
  carry the machine record; test_scheduler.py `test_stage1_machine`).
* `scripts/laptop/make_plan.py --stage1` writes the same plan statically
  (for the laptop or rented machines without the scheduler state), with
  those predictions per unit, and a PREREGISTERED file: the plan's and the
  decision table's sha256, the code commit, the totals by N' band and mode,
  the N' of the band units and the pairs of every sum in the band. The plan
  of this note is `stage1/plan-stage1.jsonl.xz` (21,298 units, 60.0 h,
  E[magic] 0.0137, 999.7 pairs of which 60.5 in the 134 band units),
  registered in `stage1/PREREGISTERED.txt` before any of its units ran.
  Its state is the laptop plan's, and its first 15,164 units (40.0 h) are
  research/laptop/plan-20261010's, unit for unit (the same arguments), so
  E is the same at 40 h (0.0109).
* `scripts/decide.py STATE` (or `STATE --plan PLAN --units-dir DIR` for a
  static plan run by scripts/laptop/run.py, or `--plan PLAN --units-dir
  DIR` alone, without a state): reports, never decides. The state must be
  fully profiled, or made with `--only` and given the same `--only` (else
  the model first profiles the whole pool, for minutes, into the state).
  * the records: by default the units launched under stage 1;
  * each plain sum (P, S) once: of its copies the one kept is a record not
    truncated, then the one with the most squares, then the first file's;
    the others are dropped with their squares, the calibration stream of a sum with a plain record is
    dropped, and a d-first sum searched in full counts only where no plain
    record of it is kept;
  * checks: observed / predicted squares, S and P traversals per N' band
    (learned and shipped class factors), and against the frozen per-unit
    predictions;
  * the f_rho posterior: the f1/p1 fit as written (scratchpad
    `p1-learn-the-sp-coup/frho_block.py`; p1 is not committed, so its fit
    is copied, not its Summary changes): per P the band's plain pairs
    Y_P against M_P, plus the d-first sums searched in full in one record
    (their est_pairs; their calibration streams left out), quasi-Poisson
    with a gamma prior, phi = the dispersion of the per-P totals
    (sum (Y_P - f M_P)^2 / (f (M - sum M_P^2 / M)), unbiased also when one P
    holds much of M) floored at 1 (`--phi-min` raises the floor); exact
    gamma quantiles. Priors: shipped (the
    calibration search's posterior, set by its median 1.06, f_rho^2 ln sd
    0.47: f_rho^2 median 1.12, mean 1.23, the search's 1.12 / 1.24; the
    first version set its mean to 1.06, median 1.041, which read P(>=1)
    0.4-0.8 percentage points low) for new records; pre (its prior, mean 1)
    to replay its own records; flat. The fit with the flat prior is always
    printed too: the shipped prior weighs as ~19 pairs (12 of them from the
    assumed extrapolation sd 0.25), so with 45 pairs a far-off truth shows
    only partly in the shipped-prior interval (section 2). `--band LO:HI`
    must be whole N' bands (edges 0, 1k, 3k, 6k, 12k, 24k, inf) when a
    state is used; anything else is an error. Pairs outside the band are a
    check only (pooling them needs a model of f_rho's N dependence, which is
    the extrapolation this stage measures), and their dispersion is shown
    as a sensitivity;
  * without a state, the fit takes the finished band units of the plan
    (plan "stage1"): Y = their pairs, M = their frozen pred_pairs x observed
    / predicted squares, so that f_rho is per square found, as with the
    state. On the verification's 10 band units both give 6 / 5.04 and
    median 1.080;
  * the predictive E and P(>=1) at 1 / 3 / 10 CPU-years: E_ship(C)
    (`stage1/ecurve.json`, forecast --shipped) x f_rho^2 draws x
    lognormal(class factors 0.115-0.097, pair 0.2, selection 0.12), with
    the flat-prior P(>=1) next to it, and the CPU for P(>=1) = 25 / 50%;
  * the decision table (`stage1/decision-table.json`, frozen with the
    plan): each row's P(>=1) against the user's threshold X (null until the
    user sets it);
  * `--simulate PLAN`: the coverage gate below.

## 2. Gates

| gate | result |
| --- | --- |
| `--stage1 off` (or none) reproduces today's plan | emit 3,000 units and a 120-hour make_plan (38,646 units, 72 of them d-first): sha256 identical to 80d15cc's. After the merge with machine calibration (integ/stage1): make_plan reproduces plan-20261010 (a2852ebe...) and the stage-1 plan (a6a4be6a...) byte for byte |
| forecast --shipped 0.120 / 0.185 / 0.283 | 0.1202 / 0.1843 / 0.2827 (the 3-year value is 0.184 to 3 digits; the target's 0.185 was rounded) |
| `--stage1` E at 1 CPU-year drops <= 1% | 0.1202 -> 0.1202 (0.0%); also 3 and 10 CPU-years unchanged to 4 digits (again after the merge: 0.0427 / 0.0718 / 0.1202 / 0.1843 / 0.2827 at 0.1-10 CPU-years with and without `--stage1`) |
| replay: the calibration search's records, 'pre' prior -> median 1.06, f_rho^2 ln sd 0.47 within 0.02 | passes only with `--band 3000:inf` (the band that search was fitted on): median 1.051, f_rho^2 ln sd 0.460 (7 pairs / 5.74: 4.85 from plain squares, 0.89 from its 7 whole d loops, which found 0); without the d loops the median is 1.099. With the default band (3-6k) it gives 1.027 / 0.487 (5 pairs / 4.21) and fails by 0.033 |
| simulation (1,000 replicates per f_rho 0.5 / 1 / 2, per-P gamma clustering, the plan's first 45 band pairs over 74 P): 90% interval covers 85-95%, posterior median unbiased within 5% | flat prior, phi 1.0 / 1.41 / 2.1: covers 86.1-91.3%, bias -2.9% to +0.2%: **passes**. phi 6.8: 77.7-86.1%, bias -1.8 to -10.1%: fails (phi_hat 4.9-6.1 runs low), but section 4 finds no such clustering |
| (same, shipped prior) | f_rho = 1: covers 93.9-95.4%, bias +1.6 to +2.1%; f_rho drawn from the prior: 87.9-91.4%; at f_rho 0.5 / 2 it covers only 35-53% at phi 1-2.1 (the prior pulls the median +32-45% / -14-20%): by design, but a far-off truth needs the flat-prior line, which decide.py prints |
| band <= x4 at 45 pairs (shipped prior, f_rho = 1) | x3.33 / 3.50 / 3.76 at phi 1 / 1.41 / 2.1: passes (x4.61 at phi 6.8, which the deduplicated records do not show) |
| ctest -R fast_ | 51/51 passed (f1/stage1); 67/67 after the merge with the machine calibration and the portable path (integ/stage1); test_scheduler.py all ok (new test_stage1, test_stage1_machine, and test_laptop_run: run.py exits on --hours, SIGINT and SIGHUP and removes partial output), test_machine_cal.py and test_calibrate.py ok |
| measurement: squares within Poisson + 20% | 5,695 / 5,820 = 0.979 [0.957, 1.000]: passes; band units 5,476 / 5,021 = 1.09 |

The simulation and the band rows are from `decide.py --simulate` on the
committed plan with the shipped prior set by its median, exact gamma
quantiles and the dispersion for unequal M_P (integ/stage1, after review);
the first version's numbers (prior set by its mean) differ by up to 6
points of coverage, the review's changes by up to 0.4.

## 3. Measurements

**The scheduler's own order** (one process, nice 10, 20 minutes, 0.32
CPU-hours): `scheduler.py run --stage1 --workers 1 --hours 0.333` on a copy
of the plan's state: 203 units, all at N' < 3k (the first 3-6k unit of the
plan comes after 10.8 reference hours), 0.297 reference hours predicted,
0.321 CPU-h measured.

| | observed | predicted (frozen at launch, f_rho = 1) | ratio [90%] | per CPU-hour (measured) | model per reference hour |
| --- | ---: | ---: | --- | ---: | ---: |
| squares | 5,695 | 5,820 | 0.979 [0.957, 1.000] | 17,700 | 19,600 |
| S traversals | 7,313 | 7,460 | 0.980 [0.962, 0.999] | | |
| P traversals | 1,732 | 1,855 | 0.934 [0.897, 0.972] | | |
| pairs (= SP traversals) | 8 | 12.16 (17.6 before the class factors) | 0.66 [0.33, 1.19] | 24.9 | 41 |

One of the 8 is an SP+P square (P = 2^12 3^7 5^4 7^2 11^2, S = 1220,
best_score 10). The 203 units are the plan's densest: the plan averages 16.7
pairs per reference hour over 60 h and 14.2 in the 3-6k band.

**Band units** (the verification of f1/stage1: 10 band units drawn at
random from the committed plan, one process, 1,323 CPU-s):

| | observed | predicted (frozen in the plan) | ratio [90%] |
| --- | ---: | ---: | --- |
| squares | 5,476 | 5,021 | 1.09 [1.07, 1.12] |
| S traversals | 4,778 | 4,579 | 1.04 [1.02, 1.07] |
| P traversals | 1,156 | 1,075 | 1.08 [1.02, 1.13] |
| pairs | 6 (3 in one P) | 4.58 | 1.31 [0.57, 2.59] |
| CPU | 0.368 h | 0.319 h (reference) | 1.15 |

decide.py on them without a state (`--plan --units-dir`): f_rho 1.080
median [0.745, 1.504] (shipped prior, 6 / 5.04 over 10 P, phi 1.61); flat
prior 1.087 [0.386, 2.355]. With the state's records (decide.py STATE
--plan --units-dir) the shipped-prior line is the same and the flat prior
gives 1.086 [0.387, 2.352].

Projected CPU for 45 band pairs: 40.3 reference hours counting every sum in
the band, 51.4 counting the band units only (what decide.py fits without a
state); x1.08-1.15 measured / reference = 44-59 CPU-h at f_rho = 1. Both
are within the 80-hour criterion (<= 1% of a 1 CPU-year budget).

## 4. The clustering of pairs, and what it means

The first version of this section read the state's raw records: at N' < 3k
166 pairs against 38 predicted (f_rho 4.3, phi 6.8), a few P with 9-20
pairs against < 1 expected, and concluded that 45 band pairs would give
x4.6 instead of x3.5 and that x4 needed ~120 band pairs. That came from
duplicated records. In the state's 1,324 files the 169 counted SP pairs are
only 33 distinct squares (x5.1; squares overall x1.21): the sums holding
them are the benchmark and regression sums, searched 8-16 times each
(the verification's `tools/dups.py`). decide.py now counts each plain sum
(P, S) once (section 1).

| records of the plan's state (1,324 files) | outside the band (flat prior) | inside the band (3-6k) |
| --- | --- | --- |
| raw records (decide.py of 3f1cb6a, `--select all --prior flat`) | 158 pairs / 26.6, phi 6.1 | 12 / 7.7, phi 2.1 |
| each plain sum once (the verification's dedupe.py) | 26 / 17.6, phi 1.12, f_rho 1.46 [1.01, 2.02] | 8 / 7.0, phi 1.0 |
| each sum once, also the calibration streams and d-first sums of a sum with a plain record (decide.py) | 25 / 18.1, phi 1.05, f_rho 1.36 [0.95, 1.88] | 8 / 7.05, phi 1.0, f_rho 1.08 [0.77, 1.47] (shipped prior) |

(decide.py dropped 12,495 repeated sums with 5,726 squares and 136 pairs,
and left out 84 d-first sums that also have a plain record. Of the copies of
a sum it keeps a record not truncated, then the one with the most squares;
keeping the first file's, as first written, dropped 300 more squares, among
them a band sum's 77 behind a truncated record with none: 0.05 model pairs.) So the
evidence is phi ~1 inside and outside the band, and the "band <= x4 at 45
pairs" criterion passes (x3.3-3.8 at phi 1-2.1). The f_rho ~1.4 outside the
band is likely biased up still, by which sums were chosen for the
benchmark and regression sets (sums known to hold squares); it is a check
only and does not enter the fit. decide.py still prints the fit at the
outside dispersion as a sensitivity line (now ~1.0, so it changes nothing),
and `--phi-min` can floor the dispersion if the band's pairs turn out
clustered after all.

What 45 band pairs are worth, as the band of E at 1 CPU-year (x5.8
before stage 1; shipped prior, f_rho = 1; the last column from the first
version's simulation):

| phi of the band's pairs | 45 pairs | pairs for x4 |
| ---: | ---: | ---: |
| 1.0 | x3.3 | 18 |
| 1.41 | x3.5 | 25 |
| 2.1 | x3.8 | 37 |
| 6.8 (not seen) | x4.6 | ~120 |

## 5. Using it

The laptop runs the static plan (research/laptop/README.md):

```
python3 scripts/laptop/run.py research/stage1/plan-stage1.jsonl.xz --out laptop-runs
python3 scripts/decide.py --plan research/stage1/plan-stage1.jsonl.xz --units-dir laptop-runs
```

With the scheduler and its state:

```
python3 scripts/scheduler.py run --stage1          # 60 h, then on as usual (resumable)
python3 scripts/decide.py data/sched               # any time: the report
# or a static plan for other machines:
python3 scripts/laptop/make_plan.py --state data/sched --stage1 --out plan.jsonl
git add plan.jsonl.PREREGISTERED.txt && git commit   # before the first unit runs
python3 scripts/laptop/run.py plan.jsonl --out runs/
python3 scripts/decide.py data/sched --plan plan.jsonl --units-dir runs/
```

Set the thresholds X in `stage1/decision-table.json` before stage 1 runs
(the plan registers its sha256). Read the flat-prior line next to the
shipped-prior one: if they disagree, the band's pairs are pulling away from
the prior and more pairs settle it. To extend stage 1, give the total:
`--stage1 120` after `--stage1 60` runs 60 more hours; for the static plan,
regenerate it with `--stage1 120` from the original state (its first 60 h
should be the same units; run.py refuses an output directory whose finished
units do not match the plan's lines) and keep running into the same output
directory.

Files: `scripts/decide.py`, `scripts/scheduler.py` (Stage1, unit_predictions),
`scripts/laptop/make_plan.py`, `research/stage1/` (plan, PREREGISTERED.txt,
decision-table.json, ecurve.json). Scratch: `scratchpad/f1/stage1/`
(fc_*.txt forecasts, run.log and state_run the measurement, decide_*.txt),
`scratchpad/f1/verify-stage1/` (bandruns/ the band units, tools/dups.py and
dedupe.py), `scratchpad/integ-stage1/` (the plans regenerated, decide_*.txt,
sim.txt, replay_*.txt).
