# Running the search on a laptop (macOS / Apple Silicon, or any Linux)

The laptop runs the **stage-1 plan**, `research/stage1/plan-stage1.jsonl.xz`
(21,298 units on plain sums, 60 predicted CPU-hours of the fast x86 build,
E[magic] 0.0137), with the tools here, without the scheduler or numpy.
It measures the SP coupling f_rho at N' 3-6k: the largest uncertainty of
the forecast (most of its x6 band) and the input to the decision whether
to pay for 1-10 CPU-years (research/stage1.md). It does so at the same E
as the greedy plan it replaces: its first 15,164 units (40 h) are exactly
those of `plan-20261010.jsonl.xz`, and the next 20 h continue in the
scheduler's order with the N' 3-6k sums kept plain.

`plan-20261010.jsonl.xz` (the greedy 40-hour plan, PREREGISTERED.txt here)
is superseded by it; it is kept because it is pre-registered. Outputs of it
count for the stage-1 plan too (the same unit numbers).

## Steps (macOS)

1. Command line tools, once: `xcode-select --install` (gives cc, git,
   python3). Check that `uname -m` prints `arm64`.
2. Get the branch: `git clone -b claude/magic-search-speedups <repo>` (or
   `git pull` in an existing clone).
3. Build and check (about 2 minutes):

       ./scripts/laptop/build.sh

   It compiles msearch, bench and fuzz_arrange into `build-laptop/`, checks
   every bench/quick.txt and bench/full.txt instance against the expected
   squares and hashes and the node counts of the fast x86 build, runs a
   brute-force differential test (also node for node), and prints the
   time of bench/quick.txt. Please report it: divided by the fast x86
   build's ~0.25-0.29 s it is the laptop's time per core against the
   reference (expected ~2.5-3.5x, see below). Stop if it says FAILED.
4. Run, plugged in, lid open (a closed lid sleeps the Mac):

       caffeinate -i python3 scripts/laptop/run.py research/stage1/plan-stage1.jsonl.xz --hours 9

   That is the plan file to pass. It uses the performance cores (`--workers
   K` to change), runs the units in the plan's order, prints progress every
   minute (CPU-hours, the reference CPU-hours they stand for, the stage-1
   band units finished and their pairs) and announces any magic square
   loudly (also written to `laptop-runs/MAGIC.txt`). Ctrl-C stops it
   cleanly; the same command the next night continues where it stopped
   (finished units are skipped, an interrupted one is rerun).
5. The decision report, any morning (needs numpy once:
   `python3 -m pip install --user numpy`):

       python3 scripts/decide.py --plan research/stage1/plan-stage1.jsonl.xz --units-dir laptop-runs

   It compares the finished units with the predictions frozen in the plan
   (squares, S and P traversals, pairs, CPU), fits f_rho on the pairs of
   the finished band units (the shipped prior, and the flat prior on the
   band's pairs alone), and prints the predictive E and P(>=1 magic square)
   at 1 / 3 / 10 CPU-years and the decision table's rows
   (`research/stage1/decision-table.json`: the thresholds X are yours to
   set; the tool only reports). Here, with the scheduler state the plan was
   made from, `python3 scripts/decide.py STATE --plan ... --units-dir ...`
   also counts the band sums of units that start below N' 3k.
6. Send the results back:

       ./scripts/laptop/collect.sh
       git add research/laptop/results-*.tar.xz && git commit -m "laptop run results" && git push

## How long

A night of ~40 CPU-hours (e.g. 8 performance cores for 5 hours) covers
roughly 40 / (2.5-3.5) ~ 11-16 reference CPU-hours of the plan. The plan
predicts 45 (square, SP traversal) pairs in its band units after 51.4
reference hours (40.3 counting every sum in the band), so the 45 band pairs
take about 4 nights, and the whole plan (60 h) 4-5.5 nights; run.py
resumes. The first band unit comes after 10.8 reference hours, so the
first night finds none. With 45 pairs the 90% band of E at 1 CPU-year
narrows from x5.8 to about x3.3-3.8 (research/stage1.md).

## What the run is for

* The measurement: the rate of diagonals with both the magic sum and
  product (SP), against the model's prediction, at N' 3-6k (the plan's
  band units sit at N' 3.0-4.2k), where the forecast for a 1-10 CPU-year
  budget extrapolates it from small N. The predictions are frozen in the
  plan (pred_squares, pred_S, pred_P, pred_pairs, pred_time per unit) and
  registered in research/stage1/PREREGISTERED.txt before any unit ran.
* A real attempt: the plan's first units are the highest expected magic
  squares per CPU-second in the pool: about 0.005-0.006 expected magic
  squares the first night (11-16 reference hours), 0.014 over the whole
  plan.
* The laptop's speed: ARM has no AVX-512, and runs the same carried search
  as the fast x86 build with plain C kernels (src/c/arrange_carry.h; the
  same nodes and squares, which build.sh checks). On x86 without AVX-512
  these kernels take, against the fast build, with clang -march=x86-64-v3
  (256-bit vectors) 2.5x on bench quick, 2.55-2.6x on plain sums and
  1.86-1.97x on d-first ones; with 128-bit vectors, NEON's width
  (x86-64-v2), 3.2-3.7x; gcc v3 2.9-3.6x (two independent verifications;
  research/ideas.md, "The carried path without AVX-512"). So the M1 is
  expected around 2.5-3.5x per core on plain sums, to be measured by the
  bench/quick.txt time of step 3 and by run.py's CPU against reference
  CPU. The plan's --time-limit is 12x the scheduler's (a factor of ~6 would
  cover 3.5x with margin; 12 also covers a much slower core), and its
  --node-limit factor of 4 does not bind.

The portable ARM build was checked here by cross-compiling and running under
qemu: bench quick/full all ok, and three plan units give the same squares,
hashes, traversal counts and vector counts as the x86 build.

## Files

* `research/stage1/plan-stage1.jsonl.xz`: the plan to run;
  `research/stage1/PREREGISTERED.txt` its registration.
* `plan-20261010.jsonl.xz`, `PREREGISTERED.txt` (here): the superseded
  greedy plan (the stage-1 plan's first 15,164 units).
* `scripts/laptop/make_plan.py`: writes such a plan from a scheduler state
  (`--stage1` for a stage-1 plan).
* `scripts/laptop/build.sh`, `run.py`, `collect.sh`: build and check, run,
  pack the results. `scripts/decide.py`: the report.
