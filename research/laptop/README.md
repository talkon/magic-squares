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
   `git pull` in an existing clone), and check that
   `research/stage1/plan-stage1.jsonl.xz` is there.
3. Build and check (about 2 minutes):

       ./scripts/laptop/build.sh

   It compiles msearch, bench and fuzz_arrange into `build-laptop/`, checks
   every bench/quick.txt and bench/full.txt instance against the expected
   squares and hashes and the node counts of the fast x86 build, runs a
   brute-force differential test (also node for node), and prints the
   time of bench/quick.txt (divided by the fast x86 build's ~0.25-0.29 s it
   is the laptop's time per core against the reference, expected ~2.5-4.5x,
   see below). It also writes that time and the search kernels into
   `build-laptop/build_info.txt`, which the results carry (step 6). Stop if
   it says FAILED.
4. Run, plugged in, lid open (a closed lid sleeps the Mac):

       caffeinate -i python3 scripts/laptop/run.py research/stage1/plan-stage1.jsonl.xz --hours 9

   That is the plan file to pass. `--hours` is wall-clock hours. It uses
   the performance cores (`--workers K` to change; 4 on an M1, 6 or 8 on an
   M1 Pro / Max), runs the units in the plan's order, prints progress every
   minute (CPU-hours, the reference CPU-hours they stand for, the stage-1
   band units finished and their pairs) and announces any magic square
   loudly (also written to `laptop-runs/MAGIC.txt`). After 9 hours, or on
   Ctrl-C, or when the terminal is closed, it stops the running units,
   removes their partial output, prints a `finished:` line and exits (so
   caffeinate lets the Mac sleep again). The same command the next night
   continues where it stopped (finished units are skipped, an interrupted
   one is rerun from its start). The "x the reference per core" in the
   progress line is the laptop's speed mixed with any other load on it and
   with the predictions' own error (the fast build itself ran at 1.32x its
   predictions on a busy shared machine).
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

   The tarball holds the finished units, the meta files (machine, plan
   hash, build), progress.log (the speed) and the build's record
   (build_info.txt with the bench/quick.txt time and the search kernels,
   check_quick.txt). Nothing else needs to be sent.

## How long

The command above runs 9 hours on K performance cores, 9 K CPU-hours, which
at 2.5-4.5x per core is 9 K / (2.5-4.5) reference CPU-hours of the plan
(if K busy cores keep the speed of one; run.py's progress line shows it):

| performance cores | a 9-hour night | 45 band pairs (51.4 h) | the whole plan (60 h) |
| ---: | ---: | ---: | ---: |
| 4 (M1) | 8.0-14.4 reference h | 3.6-6.4 nights | 4.2-7.5 nights |
| 6 (8-core M1 Pro) | 12.0-21.6 h | 2.4-4.3 nights | 2.8-5.0 nights |
| 8 (10-core M1 Pro / Max) | 16.0-28.8 h | 1.8-3.2 nights | 2.1-3.8 nights |

The plan predicts 45 (square, SP traversal) pairs in its band units after
51.4 reference hours (decide.py without the state; with it, counting every
sum in the band, after 40.3 h). run.py resumes, so nights of any length add
up. The band units start at 10.8 reference hours: a first night of 8-29
reference hours reaches 0-21 of them (1.7 predicted pairs by 12 h, 2.9 by
16 h, 6.0 by 21.6 h, 11.0 by 28.8 h). With 45 pairs the 90% band of E at
1 CPU-year narrows from x5.8 to about x3.3-3.8 (research/stage1.md).

## What the run is for

* The measurement: the rate of diagonals with both the magic sum and
  product (SP), against the model's prediction, at N' 3-6k (the plan's
  band units sit at N' 3.0-4.2k), where the forecast for a 1-10 CPU-year
  budget extrapolates it from small N. The predictions are frozen in the
  plan (pred_squares, pred_S, pred_P, pred_pairs, pred_time per unit) and
  registered in research/stage1/PREREGISTERED.txt before any unit ran.
* A real attempt: the plan's first units are the highest expected magic
  squares per CPU-second in the pool: E[magic] 0.0109 over its first 40
  reference hours and 0.0137 over the whole plan (60 h).
* The laptop's speed: ARM has no AVX-512, and runs the same carried search
  as the fast x86 build with plain C kernels (src/c/arrange_carry.h; the
  same nodes and squares, which build.sh checks). On x86 without AVX-512
  these kernels take, against the fast build (two independent
  verifications, min of 2 rounds; research/ideas.md, "The carried path
  without AVX-512"):

  | build | bench quick | plain sums | d-first sums |
  | --- | ---: | ---: | ---: |
  | clang -march=x86-64-v3 (256-bit) | 2.5-2.6x | 2.55-3.0x | 1.6-2.0x |
  | gcc -march=x86-64-v3 | 3.1-3.4x | 3.5-4.2x | 2.3-3.1x |
  | gcc -march=x86-64-v2 (128-bit, NEON's width; one verification) | 3.6x | 4.5-4.6x | 2.5-3.4x |
  | clang -march=x86-64-v2 (the stage-1 rehearsal, busy machine) | 3.2x | 3.2x (49 plan units) | - |

  The M1 builds with clang at NEON's 128-bit width, so it is expected at
  ~2.5-4.5x per core on plain sums (~3x if it follows clang on x86), to be
  measured by the bench/quick.txt time of step 3 and by run.py's CPU
  against reference CPU. The plan's --time-limit is 12x the scheduler's
  (the scheduler's limit is at least 3x a unit's predicted CPU, so 12
  leaves 2.7x at 4.5x and also covers a much slower core), and its
  --node-limit factor of 4 does not bind.

To rehearse the laptop's build on an x86 machine with AVX-512, give the
target: `CC=clang ARCH=-march=x86-64-v2 ./scripts/laptop/build.sh` (without
ARCH, build.sh adds -march=native, which picks the AVX-512 kernels).

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
