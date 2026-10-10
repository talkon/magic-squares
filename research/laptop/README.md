# Running the search on a laptop (macOS / Apple Silicon, or any Linux)

A static, pre-registered plan of the units scheduler v2 would launch first
(`plan-20261010.jsonl.xz`, 15,164 units on 10,276 P, plain search, 40
predicted CPU-hours of the fast x86 build, E[magic] 0.011), and the tools to
run it anywhere without the scheduler or numpy.

## Steps (macOS)

1. Command line tools, once: `xcode-select --install` (gives cc, git,
   python3). Check that `uname -m` prints `arm64`.
2. Get the branch: `git clone -b claude/magic-search-speedups <repo>` (or
   `git pull` in an existing clone).
3. Build and check (about 2 minutes):

       ./scripts/laptop/build.sh

   It compiles msearch, bench and fuzz_arrange into `build-laptop/`, checks
   every bench/quick.txt and bench/full.txt instance against the expected
   squares and hashes, runs a brute-force differential test, and prints the
   time of bench/quick.txt (please report it: it measures the laptop's speed
   against the fast x86 build's ~0.29 s). Stop if it says FAILED.
4. Run, plugged in, lid open (a closed lid sleeps the Mac):

       caffeinate -i python3 scripts/laptop/run.py research/laptop/plan-20261010.jsonl.xz --hours 9

   It uses the performance cores (`--workers K` to change), runs the units in
   the plan's order (best expected magic squares per CPU first), prints
   progress every minute and announces any magic square loudly (also written
   to `laptop-runs/MAGIC.txt`). Ctrl-C stops it cleanly; running the same
   command again continues where it stopped.
5. Send the results back:

       ./scripts/laptop/collect.sh
       git add research/laptop/results-*.tar.xz && git commit -m "laptop run results" && git push

## What the run is for

* A real attempt: the plan's first units are the highest expected magic
  squares per CPU-second in the pool (about 0.005 expected magic squares
  per night on an M1 Pro).
* Calibration on exactly the sums a paid run would start with: the rate of
  diagonals with both the magic sum and product (SP), the quantity the
  forecast's main remaining uncertainty rests on, plus squares and S/P
  traversals per square. The predictions are frozen in the plan
  (pred_squares, pred_magic, pred_time per unit) and in PREREGISTERED.txt.
* The laptop's speed: ARM has no AVX-512, and now runs the same carried
  search as the fast x86 build with plain C kernels (src/c/arrange_carry.h;
  the same nodes and squares, which build.sh checks), on x86 without
  AVX-512 (clang, AVX2) 2.3-2.4x the fast build's time per core, against
  4.1-4.6x for the matrix path these builds had before (research/ideas.md,
  "The carried path without AVX-512"; to be measured on the M1). The
  plan's --node-limit factor of 4 was for the matrix path's extra nodes;
  with the same nodes it does not bind.

The portable ARM build was checked here by cross-compiling and running under
qemu: bench quick/full all ok, and three plan units give the same squares,
hashes, traversal counts and vector counts as the x86 build.

## Files

* `scripts/laptop/make_plan.py`: writes such a plan from a scheduler state.
* `scripts/laptop/build.sh`, `run.py`, `collect.sh`: build and check, run,
  pack the results.
