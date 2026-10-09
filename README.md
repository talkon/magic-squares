# magic-squares

Search for 6x6 (and 5x5, 7x7) additive-multiplicative magic squares: n x n
grids of distinct positive integers whose rows, columns and both diagonals all
have the same sum S and the same product P. See `project-report.pdf` for the
background and the results of the first search (758,949 semi-magic squares,
no magic square).

The search has three steps: **enumerate** all sets of n distinct numbers with
sum S and product P ("vectors"), **arrange** them into semi-magic squares
(rows and columns only), and check whether any square has two diagonals that
can be made magic by permuting its rows and columns.

## Quick start

```
./build.sh                          # build (cmake-build-release/) and link into bin/
bin/bench bench/quick.txt           # ~1 s: arrangement benchmark + correctness check
bin/msearch --min-sum 320 --max-sum 340 10 4 3 2    # P = 2^10 3^4 5^3 7^2

# adaptive search over (P, S), using all cores, resumable (state in data/sched)
python3 scripts/scheduler.py import-legacy stats/stats_short.txt
python3 scripts/scheduler.py profile           # once: the model of 378k P (~1 CPU-hour, niced)
python3 scripts/scheduler.py plan              # what it would run next, and why
python3 scripts/scheduler.py forecast --hours 1000   # E with and without d-first
python3 scripts/scheduler.py run --hours 24    # plain or d-first per sum (--dfirst auto|off|on)
python3 scripts/scheduler.py report
```

`./build.sh t` runs the tests (`ctest -R fast_` for just the fast ones, and
`python3 scripts/test_scheduler.py` for the scheduler, `python3
scripts/test_calibrate.py` for the calibration ladder).

## Programs

| program | what it does |
| --- | --- |
| `bin/msearch` | enumerate + arrange + diagonal check for one P and a range of sums, in one process, writing JSON lines as it goes (or the legacy text format) |
| `bin/bench` | arrangement benchmark on fixed (P, S) instances, checking the squares found against expected counts and hashes |
| `bin/enumerate` | write enumeration files (drop-in replacement for `enumeration.cpp`), vector counts per sum, or the smallest sum |
| `scripts/scheduler.py` | chooses which (P, S) to search next to maximize the expected number of magic squares per CPU-second, runs `msearch` workers, refits its model as results come in (model: `scripts/amodel.py`, candidate pool: `scripts/pool.py`) |
| `submit-sc-plan.sh` | runs a static plan from `scheduler.py emit` on SuperCloud with LLsub |
| `bin/arrangement_{5,6,7}`, `bin/enumeration`, `src/py/`, `run.sh`, `submit-sc*.sh` | the original pipeline (still works) |

Details are in `src/c/README.md` (C programs) and in the docstring of
`scripts/scheduler.py`.

## What changed compared to the report

**Enumeration** (`src/c/enumerate.c`) is a depth-first search over the
divisors of P with bounds on the remaining sum, instead of materializing all
multisets of divisors with product P. Memory is proportional to the output
rather than to prod C(e_i + 5, 5), and it can be restricted to a range of sums,
which is all the scheduler ever needs. For example, all vectors of
P = 13 6 3 2 with S <= 1000 take 1 s.

**Pipelining**: `msearch` enumerates the sums it needs in-process, so there
are no enumeration files, and writes one record per sum and per square as soon
as they are done, so a killed job loses at most one sum.

**Arrangement** (`src/c/arrange.c`) finds exactly the same squares as before
(checked against the legacy program on all benchmark instances, and on all 40
squares of P = 13 6 3 2 with S <= 700), ~3x faster on the benchmark sets (see
below). The main changes are forward checking (prune as soon as some
unmatched number has no remaining candidate row/column), branching on the
most constrained number instead of the largest, label bitsets for all the
bookkeeping, and vectorized filtering/counting. It no longer requires
AVX-512.

**Scheduling** (`scripts/scheduler.py`): semi-magic squares are found much
faster near the smallest possible sum S_min(P) (about 600 to 10,000 per
CPU-hour depending on P, versus ~80 per hour on average in the first search),
and squares with a smaller S are also more likely to be magic. The scheduler models

    magic squares per CPU-second = (semi-magic squares per sum) / (CPU-seconds per sum)
                                   * P(a semi-magic square is magic)

per (P, S) with the analytic model of research/existence.md
(`scripts/amodel.py`: a max-entropy local limit theorem with exact lattice
factors and ~3 fitted constants, which predicts the number of vectors N, the
squares per sum and P(magic | square) = 5400 kappa p_pair, and was tested on
held-out data up to N = 45k), times the review corrections of that study, and
a time law in N and the predicted number of distinct entries. The candidates
are every exponent assignment of P over 2..29 within a small factor of the
sorted one (`scripts/pool.py`, 377,908 P: most of the yield is in exponents
that are not non-increasing, e.g. 13 7 4 3 0 0 1 1). The model of each P is
computed once on 24 sums; then a lazy greedy planner always runs the next
chunk of sums of the P with the best predicted yield, searching each P in
increasing order of S from ceil(6 P^(1/6)). As results come in it refits
class factors (squares, S and P traversals, by N band, number of primes,
assignment ratio and S / S_min; quasi-Poisson, since counts of one P move
together), the level of the time law (its shape stays fixed; msearch
reports CPU time per sum), and per-P factors (empirical Bayes); `report`
compares observed and predicted counts. Each sum is searched plain or
diagonal-first (`msearch --diag-first`, which finds every magic square but
not the other semi-magic squares), by the measured d-first / plain CPU
ratio (`--dfirst auto`, the default: d-first from N' ~ 4.3k on; the ratio's
level is learned online). A d-first sum gets a plain calibration stream of
~7% of its CPU (`--calib-r1-stride`), whose sampled squares keep the class
and per-P factors learning and whose plain-time estimate teaches the plain
law and the ratio, and a sum of hours is split into units of d
(`--d-range`) that are merged and resumed from msearch's checkpoint
records (research/scheduler-v2.md, "d-first units"). The previous
scheduler (a Poisson regression fitted to our runs, over non-increasing
exponents) is still available as `--model regression --pool classic`; it
runs the plain search only.

### Running on a cluster (SuperCloud)

The scheduler can also write a static plan, to run as a job array:

```
python3 scripts/scheduler.py import-legacy stats/stats_short.txt   # once
python3 scripts/scheduler.py emit --units 4800 --unit-time 600 > plan.txt
PLAN=plan.txt LLsub ./submit-sc-plan.sh [8,48,1]     # LLsub triples mode
```

Each line of `plan.txt` is one `msearch` run writing to
`data/sched/units/`. Afterwards `scheduler.py report` summarizes the results,
`scheduler.py fit` refits the model on everything found so far, and the next
`emit` continues where the plan left off.

## Performance

Arrangement benchmark (`bin/bench`, single core, search time; Intel Sapphire
Rapids, built with `-march=native` and with `-march=cascadelake` as for the
production cluster, both run on the same machine; minimum of alternating
runs, October 2026):

| | nodes | -march=native | -march=cascadelake |
| --- | ---: | ---: | ---: |
| `bench/quick.txt` (7 instances, N = 450-1700) | 1.77M | 0.27 s | 0.30 s |
| `bench/full.txt` (9 instances, N = 330-2240) | 14.96M | 3.27 s | 3.52 s |
| `bench/prod.txt` (9 production-like sums, N = 1491-2994) | 50.38M | 14.0 s | 14.8 s |

The round-2 micro-optimizations took 8-9% off the native build and 12-15%
off the cascadelake build (from 0.30 / 3.65 / 14.6 s and 0.34 / 4.02 /
16.0 s; research/ideas.md, "Integration of round 2"). For comparison, the
legacy `arrangement_6` took 3.56 s on `quick.txt` and 63.1 s on `full.txt`
(earlier measurement, Intel Xeon with AVX-512), visiting 31.7M and 598M
nodes.

The current code (cx/integrated) visits the same nodes and finds the same
squares on all three files. Its times above were measured on October 8
alternating with fdb77fc's, built the same way, on a shared machine (min of
3 rounds: fdb77fc 0.27 / 3.39 / 13.8 s native, 0.29 / 3.63 / 14.3 s
cascadelake; the mean over the rounds differs by less than 1% on prod):
the benchmark sums have at most 121 labels, below everything the
integration changed, except the pretest (0-3% here). Its gains are at
large N.

### Large sums

The benchmark sums have at most 2,994 vectors and 121 labels. The sums that
matter most for a first magic square have N ~ 3k-32k vectors (see Status),
where three changes of October 2026 (branch cx/integrated; research/ideas.md,
"Carried bitsets up to 512 labels", "The pretest of the deep children",
"Diagonal-first search in msearch" and "Measurements on the integrated
binary") speed the search up:

* the plain search runs each first row with only the 64-bit words of labels
  it needs, carries its bitsets up to 512 labels instead of switching to
  N x N matrices above 256, and pretests the deep children;
* `msearch --diag-first` searches, for every vector d of the sum, the
  squares that have d as a diagonal (the semi-magic search on the vectors
  that meet d once): it finds every magic square (twice) without
  enumerating the semi-magic squares, so it is used only where it is
  cheaper, N >= 5000 (`--diag-first-min-n`).

CPU seconds per sum (paired, r1- and d-sampled where large; relSE <= 8%;
the previous code is fdb77fc):

| P / S | N | labels | previous plain | plain | d-first | plain speedup | best speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 6 4 2 1 1 / 838 | 1,986 | 104 | 1.01 | 1.01 | 1.65 | 1.0 | 1.0 |
| 12 6 3 2 1 0 1 / 900 | 2,994 | 121 | 5.9 | 5.8 | 5.9 | 1.0 | 1.0 |
| 13 7 4 3 1 1 / 1900 | 4,111 | 137 | 35.5 | 23.5 | 29.4 | 1.5 | 1.5 |
| 12 6 3 2 1 1 / 988 | 6,671 | 159 | 260 | 231 | 149 | 1.1 | 1.7 |
| 13 7 4 3 1 1 / 2000 | 7,593 | 168 | 511 | 441 | 251 | 1.2 | 2.0 |
| 12 6 3 2 1 1 / 1200 | 11,697 | 199 | 2,023 | 1,575 | 870 | 1.3 | 2.3 |
| 13 7 4 3 1 1 / 2200 | 15,199 | 211 | 11,316 | 8,223 | 3,432 | 1.4 | 3.3 |
| 11 6 4 3 2 1 / 2174 | 16,424 | 220 | 11,503 | 8,685 | 4,056 | 1.3 | 2.8 |
| 14 7 4 4 1 0 0 1 / 3648 | 20,538 | 252 | 14,460 | 12,623 | 5,407 | 1.1 | 2.7 |
| 9 6 4 3 1 1 1 1 / 2700 | 20,896 | 259 | 53,077 | 13,835 | 4,457 | 3.8 | 11.9 |
| 13 7 4 3 1 1 / 2400 | 22,992 | 245 | 47,400 | 40,389 | 11,785 | 1.2 | 4.0 |
| 13 7 4 3 1 1 / 2500 | 26,585 | 260 | 205,729 | 54,280 | 23,366 | 3.8 | 8.8 |
| 13 7 4 3 1 1 / 2650 | 31,743 | 279 | 467,238 | 121,918 | 39,917 | 3.8 | 11.7 |

* The plain search is unchanged up to 128 labels, 1.1-1.5x faster at
  129-256 labels and 3.8x faster above 256.
* d-first costs 1.25x the plain search at N = 4.1k, 0.57-0.65x at 6.7-7.6k,
  0.42-0.55x at 11.7-21k and 0.29-0.43x at 21-32k. It crosses over at
  N ~ 4-5k (3.8k pooled over the 13 sums, 4.9k along 13 7 4 3 1 1; it
  depends on P), hence the threshold of 5000, which also keeps the plain
  search's semi-magic squares, which the models are fitted to, below it.
* Time per sum, t = a (N / 4000)^b CPU-s fitted at N = 4.1-31.7k: the
  previous plain search a = 24.8, b = 4.5 (t(32k) ~ 280k CPU-s); now, the
  better of the two modes, a = 23, b = 3.5 (t(32k) ~ 35k).
* Scheduler v2 chooses the mode per sum by the measured ratio and runs
  large d-first sums as units of d (research/scheduler-v2.md, "d-first
  units").

## Status (October 2026)

In a few CPU-hours of testing near S_min (36 values of P searched from S_min
for 4 minutes each, then 40 units chosen by the scheduler), the new pipeline
found 5,821 semi-magic squares, including 12 SP-type, one SP+S-type
(P = 2^12 3^6 5^3 7^2 11 17, S = 836) and one SP+P-type square
(P = 2^16 3^5 5^4 7^2, S = 849; the first search found one SP+P square in
total):

```
 324  120    7  128  200   70
  25   64  216  420   40   84
 300  147   50   48   16  288
 112    8  180  175   54  320
  32  360  140   18  224   75
  56  150  256   60  315   12
```

Estimating the chance that each square is magic from the observed rates of
S- and P-traversals, these runs found magic-square "probability mass" at
~2.5-3e-5 per CPU-hour, versus ~4e-7 per CPU-hour for the first search (the
same estimate applied to its 758,949 squares), i.e. ~60x more per CPU-hour,
from the faster search and from staying close to S_min.

That pace cannot be kept up: every P has only a short stretch of good sums,
so the rate falls as the best ones are used up. The figure quoted here
earlier ("40-65 CPU-years per magic square") was the starting rate, not the
time to a first magic square.

[research/existence.md](research/existence.md) estimates how many 6x6
magic squares there should be at all. It uses an analytic heuristic (exact
exponent-matrix counts and a local limit theorem), calibrated on the search
data and tested on data it was not fitted to, up to N = 45k vectors per sum.

* A 6x6 magic square very probably exists (~0.95, if nothing structural
  forbids it, as for n = 3 and 4). The smallest magic sum is ~2,800
  (10-90%: ~1,500 to 6,000-10,000), at P with 5-15k divisors and sums with
  N ~ 3k-32k vectors. The expected total is finite, ~4.6 (68%: 2.2-13).
* Expected magic squares found with the current build and an ideal choice of
  (P, S): **0.16 in 1 CPU-year, 0.34 in 10, 0.64 in 100, 1.1 in 1,000**
  (68% band x/÷2). Everything searched so far holds ~0.005.
* The scheduler as fitted forecasts only 0.015-0.018 for itself in 1
  CPU-year ([research/forecast.md](research/forecast.md)). Scored by the
  analytic model, its choices are worth ~0.05, about 3x less than the ideal
  choice at 1 CPU-year and 4x at 100. Its model of squares per sum stops
  growing at N ~ 5,500 (it predicts 88 squares where 484 were found at
  N = 5-8k, and 0.08 where 60 were found at N = 15-32k), and its pool omits
  exponent orders that are not non-increasing, which hold more than half of
  the yield.
* The analytic scheduler ([research/scheduler-v2.md](research/scheduler-v2.md),
  now the default) uses that model per (P, S) over a pool of 377,908 P
  that includes exponent orders that are not non-increasing. It forecasts
  **0.12 magic squares in 1 CPU-year, 0.25 in 10 and 0.39 in 100**
  (`forecast --shipped`; 90% band from its calibration 0.04-0.32 at 1
  CPU-year). Units it ranks highest came in at 0.7-0.9 of prediction on
  held-out data, so quote **about 0.1 and 0.2**. Scored by the same model,
  the old scheduler's choices collect 0.059 and 0.096 (x2.0 / x2.6). The
  ideal frontier above uses a slower time law and is not directly
  comparable (v2 reaches ~2/3 of it under that law).
* Speedups are best compared by the ratio of E at a fixed budget: 10x faster
  gives x2.2 at 1 CPU-year and x1.5 at 1,000; time ~ N^2 instead of ~N^4
  beyond N = 4,000 gives x1.1 and x1.7.

[research/retrospective.md](research/retrospective.md) puts every version of
the code through the same model, with times measured for 9 versions on 33
sums (N = 451-31,743):

| code | CPU to reach a given E, vs the first search's code | E(1 CPU-yr) | E(100 CPU-yr) |
|---|---:|---:|---:|
| first search (legacy `arrangement_6`) | 1 | 0.051 | 0.28 |
| df352df (new arrange.c) | 1/2.5-1/3 | 0.080 | 0.38 |
| a72cef3 (+ support filter) | 1/6.5 | 0.11 | 0.48 |
| 2d2bc6d (+ carried bitsets) | 1/11-1/13 | 0.14 | 0.58 |
| fdb77fc | 1/17-1/18 | 0.16 | 0.64 |
| current (engine 3), plain and d-first | 1/21 (E = 0.1) to 1/43 (E = 1) | 0.17 | 0.76 |

The current code (msearch engine 3: the "Large sums" table above) is
1.1-3.8x faster than fdb77fc at N >= 4k with the plain search and 1.5-12x
with the better of the plain and d-first modes. Its row assumes plain below
N = 5,000, d-first above, and a calibration stream costing 10% of each
d-first sum (retrospective.md section 7). It lowers the exponent of N
(t ~ N^3.5 against N^4.5), so the gain grows with the budget: E x1.11 at 1
CPU-year, x1.24 at 1,000, and 2.3x less CPU to reach E = 1. Scheduler v2
now launches d-first units, choosing the mode by the measured ratio
(d-first from N' ~ 4.3k). Charged under a measurement-anchored truth
(`forecast --truth anchored`), d-first adds x1.02 / 1.07 to its plan at 1
/ 10 CPU-years and needs 1.24x less CPU for the plain-only E of 10
CPU-years (the review's simulation: x1.14 and 1.78x at 100), close to the
retrospective's x1.03 / 1.07 / 1.13 for d-first on top of the plain search
(research/scheduler-v2.md, "d-first units"). The laws also have to
learn engine 3's level first: re-run `forecast` after some engine-3 units. The
17-18x is measured speed; E at a fixed budget grows only ~3x because
E(C) rises slowly with C. Search strategy multiplies with it: the first
search's expected yield (1.1 CPU-years with the legacy code) is reached in
about 6.5 CPU-hours with the current code and an ideal choice of (P, S).

[research/calibration.md](research/calibration.md) checks the per-square
part of the heuristic rung by rung (S, P, S+S, S+P, P+P, SP) on the 7,021
squares found so far: it lines up within errors.

## Directory structure

```
magic-squares
|-- bench/     -> arrangement benchmark instances (bin/bench) and legacy cross-check
|-- scripts/   -> scheduler for the search over (P, S)
|-- src/c/     -> enumeration, arrangement, msearch, bench (C), legacy arrangement
|-- src/py/    -> original python implementation and postprocessing
|-- stats/     -> summaries of the first search
|-- runs/      -> inputs of the first search's SuperCloud runs
|-- research/  -> notes
|-- data/      -> search outputs (not in git)
```

## Usage of the original pipeline

Release: `./build.sh`.

Debug: `./build.sh d`, useful for `perf` and `gdb` etc.

After building, you can do `./run.sh` to run the whole pipeline.

Individual usage is in README.md of each directory.
