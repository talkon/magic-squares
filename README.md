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
python3 scripts/scheduler.py plan              # what it would run next, and why
python3 scripts/scheduler.py forecast --hours 1000
python3 scripts/scheduler.py run --hours 24
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
| `scripts/scheduler.py` | chooses which (P, S) to search next to maximize the expected number of magic squares per CPU-second, runs `msearch` workers, refits its model as results come in |
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

    magic squares per CPU-second = (semi-magic squares per CPU-second)
                                   * 5400 * (rho * p_S * p_P)^2

where p_S, p_P are the probabilities that a traversal (a possible diagonal)
has the magic sum / product and rho ~ 2 corrects for these not being
independent. The first factor is a Poisson regression on the number of
vectors N, S / S_min, and the divisor structure of P, fitted to our runs and to
the per-P totals of the first search; p_S and p_P are fitted to the traversal
counts of the squares found. Each P also gets its own correction factors
(empirical Bayes) as its results come in, with an optimism bonus to explore
new P. Each P is searched in increasing order of S; the scheduler always runs
the next chunk of sums of the P with the best predicted yield.

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
| `bench/quick.txt` (7 instances, N = 450-1700) | 1.77M | 0.27 s | 0.29 s |
| `bench/full.txt` (9 instances, N = 330-2240) | 14.96M | 3.29 s | 3.47 s |
| `bench/prod.txt` (9 production-like sums, N = 1491-2994) | 50.38M | 13.4 s | 14.2 s |

The round-2 micro-optimizations took 8-9% off the native build and 12-15%
off the cascadelake build (from 0.30 / 3.65 / 14.6 s and 0.34 / 4.02 /
16.0 s; research/ideas.md, "Integration of round 2"). For comparison, the
legacy `arrangement_6` took 3.56 s on `quick.txt` and 63.1 s on `full.txt`
(earlier measurement, Intel Xeon with AVX-512), visiting 31.7M and 598M
nodes.

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
| fdb77fc (current) | 1/17-1/18 | 0.16 | 0.64 |

The 17-18x is measured speed; E at a fixed budget grows only ~3x because
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
