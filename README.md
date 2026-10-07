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
`python3 scripts/test_scheduler.py` for the scheduler).

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

Arrangement benchmark (single core, Intel Xeon with AVX-512):

| | legacy `arrangement_6` | new (AVX-512) | new (portable, AVX2) |
| --- | ---: | ---: | ---: |
| `bench/quick.txt` (7 instances, N = 450-1700) | 3.56 s | 1.11 s | 1.85 s |
| `bench/full.txt` (9 instances, N = 330-2240) | 63.1 s | 18.0 s | |

The new search visits 4-5x fewer nodes than the legacy search (7.5M vs 31.7M on
`quick.txt`, 126M vs 598M on `full.txt`), at a somewhat higher cost per node.

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
