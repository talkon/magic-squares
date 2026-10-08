Build with CMake from the root directory (`./build.sh`), which also links the
programs into `bin/`.

## msearch: enumerate + arrange + diagonals for one P

```
bin/msearch [options] e1 e2 ...        # P = 2^e1 3^e2 5^e3 ...
bin/msearch 10 4 3 2                   # the first 100 sums from the smallest possible
bin/msearch --min-sum 320 --max-sum 330 10 4 3 2
bin/msearch --sums 327,350 10 4 3 2
bin/msearch --vec-size 5 --time-limit 600 --out out.jsonl 10 6 3 2 1
bin/msearch --format legacy --sums 327 10 4 3 2 > out.txt   # for postprocess.py
```

Options: `--vec-size N` (5, 6, 7), `--min-sum`, `--max-sum`, `--sums`,
`--node-limit X` (per sum; the record is marked `truncated`), `--total-nodes X`
and `--time-limit T` (stop after the sum during which the limit is reached),
`--out FILE` (append), `--format legacy`, `--reduce strong`, and `--no-fc` /
`--no-mrv` / `--no-support` / `--no-cross` to disable the search heuristics.

Output is JSON lines, flushed as it goes: one `square` record per semi-magic
square (with its traversal counts, the best pair of diagonals, and the square
reordered so that pair is on the diagonals), one `sum` record per searched sum
(number of vectors, labels, nodes, squares, time), and a final `done` record
with the last sum searched. A killed run keeps everything up to the last
completed sum. A best diagonal score of 14 is a magic square, and is also
announced on stderr.

The vectors are enumerated in-process (`enumerate.c`), only for the requested
sums, so there are no intermediate files and memory stays small. Everything
but the search itself (enumeration, reduction, relabelling, output) takes
~0.3% of a unit near S_min (N up to ~1700), and ~3% of a short unit whose
sums all have N < 1000 (see research/ideas.md, "msearch overhead").

## bench: arrangement benchmark

```
bin/bench bench/quick.txt              # ~1 s
bin/bench --repeat 3 bench/full.txt    # ~5 s, closer to production sizes
bin/bench bench/prod.txt               # ~15 s, sampled from scheduler runs
bin/bench --only "S=648" bench/full.txt
bin/bench --no-fc --no-mrv bench/quick.txt   # compare search variants
bin/bench --no-support bench/quick.txt
bin/bench --no-cross bench/quick.txt
bin/bench --min-words 6 bench/quick.txt   # test the paths for more labels (carried up to 8 words)
bin/bench --no-r1-width bench/full.txt # every r1 at the width of all labels (same nodes, slower)
bin/bench_matrices bench/quick.txt     # the same with the N x N matrices (CARRY_MAX_W=0)
bin/bench_matrices --min-words 8 --gather bench/quick.txt   # matrices as for > 512 labels, N > 3072
bin/bench --node-limit 3000000 big.txt # time the first 3M nodes of larger instances
bin/bench --update bench/quick.txt     # print instances with observed values
```

Built with `-DCHILD_PROF` (e.g. `gcc -O3 -march=native -DCHILD_PROF -o
bench_prof bench.c arrange.c square.c enumerate.c -lm`), the search prints
on exit where the creation of children spends its cycles, by layer (rows
and cols placed) and phase, and which test kills them.

Each instance line is `n | exponents | S | expected squares | expected hash`.
The hash is an order-independent hash of the squares found (sum of hashes of
their canonical forms), so it checks that exactly the right squares are found,
not just how many. The expectations in `bench/*.txt` were checked against the
legacy `arrangement_*` programs with `bench/validate_legacy.py`.

## enumerate: write enumeration files

Drop-in replacement for `enumeration.cpp` / `enumeration.py` (byte-identical
output), using a depth-first search over divisors instead of materializing
every multiset of divisors with product P:

```
bin/enumerate --file out.txt 10 4 3 2
bin/enumerate --max-sum 1000 --file out.txt 13 6 3 2
bin/enumerate --counts --reduce none --min-sum 500 --max-sum 600 --file counts.txt 13 6 3 2
bin/enumerate --print-min-sum 13 6 3 2
```

## Search algorithm (arrange.c)

Same overall approach as the legacy search (the first row contains the largest
label of the square; every step adds a row or column through an "unmatched"
cell, i.e. a number that is in a placed row but no placed column, or vice
versa), plus:

- vectors are bitsets over labels (only ~50-200 distinct numbers per (P, S)),
  so the placed cells, unmatched cells and candidate counts are a few bitwise
  operations;
- the search is split by the largest label x of the square (r1 and c1 go
  through x, every other vector has only labels < x), and the labels are a
  degeneracy order: the largest label goes to the rarest number, the next one
  to the rarest number among the vectors left without it, and so on, so that
  each x is rare in its own part of the search (3% fewer nodes than labels by
  plain frequency on `bench/full.txt`, 6% on `bench/quick.txt`);
- forward checking: a node is pruned as soon as some unmatched cell has no
  remaining candidate through it;
- support filter: every cell of a candidate column must lie in a placed row
  or in some candidate row (the row of the square through it), and vice
  versa, so candidates using a number that no candidate of the other axis
  covers are dropped, repeatedly until both lists are closed. This prunes most
  nodes with two rows and two columns placed, which forward checking lets
  through: 4x fewer nodes on `bench/full.txt`, and half the time with the
  carried bitsets below (30% less with the matrices);
- cross support (n <= 6, when the children with two rows and two columns
  placed are created): the column through an unmatched cell y of a placed
  row is one of the candidate columns through y, and every remaining row
  meets it, so a candidate row that misses the union of the candidate
  columns through some row-unmatched cell is dropped, and vice versa. One
  pass per axis, with the support filter's test folded in, before the
  support filter: it kills most of these children, and nearly all of the
  ones that the support filter let through and that were searched (with a
  few children each, all dead). Where the support filter alone kills
  nearly all of these children, cross support is cheaper after it, on the
  few it lets through, so the search estimates that kill rate as it goes
  (on one child in 16) and switches. 41% fewer nodes and ~10% less time on
  `bench/full.txt`, 61% fewer nodes and 21-23% less time on the
  production-like sums of `bench/prod.txt` (built with -march=native or
  -march=cascadelake; 24-27% on the three with N = 2161-2994). With the
  carried bitsets, a pass builds the unions 8 cells at a time (the word
  and bit of each cell from a pdep per word, without a branch; 8 entries
  of the other list: a test per cell and a masked or per word, then one
  transposing reduction) and tests 8 candidates at a time against all of
  them; with the matrices it runs
  only with AVX-512 and byte counters (39% fewer nodes, 22% less time on
  `bench/full.txt`), elsewhere it cost more than it saved;
- branching on the unmatched cell with the fewest candidates, ties broken by
  the smallest label; the counts are exact (saturating byte counters) with
  AVX-512BW and up to 512 labels (256 on the matrix path), where the cell is
  found with one masked
  minimum over the counters instead of a scan of the counts class by
  class (0-3.5% less time on `bench/prod.txt`, depending on the rest of
  the code; `fuzz_arrange --mode 6` tests the saturated counts, >= 255
  candidates through every unmatched cell), and
  bit-sliced otherwise (exact up to 7 with AVX-512, 3 in plain C);
- the axis whose cells usually lose all candidates is filtered first, so
  most dead children are discarded after one pass; the forward check and the
  support filter only need the union of each candidate list, so the full
  counts are computed only for the children that survive them (with the
  N x N matrices below, from depth 4 on: above it, counting while filtering
  is cheaper);
- the children of a node (the candidates through the branching cell) are
  selected with a vectorized filter (from the bitsets the lists carry, or
  from a label -> vectors bit matrix, see below), rather than tested one by
  one;
- with AVX-512BW and up to 512 labels, the candidate lists carry the label
  bitsets of their entries (one array per 64-bit word), so filtering a list
  by the vector v just placed computes |u & v| for 8 entries at a time from
  contiguous loads (vpopcntq), and compresses the kept entries with their
  bitsets. The support filter is a test of the same carried words against
  the uncovered numbers, 8 entries at a time, and its first pass over the
  axis disjoint from v is fused into that axis' filter (one test against
  v | uncovered). The labels are counted (8 entries at a time: byte transpose
  with vpermb, bit transpose with gf2p8affineqb, vpopcntb) only for the
  children that survive all the filters. No N x N intersection matrices:
  setup is 4-8x faster and memory is O(N) instead of O(N^2) (7.6 MB of
  matrices at N = 5500). With the support filter this halves the search
  time (`bench/full.txt` 8.2 s -> 4.4 s on Sapphire Rapids). AVX-512 CPUs
  without VPOPCNTDQ or VBMI/GFNI/BITALG (Skylake-X, Cascade Lake) use
  fallbacks, and are also much faster than with the matrices (built with
  -march=cascadelake, run on Sapphire Rapids: 8.8 s -> 5.2 s); build with
  -DCARRY_MAX_W=0 to use the matrices instead. The filter that creates the
  children's lists needs no lane masks: each node pads its lists with 8
  full sets, which the filter drops, before creating its children; the
  label counts and the selection of the children run on the full groups
  of 8 entries without a lane mask, then on the last group (together ~5%
  less time than a lane mask in every iteration). |u & v| = 1 is tested
  on one word: the words of u & v rotated into one, with rotations that
  keep the words of v disjoint (one vpopcntq, or without VPOPCNTDQ one
  exactly-one-bit test instead of a test per word and across words: 5%
  less time with -march=cascadelake);
- each first row r1 is searched with only the 64-bit words it needs: its
  subproblem has the vectors after it, whose labels are all at most its
  largest label x, so the r1 are searched in runs of equal width
  ceil((x + 1) / 64) (at least 2), each run filling the depth-0 lists with
  that many words (the same nodes as at a single width; `--no-r1-width`
  turns it off). Each word costs ~x1.3 time at the same nodes, so this is
  1.2-1.3x less time where the r1 with the most work need one word less
  than all the labels (13 7 4 3 1 1, S = 2200, N = 15k: 1.23x; 12 6 3 2 1 1,
  S = 1200: 1.31x), and nothing where they need them all. Above 256 labels
  it replaces the matrices, which also lacked cross support and exact MRV:
  2.8-3.1x less time at N = 21-32k with 259-279 labels (5 words only for
  the first r1), 2.2x inside the d-first V_d searches (research/ideas.md,
  "Carried bitsets up to 512 labels");
- otherwise (more than 512 labels, or no AVX-512BW), the candidate lists are
  indices filtered with N x N intersection bit matrices: with AVX-512, 16 at
  a time against a row of a matrix, looking the bits up with permutes from
  the row held in registers when N <= 3072 (gathers otherwise), and
  compressed; vectorized counting; plain C otherwise.

## Legacy arrangement program

```
arrangement_6 --file {file}
arrangement_6 --file {file} --sum 327
arrangement_6 --file {file} --min-sum 320 --max-sum 330 --count 100000000000
```

short arguments:
- `-f: --file`
- `-s: --sum`
- `-n: --min`
- `-x: --max`

It requires AVX-512, so it is only built where the compiler supports it.
