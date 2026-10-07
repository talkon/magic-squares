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
`--no-mrv` / `--no-support` to disable the search heuristics.

Output is JSON lines, flushed as it goes: one `square` record per semi-magic
square (with its traversal counts, the best pair of diagonals, and the square
reordered so that pair is on the diagonals), one `sum` record per searched sum
(number of vectors, labels, nodes, squares, time), and a final `done` record
with the last sum searched. A killed run keeps everything up to the last
completed sum. A best diagonal score of 14 is a magic square, and is also
announced on stderr.

The vectors are enumerated in-process (`enumerate.c`), only for the requested
sums, so there are no intermediate files and memory stays small.

## bench: arrangement benchmark

```
bin/bench bench/quick.txt              # ~1 s
bin/bench --repeat 3 bench/full.txt    # ~20 s, closer to production sizes
bin/bench --only "S=648" bench/full.txt
bin/bench --no-fc --no-mrv bench/quick.txt   # compare search variants
bin/bench --no-support bench/quick.txt
bin/bench --min-words 8 --gather bench/quick.txt   # test the paths for more labels / larger N
bin/bench --update bench/quick.txt     # print instances with observed values
```

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
- forward checking: a node is pruned as soon as some unmatched cell has no
  remaining candidate through it;
- support filter: every cell of a candidate column must lie in a placed row
  or in some candidate row (the row of the square through it), and vice
  versa, so candidates using a number that no candidate of the other axis
  covers are dropped, repeatedly until both lists are closed. This prunes most
  nodes with two rows and two columns placed, which forward checking lets
  through: 4x fewer nodes and ~30% less time on `bench/full.txt`;
- branching on the unmatched cell with the fewest candidates, ties broken by
  the largest label; the counts are exact (saturating byte counters, one
  masked add per 64 labels) with AVX-512BW and up to 256 labels, and
  bit-sliced otherwise (exact up to 7 with AVX-512, 3 in plain C);
- the counts are computed while filtering the candidate lists, filtering first
  the axis whose cells usually lose all candidates, so most dead children are
  discarded after one pass; from depth 4 on, where most children are pruned
  by the support filter, only the unions are computed while filtering, and
  the full counts only for the nodes that are searched;
- the children of a node (the candidates through the branching cell) are
  selected with the same filter as the candidate lists, from a label ->
  vectors bit matrix, rather than tested one by one;
- with AVX-512, the candidate lists are filtered 16 at a time against a row
  of a bit matrix, looking the bits up with permutes from the row held in
  registers when N <= 3072 (gathers otherwise), and compressed; vectorized
  counting; plain C otherwise.

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
