# Ideas

General approach: for each $(P, S)$, **enumerate** all possible "rows" of $n = 6$ numbers with product $P$ and sum $S$, and then search through ways to **arrange** these rows/columns into a magic square.

## Enumeration

Enumerating all possible rows with product $P$ can be done in a reasonable amount of time. 

## Arrangement

* Find *semi-magic* squares first, because they have more symmetry, and then hope that with enough of them we can find one that is also magic (possibly after reordering rows/columns).
* Frame this as a subgraph isomorphism problem: let $R$ be the set of $(P, S)$-rows, and let $N$ be the set of numbers appearing in some rows in $R$. Construct a bipartite graph $G$ on vertex set $N\sqcup R$, with an edge $(n, r)$ if row $r$ contains $n$. Then we're looking for a subgraph of $G$ with $n^2$ $N$-vertices and $2n$ $r$-vertices with the correct edges.
  * The state of the art code for subgraph isomorphism is [Glasgow subgraph solver (2020)](https://github.com/ciaranm/glasgow-subgraph-solver) and [VF3 (2017)](https://github.com/MiviaLab/vf3lib)

* most of the work is in filtering valid_rows and valid_cols
	* like the ways we use valid_rows are:
		1. take subset of valid_rows that have intersection 0 or 1 with last_row
		2. take subset of valid_rows that have value at least minvec
		3. iterate over valid_rows
	* lists are good at 2 and 3 but not great at 1?
	* bitsets are good at 1 and 2 but bad at 3 because our sets are sparse?
	* if valid_rows/valid_cols is dense enough we can efficiently use bitsets instead of gather/scatter to get the intersections
* valid_rows and valid_cols have size around 10 each
* for 10-4-3-2, the number of test_rows is around 200, and numvecs around 450
* for 10-6-3-1-0-1, the number of test_rows is either around 1000 or 500, and numvecs around 1700
* numvecs gets larger for larger products

* idea for improving efficiency: after adding the first row and col, force each of the next rows and cols to intersect max_unmatched (without caring whether the set added is earlier than the remaining sets)

* vectorization
	* also for later: https://quickwit.io/blog/filtering%20a%20vector%20with%20simd%20instructions%20avx-2%20and%20avx-512
  * discussion about the article: https://news.ycombinator.com/item?id=32674040

## Notes from the 2026 rewrite (src/c/arrange.c, msearch, scheduler)

Where the time goes: in scheduler runs near S_min, ~99% of the CPU time is
the arrangement search itself (enumeration ~0.05%, relabelling and the
intersection tables ~1%), almost all of it on sums with N = 1000-5000
vectors. Sums with N < 1000 rarely produce squares.

What helped (see README, "Performance"): forward checking, branching on the
most constrained unmatched number (MRV), label bitsets for all bookkeeping,
fused filter + count, AVX-512 gather/compress filtering, root symmetry
breaking. Together 4-5x fewer nodes and ~3.2x less time than the legacy
search.

Things that were tried and did **not** pay off (numbers on `bench/quick.txt`
or `bench/full.txt`):

* Bitsets over vector indices (the N x N `inters0`/`inters1` rows) for the
  valid row/col sets instead of index lists: the sets stay sparse (tens to a
  few hundred out of thousands), so iterating and counting is slower than
  filtering short lists.
* Lazy candidate counting (count only when a cell is about to be chosen):
  slower, since the counts are what forward checking and MRV feed on.
* Arc consistency on the cell/candidate graph (repeatedly removing candidate
  vectors that would leave some cell uncoverable): 28% fewer nodes but much
  more work per node, slower overall.
* A pairwise Hall check (two unmatched cells whose only candidates are the
  same single vector): 7% fewer nodes, no gain in time.
* Pre-checking the cells a new vector would add before placing it: no
  effect (forward checking at the child catches the same cases).
* Pruning with the diagonals: keep, alongside the valid rows/cols, the list
  of vectors that meet every placed row and col exactly once (the possible
  diagonals), and prune when fewer than n remain or no two of them are
  disjoint. Correct (checked by planting the diagonals of found squares),
  but only ~7% fewer nodes: the diagonal list stays about as large as the
  valid col list, because the vectors are all built from the same few
  hundred numbers and are strongly correlated.
* Searching for magic squares directly, diagonals first: fix two disjoint
  vectors d1 < d2 as the diagonals, then search rows/cols among the vectors
  meeting both exactly once (a bitset AND of two precomputed rows), with
  forward checking. Also correct (finds every planted diagonal pair, each
  square once per transposition), but 15x slower than the semi-magic
  search on `bench/quick.txt` instances (e.g. P = 13 6 3 2, S = 517:
  0.18 s semi-magic vs 2.8 s): about 90% of the N^2/4 disjoint pairs
  survive the "every diagonal number has a row and a col candidate" filter,
  and each costs a small search. So for exhaustive search per (P, S), the
  semi-magic search followed by the diagonal check is the cheaper route.

## End-game shortcuts and the support filter (October 2026)

Where the nodes are (search as above, `bench/full.txt`, 126M nodes, by rows
and cols placed after the step): (2,2) 13%, (2,3) 36%, (3,2) 19%, (3,3) 25%,
(1,2) 2%. Most nodes four to six vectors deep die as soon as they are
created. The literal end-game (one axis complete, or one row/col missing) is
negligible: about 250 nodes have all rows or all cols placed and ~5000 have
8 or more vectors placed, out of 126M. So shortcuts there (looking up the
last row/col, which is determined by the unmatched cells; exact covers when
one axis is complete) cannot save anything measurable.

* One-axis feasibility check (did not pay off): the remaining rows must be
  n - r pairwise disjoint candidates, each taking one unmatched cell of every
  placed col, so a node is dead if no such set exists. Exact, and it kills
  ~83% of the (2,3)/(3,2) nodes that pass forward checking (12.6% fewer nodes
  on `quick.txt`, 20% on `full.txt`), but those nodes would have died one
  level down anyway, and the check costs ~300 cycles however it is done
  (recursive filtering of the list with the intersection matrices; a clique
  search on 64-bit masks over the local candidate list; one class of
  candidates per unmatched cell of a placed vector, with AVX-512 disjointness
  masks): 1.13 s -> 1.14-1.15 s on `quick.txt`.
* Support filter (kept, the big win): every cell of a candidate col must be
  in a placed row or in some candidate row, and vice versa, iterated to a
  fixpoint. It kills ~85% of the (2,2) nodes that pass forward checking and
  virtually all deeper ones. Nodes 7.5M -> 2.3M (`quick.txt`), 126M -> 31M
  (`full.txt`); time 1.12 s -> 0.67 s and 18.3 s -> 9.2 s. Cheap because the
  union of a list (slice 0 of the counts) says whether filtering is needed at
  all, and the first pass is fused into the filtering of the second axis.
  Variants: stopping after 1-6 passes instead of at the fixpoint is slower;
  filtering the disjoint axis first in deep children is 5% slower; computing
  only the unions while creating children (full counts only for searched
  nodes) helps from depth 4 on (3%) but not above; SIMD loops for the unions
  and the in-place filtering: 8%.
* Pairwise support (not kept): a col is kept only if every cell of it is
  covered by a row candidate meeting it exactly once. On top of the support
  filter it kills 39% of the searched (1,2) nodes and 99% of the searched
  (2,2) nodes (nodes 2.3M -> 1.6M on `quick.txt`), but done naively (one
  filtered list per candidate) the search is 8x slower; it would need a much
  cheaper formulation. Similarly, among the row completions of a (2,2) node
  (sets of 4 disjoint candidate rows, ~90 on average), only ~0.6 of its ~30
  candidate cols meet every row of some completion exactly once.
* The cheap special case of pairwise support (not kept): if a number is in
  exactly one row candidate r (count 1 in the bit-sliced counts), every col
  candidate containing it must meet r exactly once. Iterated together with
  the support filter it kills 28% of the searched (1,2) nodes and 97% of the
  searched (2,2) nodes, but a single pass kills only 4% / 29%: the power is
  in the cascade, which needs recounting both lists several times. 2.8x
  slower at the (1,2) nodes, 13% slower at the (2,2) nodes only.

## Candidate lists that carry their bitsets (October 2026)

With AVX-512BW and at most 256 labels, the candidate lists hold the label
bitsets of their entries (one array per 64-bit word) instead of vector
indices, and the filters compute |u & v| from them rather than looking up
bit u of row v of the N x N intersection matrices (see README). Before the
support filter this was 13-15% faster than the matrices; with it, where the
support filter's per-candidate test becomes a contiguous vectorized test of
the carried words, the search takes half the time (min of alternating runs,
Sapphire Rapids, `-march=native`: `quick.txt` 0.62 s -> 0.34 s, `full.txt`
8.2 s -> 4.4 s, P = 13 6 3 2 1 1, S = 905: 4.1 s -> 2.1 s, P = 12 7 4 2 1,
S = 900: 5.1 s -> 2.6 s; same nodes). Built with `-march=cascadelake` (no
VPOPCNTDQ/VBMI/GFNI/BITALG, so the fallbacks; run on the same machine):
`quick.txt` 0.64 s -> 0.39 s, `full.txt` 8.8 s -> 5.2 s; without the support
filter the fallbacks had been 5% slower than the matrices.

Where the time goes now (rdtsc around each phase, P = 14 7 5 3, S = 1460):
for the children with four vectors placed, the exactly-once filter 23%, the
disjoint filter 13%, the support filter 27% (3.4 passes per child that gets
there, two thirds of those children die in it), counting 4%; the children
one level deeper (nearly all dead after the disjoint filter) 17%.

Variants that did not pay off (min of alternating runs, `bench/quick.txt`):

* Before the support filter (matrices 0.67 s): gathering the bitsets of the
  list entries (2 qword gathers per word per 16 entries): 0.97 s, and still
  19% slower than the matrices at N = 2900-5500; rebuilding the matrix row
  of v in registers from the 6 label -> vectors rows of v's labels: 0.81 s.
  Carrying the bitsets in the lists, so that the loads are contiguous, is
  what made it faster (0.58 s). Also: counting the labels in the same pass
  as the filter (4-5% slower than counting only for the children that
  survive the forward check); deferring the compression of the lists until
  the child survives (saving the kept-entry masks, then compressing): 18%
  slower; running the first forward check for all the kids of a node before
  searching any: 14% slower; unrolling the filter loop by 2: 6% slower;
  skipping the words of v without labels when W = 2: 2% slower.
* With the support filter (carried bitsets 0.333-0.342 s): the support
  filter's first pass over the disjoint axis as a separate pass after its
  filter, rather than fused into it: 0.394 s (+17%).
* Filtering the disjoint axis first, with the support pass fused into the
  exactly-once filter instead: from depth 3 on 0.385 s (+13%), from depth 4
  on 0.383 s (+12%), from depth 5 on 0.353 s (+3%); `full.txt` +7% / +1%.
* Counting while filtering at shallow depths (the vpermb / gf2p8affineqb
  transpose on the compressed kept entries, recounting a list only if the
  support filter then drops some): below depth 3: 0.349 s (+5%), below 4:
  0.367 s (+10%), below 5: 0.417 s (+25%); the extra code alone (never
  counting early) costs 4%. So counting only the children that survive all
  the filters replaces LAZY_DEPTH on this path.
* One OR reduction for both words of the unions (W = 2: unpack, then one
  512 -> 128 bit reduction) instead of one per word: no change (0.341 s).
* LAZY_DEPTH, still used with the matrices (more than 256 labels, no
  AVX-512BW): byte counters (`-DCARRY_MAX_W=0`): 1: 0.633 s, 2: 0.632-0.655,
  3: 0.628-0.637, 4: 0.610-0.613, 5: 0.69-0.71, never lazy: 0.717;
  bit-sliced W = 8 (`--min-words 8`): 1: 0.733, 2: 0.737, 3: 0.713-0.732,
  4: 0.719-0.746, 5: 0.879, never: 0.916; AVX2 build (plain C, 4 slices):
  1: 1.356, 2: 1.302, 3: 1.334, 4: 1.347, 5: 1.607, never: 1.579. 2-4 are
  within noise, so it stays 4.

## Label order and root strategy (October 2026)

The search finds each square from its largest label x (r1 and c1 go through
x, every other vector has only labels < x), so the label order decides how
the search is split into one subproblem per number. On `bench/full.txt` the
work is spread over many r1 (the top 100 r1 hold 20-40% of the nodes) and
over the 20-30 rarest anchors, whose subproblems keep most of the vectors
(85-95% on the largest instance); nearly all nodes are 4-5 vectors deep.

* Degeneracy order (kept): the largest label goes to the rarest number, the
  next to the rarest among the vectors left without it, etc. 2-8% fewer
  nodes on every instance (quick 2.07M -> 1.94M, full 26.0M -> 25.2M, the
  largest instance only 2%), ~2% less time on `full.txt`.
* Other label orders, nodes on `quick.txt` relative to plain frequency:
  ties by larger value -0.4%; frequency ascending +43%; random +29%; number
  of distinct co-occurring numbers +6% (degeneracy version +3%); pairs of
  vectors meeting exactly at x +0.3% (degeneracy version -5%). Variants of
  the degeneracy order are all within 0.5% of it (-6.3%): other tie-breaks,
  counting only the vectors with a partner meeting them exactly at x, or
  also removing, after each step, the vectors that can no longer be in a
  square (a cell with no vector meeting it exactly there, or fewer than
  n - 1 disjoint vectors).
* Upper bound for orders: a local search over label orders (moving a number
  up to 12 places, keeping improvements, one search per step) found only
  2.7% (S = 506) and 3.8% (S = 561) fewer nodes than the degeneracy order.
* Vector orders that are not by label (r1 = first vector of the square in a
  greedy order, e.g. fewest partners through its most constrained cell):
  -3.6% at best, worse than the label blocks, where c1 shares the rare
  anchor with r1.
* Depth 1: always branching on the anchor x (choosing c1 right after r1):
  +1.7% nodes with labels by frequency, +0.1% with the degeneracy order (MRV
  picks x for ~75% / ~95% of the r1 anyway).
* Shrinking each anchor's universe by the reduction rules before searching
  it: -0.04% nodes (forward checking and the support filter catch the same).
* One-step lookahead at shallow depths (the cell with the fewest children
  that survive their creation, instead of the fewest candidates): more
  nodes (+0.1% at depth 1, +3% to depth 2, 3x to depth 3).
* The one-axis feasibility check (see above) near the root, where it could
  afford to be slow: it kills 3-13% / 0.4-4% / 2-20% of the searched nodes
  at depths 1 / 2 / 3 (all but the smallest instance of `quick.txt`), but
  those die within ~2-5 nodes anyway (2% fewer nodes on S = 517). The nodes that cost are the ones that look alive for a few
  more levels.
* Pairwise support (see above) only near the root: -2% nodes when applied
  up to depth 2, -26% up to depth 3. The naive version (a pass over the other
  list per candidate) is 10-75x slower; even vectorized (filter the other
  list by the candidate's inters1 row, then a union) it would take a few
  microseconds per depth-3 node, against ~0.4 us of nodes saved per node.

## Per-node overhead on the matrix path (October 2026, not merged)

Measured on the matrix path before the support filter and the carried
bitsets (branch `opt/node-overhead`, 8-13% faster there, same nodes):
children created in batches of 16 (filter all their other-axis lists, then
count, then forward-check, with no branch per child), 16-bit vector indices
with AVX-512 VBMI2 (32 entries per register in the filter; ~10% alone),
prefetching the batch's matrix rows (~2.5%). Profile of that code: counting
~35% (~2 cycles per vector: mask loads and a single-port saturating add),
filtering ~30%, choosing the cell and the children ~10%, the rest
bookkeeping and the mispredicted "dead?" branch after a child's first pass.
Not ported: since the lists carry their bitsets, the matrix path only runs
for more than 256 labels or without AVX-512BW, and batching the children's
first pass on the carried lists was 14% slower (see above).

## msearch overhead outside the search (October 2026)

Time split of `msearch` units near S_min (instrumented build, one core,
noisy machine; "overhead" = enumeration + reduction + per-sum setup +
output, everything but `search_root`):

| unit | sums | search | overhead before | overhead after |
| --- | ---: | ---: | ---: | ---: |
| P = 13 6 5 2, S = 595-760 (N <= 890) | 157 | 0.57 s | 0.100 s (15%) | 0.022 s (3.6%) |
| P = 9 5 4 2 1, S = 356-470 (N <= 1135) | 111 | 1.05 s | 0.091 s (8%) | 0.022 s (2.1%) |
| P = 13 6 5 2, S = 595-889 (N <= 1700) | 286 | 26 s | 0.38 s (1.4%) | 0.077 s (0.3%) |

Before, the overhead was (first unit) the reduction 0.035 s, the search
setup 0.046 s and the enumeration 0.018 s; in the reduction and the setup,
nearly all of it was sorting: a qsort of all 6N elements to find the ~20-200
distinct values (then a binary search per element), and a qsort of the N
label arrays with a comparator. Now the values get dense ids from a small
hash table (`dense_ids`, enumerate.c), the label arrays are packed into
64-bit keys (n labels of ceil(log2 L) bits, when that fits: always for
n = 6 and L <= 1024) and radix sorted, the weak reduction no longer builds
the vector lists only the strong one uses, and msearch enumerates with
`enum_vectors_grouped`, which skips the per-sum sort that only the
`enumerate` file output needs (a third of the enumeration). The output is
identical (same squares in the same order, same node counts), since the
labels and the vector order of the search are the same.

What is left: enumeration 0.010 s (the divisor DFS and first-touch page
faults of the vector lists), relabelling + key sort 0.005 s, allocations
0.001 s, output 0.001 s. Wall time of the first unit 0.687 s -> 0.606 s,
second 1.148 s -> 1.118 s (min of 9 alternating runs; the second is within
the noise of the shared machine). For the units the scheduler runs (about
two minutes, dominated by the sums with N = 1000-2500) the gain is ~1%: the
search is 99.7% of the time, and the sums with small N, where the overhead
mattered, are cheap anyway. Not done, as they would gain < 0.3% of a real
unit: reusing the search allocations across sums, a faster enumeration DFS
(pow() per node), skipping sums with small N (not provably square-free, and
they cost nothing). The weak reduction removes almost nothing near S_min
(184 of 222,070 vectors on the third unit), and the strong one 3% of the
vectors but only 0.2% of the nodes (the search discards those vectors at
the first level anyway) and costs more than it saves (first unit: 0.30 s
of setup + reduction instead of 0.08 s before this change), so it stays
off.

## Stronger pruning on top of the support filter (October 2026)

(Measured on the matrix path, before the candidate lists carried their
bitsets, on branch `opt/stronger-pruning`; for the carried lists see the
next section.)

Where the work is with the support filter (`bench/full.txt`, 26.0M nodes;
"(r,c)" = rows and cols placed): creating the (2,2) children of the (1,2)
and (2,1) nodes is ~60% of the time (13.2M children, 80% pruned, mostly by
the support cascade, ~900 cycles each), and the (2,3)/(3,2) children of
the surviving (2,2) nodes ~15% (10.4M children, all pruned, ~250 cycles
each). Kill rates below are of the searched (2,2) nodes on P = 13 6 3 2,
S = 561 or 699 unless noted, each followed by the support filter.

* Pair support (PS): every cell of a candidate outside the placed vectors
  of the other axis lies in a candidate of the other axis meeting it
  exactly once. Iterated: 98% of (2,2) killed, -29% nodes; one pass on
  both axes 93%, on one axis 80% (-24% nodes); at (1,2)/(2,1) -16..20%
  nodes. Even with label -> candidate-mask tables (T[x] = 64-bit mask of
  the other list) it costs ~1500 cycles per node: 20% slower overall.
* Count-1 rule (a number in exactly one candidate u of the other axis:
  candidates through it must meet u once): iterated 96% at (2,2), one
  pass 83%; -10% nodes at (1,2)/(2,1). Needs per-label lookups like PS.
* One-axis exact cover (the n - r rows must be disjoint): 17% at (2,2),
  -2.5% nodes. Same-axis pair support (candidates disjoint from v must
  cover the cells v does not): 73% at (2,2), -17% nodes, pairwise. A count
  of free numbers covered by both axes (>= (n-r)(n-c)): never fires.
* Static pair consistency (drop pairs of vectors whose 2-vector node dies
  under forward checking + support, with all vectors as candidates):
  removes thousands of pairs but only 0.2% of the nodes.
* Branch support: exactly one child through the branching cell is in any
  completion, so the other axis must meet some child once and the same
  axis must be disjoint from some child (unions of inters rows), then the
  support filter. 58-70% of (2,2) killed, -25..36% nodes, but ~800 cycles
  per node: no gain in time. At (1,2) only 4-10% killed.
* Cell support: for every unmatched cell y, the other axis must meet some
  candidate through y exactly once: 99.96% of (2,2) killed (93% even
  without the support filter after it); at (1,2)/(2,1) 20%/29%.
* Cross support (kept, `CROSS` in arrange_core.h): the >= 1 relaxation of
  cell support, "a candidate must meet the union of the candidates of the
  other axis through each unmatched cell", needs only bitset ANDs: 99.7%
  of (2,2) killed (85% without the support filter). Applied at searched
  (2,2) nodes it removes the (2,3)/(3,2) layer (-39% nodes, -9% time);
  applied to the (2,2) children before the support filter, where it kills
  70-90% in one pass per axis, the support cascade runs 7-10x less often:
  -19% time on `full.txt`, -7% on `quick.txt` (fewer (2,3)/(3,2) nodes
  there). Cost ~350-500 cycles per call with AVX-512 (pext to the cells'
  lanes, masked ORs, one test per candidate). Not kept: at (1,2)/(2,1)
  nodes (kills 12-21%, costs ~1600 cycles: slower), at the (1,3)/(3,1)/
  (2,3)/(3,2) children (no gain), with 8 words of labels (bit-sliced
  counters: 10% slower), in plain C / AVX2 (10% slower), and for n = 7
  (P = 12 6 3 2 1, S = 290..380: 6% fewer nodes but 25% slower; the 7x7
  search dies further down). Gains vary with the instance: -64% nodes and
  -30% time on P = 12 6 3 2 1 0 1, S = 900 (N = 2994, 18 squares), but
  only -5% nodes and no gain on P = 16 5 4 2, S = 1200, where the support
  filter already kills 98% of the (2,2) children.
* Cheaper cross support (none kept): fusing the two cross passes into the
  child's filtering passes (unions of the new o candidates built with their
  union, b candidates tested while filtered) was slower, as both halves
  then always run (the first one alone kills about half); a byte table
  (label -> cells) looked up with vpermi2b for 8 candidates at a time made
  the test ~3% cheaper per call, not measurable overall; iterating cross
  and support, cross after the support filter, or rows first: no gain.
  Unions over the parent's list (shared by the ~6 siblings): only 6% of
  the children killed (31% when that list is first filtered by
  disjointness from the new vector, which is the child's own work).

## Cross support on the carried bitsets (October 2026)

Cross support (above) ported to the candidate lists that carry their
bitsets (AVX-512BW, up to 256 labels: production), and re-tuned there,
since not all of the conclusions on the matrix path carried over: on the
carried lists the support filter is a cheap vectorized test of
contiguous words, and so are the (2,3)/(3,2) children that cross support
removes. Min of alternating runs on the shared Sapphire Rapids machine
(3-5% noise, often thread CPU time to cut it); "prod" = P = 13 6 3 2 1 1,
S = 905, P = 12 7 4 2 1, S = 900 and P = 12 6 3 2 1 0 1, S = 900 (N =
2161-2994, 106-121 labels).

Implementation (`CROSS_AXIS`, `CROSS_TEST`): for each group of 8 unmatched
cells, 8 entries of the other list at a time, a test per cell (the entries
through it) and a masked or per word into that cell's accumulator, then one
transposing or-reduction of the 8 accumulators of a word (`or_lanes8`: lane
j = U_y of cell j). Then 8 candidates at a time, t_y = u & U_y (an and per
word) for every cell, kept if the minimum of the t_y is nonzero (a tree of
vpminuq; a chain of masked tests was 2-5% slower per call), and the
support filter's test folded in for free; compressed in place. The lists
are padded with empty sets, so there are no lane masks. In a replay
harness over dumped (2,2) states: ~250-300 cycles per call (two passes,
unless the first one kills) on `full.txt` (lists of ~35 entries), ~400 on
prod (~47): building the unions ~60%, the test ~40%. The setup was not
unrolled at first (accumulators spilled and reloaded): -4% per call;
padding instead of lane masks: -1..3%; the cell list with pdep instead of
a loop: no change.

Where it runs, at the (2,2) children (`full.txt` with e9a0540: 25.2M
nodes, 4.3 s; prod 91.6M nodes, 13.6 s):

* Room: skipping (wrongly) every (2,2) child that survives the support
  filter saves 14% of the time on `full.txt` and 37% on prod; no filter of
  that layer can save more.
* The branch's placement (before the support filter, cols first, no
  support test in the passes): -40% nodes but only -1% time on `full.txt`
  (10.4M calls at ~570 cycles). After the support filter, on its 2.6M
  survivors: -34% nodes, -3%; and then the support filter again when cross
  support dropped something: -41% nodes, -3..5% (prod -9..16%). Before the
  support filter with its test folded into each pass: -5% / -17%; with the
  pass on the axis that the support filter checks first (o, not the axis
  of the vector just placed) first: -6% / -20%; with the micro-
  optimizations above: -9% / -25..30%.
* But where the support filter alone kills nearly all the (2,2) children
  (P = 16 5 4 2, S = 1200, 133 labels: 98%), running cross support on all
  of them first costs more (+10-14%, part of it the inlining below) than the
  few searched (2,2) nodes it saves, while after the support filter it is
  nearly free. Before is better where the support filter kills about half of
  them (prod: 46%; ~8% faster than after), and the two are equal on
  `full.txt` (75%). So the mode adapts (`cross_after`): one (2,2) child in
  16 runs the support filter first, and once it has killed more than 3/4 of
  the last ~256 of those, all of them do (and measure it). Same time as
  "before" on `full.txt` and prod (it stays mostly before there), 4-8% less
  than "before" on P = 16 5 4 2, 1-4% less on `quick.txt`. Thresholds 5/8,
  3/4 and 7/8: within noise on `full.txt`.
* Not kept: also at the (1,2)/(2,1) children (12% fewer (2,2) children,
  but 10-15% more time, before or after the support filter); at the
  (2,3)/(3,2) or (1,3)/(3,1) children (no gain); one pass only, on either
  axis (+5-7%: with the support filter after it, 85% of the children
  killed instead of 98%); the support cascade between the two passes (more
  kills, ~3% slower); rows first (same); n = 7 (P = 12 6 3 2 1, S =
  352-355: -7% nodes, +33% time, so n <= 6 by default as on the matrix
  path).
* Kill structure (dumped states): about half of the candidates tested in a
  pass miss some U_y, spread over cells with 1 to ~20 candidates of the
  other axis through them, so testing only the cells with few candidates
  would lose most of the kills (and building the unions is the larger
  part anyway). 55% of the children that a pass kills die by count (fewer
  candidates than vectors still to place), the rest by forward checking.
* Inlining: with cross support, GCC no longer inlined TRY_CHILD into
  SEARCH_REC for W = 3 (+5% on P = 16 5 4 2, also with --no-cross). Forcing
  it on the carried path is faster for W = 2 too, also on e9a0540 alone
  (`quick.txt` -6%, `full.txt` -3%, prod no change), but 7% slower with the
  matrices (-DCARRY_MAX_W=0), where it is not forced.
* W = 3 / 4 (`--min-words`): -2..5% / no change on `full.txt`.

Result (bench wall time, min of alternating runs; e9a0540 -> this):

| | nodes | -march=native | -march=cascadelake |
| --- | ---: | ---: | ---: |
| `quick.txt` | 1.94M -> 1.77M | 0.324 -> 0.303 s | 0.363 -> 0.340 s |
| `full.txt` | 25.2M -> 15.0M | 4.33 -> 3.71 s | 4.98 -> 4.16 s |
| P = 13 6 3 2 1 1, S = 905 | 14.4M -> 5.6M | 2.14 -> 1.62 s | 2.33 -> 1.78 s |
| P = 12 7 4 2 1, S = 900 | 18.2M -> 7.0M | 2.52 -> 1.92 s | 2.86 -> 2.16 s |
| P = 12 6 3 2 1 0 1, S = 900 | 59.1M -> 19.8M | 9.05 -> 6.93 s | 10.08 -> 7.35 s |
| P = 16 5 4 2, S = 1200 | 3.62M -> 3.44M | 0.93 -> 0.95 s | 1.09 -> 1.07 s |

On the time-weighted sample of production sums (`bench/prod.txt`, added
on the integration branch in 10df892: nine sums, N = 1491-2994): 128.9M ->
50.4M nodes, 18.9 -> 14.9 s native, 21.5 -> 16.5 s cascadelake (-21% /
-23%; per sum -6% (P = 10 4 2 2 1 1, S = 460) to -28%).

Of this, the forced inlining is ~6% on `quick.txt` and ~3% on `full.txt`
(e9a0540 with it: 0.302 s, 4.15 s), nothing on prod. On the matrix path
(-DCARRY_MAX_W=0, byte counters), the branch's code as it was: `quick.txt`
-6%, `full.txt` -22% (8.15 -> 6.37 s CPU).
## Diagonal-first search (October 2026, second look)

A square has a vector d as a traversal (a possible diagonal) iff every row
and col meets d in exactly one number, so the squares with d as a traversal
are exactly the semi-magic squares that can be built from
V_d = { v : |v & d| = 1 }. Running the semi-magic search once per d on V_d
("d-first", `bin/dsearch`, src/c/dsearch.c) finds every (square, SP
traversal) pair exactly once (checked on the three SP-type squares of the
October runs: P = 13 5 3 2 0 1 / S = 632, 16 5 4 2 / 849, 12 6 3 2 1 0 1 /
836), and only those, which is what a magic square needs twice. The
restricted searches are cheap (|V_d| ~ 0.25 N, ~1000-50000 nodes each),
but there are N of them, and the same (r1, c1) is re-explored under every d
meeting both: nodes relative to the plain search (same P, S; the plain
search then still has to check the diagonals, which is negligible):

| P | S | N | plain nodes | d-first nodes | ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| 10 4 3 2 | 327 | 452 | 8.7K | 118K | 13.6 |
| 13 6 3 2 | 517 | 976 | 409K | 1.72M | 4.2 |
| 12 6 4 2 1 1 | 994 | 1153 | 663K | 2.41M | 3.6 |
| 10 6 3 1 0 1 | 648 | 1682 | 1.02M | 3.78M | 3.7 |
| 12 6 4 2 1 1 | 1026 | 1996 | 10.9M | 16.0M | 1.46 |
| 14 7 5 3 | 1488 | 2450 | 22.5M | 26.3M | 1.17 |
| 14 7 5 3 | 1560 | 3018 | 59.2M | 52.0M | 0.88 |
| 14 7 5 3 | 1656 | 3787 | 161M | 106M | 0.66 |
| 14 7 5 3 | 1760 | 4510 | 289M | 169M | 0.58 |
| 14 7 5 3 | 1900 | 5251 | 371M | 217M | 0.59 |
| 14 7 5 3 | 2050 | 6206 | 608M | 324M | 0.53 |

(time ratios are within 10% of the node ratios, with per-d setup included
in the d-first time). So d-first only wins beyond N ~ 3000, and then by at
most ~2x: the ratio flattens at 0.5-0.6, because past N ~ 4000 the plain
search's nodes grow only like N^2.7 at this P (S/S_min = 1.5-1.9), not the
N^5 seen near S_min. Large N is also where the squares are least likely to
be magic: for a fixed P, squares per sum grow like N^4 near S_min
(S/S_min < 1.2, big-tau P) but only like N^1.5 at S/S_min ~ 1.5-1.9, while
p_S and p_P fall with S (e.g. P = 14 8 5 3 1 1: p_P 4.6e-4 at N ~ 1200 ->
1.7e-4 at N ~ 3500, p_S 1.6e-3 -> 6.9e-4), so the predicted magic squares
per CPU-hour at N > 4000 are 5-100x below the near-S_min values. Fixing
both diagonals first (N^2/4 roots, V_{d1,d2} ~ 0.04 N) has the same
redundancy with a worse constant (the 15x of the earlier note at N ~ 900),
and pruning the plain search by "some vector disjoint from d1 still meets
every placed vector once" only bites below the depths where the nodes are.
Conclusion: not a lever for the search as scheduled (N ~ 1000-2500 near
S_min); kept as an experiment driver only.

Traversal counts per square (5024 squares of the October runs, N ~
1000-3500): s_count is close to Poisson(0.98) (0: 1990, 1: 1759, 2: 830,
3: 314, 4: 105, 5: 20, 6: 5, 7: 1), p_count close to Poisson(0.23)
(0: 4021, 1: 860, 2: 126, 3: 13, 4: 4), and sum s*p/720 = 1.46 against 5
SP traversals observed (rho ~ 2-3 with large error bars). The legacy
hSS = sum C(s, 2) per P is also within ~20% of the Poisson value. So there
is no visible class of "structured" squares with many S- or P-traversals
to target: a square's chance of being magic is just 5400 (rho p_S p_P)^2.

3-prime P (2^a 3^b 5^c, tau 765-2275) have fewer than 600 vectors up to
1.5 S_min (S_min 486-4332), where squares essentially never occur, so the
higher p_P one would expect with fewer primes is not reachable there.

Model check on fresh P (October 2026): with the model refit on the legacy
totals plus ~3 CPU-hours of new runs, `plan` ranked P = 12 9 4 2, 14 6 5 2
and 11 6 5 3 at the top with 10,600-13,300 predicted squares per CPU-hour
(optimistic scores, i.e. with the exploration bonus). Running each for 15
minutes from S_min gave 1700-2000 squares per CPU-hour (0.8-1.0 per sum
over ~500 sums), so the top of the ranking over-predicts the semi-magic
rate by ~6x: the ranking selects the P with the largest positive model
error (winner's curse), on top of the bonus. The refit model's 600-hour
forecast (19 CPU-years per magic square, 2e-5 magic/CPU-hour at the start
decaying ~3x by 600 hours) should be read with that in mind; the direct
estimate from the new runs (25 SP-type squares in 20 CPU-hours, each with
P(magic) ~ 15 rho p_S p_P ~ 1e-5) is ~1.5e-5 magic/CPU-hour, i.e. ~8
CPU-years at a rate that cannot be sustained once the best sums of the best
P are used up. Wider candidate families (no factor 7, 11^2 / 13^2, prime
23, exponents beyond gen_candidates' ranges: scripts/cand_explore.py,
samples of 1500 and 1200 P) put no new P in the model's top 200 and ~1% of
the score mass in the top 1000, so the candidate space is not the limit.


## The support / cross cascade at the (2,2) children (October 2026, no gain)

(Branch `opt2/cascade`, commit c6a487c, from f61e719; not merged. It has
the measurement switches `-DCASCADE_STATS` and `-DCASCADE_DUMP=K` in
arrange.c, which compile to the same code when off, and
`research/cascade_replay.c`. No change to the search: nothing below beat
the noise of the shared machine.)

How it was measured. `-DCASCADE_STATS` counts, per class of child (rows,
cols placed) and step of TRY_CHILD on the carried path, the calls,
candidates in and out, children killed and cycles (rdtsc; ~25 cycles of
overhead per step, prod 15 s -> 25 s, so only the shares are meaningful).
`-DCASCADE_DUMP=256` writes one (2,2) child in 256 as it enters cross
support (143,587 states on `bench/prod.txt`, 210 MB), and
`research/cascade_replay.c` replays them through variants of the cascade,
timing each state alone (min of 3 rounds, ~50 cycles of timing overhead
subtracted below), which resolves differences of 1-2% that whole runs on
this machine cannot. In-place timings: per-instance minima over 2-5
alternating rounds, summed (native and `-march=cascadelake`, both on
Sapphire Rapids); even so, two builds of the same code differ by up to
~1.5% (f61e719 vs this branch, whose search code is identical: prod 15.00
vs 14.79 s native, 16.48 vs 16.51 s cascadelake; full 3.65 vs 3.70 s,
4.10 vs 4.09 s).

Where the time goes (prod, native, default adaptive mode; the instrumented
cycles): the (2,2) children are 85% of all child creation (83% with
`-march=cascadelake`), the (1,2)/(2,1) children 11%, (2,3)/(3,2) 2%. Per
(2,2) child (39.2M of them, 96.6% killed, 1.32M searched), by step:

| step | runs on | in -> out | kills | cycles/call |
| --- | ---: | ---: | ---: | ---: |
| filter o (exactly once) | 100% | 187 -> 44 | 1.3% | 240 |
| filter b (disjoint, 1st support pass) | 99% | 89 -> 36 | 5.1% | 121 |
| cross o (before mode) | 73% | 48 -> 21 | 27% | 264 |
| cross b | 53% | 46 -> 10 | 58% | 234 |
| support passes 1, 2, 3, ... | 43%, 26%, 15%, ... | 35 -> 16, 27 -> 17, ... | 37%, 33%, 29%, ... | 83-93 |
| after mode (22% of the children): support, cross o, cross b, support | | | | |
| count (survivors) | 3.4% | | | 129 |

So filters 34%, cross passes 32%, support passes 10% (2.5 passes per child
that reaches them), counting 0.4%, the rest rdtsc and bookkeeping. With the
before mode forced: cross o kills 36%, cross b 60% of the rest, the
support loop 86% of the rest (2.4 passes). With the after mode forced, the
support loop alone kills only 50%, after 3.4 passes on average (10% of
the children need 7 or more), which is why cross support goes first.
With `-march=cascadelake`: filter o 284 cycles (no VPOPCNTDQ), counting
207, the cascade the same.

Replay (cycles per state, net): before mode 369, after mode 446, support
filter alone 233 (kills 50%), the two cross passes alone 331 (74%), one
cross pass 215 (36%). One cross pass (W = 2, 8 cells, lists of ~36 and
~45): call and setup ~17 + 30, building the 8 unions 63, the transposing
reduction 14, the test 100. Both loops are bound by ports 0/5 (24 and 32
vector uops per 8 entries); so is the support pass (~45 cycles).

Tried, none kept:

* The cells of one placed vector only (4 of the 8; every candidate of the
  other axis passes through exactly one of them, so the test stays valid):
  no cheaper, as the pass works on groups of 8 cells, and much weaker
  (prod nodes 50.4M -> 64.9M; first pass kills 16% instead of 27%).
* Cross b first: the first pass kills 14% instead of 36% (cross b is strong
  only on the o list that cross o has shrunk); 52.0M nodes, replay 407 vs
  369. Choosing the order per child: b first is 1-2% cheaper only when
  ko < kb (17% of the states).
* The support loop between the two cross passes (to spare cross b on the
  children it kills): replay 397 (one support pass: 379) vs 370.
* One support pass as a probe, then after mode if it dropped more than a
  fraction of the o list: the probe costs ~33 cycles, best 389 vs 369.
* Choosing before / after per child: a per-child oracle would save 7.7% of
  the cascade (11% with b-first as a third order), the best mode per
  instance and list-size bucket 2.6% (~1% of the time). Rules tried in
  place: by list size (ko + kb < 50-70: within noise on full, +4-10% on
  P = 16 5 4 2, S = 1200, where the sampled switch is right), by the
  number of o labels the support test would reject (>= 11, alone or with
  the sampled switch: replay -2.7% of the cascade; in place prod -0.6 to
  -1.9%, P = 16 5 4 2 -0.7 to -2.9%, full +1.0 to +2.4%: noise). It does
  match the sampled switch without sampling, if that switch is ever to be
  simplified.
* The sampled switch's threshold: all prod sums prefer the before mode
  (S = 460, where the support filter kills 82%, is a tie), and the replay
  puts the break-even at a support kill rate of 0.80-0.85 per sum and
  0.86-0.92 per size bucket, so 3/4 looked low; but forced before mode is
  only -1.1..-1.5% on prod (and +7-14% on P = 16 5 4 2), and 13/16, 7/8
  and 7/8 with sampling 1/32 instead of 1/16 are all within noise (7/8,
  1/32: prod -1.2% native, +0.9% cascadelake; P = 16 5 4 2 +1.8/+3.0%).
* Iterating cross support and the support filter to their common fixpoint
  (a cross pass reruns when the other list has shrunk since its unions
  were built): the (2,2) children searched drop from 1.32M to 21k, their
  (2,3)/(3,2) children from 6.0M to 0.2M, nodes 50.4M -> 44.6M (-12%); but
  time: native +1.4% and -0.6%, cascadelake -0.3% and -1.1% (two sets of
  runs): neutral. One extra round only: 46.6M nodes, +0.7..+2.1%. A
  searched (2,2) node costs only ~700 cycles (counting, choosing the cell,
  ~4.5 children at ~120 cycles that nearly all die in their first filter),
  about what the second round costs per node it kills.
* Hall's condition (stronger than cross support): a candidate u of axis f
  meets the n - 2 remaining vectors of axis h at its n - 2 free cells, and
  each of those passes through one unmatched cell of each placed f vector
  p, so the sets t_y = u & U_y (y: the unmatched cells of p) must have a
  system of distinct representatives (cross support checks |t_y| >= 1).
  Replay, scalar: cross o + Hall kills 44% instead of 36%, before mode +
  Hall 99.75% instead of 96.39%. But the extra kills are mostly nodes
  worth ~700 cycles (as above), and a vectorized test (pairs >= 2 bits,
  triples >= 3, per placed vector) would cost several times the cross test.
* Cross support at the (1,2)/(2,1) children, iterated: kills 12% of them
  (lists of 81 and 167), 6% fewer (2,2) children, ~2000 cycles more per
  child: slower, as on the matrix path.
* The cross pass itself (replay, net): a single-group version with
  branchless cell extraction (pdep per word, or a select chain, which GCC
  turned into branches and spills): 219-234 vs 216 cycles per pass; the
  cells extracted in SIMD lanes (clear the lowest bit R times): setup 66 vs
  47; forcing CROSS_AXIS inline: 375 vs 369 per state. Starting the
  support filter with the other axis known closed (it always is): no
  change (369.1 vs 369.6).
* By analysis, not built: an early exit when a list falls below the needed
  count can skip at most the last block of a pass (lists of ~40 losing at
  most 8 per block, 4 needed); incremental unions (per-label counts) save
  only the 2 ORs per block of a support pass, whose floor is the scan,
  and re-testing only the entries with newly uncovered labels needs an
  inverted index per child (a transpose, ~100+ cycles), more than the
  whole support loop costs where it runs (~100 cycles on 26% of the
  children).

What is left at the (2,2) children is the filters (34%, mostly filter o,
the largest single step: the exactly-once test over the ~187 entries of
the parent's list, of which 23% are kept), and the first cross pass on
nearly every child. An idea for filter o, outside this direction: an
inverted index of the parent's list (label -> mask of entries, 192 bits),
from which "meets v exactly once" is a few mask operations over v's 4
labels, leaving only the compress work per child.


## Stronger pruning above the (2,2) children, carried path (October 2026, not merged)

Question (branch `opt2/shallow`, from f61e719): now that the candidate
lists carry their bitsets, do the stronger rules of "Stronger pruning on
top of the support filter" pay off at the searched nodes with one or two
rows and cols placed, where a kill saves the whole subtree? No: at
production sizes the nodes that hold the time are not killed by any of
them, even at zero cost. Nothing kept; the instrumentation (rdtsc per
subtree, reference implementations of the rules) is not committed.

Where the time is (rdtsc around each searched node's subtree,
`bench/prod.txt`, -march=native): the subtrees of the searched (1,2)/(2,1)
nodes are 85-92% of the search (S = 900: 12.8 of 14.0 Gcyc), 3400-8600
cycles each; a subtree is the node's 6-12 children with two rows and two
cols, nearly all of them killed as they are created (by cross support and
the support filter, 620-810 cycles each). The searched (2,2) nodes'
subtrees are 0.1-5%; the (1,1) nodes' own work, creating the (1,2)/(2,1)
children, most of the rest.

Method: reference implementations in plain C, run on copies of the lists
of 1 in 8 searched nodes, each rule iterated with the support filter and
forward checking to a fixpoint; the cycles of the node's subtree (the
evaluation excluded) give what a kill would save at zero cost. Rules, for
a candidate u of axis f (h the other axis, y an unmatched cell of a placed
f vector, i.e. one that needs an h vector):

* cross: u meets the union of the h candidates through y, for every y;
* cell: u meets some h candidate through y exactly once, for every y;
* pair: every cell of u outside the placed h vectors is in an h candidate
  meeting u exactly once;
* p1 (count-1 rule): a number in exactly one h candidate c: the f
  candidates through it meet c exactly once;
* also same-axis support, forced numbers (the intersection of the
  candidates through a cell), Hall (with two placed h vectors, the f
  candidates are edges between their unmatched cells and must lie in a
  perfect matching), and exact one-axis feasibility (n - r pairwise
  disjoint candidates, a DFS).

No rule ever killed a node with a square in its subtree. Share of the
(1,2) subtree time that the kills would save (the (2,1) nodes: 1-2x
these), on a machine shared with other benchmarks:

| sum | N | cross | cross + p1 | cell | cell + pair | all rules |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| quick, P = 10 6 3 1 0 1, S = 648 | 1682 | 38% | 58% | 55% | 77% | 77% |
| full, P = 13 6 3 2, S = 561 | 1195 | 10% | 20% | 23% | 40% | |
| full, P = 14 7 5 3, S = 1460 | 2237 | 2.6% | 5.3% | 5.9% | 10% | |
| prod, P = 14 7 5 3, S = 1361 | 1491 | 4.7% | 9.5% | | 27% | |
| prod, P = 10 4 2 2 1 1, S = 460 | 1932 | 7.9% | 16% | | 29% | |
| prod, P = 10 6 4 2 1 1, S = 838 | 1986 | 2.4% | 5.0% | 6.0% | 11% | 12% |
| prod, P = 14 7 4 2 1, S = 1085 | 2082 | 1.3% | 3.0% | | 7.2% | |
| prod, P = 12 6 3 2 1 0 1, S = 900 | 2994 | 0.9% | 1.8% | 2.1% | 3.7% | |

Weighted by the time of each sum, cell + pair would save ~8% of prod at no
cost, cross + p1 ~4%; one pass of each rule instead of the fixpoint, a
sixth to a third of that. The exact one-axis check kills 0.6-1.2% of the
nodes (0.1% of the time, S = 838) and nothing on top of the others, Hall
0-1% on top of cross.

Why: the rules kill the cheap nodes. By the number k of children of the
node (the MRV count), S = 900: the nodes with k >= 9 hold 81% of the
(1,2)/(2,1) subtree time and cell + pair kills 0.0-0.4% of it; k = 5-8,
17% of the time, 17% killed; k <= 4, 2% of the time, 60-91% killed. These
nodes are consistent under every local rule: their children die only
when the second row (or col) is placed. The rules do not shrink the lists
of the nodes they leave either (their MRV count drops by 1-4%), and at the
(1,1) nodes cell + pair kills 5% of them, with 0.7% of the nodes below
them, and leaves the lists of the others the size they were (S = 1361).
Singleton consistency (place each candidate, run the filters that create
the child, drop the candidate if the child is pruned) finds the same: at
the (1,1) nodes 4% of the candidates go (11% of the nodes die), -2.3%
nodes, for probes costing twice the subtree; at the (1,0) nodes 10% of
the candidates go, -0.1% nodes (S = 1361). So the (2,1)/(1,2) nodes
under a (1,1) node are nearly all consistent too, and a cache of dead
(2,1) nodes to skip the (2,2) children that contain them would not help.

Cost, against the ~300-900 cycles per node that the best rule could save
on prod (S = 900: 3.7% of ~8600): cross support iterated with the support
filter at the searched (1,2)/(2,1) nodes, with the tuned CROSS of the
(2,2) children (14 cells, lists of 60-220): `full.txt` 14.96M -> 12.36M
nodes but 3.78 -> 5.51 s (~3000 cycles per node). The building block of
the pairwise rules, label -> candidate-mask tables for both lists (scalar,
masks of up to 512 bits), alone costs 3300-4700 cycles per node at prod
sizes (167-232 entries in the two lists, S = 1361 / 838 / 1085): the
lists are in vector order, so consecutive entries share their largest
numbers and the read-modify-writes of the table rows form chains; a
64 x 64 bit transpose might be ~4x cheaper (estimate), still more than
the saving.

At the (2,2) nodes:

* The searched (2,2) nodes (which survived cross support and the support
  filter as they were created) die under cross support iterated with the
  support filter (98-99.6%): one round of cross support is not a fixpoint,
  as the second pass changes the unions the first one used. Running
  CROSS and SUPPORT in turn at the (2,2) children until neither drops
  anything, when the support filter dropped something after the cross
  passes (198k calls on P = 13 6 3 2 1 1, S = 899, 98% killed, ~590
  cycles each): prod 50.4M -> 44.6M nodes, `full.txt` 14.96M -> 14.75M,
  but no gain in time: prod native 14.78 -> 14.61 s (sum of the per-sum
  minima of 3 alternating rounds), cascadelake 16.11 -> 16.00 s and, in 4
  more rounds, 16.14 -> 16.10 s (min of totals 16.19 -> 16.30 s);
  `full.txt` native +2.5%, cascadelake -0.6%; `quick.txt` +-0. A searched
  (2,2) node costs about what the extra round costs (counting two lists of
  ~30, ~3 children dying at their first filters). Unconditionally (on
  every surviving child): the same nodes, no gain; a single cross pass on
  either axis instead of the round: +1% (three largest prod sums).
* The (2,2) children: where they die (S = 900 / S = 1361): first filter
  0.6 / 2.9%, second 2.8 / 11.2%, first cross pass 15.6 / 27.9%, second
  34.1 / 21.5%, support filter after them 43.1 / 36.2%, alive 3.9 / 0.3%.
  Kill rates right after the two filters (S = 1361 / S = 838), one pass on
  each axis: support 29 / 15%, Hall 2 / 1%, p1 19 / 10%, pair 25 / 13%,
  cross 78 / 57%, cross with the support test folded in (the current
  pass) 91 / 74%; to a fixpoint with the support filter: support 72 / 49%,
  Hall 81 / 60%, p1 98 / 91%, cross 100%. So cross support stays the
  cheapest strong test there. Its pass on the axis of the vector just
  placed first: +5% (prod). Choosing per child to run the support filter
  first where the other axis is not closed under it: not selective (that
  is the case for 98.5-99.6% of the children).

So above the (2,2) children the search is locally consistent, as far as
these rules can tell at a cost below the subtree they would save; what is
left is the cost of creating those children, not missing pruning.

## Micro-optimizations and the build (October 2026, round 2)

(Branch `opt2/micro-pgo`, commit 99d2aaf, from f61e719, where everything
below was measured. Integrated with a fix of the masked minimum (it could
drop subtrees when every unmatched cell had saturated counts), and with
opt2/children's padding instead of the loop split in FILTER_CARRY and the
lane mask kept in KEEP_CARRY: see "Integration of round 2" below.)

Per-node overheads on the carried path, and the build (pins, PGO, clang),
re-measured on f61e719 with `bench/prod.txt` as the main target.

How it was measured: the machine is shared (load 3-4.6 from other jobs),
and the minimum of alternating runs drifted by 5-10% between rounds, so
the comparisons below are paired: A and B started at the same time (so
that they see the same load), user CPU time from `wait4`, the ratio B/A
per round, median of 4-8 rounds (`bench/prod.txt`: 15 s per run). Rounds
agree within 1-2%, and the binary started first gets ~1% (both orders
were run where it mattered). There is no `perf` in the VM; profiles come
from a sampler (a `timer_create` signal at 2.5 kHz recording the RIP,
mapped with `addr2line -i`), which shows where the instructions are, not
mispredicts or latency. All variants visit the same nodes (checked on
`quick.txt` for every option set and build of the CMake tests) and find
the same squares.

Profile of f61e719 on `prod.txt` (native): FILTER_CARRY (inlined) 45%,
CROSS_TEST 17%, CROSS_AXIS with or_lanes8 15%, KEEP_CARRY 7%, TRY_CHILD's
own code 4%, COUNT_CARRY 3% (6% with -march=cascadelake), SUPPORT 3%,
SEARCH_REC's own code (choice of the cell, selection of the children)
3%, the root loop 0.01%.

Kept (median paired ratios on `prod.txt`, native / -march=cascadelake):

* The branching cell by a masked minimum over the byte counters (t =
  count - 1 for the unmatched cells, 255 elsewhere; vpminub across the
  words, a 64-byte horizontal minimum with phminposuw; then the first or,
  for saturated counts, the last label with that t) instead of scanning
  up to 254 classes (4 loads and compares per class, ~7 classes at the
  searched nodes of `full.txt`, where the minimum count is 1-20). In the
  final code, going back to the class scan costs +2.7% / +3.2% on
  `prod.txt` and +5% on `full.txt`; on f61e719 alone it was -1.6% / 0%.
  The sampler puts only ~1% of the time in the class scan, so most of the
  gain is the hard-to-predict exit of that loop on the way to the
  children.
* Splitting the filter loops (FILTER_CARRY, KEEP_CARRY, COUNT_CARRY with
  VBMI/GFNI, the selection of the children) into the full groups of 8
  entries without a lane mask and a last masked group: the per-iteration
  lane mask was ~8 scalar instructions plus masked loads in a ~35
  instruction loop. This is what gcc's profile-guided build did to the
  hottest loop (its FILTER_CARRY took 14% fewer samples). -4.3% / -5.2%
  (on top of the previous item). Padding the lists with empty sets
  instead (COUNT_CARRY, children) measured the same; a 64-byte padding
  store just before loads that partly overlap it costs a store-forwarding
  stall, so those loops are split too.
* |u & v| = 1 on one word: word w of v is rotated by r_w so that the
  rotated words are disjoint (always possible: n <= 8 labels in all, at
  most 16 of the 64 rotations collide), and c = OR_w rot(u_w & v_w, r_w)
  has |c| = |u & v|. For W = 2: and + vprolvq + vpternlogq + one test,
  instead of a popcount per word and an add (VPOPCNTDQ: 5 instead of 6
  vector instructions per 8 entries) or, without VPOPCNTDQ (Cascade Lake),
  instead of a power-of-two test per word plus a min test across words (6
  instead of 11). 1.003 native (no change) / -5.2% cascadelake.

Total, f61e719 -> this (same nodes: quick 1.77M, full 15.0M, prod 50.4M):

| | native, paired ratio | native, min wall | cascadelake, paired | cascadelake, min wall |
| --- | ---: | ---: | ---: | ---: |
| `quick.txt` | 0.920 | 0.302 -> 0.283 s | 0.887 | 0.339 -> 0.304 s |
| `full.txt` | 0.919 | 3.64 -> 3.39 s | 0.903 | 4.17 -> 3.72 s |
| `prod.txt` | 0.941 | 14.70 -> 13.58 s | 0.900 | 16.07 -> 14.66 s |

Profile after (native, `prod.txt`): FILTER_CARRY 46%, CROSS_TEST 18%,
CROSS_AXIS with or_lanes8 17%, TRY_CHILD 4%, CROSS 3%, KEEP_CARRY 3%,
SUPPORT 3%, SEARCH_REC 3% (the new choice of the cell 0.3%), COUNT_CARRY
2%. The filter loop is now bound
by the vector ports: per 8 entries 2 vpcompressq (2 uops each on port 5),
the test, 2 ors for the union.

Not kept:

* gcc PGO (gcc 13.3, -fprofile-use -fprofile-partial-training, LTO;
  trained on sums disjoint from quick/full/prod: `6 | 11 6 4 2 1 | 688`,
  `6 | 12 5 3 2 1 1 | 687`, `6 | 13 7 3 2 1 | 767`, `6 | 15 6 3 2 1 | 816`,
  `6 | 11 7 3 2 1 0 1 | 897`, `6 | 12 6 4 2 0 1 1 | 1104`,
  `6 | 14 6 3 2 1 0 1 | 1042`, `6 | 12 7 4 2 1 | 866` (N = 1636-1828,
  96-104 labels, 1-5 squares each) plus `6 | 16 5 4 2 | 1200` (133
  labels) and a 5x5 sum; expectations checked with the matrix build
  without support or cross filters). On f61e719 it gave -4.8% native
  on `prod.txt` (held out), mostly the loop split above. On the final code:
  native 0.997, cascadelake 0.984 on `prod.txt` (both orders, corrected
  for the start-order bias), 0.99-1.00 on `full.txt`, and the sampler sees
  the same total. The default build of this code is also 1% (native) and
  4.6% (cascadelake) faster than f61e719 with PGO. Not worth a two-pass
  build on the cluster, so no MAGIC_PGO option; the build there
  (build-sc.sh, -march=native) is unchanged.
* Pins: SEARCH_REC noinline gives the same binary as without it, with
  and without PGO: gcc 13 already keeps search_rec out of line with
  TRY_CHILD inlined (always_inline on the carried path since e9a0540),
  and the PGO build kept that layout too (it out-of-lines FILTER_CARRY and
  COUNT_CARRY at their cold call sites and inlines SUPPORT; forcing those
  back changed nothing). Forcing SUPPORT and KEEP_CARRY inline into
  TRY_CHILD in the default build: 1.0005 / 1.002.
* clang 18 (-O3 -flto): 6% / 7% slower than gcc 13 on the final code.
* -fno-stack-protector -fcf-protection=none (Ubuntu's hardening
  defaults, which put a canary check in search_rec, cross_axis, support):
  0.989 / 0.989, within the noise.
* The union of a filter pass from the loaded entries (masked or with
  keep) instead of from the compressed ones, in FILTER_STEP and
  CROSS_TEST: -1..-2% native, +0.4% cascadelake. Removing the two
  register copies of the accumulators that gcc emits in the split loop
  (same, or a precomputed loop bound): no change.
* Not retried: one OR reduction for both words of the unions (no change
  before, see above). The root loop (N^2 / 2 entries filtered per sum, a
  few ms at N = 3000), copying the cells and the cross-support sampling
  in TRY_CHILD are too small to measure.

## The cost of creating children (October 2026, round 2)

(Branch `opt2/children`, commits 9fa7162..8f5cb3a, from f61e719, where
everything below was measured. Integrated (see "Integration of round 2"
below): the padding of the lists with full sets for the exactly-once
filter, with the support passes keeping their lane mask as on this
branch (the loop split of the section above measured the same there,
within ~1%), the branch-free cells of cross support, and the
`-DCHILD_PROF` profile; the folded exactly-once test is the same as in
the section above.)

Where the time goes (f61e719, `bench/prod.txt`, `-march=native`, 15.0 s,
50.4M nodes; profile build `-DCHILD_PROF`, see arrange.c, by the layer of
the child, rows and cols placed): the (2,2) children of the (1,2)/(2,1)
nodes are 39.2M of the 50.4M nodes and two thirds of the cycles; each
(1,2)/(2,1) node has 8.2 of them (rows through its branching cell, or
cols), and 3.4% survive. Their lists come from parent lists of 187 (o, the
axis that must meet the new vector once) and 88 (b) entries, and are 43.5
and 36.2 entries after the filters. Of the (2,2) children, 1.3% die after
the o filter, 5% after the b filter, 50% in cross support, 36% in the
support filter, 4% in the "after" cross (see `cross_after`), 3.4% live.
Cross pass 0 (the o candidates against the unions of the b candidates
through the unmatched cells of the placed o vectors) shrinks the o list
from 47.6 to 20.5 entries and kills 27% of the children that reach it;
pass 1 shrinks the b list from 46 to 10 and kills 58% of the rest, mostly
by count. The (1,2)/(2,1) children (4.85M, from (1,1) nodes with 26.6 kids
each, lists of ~300) are ~11% of the time, the (2,3)/(3,2) children (6M,
nearly all dead after their b filter) ~2%.

The rdtsc pairs of the profile build cost ~50% more time and inflate the
short phases, so the shares below come from builds that stop the (2,2)
children after a stage (wrong results, same work up to that stage; cross
support always before the support filter), min of 3 alternating runs in
one session:

| (2,2) children stopped | f61e719 | this branch |
| --- | ---: | ---: |
| before creating them (the other layers) | 2.10 s | 1.85 s |
| at once (the kids loop) | 2.27 s (+0.17) | 2.05 s (+0.20) |
| after the o filter | 6.04 s (+3.77) | 5.27 s (+3.22) |
| after the b filter | 7.58 s (+1.54) | 6.58 s (+1.31) |
| after cross support | 13.58 s (+6.00) | 12.65 s (+6.07) |
| not stopped | 14.90 s (+1.32) | 13.50 s (+0.85) |

So per (2,2) child: o filter ~95 ns (187 entries; 82 ns now), b filter
~40 ns, cross support ~165 ns per call (~430 cycles for both passes: four
loops of 3-6 blocks of 8 entries, ~100 cycles of setup and reductions),
support filter and the survivors' subtrees ~30 ns. A replay harness (the
(1,2)/(2,1) nodes of prod dumped with their kids, 1 in 150, and TRY_CHILD
replayed on them without the recursion) and a cross harness (the (2,2)
states at the cross call, 1 in 200) gave per-phase cycles within ~2%,
which the end-to-end runs (3-5% noise on this shared machine) cannot; but
they replay states out of order (siblings are not consecutive), so they
overstate the cost of branches that siblings predict.

Kept:

* The exactly-once filter on folded words (FILTER_CARRY): the words of
  u & v are or-ed into one, word w rotated left by r_w with the rotations
  chosen so that the rotated words of v are disjoint (v has at most 8
  labels, so they exist; nearly always r = 0 or 1), which keeps the count:
  one vpopcntq and compare, or without VPOPCNTDQ one "single bit" test
  (a != 0, a & (a - 1) == 0), instead of a popcount per word and their
  sum, or the and/sub/ternlog/min chain per word. And no lane masks: a
  node pads its two lists with 8 full sets (PAD_LIST) before creating its
  children, which both filters drop (n >= 3). Microbenchmark, 187 entries
  of 6 of 100 labels, v with 4 labels, TSC cycles per call: Sapphire
  Rapids 171 -> 142 (padding) -> 136 (fold), Cascade Lake build 218 ->
  179 -> 137.
* Cross support: the word and bit of each cell from one pdep per word,
  without a branch per cell on the split of the cells between the two
  words: 421 -> 395 cycles per call in the cross harness (native), 421 ->
  389 (cascadelake build), but in place, where the siblings share most of
  their cells and the branch was mostly predicted, prod 13.86 -> 13.86 s
  native, 14.39 -> 14.15 s cascadelake (min of 4).

Together (min of 3 alternating runs; same nodes): `bench/prod.txt`
14.69 -> 13.72 s native (-6.6%), 16.41 -> 14.18 s cascadelake build
(-13.6%); `full.txt` 3.70 -> 3.41 s and 4.13 -> 3.54 s; `quick.txt`
0.307 -> 0.294 s and 0.355 -> 0.303 s. The cascadelake build gains more:
without VPOPCNTDQ the per-word test was 11 vector ops per 8 entries, now 6.

Not kept (replay = cycles per (2,2) child in the replay harness, before /
after; prod in seconds, min of alternating runs):

* Node-local renumbering (the children's lists as bitmasks over the
  parent's list positions, from a kids x list relation computed once per
  node). Per 8 entries the o filter is ~5 ops of test and ~8 of compress,
  store and union; the relation needs the same tests (batched over the
  kids it shares only the loads), and the masks must still be compressed
  into lists for cross support and the support filter, which read the
  carried words. How many kids of a (1,2) node keep an entry of its o
  list: none 12%, 1: 25%, 2: 26%, 3: 19%, 4: 10%, 5 or more 7% (each kid
  23%), so the kids' lists hold ~1.9x the parent's entries and a scatter
  would write about as much as the compresses do. One level up, a full
  local relation at the (1,1) nodes (rows x rows, cols x cols, rows x cols:
  ~180K pair tests per node, ~20% of a (1,1) subtree's time) would turn
  both filters into ANDs, but cross support on position masks tests one
  entry per instruction (a 512-bit vptestmq per cell) instead of 8, 3-4x
  the cost of the current test, and the unions the forward check and the
  support filter need take a test per label. Not implemented.
* Cheaper tests before building a child's lists: cross support of the kids
  themselves at their parent (unions of the parent's other list through
  the unmatched cells of the placed vectors of the kids' axis): 0.4% of the
  kids. Cross support of the (1,2)/(2,1) node's own lists: 2.9% of the
  nodes, o list 165 -> 139 entries, b list 79 -> 78.5 (less than it costs,
  as found on the matrix path). A perfect matching of the b candidates on
  the free cells of the two placed vectors of the other axis, and of the o
  candidates likewise (4 x 4, from co-occurrence): 0.7% of the (2,2)
  children at the cross call.
* Reordering: one support pass on the o list before cross support
  (replay 643 -> 660); the first cross pass on the longer list (627), the
  shorter one (667) or always b (656) instead of always o (619), with 2%,
  28% and 30% more children surviving; cross support with half the cells (the cells of one
  placed vector) in either or both passes: 1.7-3.5x more children survive
  cross support (8.6K -> 14.8-30K of 263K kids) and 5-24% slower. In pass 1,
  the cell with the fewest b candidates loses all of them in 32% of the
  calls (54% of its kills), so testing that class first would save ~12
  cycles per call: not tried.
* vpermq with indices from a 256-entry table instead of vpcompressq (2 uops
  on port 5): -12% in the filter microbenchmark, but in place: replay
  620.6 -> 616.5 (native), 625.7 -> 634.5 (cascadelake), prod (cascadelake)
  14.17 -> 14.23 s; with an 8-byte table and vpmovzxbq: slower; in the
  cross test: 400 -> 393 / 390 -> 393 (cross harness).
* No lane masks in KEEP_CARRY (pad with full sets at entry): replay 626 ->
  643 / 632 -> 635. Choosing the fold rotation without a branch (the
  smallest of 0-3 by cmov, loop beyond): replay 643 -> 657 / 630 -> 651.
* MRV by a minimum over the byte counters (count - 1, wrapping, 255 off
  the unmatched cells; then one compare per word) instead of scanning the
  classes 1, 2, ... (one compare per word per class, 5-20 classes at the
  (1,2) nodes): same nodes, `full.txt` 3.42 -> 3.41 s user, no measurable
  gain.
* Without GFNI, counting the labels of long lists with few labels to count
  by one test per label per 8 entries instead of a masked add per entry:
  167 cols and 4 labels (the free cells of the placed row at a (1,2)
  child) 249 -> 141 cycles in a microbenchmark, but 81 rows and 10 labels
  125 -> 164, and much slower on short lists; for lists of 64 or more and
  at most 4 labels per word, prod (cascadelake) 14.47 -> 14.47 s.
* Forcing SUPPORT inline into TRY_CHILD: prod -2.5% native, +1.1%
  cascadelake (noise ~2%); forcing CROSS inline: `full.txt` +1%.
* Store forwarding: the lists are read right after the filters write them
  with unaligned 8-lane stores; rewriting them that way or with aligned
  stores before timing cross support in the cross harness: no difference.

## A residual solver near the leaves (October 2026, branch opt2/residual)

(Commits d574254 and 4690dc0, from f61e719, where everything below was
measured. The label -> position masks kept on that branch are not
integrated: they replace the exactly-once filter, which the folded test
of the two sections above made about twice as cheap without VPOPCNTDQ,
and on top of it they cost time with -march=cascadelake (paired runs,
median ratio of user CPU time with / without the masks):
`bench/prod.txt` 1.022 from 8 children (the branch's threshold), 1.002
from 16, 0.995 from 32; `full.txt` 1.030 from 8, 1.010 from 16. See
"Integration of round 2" below.)

The idea: at a node with r rows and c cols placed, the rest is a small
exact-cover problem (the n - r rows are pairwise disjoint candidate rows,
one through each col-unmatched cell of a placed col, the n - c cols
likewise, every remaining row meets every remaining col exactly once), so a
specialized solver with node-local bitmasks over list positions might
finish it faster than the generic search. It does not: below the nodes
where it applies there is almost nothing left to save, and where the time
is, its setup costs more than the generic filters. What came out of it is a
use of those masks one level up (label -> position masks, kept, see the
end).

Where the nodes and the time are (`bench/prod.txt`, f61e719, native build,
rdtsc at every change of layer, so the cycle shares include ~15% of timer
overhead; "(r,c)" = rows and cols placed in the child):

| children | created | searched | cycles creating them | lists at the searched nodes |
| --- | ---: | ---: | ---: | --- |
| (1,1) | 182k | 182k | 0.8% | 250 rows, 234 cols, 26.8 kids |
| (1,2) | 3.41M | 3.37M | 7.1% | 82 rows, 168 cols, 8.5 kids |
| (2,1) | 1.44M | 1.41M | 2.8% | 159 rows, 74 cols, 7.6 kids |
| (2,2) | 39.24M | 1.32M | 75% | 34 rows, 34 cols, 4.5 kids |
| (2,3), (3,2) | 6.01M | 1637 | 2.3% | 3-4 entries |
| (3,3) and deeper | ~2k | ~600 | 0.0% | |

The searched nodes themselves (choosing the cell and the kids) take another
~12% (6.5% at (1,2)). So 78% of the nodes are the (2,2) children and
creating them is three quarters of the time; per child (cycles with the
timers, and the share killed at that stage): exactly-once filter of the
other axis 245 (1.3%), disjoint filter of the same axis 129 (5.1%), cross
support 331 (50.3%), support filter 124 (35.7%), cross support after it 57
(4.3%); 3.4% survive. On `bench/full.txt` 0.4% survive (51k of 12.6M).

* Room below the searched (2,2) nodes: dropping them (wrongly, after all
  their filters, before their counts) takes prod from 14.77 to 14.24 s and
  15.02 to 14.58 s (3.0-3.6%); on `full.txt` ~0.5%. Their subtree costs
  714-719 TSC cycles per searched (2,2) node in the native build, 933-1025
  built with -march=cascadelake (rdtsc around counts + SEARCH_REC). Of the
  1.32M, at most 55 (the squares) have a completion; the rest die one
  level down (6.0M children, 1637 searched), 40% by the forward check after
  the first filter, 57% after the second (with the support filter's first
  pass), 3% by the support filter.
* Residual solvers at the searched (2,2) nodes (lists of at most 64
  entries, masks over the positions, a pass over a list per chosen vector):
  (a) rows first, an exact cover of the col-unmatched cells by disjoint
  rows (branching on the uncovered cell with the fewest rows left), keeping
  the mask of the cols that meet every chosen row once, then the cols as an
  exact cover of the row-unmatched cells by disjoint cols of that mask:
  24.8M solver nodes instead of 6.0M children, prod ~10% slower. (b) Both
  axes, branching on the uncovered cell of the 16 unmatched cells of the
  node with the fewest candidates, forward checking on those cells only:
  11.5M solver nodes, 4-6% slower. (c) As (b), with the forward check of
  TRY_CHILD on all unmatched cells (the unions of the kept entries) and the
  support filter's first pass (the same axis within the other axis' cells
  and union): 6.15M solver nodes (the generic search: 6.0M + 1637), and the
  same time: prod CPU time 14.95 vs 15.05 s (min of 3), 15.08 vs 15.28
  (median), 730-745 vs 714-719 TSC cycles per searched (2,2) node native,
  1023-1068 vs 933-1025 cascadelake; `full.txt` 3.72 vs 3.68 s. The setup
  (masks of the rows / cols through each of the 16 unmatched cells: a test
  per cell and 8 entries) and two passes per solver node cost what the
  generic search spends on its counts and children, whose lists are just
  as short. With pairwise tables instead of passes (rows disjoint from each
  row, cols meeting each row once, ~34 x 68 pair tests) the setup alone
  would be ~1000 cycles, more than the whole generic subtree.
* At the (2,2) children (75% of the time), an exact solve would need those
  pairwise tables on lists of ~47 entries (~2200 pair tests, over 1000
  cycles) where cross support and the support filter cost ~430 and leave
  3.4%. At the (1,2) nodes, the 5 remaining rows from ~82 candidates (5
  classes of ~16, pairwise disjoint with p ~ 0.5) give ~16^5 p^10 ~ 1000
  row completions, each still to be matched with the cols, against ~6500
  cycles for all the children of such a node in the generic search. Not
  tried.

Label -> position masks (kept, `LABEL_MASKS` in arrange_core.h). The
exactly-once filter of the other axis' list, done for every child (at a
(1,2) node: 8.5 children x 168 cols, a popcount test per entry and word),
only depends on the child's labels outside the placed vectors of that
axis (4 of them), so with the mask of the list entries containing each such
label (built once per node) the entries meeting a child exactly once are
those in exactly one of its 4 masks, and the child only compresses them
(FILTER_KEEP). In a microbenchmark (168 entries, 2 words) the exactly-once
filter takes ~150-230 TSC cycles native and ~200 built with
-march=cascadelake, the compress with a known mask ~90-95. Building the
masks is a 64 x 64 bit transpose per 64 entries and word of labels.

* First version: byte transpose of each 8 entries (vpermb, or vpshufb +
  vpermw without VBMI), an 8 x 8 qword transpose across the registers, a
  vptestmb per label, per-child scalar loops over the 64-entry groups. On
  the cascadelake build: `full.txt` -3% from 4 children, -5% from 12; prod
  16.57 -> 16.39 (8), 16.04 (12), 16.18 s (20). In situ (prod, with
  timers), the filter went from 354 to 262 cycles per (2,2) child but the
  masks cost ~650 cycles per (1,2) node. Native: no gain from 12
  children, 3-5% slower from 4-8.
* Final version: vpshufb (byte j of the 2 entries of a lane into word j)
  then an 8 x 8 transpose of words with three rounds of unpacks, no vpermw
  (32 single-uop shuffles per 64 entries and word); the bits come out
  permuted (entry 8 b + 2 L + t at bit 16 L + 2 b + t), which a pext per 8
  entries undoes in FILTER_KEEP; the per-child exactly-one-of with the 8
  words of a label's masks in one register (2 operations per label).
  `LABEL_MASKS` costs 515-520 TSC cycles per call in situ (307 in a loop),
  on prod (cascadelake) 2.76M calls with 11.9 children and 206 entries on
  average (1.90M of the 3.37M (1,2) nodes, 0.65M of the 1.41M (2,1), 0.17M
  (1,1)): ~4% of the time to build masks, ~9% saved in the filter.
  Applying the masks with a loop over the bytes of the needed labels and a
  doubling chain + vpmovb2m instead of vptestmb (to take the tests off
  port 5): 3.5x slower in a microbenchmark (branches).
* Threshold (children per node; lists of 64-512 entries), cascadelake,
  CPU time, min of 3: prod 16.13 s without, 15.57 (4), 15.66 (6), 15.53
  (8); `full.txt` 4.11 without, 3.94 (4), 3.97 (6), 4.00 (8), 3.98 (12).
  Native (VPOPCNTDQ): prod 14.74 without, 15.29 (4), 14.91 (8): the masks
  are only used without VPOPCNTDQ (Skylake-X, Cascade Lake), from 8
  children.
* With the masks' child loop inlined into SEARCH_REC (two copies of
  TRY_CHILD there) the native build was 1.7% slower on `full.txt` even with
  the masks unused, so that loop is a separate function.

Result (f61e719 -> this, same nodes; wall time, min / median of alternating
runs on the shared machine, `bench` built with -march=cascadelake -flto as
for production, and the native build):

| | cascadelake | native |
| --- | ---: | ---: |
| `quick.txt` | 0.347 / 0.349 -> 0.339 / 0.344 s | 0.306 / 0.314 -> 0.301 / 0.306 s |
| `full.txt` | 4.14 / 4.21 -> 3.95 / 4.05 s | 3.70 / 3.75 -> 3.67 / 3.71 s |
| `prod.txt` | 16.19 / 16.36 -> 15.64 / 15.82 s | 14.98 / 15.09 -> 14.89 / 14.89 s |

In CPU time (thread clock, 3 alternating runs), cascadelake: prod 16.45 /
16.49 -> 15.70 / 15.73 s (-4.5%), `full.txt` 4.13 / 4.26 -> 4.00 / 4.05 s;
n = 7 (P = 12 6 3 2 1, S = 352-370) 8.17 / 8.37 -> 8.13 / 8.18 s; 4 words
of labels (`--min-words 4`, P = 14 7 5 3, S = 1460) 4.23 / 4.35 -> 4.09 /
4.10 s.

## Integration of round 2 (October 2026)

Branch `opt2/integrated`, from 23a34f7 (whose search code is f61e719's):
the three round-2 branches above combined piece by piece, where they
overlap by measurement. Method as in "Micro-optimizations and the build":
paired runs (both binaries started together, user CPU time from `wait4`,
median of the ratios of 6-8 rounds, start order alternating) on
`bench/prod.txt`, plus minima of alternating runs; the native build and
-march=cascadelake (the cluster), both on the shared Sapphire Rapids
(load 0-6 from other jobs; the machine restarted twice during the work, so
some comparisons were repeated in quieter sessions). Every build visits
the same nodes (quick 1,770,779; full 14,958,507; prod 50,375,738) and
finds the same squares, also with every option set of the CMake tests.

What went in:

* opt2/micro-pgo: the branching cell by one masked minimum over the byte
  counters, the exactly-once test on folded words, and the loop split
  (full groups of 8 without a lane mask, then a masked last group) for the
  label counts and the selection of the children. Its masked minimum had a
  soundness bug, fixed here: the masked subtraction took the running
  minimum as its source, so for words w >= 1 the cells that are not
  unmatched carried earlier words' values, and where every unmatched cell
  had >= 255 candidates (saturated counters) the scan for the largest
  saturated label could pick such a cell, which has no children: the
  subtree was dropped (squares lost). Not seen on the benchmarks or in
  production (the top labels are rare there), but possible for 6x6 sums
  with N ~ 11k and 4 words of labels, and for 7x7. `fuzz_arrange --mode 6`
  builds such families (an r1 whose numbers are the rarest, each in >= 262
  vectors meeting r1 only there, a planted square through r1): the buggy
  code fails every seed tried (the default, w3 and w4 variants find no
  square), the fix and f61e719 pass; ctest runs it with the carried lists
  and with the matrices.
* opt2/children: the exactly-once filter (FILTER_CARRY) on lists padded
  with 8 full sets (a precondition now in its contract, checked with
  -DARRANGE_DEBUG) instead of the split, the support passes (KEEP_CARRY)
  with their lane mask per iteration as on that branch, the branch-free
  cells of cross support, and the `-DCHILD_PROF` profile.
* Not opt2/residual's label -> position masks (slower on top of the folded
  test, below), nor PGO (see "Micro-optimizations and the build").

By piece (`bench/prod.txt`, paired median ratio B/A of user CPU time, B
the variant; native / cascadelake):

| A -> B | native | cascadelake |
| --- | ---: | ---: |
| 23a34f7 -> opt2/micro-pgo with the fix | 0.929 | 0.910 |
| that -> + branch-free cross cells | 0.981 | 0.985 |
| that -> FILTER_CARRY on padded lists instead of the split | 0.982 | 0.997 |
| split everywhere -> KEEP_CARRY with a lane mask | | 0.994 |
| split everywhere -> the selection of the children with a lane mask | | 1.011 |
| split everywhere -> the final code (padding, KEEP_CARRY masked) | 0.997, 1.000 | 0.979, 0.993 |
| split everywhere -> the cell by scanning the classes | 1.039 | 1.034 |
| the final code -> the cell by scanning the classes | 0.995 | 1.007 |
| the final code -> + label masks from 8 / 16 / 32 children | (not used) | 1.022 / 1.002 / 0.995 |
| opt2/children alone -> the final code | 1.000, 0.995 | 1.027, 1.000, 1.009 |
| opt2/children + the masked minimum -> the final code | 1.015 | 1.006 |

(two or three values: separate sessions; label masks on `full.txt`: 1.030
from 8 children, 1.010 from 16). So the variants of the lane-mask handling
and the cell choice are all within 1-3% of each other once the folded test
is in, at the noise floor of these runs: the masked minimum gained 3.5%
where the filters used the split, ~0 where FILTER_CARRY uses the padding;
opt2/children alone is as fast as the final code natively and 0-3% faster
with -march=cascadelake in different sessions (opt2/children plus the
masked minimum: the same as the final code). The label masks lost their
use: they replace the exactly-once filter, which the folded test made
about twice as cheap without VPOPCNTDQ.

Result (23a34f7 -> opt2/integrated; minimum of 3-4 alternating runs, bench
search time, and the paired ratio):

| | nodes | -march=native | paired | -march=cascadelake | paired |
| --- | ---: | ---: | ---: | ---: | ---: |
| `quick.txt` (`--repeat 10`) | 1.77M | 0.296 -> 0.269 s | 0.908 | 0.337 -> 0.293 s | 0.846 |
| `full.txt` (`--repeat 3`) | 14.96M | 3.65 -> 3.29 s | 0.918 | 4.02 -> 3.47 s | 0.864 |
| `prod.txt` | 50.38M | 14.65 -> 13.38 s | 0.920 | 16.02 -> 14.23 s | 0.882 |

AddressSanitizer: GCC 13 with AVX-512 crashes in sanitizer builds with
use-after-return detection on (its default at run time): a function using
zmm registers realigns its frame to 64 bytes and GCC then accesses local
arrays at 64-aligned offsets with aligned moves (vmovdqa64, also for
`_mm512_loadu_si512`), while the fake stack of that detection keeps them
only 32-aligned. A ten-line function with a local `uint64_t[8]` passed to
`_mm512_loadu_si512` reproduces it, and no attribute on the array helps, so
it is not the code's over-alignment: the Debug build (CMakeLists.txt)
leaves that detection out of GCC's instrumentation
(`--param=asan-use-after-return=0`, the same as
`ASAN_OPTIONS=detect_stack_use_after_return=0`).

## Carried bitsets up to 512 labels, each r1 at its own width (October 2026, branch cx/wide)

Part of the "cost per square at large N" study. The profile of that study
(cx/profile) found two implementation steps on top of the algorithm: above
256 labels the search fell back to the N x N matrices (8-slice MRV, no
cross support: x2.7 time and x2.8 nodes on the same r1 at N = 23k), and
below, each extra 64-bit word of label bitset cost x1.38 time at the same
nodes. 13 7 4 3 1 1 passes 256 labels at N ~ 25k, and the existence study's
"local slope 5.2 above N = 2x10^4" was that switch.

What changed (src/c/arrange.c, arrange_core.h):

* The carried path runs up to CARRY_MAX_W = 8 words (512 labels), with W =
  ceil(L / 64) exactly (instantiations for 5, 6 and 7 words; the matrix
  path keeps 2, 3, 4, 8, 16, 64 and its own counters, so `-DCARRY_MAX_W=0`
  builds visit the same nodes as before). The carried code was already
  generic in W (the rotations of the folded exactly-once test exist for any
  W since v has <= 8 labels; the bit transpose of COUNT_CARRY, KEEP_CARRY,
  PAD_LIST and the label reconstruction from placedw work per word); the
  unions and placed bitsets now hold 8 words. The one W-specific cost was
  cross support's 8 W accumulators (40 zmm at W = 5): it now builds the U_y
  4 words at a time, two passes over the other list for 5-8 words (1.6%
  less time on the min, 7% on the median, on two 5-word r1; the same code
  for W <= 4).
* Each r1 is searched with only the words it needs (`opts.r1_width`, on by
  default; `bench --no-r1-width` turns it off). The vectors are sorted by
  their labels, descending, and the subproblem of r1 has only the vectors
  after it, whose labels are all <= x(r1), so the root loop runs in runs of
  equal width ceil((x + 1) / 64) (at least 2, at least --min-words); each
  run fills the depth-0 lists from its first r1 on with the first W words
  of the bitsets (no repacking: `bits` keeps the widest stride). The state
  of the adaptive cross support carries over from run to run, so the nodes
  are exactly those of a single width (no node difference anywhere, not
  only up to the cross adaptivity). The runs of the large sums: S = 2650
  (279 labels) has r1 0-2233 at 5 words, then 4 words up to 20113, 3 up to
  29591, 2 to the end; S = 2700 (259 labels) has only 34 r1 at 5 words;
  S = 2200 (211 labels) only r1 < 1273 at 4 words.
* `-DR1_SAMPLE`: the root loop runs only the r1 of env SAMPLE_STRIDE /
  SAMPLE_OFFSET, or SAMPLE_LIST, by global index (whatever their width), and
  SAMPLE_LOG gets one line per r1 (r1, x, W, squares, nodes, thread CPU
  seconds) and the runs. The same patch on fdb77fc gave the base binary.

**Cost per word** (`full.txt`, `--min-words`, same 14,958,507 nodes,
thread CPU, min of 2 alternating runs): 2-8 words cost 3.29, 4.78, 6.51,
8.01, 9.65, 11.29, 13.09 s, i.e. 1.00, 1.45, 1.98, 2.43, 2.93, 3.43, 3.98:
the time is within 3% of proportional to the words (+0.5 of the 2-word
time per word). Nearly all of the per-node work is moving and testing
carried words, which is why the profile's per-word factor was x1.38.

**Paired measurements** (base fdb77fc and new built alike, -O3
-march=native -flto; the same r1 in both, thread CPU per r1 from the logs,
min of 2 alternating runs per r1; shared 4-core Xeon at load 3-10, run to
run spread up to 5-12%). Squares (and their hashes, rechecked on the r1
with squares) identical everywhere; nodes identical wherever both binaries
are on the carried path.

| sum | N | labels | base path | sample | base CPU-s | new CPU-s | base / new | nodes base / new |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: |
| 13 7 4 3 1 1, S = 2650 | 31,743 | 279 | matrices (8 words) | 11 r1: 364, 1564, 2199, 3364 (heavy, < 5,164) + 7 at stride 4096 | 64.66 | 22.80 | **2.84** | 3.33 |
| (its 3 r1 at 5 words / 5 at 4 words) | | | | | 26.51 / 36.85 | 11.78 / 10.56 | 2.25 / 3.49 | |
| 13 7 4 3 1 1, S = 2500 | 26,585 | 260 | matrices | stride 5000 (6 r1) | 37.08 | 11.84 | **3.13** | 2.87 |
| 9 6 4 3 1 1 1 1, S = 2700 | 20,896 | 259 | matrices | stride 800 (26 r1) | 53.28 | 18.58 | **2.87** | 5.14 |
| d-first V_d at 13 7 4 3 1 1, S = 2650 (`dsample`, \|V_d\| ~ 4.5k) | 4,410-4,557 | 274-279 | matrices | 16 d at stride 2000 | 54.78 | 25.18 | **2.18** (1.93-2.44 per d) | 1.21 |
| 13 7 4 3 1 1, S = 2400 | 22,992 | 245 | carried, 4 | stride 1600 (15) | 23.06 | 23.04 | 1.00 | = |
| 12 9 6 2 1 1, S = 3500 | 17,142 | 256 | 4 | stride 1000 (17) | 6.36 | 6.20 | 1.02 | = |
| 13 7 4 3 1 1, S = 2200 | 15,199 | 211 | 4 | stride 400 (38) | 27.46 | 22.41 | **1.23** | = |
| 12 6 3 2 1 1, S = 1200 | 11,697 | 199 | 4 | stride 100 (117) | 22.03 | 16.79 | **1.31** | = |
| 12 6 3 2 1 1, S = 1080 | 8,891 | 179 | 3 | stride 40 (223) | 15.28 | 15.20 | 1.005 | = |
| 14 7 4 4 1 0 0 1, S = 3231 | 7,901 | 183 | 3 | stride 100 (79) | 3.95 | 3.73 | 1.06 | = |
| 12 6 3 2 1 1, S = 988 | 6,671 | 159 | 3 | stride 100 (67) | 3.04 | 2.79 | 1.09 | = |
| `bench/prod.txt` | 1,491-2,994 | 95-121 | 2 | all, 3 alternations | 13.64 | 13.59 | 1.004 | = |

* Above 256 labels the gain is 2.8-3.1x on the plain search (the r1 at 5
  words 2.2-2.3x, the r1 at 4 words or less 2.8-4.2x) and 2.2x inside the
  V_d searches, where more of the work is in 5-word r1. Most of it is the
  algorithm, not the width: 2.9-5.1x fewer nodes in the plain search
  (1.2x in V_d) with cross support and exact MRV.
* Below 256 labels the gain is the width ratio of the r1 that hold the
  work: 4/3 where the heavy r1 need one word less than all the labels
  (S = 2200: the r1 at 3 words 1.31x; b1200 1.32x), 3/2 for the r1 that
  drop to 2 words (1.4-1.8x, but they hold little time), and nothing where
  the heavy r1 need all the words (S = 2400: x(r1) >= 192 for r1 < 8,898).
* Nothing got slower: at the same width the binaries are within noise
  (the 4-word r1 of S = 2400 one r1 per process, 5 interleaved rounds:
  base / new 0.996 in total, -2.4% to +4% per r1 with no consistent sign;
  prod.txt base / new 1.004 in total, within +-3% per instance).

**Time per sum and cost per square** (stride x sampled CPU; squares from
the cx/profile study and the existence study, model.md):

| sum | N | base CPU-s | new CPU-s | squares | base s / square | new s / square |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 12 6 3 2 1 1, S = 988 | 6,671 | 304 | 279 | 245 | 1.24 | 1.14 |
| 14 7 4 4 1 0 0 1, S = 3231 | 7,901 | 395 | 373 | 128 | 3.1 | 2.9 |
| 12 6 3 2 1 1, S = 1080 | 8,891 | 611 | 608 | 275 | 2.2 | 2.2 |
| 12 6 3 2 1 1, S = 1200 | 11,697 | 2,200 | 1,680 | 384 | 5.7 | 4.4 |
| 13 7 4 3 1 1, S = 2200 | 15,199 | 11,000 | 9,000 | 2,533 | 4.3 | 3.5 |
| 9 6 4 3 1 1 1 1, S = 2700 | 20,896 | 42,600 | 14,900 | 1,600 +- 780 | 27 | 9.3 |
| 13 7 4 3 1 1, S = 2400 | 22,992 | 36,900 | 36,900 | 2,745 | 13.4 | 13.4 |
| 13 7 4 3 1 1, S = 2500 | 26,585 | 185,000 (6 r1) | 59,000 | | | |
| 13 7 4 3 1 1, S = 2650 | 31,743 | ~310,000-340,000 | ~110,000-120,000 | 4,600 +- 3,400 | ~70 | ~25 |

(S = 2650: the existence study's stratified estimate, 342-375k CPU-s on
its older copy of the code, times 0.9, the fdb77fc / older ratio on the
three r1 both timed; new = that / 2.84.) Along 13 7 4 3 1 1 the time per
sum now grows roughly like N^3-3.4 from S = 2400 to 2500 and 2650 (rough:
those samples have 6-15 r1, and the profile's 29-r1 estimate of S = 2400 is
42k CPU-s against 37k here), where it was N^5.2 above 2x10^4 with the
switch to the matrices; so the cost per square at N = 25-32k is ~2.8x below
the existence study's figures. The
cost per square still grows with N: this removes an implementation step,
not the N^1.5-per-(1,2)-node growth of the algorithm (see cx/profile).

What was not done, and why:

* Narrowing within an r1 (a subtree whose lists and unmatched cells no
  longer touch the top word could run at one word less; the word-major
  layout makes the narrower lists free): in the heavy r1 of S = 2400 ~20%
  of the root's candidates have a label in the top word, so it would only
  apply near the leaves. Relabeling the labels of a node's lists into fewer
  words would apply at the (2,2) nodes (lists of ~100-150 entries over
  ~150-200 labels: still 3-4 words), not enough for the relabeling cost.
* The matrix path stays for more than 512 labels (N well above 30k for the
  P studied) and for CPUs without AVX-512BW.

Correctness: quick/full/prod identical nodes and hashes with per-r1 width
on and off and with `--min-words` 4-8 (the carried search is width-
agnostic); the matrices build (`bench_matrices`) node-identical to fdb77fc's
on quick/full, also with `--min-words 8` and `--no-cross`; `fuzz_arrange`
has a mode 7 (planted grids over 257-500 labels; squares land in r1 runs of
every width 2-8) and variants w5, w6, w7, fixw, w5_fixw, w6_nomrv,
xcross_w5: 0 fails on 1,000 fresh default seeds and 1,000 fresh mode-7
seeds (all variants), 300 mode-7 seeds each with the matrices and without
GFNI / VPOPCNTDQ, `--mode 6` (saturated counters) with w3-w8 and fixw, and
the ASan debug build on 120 seeds; `ctest -R fast_` passes (new tests:
`fast_bench_quick_wide5/7`, `wide6_fixed`, `matrices_wide8`,
`fast_fuzz_arrange_wide`, `fast_fuzz_arrange_matrices_wide`; the gather and
plain-C cross tests at 8 words now run on `bench_matrices`, since the
native build carries 8 words).
## The pretest of the deep children (October 2026, branch cx/pretest)

At N >= 11k the (2,2) children that survive cross support and the support
filter are 30-40% of them, each creates ~15-18 (2,3)/(3,2) children, and
nearly all of those die at the count or forward check right after their
two filters (85% of all nodes at N = 15-23k, 26-32% of the cycles, the
fastest-growing term of the plain search; see the cx/profile study). The
pretest (`opts.pretest_min`, `--pretest-min K` in bench and msearch, 0 =
off; default 5, i.e. the (2,3)/(3,2) children on; carried lists only) runs
those two filters first as FILTER_COUNT: the same tests, but per 8 entries
a masked or into the union and a byte of keep mask instead of a compress,
a store and an or per word. The count and forward checks follow on those
counts and unions, so exactly the same children die; only the ones that
pass get their lists, compressed from the keep masks (PT_APPLY, not
inlined), and go on with cross support / the support filter as before.
Nodes and squares are unchanged.

The cx/prune version (3dc9b05, `-DPRETEST`) let the survivors run the two
filters again; compressing from the keep masks instead was 3.7% faster at
a2400 (27.5 vs 28.5 s, same session) and makes the pretest nearly free
where it kills little.

Measurements (fdb77fc vs this branch, `-O3 -march=native -flto`, both
with the R1_SAMPLE patch of cx/wide: the same r1 sample, two runs of each
alternating, per-r1 minimum of thread CPU time, summed; the machine was
shared, load 2-7 on 4 cores):

| P / S | N | r1 sampled | base CPU-s | pretest CPU-s | speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| 13 7 4 3 1 1 / 1850 | 2,523 | all | 2.91 | 2.88 | 1.009 |
| 13 7 4 3 1 1 / 2000 | 7,593 | 190 (stride 40) | 13.29 | 12.50 | 1.063 |
| 12 6 3 2 1 1 / 1080 | 8,891 | 111 (stride 80) | 9.10 | 8.60 | 1.058 |
| 13 7 4 3 1 1 / 2200 | 15,199 | 19 (stride 800) | 14.37 | 12.98 | 1.107 |
| 11 6 4 3 2 1 / 2174 | 16,424 | 19 (stride 900) | 14.49 | 12.93 | 1.120 |
| 14 7 4 4 1 0 0 1 / 3648 | 20,538 | 12 (stride 1800) | 9.19 | 8.96 | 1.026 |
| 13 7 4 3 1 1 / 2400 | 22,992 | 15 (stride 1600) | 32.14 | 27.46 | 1.170 |

* Geometric mean over the four sums with N >= 15k: 1.105. The gain
  follows the (2,2) survival: 34-37% at S = 2200-2400, 28% at 2174, 11%
  at 3648 (where the dead (2,3)/(3,2) layer is 64% of the nodes, not 85%).
  The cx/prune port measured 1.05-1.19 on the same samples (one loaded
  session gave 1.10 at a2400, another 1.19: ~5% session noise).
* Every pair had the same squares and the same nodes at every sampled r1.
  So the cost per square falls by the same factors: with the cx/profile
  estimates, 4.0 -> 3.6 CPU-s per square at S = 2200, 10.8 -> 9.6 at 2174,
  15.5 -> 13.2 at 2400. It is a constant factor that grows slowly with N,
  not a change of the exponent.
* `bench/prod.txt` (min of 4 alternating runs) 13.85 -> 13.49 s (-2.6%),
  `full.txt` (`--repeat 3`) 3.359 -> 3.333 s (-0.8%), `quick.txt`
  (`--repeat 10`, 5 rounds) 0.2719 -> 0.2753 s (+1.3%: its two smallest
  instances, 12-50 ms, by 3-6%; the same binary with `--pretest-min 0`
  0.2733 s).
* `--pretest-min 4` (also the (2,2) children) in the plain search: 0.97x
  base at S = 2200. Few (2,2) children die at those checks there.

The diagonal-first search (V_d, cx/dfirst `--diag-first`). Its (2,2)
children die at those checks 88% of the time, so a lower threshold could
pay there. Death stages (S = 2200, 60 d, `-DCHILD_PROF`, pretest off):

| children | made | lists in (o / b) | count o | fc o | count b | fc b | cross | support | live |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| (1,2) | 2.77M | 429 / 302 | 0% | 6.4% | 0% | 0.4% | | 69.1% | 24.2% |
| (2,1) | 6.02M | 304 / 471 | 0% | 11.9% | 0% | 1.3% | | 49.1% | 37.6% |
| (2,2) | 8.41M | 146 / 81 | 0% | 29.6% | 22.7% | 35.6% | 1.9% | 10.1% | 0.0% |
| (1,3), (3,1) | 0.42M | 70-75 / 100-130 | 0.7-0.8% | 72-74% | 25-27% | 0.5-1% | | 0.1% | 0% |

Creating the (1,2)/(2,1) children is 61% of the rdtsc cycles (the support
filter 24%), and they mostly die in the support filter, which the pretest
does not replicate: a threshold of 3 only adds work there. Timing
(cx/dfirst + this patch, d-sampled with `--sample-seed 11`, per-d minimum
of two alternating runs, against cx/dfirst's own binary):

| threshold | off | 1 | 3 | 4 | 5 (default) |
| --- | ---: | ---: | ---: | ---: | ---: |
| S = 2200 (N 15.2k, 60 d) | 1.020 | 0.952 | 0.987 | 1.053 | 1.015 |
| S = 2400 (N 23k, 30 d) | 0.980 | | 0.941 | 1.062 | 1.020 |

So 4 is the best value for V_d, but at 1.05-1.06x (1.03-1.08x against
the same binary with the pretest off) it is below the 1.08x on both sums
that would justify a separate setting. V_d keeps the default, which costs
nothing there (1.015, 1.020). The d-sampled runs also had the same nodes
at every d.

Not kept: one more round of cross support and the support filter at the
surviving (2,2) children (cx/prune's SURV_FIX), gated by their survival
(on when more than a fifth of the last ~256 survived, in the style of
the cross_after switch), on top of the pretest, same binary:

| P / S | nodes | speedup |
| --- | --- | ---: |
| 13 7 4 3 1 1 / 2200 | 28.0M -> 15.9M | 1.008 |
| 11 6 4 3 2 1 / 2174 | 24.2M -> 13.9M | 1.021 |
| 14 7 4 4 1 0 0 1 / 3648 | 7.96M -> 6.09M | 1.019 |
| 13 7 4 3 1 1 / 2400 | 52.5M -> 28.6M | 1.039 (always on: 1.050) |

It halves the nodes, but with the pretest the dead children it spares
are cheap, so it saves 1-4%, less than the 3% on two of these sums that
was the bar for keeping it. (On `full.txt` the gate is nearly always off:
14,958,507 -> 14,957,941 nodes.)

Correctness: `ctest -R fast_` passes; `fuzz_arrange` has variants with
the pretest from 1, 3, 4 vectors placed (also with 4 words, without the
support filter, without MRV, with cross support everywhere). For each of
them it also checks that the search with the pretest off visits the same
nodes. It passes on 1000 fresh seeds (7300000-7300999), and the
d-first fuzz (`fuzz_arrange --dfirst` on cx/dfirst, with the pretest from
1-5 vectors placed per seed) on 997 more. Three mutants fail on 300
seeds: the count b check of the pretest off by one (38 failures, all by
the node check alone: the children it lets through die later), one word
missing from the union (359 failures, squares lost), and one lane
dropped from the applied lists (from the first seed on).

CPU used: ~1.3 CPU-hours (~10 minutes of it a runaway mutant).
## Diagonal-first search in msearch (October 2026, branch cx/dfirst)

`msearch --diag-first` (src/c/dfirst.c) searches the sums with at least
`--diag-first-min-n` vectors (after reduction; default 3000 here, 5000 since
the integration) diagonal-first,
as `bin/dsearch` did (see "Diagonal-first search (second look)" above): for
every vector d of the sum's unreduced list, the semi-magic search on
V_d = {v in the reduced list : |v & d| = 1} finds the squares with d as an
SP traversal. So the loop over d finds every (square, SP traversal) pair
once, and every magic square at least twice, once per diagonal (dedupe by
`hash`). It does not enumerate the semi-magic squares, so `--calib-r1-stride`
adds an r1-sampled plain search for the models.

* V_d comes from a number -> reduced-vector posting index: a hit counter per
  vector is bumped along the postings of d's 6 numbers, the vectors hit once
  are kept (in list order, through a bitmap of the touched words), and only
  the touched counters are cleared. That is O(sum of the posting lengths)
  per d instead of dsearch's 36 compares per vector: 0.15 s instead of
  0.24 s for the whole d loop on P = 10 6 3 1 0 1, S = 391, with the same
  nodes. V_d plus the search's setup are 1.2% of the time at N = 7.6k and
  0.2-0.3% at N = 15-23k.
* d is processed in chunks of `--d-chunk` (256) indices, each followed by a
  "dchunk" record, and `--d-range lo:hi` searches one range of d (units of a
  sum, or a restart after the last completed chunk). `--d-stride k` with a
  random `--d-offset` samples d; the "dsum" record then holds unbiased
  estimates (`est_pairs`, `est_time`) and simple-random-sampling standard
  errors. `--d-log` writes one line per d: |V_d|, labels, nodes, pairs and
  the V_d, setup and search seconds.
* The r1 sampling of cx/profile is now in search_opts_t (`r1_stride`,
  `r1_offset`, and stratified `r1_nstrata` / `r1_sstride` / `r1_soffset`,
  plus a per-r1 log). msearch exposes it as `--r1-stride`, `--r1-strata`,
  `--r1-log` and the `SAMPLE_*` environment variables (msearch no longer
  reads those since the integration's fixes), and bench reads the
  environment variables. A sampled plain sum writes a "csum" record with
  `est_squares`, `se_squares`, `est_time` and `se_time`. The default root
  loop is unchanged: the same nodes and hashes on quick/full/prod, and on
  `bench/prod.txt` 13.70-13.87 s against 13.86-14.04 s for fdb77fc's
  bench, in alternating runs.

**The root of the V_d searches (new).** Every square of V_d contains all 6
numbers of d: its rows meet d once each and are disjoint. So d's numbers get
the top labels, with the one in the fewest vectors of V_d on top. Then every
square's largest label is that number, and only the first rows through it
(its class, ~1/6 of V_d) are roots, instead of all of V_d. On the same d,
this cuts nodes to 0.51-0.57x and time to 0.64-0.70x at every N measured
(4.1k-23k). The pairs are identical: gate (a) below on all five sums, and
gates (b) and (c). `--d-plain-root` restores the plain root, which
reproduces dsearch's nodes. The per-d cost varies more with the new root
(CV 0.54-0.61 against 0.37-0.41: the size of the rarest class varies), but
the top 10% of d still hold only 19-23% of the time.

Correctness:

* (a) With `--d-plain-root` and d-stride 1, msearch reproduces dsearch
  exactly: 622,221 nodes on 10 6 3 1 0 1 / 391 and 1,723,937 on
  13 6 3 2 / 517. It also gives the same nodes and the single SP pair on
  13 5 3 2 0 1 / 632, 16 5 4 2 / 849 and 12 6 3 2 1 0 1 / 836. The new root
  finds the same pairs.
* (b) Differential test against the plain search on 39 sums with N = 450-3600
  (bench quick/full/prod plus the 14 sums of the scheduler runs that have SP
  squares): 125 squares and 15 SP pairs, the same multiset with both roots.
  The expected multiset is {(hash, t)} over the plain squares' traversals t
  with sum S and product P, computed in Python from their grids. On every
  pair, `set_count == sp_count` and `partner == magic`. d-stride 3 (offsets
  0, 1, 2), `--d-range` splits and chunk sizes add up to the same pairs and
  nodes.
* (c) `fuzz_arrange --dfirst`: random families with random candidate
  diagonals, some outside the family, the oracle squares' random
  traversals, sometimes the family's own vectors, and a planted pair (a
  square's traversals on the diagonal and anti-diagonal together). Brute
  force lists all (square, traversal in the set) pairs, with the number of
  traversals in the set and the partner flag computed from all pairs of
  them. The planted square must come twice, flagged both times. Both roots,
  the default options plus one other option set per seed, and a d-stride 3
  split. 1000 fresh seeds: 0 fails, with 53,639 pairs, 29,433 flagged and
  1,426 planted squares. ctest runs it as `fast_fuzz_dfirst` and
  `fast_fuzz_dfirst_matrices`. Mutations of the partner test and of the
  bitmap scan fail 33 and 84 of 120 checks (60 seeds). One mutation,
  keeping the vectors
  that meet d more than once, is not caught, but it is equivalent: a square
  of such vectors still has rows that meet d once each.
* (d) `ctest -R fast_` passes (33 tests), the bench quick/full/prod hashes
  and nodes are unchanged, and `fuzz_arrange` passes 1000 fresh seeds.

Measurements (fdb77fc search code; thread CPU seconds per sum, estimated).
Method:

* Plain search: r1-sampled with 4 strata (r1-index quartiles, strides
  k : 2k : 6k : 40k, after the per-quartile costs and CVs of the cx/profile
  logs).
* d-first: d-sampled; the plain root on the same d as the new root's first
  replicate.
* Two independent replicates of each, run alternately (d-first A, plain A,
  plain-root A, plain B, d-first B) in one session, with one binary.
* The machine was shared (load 4-8 on 4 cores) for both arms alike.
* relSE: standard error / estimate, by the stratified (plain) or
  simple-random-sampling (d) formula.

| P | S | x | N | labels | plain CPU-s (relSE) | d-first (relSE) | ratio (relSE) | plain root | ratio | new / plain root, same d |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 7 4 3 1 1 | 1900 | 0.10 | 4,111 | 137 | 36 (2%) | 31 (3%) | 0.86 (4%) | 46 | 1.28 | 0.68 |
| 13 7 4 3 1 1 | 1950 | 0.13 | 5,896 | 153 | 187 (4%) | 115 (3%) | 0.61 (5%) | 171 | 0.91 | 0.65 |
| 12 6 3 2 1 1 | 988 | 0.34 | 6,671 | 159 | 274 (4%) | 156 (3%) | 0.57 (5%) | 229 | 0.84 | 0.69 |
| 13 7 4 3 1 1 | 2000 | 0.15 | 7,593 | 168 | 509 (4%) | 282 (3%) | 0.56 (5%) | 416 | 0.82 | 0.68 |
| 13 7 4 3 1 1 | 2100 | 0.20 | 11,306 | 192 | 2,978 (6%) | 1,019 (5%) | 0.34 (8%) | 1,513 | 0.51 | 0.64 |
| 13 7 4 3 1 1 | 2200 | 0.25 | 15,199 | 211 | 9,395 (7%) | 3,242 (6%) | 0.35 (9%) | 4,659 | 0.50 | 0.66 |
| 11 6 4 3 2 1 | 2174 | 0.25 | 16,424 | 220 | 11,598 (9%) | 4,262 (6%) | 0.37 (10%) | 6,280 | 0.54 | 0.70 |
| 14 7 4 4 1 0 0 1 | 3648 | 0.25 | 20,538 | 252 | 18,031 (7%) | 4,825 (5%) | 0.27 (9%) | 7,442 | 0.41 | 0.64 |
| 13 7 4 3 1 1 | 2400 | 0.33 | 22,992 | 245 | 45,614 (9%) | 12,672 (5%) | 0.28 (10%) | 20,476 | 0.45 | 0.64 |
| 12 6 3 2 1 1 | 1200 | 0.53 | 11,697 | 199 | 2,188 (7%) | 1,012 (8%) | 0.46 | | | |
| 12 6 3 2 1 1 | 1480 | 0.74 | 15,074 | 236 | 3,347 (7%) | 1,224 (7%) | 0.37 | | | |

* **At N >= 15k, d-first costs 0.27-0.37x the plain search per sum**
  (geometric mean 0.31). The design without the new root costs 0.41-0.54x
  (geometric mean 0.47).
* The cost per expected magic square is the ratio of the CPU per sum, since
  both find every magic square. So it falls 2.7-3.7x at N = 15-23k, and
  ~1.8x at N = 6-8k.
* Measured against the cost per semi-magic square of the plain search
  (4.2-33 CPU-s at N = 15-23k), d-first spends the equivalent of 1.4-8.8
  CPU-s per square.
* The ratio falls like N^-0.65. Fitting the new root's ratios at x <= 0.4
  gives a crossover at **N ~ 3,000**. These are the table's first nine
  rows, the d-first-only full runs of gate (b) at N = 2.2-3.6k (0.93-1.28)
  and 13 7 4 3 1 1 / 1850 (N 2.5k, 1.2).
* The plain root crosses over at N ~ 5,500.
* `--diag-first-min-n` defaults to 3000. Small N with large x is the only
  place where d-first loses badly (10 6 3 2 1 / 1360, N 3.1k: 2.6x), and
  those sums essentially never have squares. At N = 11.7-15k with
  x = 0.53-0.74 it still wins (0.46, 0.37).
* The plain estimates agree with earlier runs: 1900 is exact (35.9 s), and
  988 gives 274 s against 273-290 s. 2200 and 2400 give 9.4k and 45.6k,
  against the cx/profile estimates of 10.2k and 42.4k (SE 15-34%).

Caveats:

* All of these sums have <= 256 labels (the carried path).
* The search improvements in progress elsewhere (wider carried path,
  cheaper (2,2) children) speed up the plain search and the V_d searches
  alike, but not necessarily by the same factor: V_d's lists are ~6x
  shorter. The ratios should be re-measured on the integrated binary
  (done: they did not speed up alike, see "Measurements on the integrated
  binary" below).
* d-first sums write "dsquare" / "dchunk" / "dsum" records, which the
  scheduler does not read yet. It would need to mark them covered and to fit
  its squares model to the calibration stream. (Both schedulers read them
  since the integration, and keep them out of their fits.)

Tried on top of the new root, not kept:

* **Forward checking on d's numbers.** Every number of d that is in no
  placed vector of an axis must be in a candidate of that axis. This is two
  more masked tests per child, using the unions the forward check already
  has. It pruned 1% of the nodes (P = 10 6 4 2 1 1 / 855: 15.27M -> 15.19M;
  S = 2200, 20 d: 7.14M -> 7.06M), with time within noise (0.96x). The
  generic forward check already catches nearly all of it, so the core was
  left unchanged.
* **Reducing V_d** (reduce_vectors on each V_d before its search). It
  removes 0-0.05% of V_d with the default reduction and 1-2% with the
  strong one, and nodes fall by 0.6% (S = 2200, 10 d). Slower.

CPU used: ~31 min (gates ~6, experiments ~3, measurements ~19, the check of
the plain path ~1.5).

## Integration of cx/wide, cx/pretest and cx/dfirst (October 2026, branch cx/integrated)

The three prototypes, each one commit on fdb77fc, cherry-picked in that
order onto a3812f0, then the verifiers' required fixes in separate commits,
then rebased onto scheduler v2 (df6b62b).

**The root loop.** All three changed it: per-r1 width runs (wide), the
stratified r1 sampling (dfirst, in search_opts_t), the top-label root
(dfirst, `top_root_only`), and wide's env-driven `-DR1_SAMPLE`. Combined:

* the roots are [0, root_limit) (all of N, or with `top_root_only` the
  prefix of vectors through the top label), and search_root partitions
  them into the width runs (each run at ceil((x + 1) / 64) words, at least
  2 and at least `--min-words`; one run at the full width with
  `--no-r1-width` or on the matrix path);
* r1 sampling selects on the global index into the root list: strata are
  equal index ranges of [0, root_limit), each run asks `r1_next` for the
  first sampled r1 at or after its start, so every sampled r1 is searched
  exactly once, at the width of its run, and the sample does not depend on
  the widths (checked: the same r1, strata, strides, squares and nodes per
  r1 as cx/dfirst's binary, with and without the pretest);
* `-DR1_SAMPLE` is gone: its stride / offset / log are the opts sampling
  (`SAMPLE_*` in bench, and at first in msearch: see the second
  verification below), and its `SAMPLE_LIST` became
  `opts.r1_list` (bench reads `SAMPLE_LIST`; the estimates are then the
  sampled totals, no standard error). The r1 log keeps cx/dfirst's six
  columns and adds wide's largest label and width, plus wide's "# run"
  lines.

The default path is node-identical to fdb77fc on bench quick / full / prod
(1,770,779 / 14,958,507 / 50,375,738 nodes, every instance's squares and
hash), also with `--no-r1-width`, `--pretest-min 0` and `--min-words 6`.

**Verifier fixes.** The pretest's `filtered:` label is followed by an
empty statement (a declaration after a label is C23 only); its early
return counts the child in the `-DPROFILE` counters (the per-depth profile
is now the same with the pretest on and off); `fuzz_arrange --dfirst`
applies each variant's pretest and fixw, so the pt_* variants run in the
V_d searches; new variants run the pretest at 5-8 words (pt_1_w5, pt_3_w6,
pt_4_w7_fixw, pt_2_w8_xcross, pt_5_w8_fixw). For d-first: "dsum" records
carry `"mode":"dfirst"` and `complete`, the "done" record of a
`--diag-first` run `"mode":"dfirst"`, and scheduler.py's Results keeps
d-first sums out of its fits (`results.dsums`, `results.dmagic`), covering
a complete one and cutting incomplete ones out of the file's done ranges;
with `--d-range lo:hi`, lo > 0, plain sums are left to the unit from d 0
("skip" records); dfirst.c asserts at compile time that inv[] holds the
involutions of n <= SQ_MAX_N; bench closes the SAMPLE_LOG file; strata
with fewer than 2 sampled r1 are counted (`r1_strata_nose`,
"se_strata_missing"), since the standard errors cannot include them.

**The crossover moved to N ~ 5000.** Full searches on the integrated
binary, process CPU seconds, min of 2 alternating runs, one at a time
(`nice`, shared machine at load 1-2):

| P | S | N | labels | plain | d-first | ratio | cx/dfirst ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 13 7 4 3 1 1 | 1850 | 2,523 | 115 | 3.00 | 3.63 | 1.21 | 1.2 |
| 13 7 4 3 1 1 | 1900 | 4,111 | 137 | 23.79 | 31.30 | 1.32 | 0.86 |
| 13 7 4 3 1 1 | 1950 | 5,896 | 153 | 149.58 | 116.11 | 0.78 | 0.61 |
| 12 6 3 2 1 1 | 988 | 6,671 | 159 | 245.17 | 155.10 | 0.63 | 0.57 |

The d-first times are those of cx/dfirst (31 s, 115 s, 156 s): the V_d
searches have ~2 words of labels already. The plain search got 1.5x
faster at S = 1900 (36 -> 24 s) and 1.25x at 1950 (187 -> 150 s): with
137-159 labels fdb77fc ran every r1 at 3 words, the per-r1 widths run all
but the first few at 2. Interpolated, the ratio crosses 1 at N ~ 4,960, so
`--diag-first-min-n` defaults to 5000 (was 3000), which also keeps the
plain search's semi-magic squares, for the models, for the sums below.
The ratios at N >= 11k were re-measured on this binary after the
verification: see "Measurements on the integrated binary" below (5000
stays the default).

Gates (before the rebase): `ctest -R fast_` (40 tests) passes;
`fuzz_arrange` on fresh seeds (52000000-52700199): 500 default, 500 mode 7,
500 `--dfirst` (every variant, pt_* included), 200 `--dfirst --mode 7`, 200
`fuzz_arrange_matrices --dfirst`, 100 matrices mode 7, 200 portable, 20
mode 6 with the wide and pretest variants: 0 fails; `scripts/
test_scheduler.py` has a d-first test (synthetic records, and msearch
output with complete, sampled and d-range units).

**Rebased onto scheduler v2.** Scheduler v2 (`scripts/amodel.py`, the
Summary of scheduler.py) reads msearch output through its own incremental
Summary, not v1's Results, and fits its time law to a new "cpu" field of the
sum records. The rebase kept both:

* msearch: the dsum records get "cpu" too (the d loop, the reduction and
  the sum's share of the enumeration; process CPU, as in "sum"), and so do
  the csum records of `--r1-*` runs; with `--calib-r1-stride` the csum's
  cpu is the sampled search only (the dsum holds the rest), so the cpu
  fields of a file add up without double counting. cx/dfirst's own
  `cpu_time` (thread CPU) gave way to v2's (process CPU; msearch is
  single-threaded).
* A run that r1-samples its plain sums wrote an ordinary "done" record
  over sums with only "csum" records, the same trap as the d-first sums:
  its done record now has `"r1_sample":1` (and `"mode":"sampled"` without
  `--diag-first`), one of v2's SAMPLED_KEYS, so both schedulers skip it.
* v2's Summary: an incomplete dsum, a "skip" and a sampled csum are holes
  in the file's done range, kept in the file's state so they hold across
  incremental reads and a reload; a complete dsum covers its sum. d-first
  records stay out of the cells (squares and traversal fits) and the time
  law; their sums, parts and CPU are counted apart (`totals["dfirst_*"]`),
  and a magic square found d-first joins the notable squares once (flagged
  `dfirst`), so `run` announces it and `report` lists it. SUMMARY_VERSION
  in the summary's tag rebuilds older summaries. v1's Results also skips
  sampled records now. `test_dfirst_summary` checks all of it against the
  plain sums alone, and on msearch output with both schedulers.
* A conflict resolution of the dfirst cherry-pick had dropped the `break`
  after msearch's `--pretest-min`, which then also switched on
  `--diag-first` (fixed in that commit; the new test runs `--pretest-min`
  and checks for a plain done record).

Gates after the rebase (build-integ, Release): `ctest -R fast_` passes (40
tests, fast_scheduler and fast_calibrate included); `fuzz_arrange` on
fresh seeds 61000000-61600099: 600 default (1 skipped by the oracle
budget), 500 mode 7, 500 `--dfirst` (every variant, pt_* included), 200
`--dfirst --mode 7`, 200 `fuzz_arrange_matrices --dfirst`, 100 matrices
mode 7, 200 portable: 0 fails; bench quick / full / prod: 1,770,779 /
14,958,507 / 50,375,738 nodes, every instance's nodes, squares and hash
equal to fdb77fc's bench.

### Measurements on the integrated binary

Paired measurements of the plain search of fdb77fc (base) and of this
binary (new), and of this binary's d-first search (process CPU seconds per
sum; both CMake Release, `-O3 -march=native` with LTO; fdb77fc's bench with
cx/wide's R1_SAMPLE patch):

* plain: bench on identical r1 lists in both builds, 7 strata at r1-index
  fractions 0, 1/16, 1/8, 1/4, 3/8, 1/2, 3/4, 1 with random offsets; base
  runs a stratified subset of new's sample; runs alternate (base, new,
  d-first, base, new, d-first), per-r1 min of 2; sums of <= 1,200 CPU-s run
  in full, twice;
* d-first: `msearch --diag-first --diag-first-min-n 0 --d-stride k` with an
  explicit offset (113-1,370 d per sum), min of 2 per d where repeated;
* on identical lists base and new have the same squares and the same nodes
  per r1 up to 252 labels (above 256, base is on the matrices: 2.7-4.5x
  the nodes); relSE <= 8% (plain), <= 6% (d-first); load 3-6.6 from another
  workflow; ~63 CPU-min.

| P / S | N | labels | base plain | new plain | new d-first | new / base | d-first / base | d-first / new |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 6 4 2 1 1 / 838 | 1,986 | 104 | 1.01 | 1.01 | 1.65 | 1.007 | 1.6 | 1.60 |
| 12 6 3 2 1 0 1 / 900 | 2,994 | 121 | 5.9 | 5.8 | 5.9 | 0.984 | 0.97 | 0.99 |
| 13 7 4 3 1 1 / 1900 | 4,111 | 137 | 35.5 | 23.5 | 29.4 | 0.663 | 0.83 | 1.25 |
| 12 6 3 2 1 1 / 988 | 6,671 | 159 | 260 | 231 | 149 | 0.892 | 0.575 | 0.645 |
| 13 7 4 3 1 1 / 2000 | 7,593 | 168 | 511 | 441 | 251 | 0.863 | 0.491 | 0.569 |
| 12 6 3 2 1 1 / 1200 | 11,697 | 199 | 2,023 | 1,575 | 870 | 0.778 | 0.430 | 0.552 |
| 13 7 4 3 1 1 / 2200 | 15,199 | 211 | 11,316 | 8,223 | 3,432 | 0.727 | 0.303 | 0.417 |
| 11 6 4 3 2 1 / 2174 | 16,424 | 220 | 11,503 | 8,685 | 4,056 | 0.755 | 0.353 | 0.467 |
| 14 7 4 4 1 0 0 1 / 3648 | 20,538 | 252 | 14,460 | 12,623 | 5,407 | 0.873 | 0.374 | 0.428 |
| 9 6 4 3 1 1 1 1 / 2700 | 20,896 | 259 | 53,077 | 13,835 | 4,457 | 0.261 | 0.084 | 0.322 |
| 13 7 4 3 1 1 / 2400 | 22,992 | 245 | 47,400 | 40,389 | 11,785 | 0.852 | 0.249 | 0.292 |
| 13 7 4 3 1 1 / 2500 | 26,585 | 260 | 205,729 | 54,280 | 23,366 | 0.264 | 0.114 | 0.430 |
| 13 7 4 3 1 1 / 2650 | 31,743 | 279 | 467,238 | 121,918 | 39,917 | 0.261 | 0.085 | 0.327 |

* **Plain search.** Unchanged up to 128 labels (1.007, 0.984). Above, by the
  width each r1 runs at: an r1 at the same width 0.84-0.96x (the pretest),
  one word fewer 0.66-0.81x (per-r1 widths; at 137 labels nearly every r1
  drops from 3 words to 2, hence 0.663), and the matrices replaced by the
  carried path above 256 labels 0.26x. Per sum 0.73-0.89x at N 7.6-23k
  with <= 256 labels, 0.26x above.
* **d-first did not speed up alike.** Its CPU per sum is 0.86-1.12x
  cx/dfirst's fdb77fc estimates at <= 256 labels: with the top-label root
  every root of V_d contains the top label, so the per-r1 widths never
  narrow a V_d search, and the pretest does ~nothing there (see its
  section). Above 256 labels every V_d root runs at 5 words (S = 2500: 79%
  of the V_d have 257-260 labels, 1.12 us/node against 0.88 at 4 words;
  S = 2700: 46%), which is why S = 2500 is at 0.43 against 0.29 at S =
  2400. So d-first / plain is 0.55 at 11.7k, 0.42-0.47 at 15-21k, 0.29 at
  23k (geometric mean 0.38 at N >= 15k, against 0.31 on fdb77fc), and
  0.32-0.43 at 21-32k with more than 256 labels. Keeping V_d at or under
  256 labels (or narrowing its roots) is a possible later optimization.
* **Best mode against fdb77fc's plain search:** 0.66 at 4.1k (plain),
  0.49-0.58 at 6.7-7.6k, 0.43 at 11.7k, 0.25-0.37 at 15-23k (<= 256
  labels) and 0.085-0.114 at 21-32k (> 256 labels): d-first stays worth
  1.5-3.4x on top of the plain search's gain at N >= 6.7k, which scheduler
  v2 does not use yet (it never launches d-first units; a follow-up: launch
  `msearch --diag-first --calib-r1-stride k` at N' >= 5000, with a d-first
  time law such as the fit below).
* **Crossover.** Along 13 7 4 3 1 1 (with the full runs above, 1850 = 1.21,
  1950 = 0.78) d-first / plain crosses 1 at N0 ~ 4.9k; pooled over the 13
  sums, ln r = -0.028 - 0.566 ln(N / 4000) (resid sd 0.17), N0 ~ 3.8k. It
  depends on P: 12 6 3 2 1 0 1 is at 0.99 at N 3.0k, 13 7 4 3 1 1 at 1.25
  at 4.1k. `--diag-first-min-n` stays 5000: sums at 3.8-5k would save at
  most ~10%, and below 5000 the plain search's semi-magic squares feed the
  models. Scheduler v2 chooses the mode by this pooled ratio (with its 7%
  calibration stream: d-first from N' ~ 4.3k), not by the quotient of its
  two time laws, which put the switch at ~9-10k (research/scheduler-v2.md,
  "d-first units"); six more pool sums at N ~6.1k measured 0.69-0.82, on
  the pooled law.
* **Time per sum**, t = a (N / 4000)^b CPU-s, least squares on ln t over the
  11 sums with N >= 3k (4.1-31.7k, the 3 above 256 labels included):

  | search | a | b | resid sd | t(16k) | t(32k) |
  | --- | ---: | ---: | ---: | ---: | ---: |
  | fdb77fc plain | 24.8 | 4.49 +- 0.23 | 0.47 | 12.5k | 280k |
  | new plain | 25.9 | 4.04 +- 0.14 | 0.29 | | |
  | new d-first | 25.6 | 3.46 +- 0.13 | 0.27 | | |
  | new, best of the two | 23.0 | 3.53 +- 0.14 | 0.27 | 3.1k | 35k |

  (13 7 4 3 1 1 alone: new best b = 3.63, base 4.56.)
* **CPU per expected (square, SP traversal) pair** (E = amodel E_semi x
  SQ12 x 720 p_SP, 0.01-0.66 per sum; no pair occurred in the samples, and
  the full runs at N 2-3k found 0 against 0.01-0.02 expected): best new
  ~300 CPU-s at N 3-4k, 1.0-1.1k at 6.7-7.6k, 5.4-10.6k at 11.7-16.4k,
  17.8k at 23k, 38k at 26.6k, 80k at 31.7k; fdb77fc 450, 1.9-2.0k,
  17.7-25k, 71k, 338k, 939k. Along 13 7 4 3 1 1 it grows like N^2.7 (new
  best) against N^3.7 (fdb77fc).
* Scheduler v2's engine-2 time law against these sums (obs / pred): fdb77fc
  0.39-2.6; new plain 0.24-0.72 at N >= 11k. Hence engine 3 below.

Also measured here, the top-label root on the matrix path (CARRY_MAX_W=0
build, same d, min of 2): 12 6 3 2 1 1 / 988 (159 labels, 167 d at stride
40) 7.33 s against 10.29-10.97 s for the plain root (0.68x; the carried
build 3.7 against 4.8 s), but 13 7 4 3 1 1 / 2500 (260 labels, 18 d at
stride 1500, the 8-word matrices) 48.5 s / 30.3M nodes against 37.6 s /
29.2M (1.27x): on the 8-word matrices the top labels cost more nodes than
on the carried path (13.2M). dfirst.c now uses the top root only where it
pays (`top_root_pays`: the carried path, or the matrices of up to 4
words), deciding per V_d by its own labels when the sum's do not decide
it: 35.1 s on that sample (some V_d have <= 256 labels). Native builds,
which carry up to 512 labels, are unaffected.

### Fixes after the second verification

Three verifiers (correctness, performance, adversarial) passed the
integration with issues; these were fixed:

* **Engine 3.** The integration changes the plain search's CPU per sum
  (above), but msearch still wrote `"engine":2`, so scheduler v2 would have
  fitted one law to both builds' timings, and its +0.39 step at > 128
  labels (the third word, x1.48) now describes a cost that has mostly gone
  (12 6 3 2 1 1 / 890, 131 labels: 31.9 -> 21.8 CPU-s, 1.47x, both engine 2).
  msearch writes `"engine":3` and scheduler.py has ENGINE = 3. Each engine's
  law starts from the previous engine's posterior; `amodel.time_prior` now
  also shifts it by `ENGINE_TIME_SHIFT`: engine 3's > 128-label step starts
  0.227 lower (0.39 -> 0.17), the geometric mean of the paired new / base
  ratios at 137-252 labels (0.80). No feature of the law describes the
  0.26x above 256 labels (3 sums at N 21-32k); its shape (N slope, hinge,
  W3 step) should be refit for engine 3 once large-N data arrives, with the
  13 sums above as anchors. v1 falls back to its default model until it
  has 20 engine-3 sums.
* **No SAMPLE_* in msearch.** msearch had taken bench's `SAMPLE_STRIDE`,
  `SAMPLE_OFFSET` and `SAMPLE_LOG` from the environment, so a variable left
  exported from a bench session turned every scheduler unit into an
  r1-sampled run: only csum / csquare records and a sampled done record,
  which neither scheduler counts, so v2 re-planned the same unit forever.
  msearch reads only its `--r1-*` options now, and both schedulers start
  msearch without `SAMPLE_*` (`msearch_env`). ctest `fast_msearch_no_sample_env`.
* **Option checks.** Numeric values, `--d-range lo:hi` (lo < hi), the
  d-first options without `--diag-first`, more than 8 `--r1-strata` and
  log files that cannot be opened are refused (exit 2); the usage text
  lists every option (ctest `fast_msearch_bad_args`). With `--d-stride k
  --d-offset o`, o mod k != 0, the plain sums get "skip" records as with
  `--d-range` lo > 0. `--time-limit` is checked between the chunks of a
  d-first sum too (the sum is then not complete); `--node-limit` applies
  per V_d there (documented).
* **Strict C17.** `#define _POSIX_C_SOURCE 200809L` in arrange.c, dfirst.c,
  msearch.c and fuzz_arrange.c: clock_gettime with the CPU-time clocks is
  declared under `-std=c17` / `-std=c11 -pedantic-errors` again (the
  earlier claim that `-std=c17 -pedantic` was clean was wrong for this
  tree; the CMake build uses gnu17 and was not affected).
* **`--calib-r1-stride 1`** wrote est_* = 0: a stride of 1 is now a
  sampling plan over every r1 (exact estimates); the default 0 stays no
  plan.
* The dsum `se_time` uses the sum of the per-d CPU (it mixed the whole
  loop's thread CPU with the per-d squares); the README and msearch's header
  say which dsum time fields are CPU and which wall.
* The csum `est_nodes` / `est_time` are unbiased only up to the adaptive
  cross-support state, which carries over from r1 to r1 (a sampled r1 can
  differ by +-1 node from the full run; mean over all offsets 5,009,338
  against 5,009,395 nodes at 13 6 3 2 1 1 / 899, exact with `--no-cross`);
  `est_squares` is exactly unbiased (18.000000 over all offsets). Wording
  fixed in arrange.h, msearch.c and the README.
* ARRANGE_DEBUG: FILTER_COUNT checks the padding precondition it shares
  with FILTER_CARRY (the pretest replaced FILTER_CARRY at depth >= 5).
* Summary bookkeeping (SUMMARY_VERSION 3): `dfirst_sums` counts distinct
  (P, S); a notable square is kept once per (P, S, hash) also on the plain
  path (a sum searched twice, or a square found d-first and then plain);
  the CPU of sampled csum records is totalled (`sampled_cpu`); v1's report
  lists the magic squares found d-first, not only their number.
* New ctests for the paths that had none: the no-GFNI / no-VPOPCNTDQ build
  on the carried path with 7 words and without the pretest
  (`fast_bench_quick_no_gfni_wide7`, `_no_pretest`), its fuzz at 257-500
  labels and d-first (`fast_fuzz_arrange_no_gfni_wide`,
  `fast_fuzz_dfirst_no_gfni`), and d-first in the portable build
  (`fast_fuzz_dfirst_portable`).

Rejected or left open (the first two since done in scheduler v2, see
research/scheduler-v2.md, "d-first units"): the d-first parts of one sum (`--d-range` units,
`--d-offset` units) are still never merged into coverage, and nothing reads
the "dchunk" records (no resume): a split sum counts as unsearched and
would be planned again (documented in src/c/README.md). V_d searches start
with a fresh adaptive cross-support state (a possible small gain, not
measured). d-first
SP squares that are not magic are not notable (every dsquare has an SP
traversal, so they would flood the list).

Gates after these fixes (build-integ, Release): `ctest -R fast_` passes (47
tests); every commit compiles with `-Wall -Werror`, and arrange.c / dfirst.c
with `-std=c17 -pedantic-errors`. `fuzz_arrange` on fresh seeds
(91000000-95000004), 500 seeds per mode, 0 fails anywhere: default (1
skipped by the oracle budget), `--mode 1` (forced widths, every variant),
`--mode 5` (hubs: 442 ran, 58 skipped), `--mode 6` (saturated counters,
variants default, nomrv, w5, w8, pt_1_w5), `--mode 7`; `--dfirst` with
modes default (1 skipped), 1, 5 (443 ran, 57 skipped) and 7; the matrices
build default, mode 7, `--dfirst` (2 skipped) and `--dfirst --mode 7`
(which takes both roots: 76 of 13,313 V_d of 40 of its seeds had > 256
labels); the portable build default and `--dfirst` (1 skipped); the
no-GFNI build default, mode 7 and `--dfirst`. The d-first runs checked
623,403 (square, d) pairs. bench quick / full / prod: 1,770,779 /
14,958,507 / 50,375,738 nodes, every instance's nodes, squares and hash
equal to fdb77fc's, also with `--no-r1-width`, `--pretest-min 0` and
`--min-words 6`. Mode 4 (large sparse families, limited by its brute-force
oracle) was not rerun.

## Class support in the V_d searches (October 2026, branch c2/classsup)

**The rule** comes from the prototype c2/vdclass 65c6495: its CLASS_VEC
level 2 at the (1,1) children, i.e. CLASS_SUP=6 CLASS_LVL=2 CLASS_VEC=1.
Every vector of V_d holds exactly one of d's numbers, its class, and these
numbers have the top labels. So the rows (cols) of a square of V_d have
the n classes: the row of class j meets the col of class j at d_j, and
every other col at a number outside d. At a node, take a candidate x of
axis a and class j, with oa the other axis and C the cells of oa's placed
vectors:

* (bad) d_j must lie in C or in a candidate of oa of class j, and every
  other cell of x in C or in a candidate of oa of another class. This is
  the support filter, with the classes.
* (meet) for every class c != j not placed on oa, x must meet a candidate
  of oa of class c outside C.

The passes run on both axes, alternately, to the fixpoint. The node dies
when an axis has too few candidates or an unmatched cell without one.
This is the greatest fixpoint (a test only gets stricter as the lists
shrink), so the nodes do not depend on how the passes are organized.

The candidate lists are sorted by descending labels, and the class label
is a vector's largest. So each class is a block of each list, and the
tests of block j use only the unions of the other axis' blocks.

**Production version** (`search_opts_t class_support`; `CLASS_SUP` in
arrange_core.h).

* It is on when all of these hold: `top_numbers` holds n numbers
  (`n_top == n`); the carried path is in use (W <= CARRY_MAX_W with
  AVX-512BW); the support filter and forward checking are on; and every
  vector has exactly one top label. The last is checked once at setup: a
  vector's largest label is a top one and its second largest is not, which
  V_d guarantees.
* Otherwise it is a no-op. The plain search and the matrix path are
  unchanged.
* It runs right after the support filter at the (1,1) children (the
  prototype's CLASS_SUP=6 hook).
* `msearch --diag-first` sets it by default. `--no-class-support` turns it
  off for A/B runs, and is refused without `--diag-first`.
* `fuzz_arrange --dfirst` runs every seed three ways: without it, with
  it, and with a rotating variant. The `cls_*` variants cover forced
  widths 3-8, `--no-r1-width`, the pretest from depths 1-3 and off, cross
  support everywhere or nowhere, no MRV, and runs without the support
  filter or forward checking (where it is off).
* With `--noracle` / `--diff`, the plain search's squares stand in for
  the oracle's (modes 4-6, whose oracle is slow or over budget).
* Through the API it also works with the plain root (`top_root_only` 0,
  which dfirst.c never combines with top numbers).
  * There an r1 without the top label still has a class label as its
    largest, x, and runs at its own width. That width is one word short of
    the class labels when they straddle a word boundary (130 labels: the
    classes are labels 124..129, and the r1 with x = 124..127 run at 2
    words).
  * `CLS_AT` treats the words above the run's width as empty: the classes
    above x are absent from that subproblem.
  * The first version read those words, i.e. stale lists of the wider
    runs. No square can be lost there (every square of V_d holds all the
    top numbers, so its r1 is in the full-width run), but ARRANGE_DEBUG's
    class-block check failed on 1 family in 1,000.
  * Debug builds now poison the words a narrower run must not read (all
    ones), which fails that check within 50 families on the old code.
  * `class_support_test` (ctest `fast_class_support_plain_root`, an
    ARRANGE_DEBUG build) checks 400 such families: the same squares under
    every combination of class support, root and per-r1 widths, and the
    same nodes with and without the per-r1 widths.
* `-DCLASS_PROF` (with CMake: `-DMAGIC_DEFS=CLASS_PROF`, a new option for
  instrumentation builds) prints its cycles, calls, passes and drops.

**Removing the prototype's overhead.** First a faithful port, then a
profile on a2200 = 13 7 4 3 1 1 / 2200 (31 d, stride 500). The table
gives TSC cycles per call of CLASS_SUP and its share of the search's
cycles. Every version has the same nodes and calls (189,567, 13.3% of
them dead) and, except the Gauss-Seidel one, 7.14 passes per call.

| version | cycles / call | share | what changed |
| --- | ---: | ---: | --- |
| faithful port | 12.1k | 18.0% | the prototype's per-call bounds (1.9%), union scans (0.55%), full passes until quiet (15.3%) |
| lazy keep masks, incremental passes, a branch per group | 19.2k | 26.3% | worse: mispredictions, and 2.8M union rescans (255M entries) |
| branch-free passes, testing only what changed | 13.6k | 20.4% | |
| Gauss-Seidel block order (unions updated block by block) | 14.2k | 21.1% | 5.09 passes instead of 7.14, but costlier passes: rejected |
| vectorized pass setup, loops specialized by the number of classes tested | 10.7k | 17.1% | |
| bounds by a branchless binary search on the class words (a histogram before) | 10.9k | 17.1% | bounds 1.44% -> 0.92% |
| final layout | 9.05k | 14.65% | see the list below |
| final, a generic loop over the classes | 10.1k | 15.7% | +9%: rejected |
| final, starting on the other axis | 9.6k | 14.8% | +4% (7.25 passes): rejected |

The final layout:

* Per class, the union G and the labels it lost since the last pass over
  the other axis (Lost) are one zmm each (a lane per word), with prefix and
  suffix ORs for "the other classes".
* A dropped entry is cleared in its group's keep mask (a byte per 8
  entries). The lists are compacted once, only if the node lives, so the
  blocks keep their bounds through the passes.
* After an axis' first pass, a pass over block j tests only two things:
  the cells that became bad, and the classes whose unions lost labels
  inside the block's own union. An alive entry passed every older test, so
  these are the only tests it can fail now.
* The class tests are or-ed per word and combined with an unsigned minimum
  (nonzero iff every class is met), 8 entries at a time. The loops are
  specialized (switch on the number of classes tested, 0-7, and on whether
  there are bad cells).
* A block's new union comes out of the pass that scans it (or_lanes8 over
  the per-word accumulators), and so do the lists' unions at the end. A
  pass that drops nothing skips the reduction.
* The fixpoint is kept, and so are the early returns: the death check
  (count, then unmatched cells against the unions in one ternlog) runs
  after every pass that drops entries.

Two items of the brief were done differently:

* **Bounds.** A per-call search remains, but it is branchless on the 1-2
  words that hold the class labels (clz), at 0.9% of the search's cycles.
  A histogram was slower. Carrying the block bounds through the filters
  that create the children would save at most that 0.9%, and needs
  class-aware variants of the hot filter loops.
* **Initial unions.** The other axis' unions at entry still come from one
  masked-OR scan (0.55-0.6%). Taking them as a by-product of the last
  support pass would need per-class unions in that hot loop.

On the other sums, the cost is 6.7k (a1950) to 14.8k (a2650) cycles per
call, at 7.0-9.5 passes. It is 11.6-35.1% of the search's cycles: the
lowest at 13 7 4 3 1 1 / 2650 (11.6%), the highest where it kills the most
(b1200 35.1%, c2700 32.4%: 40-48% of the calls dead).

**Gates** (all before any timing):

* **(a) Nodes equal to the prototype.** The prototype is the 65c6495
  build with CLASS_SUP=6 CLASS_LVL=2 CLASS_VEC=1. Equal on:

  | sample | stride | offset | nodes |
  | --- | ---: | ---: | ---: |
  | a2000 | 200 | 13 | 2,420,190 |
  | b1200 = 12 6 3 2 1 1 / 1200 | 400 | 13 | 1,472,383 |
  | a2200 | 500 | 13 | 7,205,201 |
  | c3648 = 14 7 4 4 1 0 0 1 / 3648 | 1000 | 13 | 2,941,370 |
  | c3648 | 1000 | 517 | 2,027,116 |
  | a2400 | 1200 | 13 | 9,367,374 |
  | a2400 | 1200 | 611 | 5,563,046 |
  | a2650 | 3000 | 13 | 5,634,085 |

  The other 7 measured samples are also equal, and all 15 again on the
  final binary.
* **Mutants.** Bounds off by one, a keep-mask update that clobbers other
  blocks, two wrong bad-cell formulas, and the count death check off by
  one each fail 70-76 of 100 `--dfirst` fuzz seeds. Two mutants
  under-prune: a wrong Lost update, and no incremental class tests. They
  find the same squares but take 2,718,161 / 2,644,286 nodes on a2000
  instead of 2,420,190, which gate (a) catches.
* **(b) Fuzz.** `fuzz_arrange --dfirst` on the final binaries, with 0
  mismatches anywhere. Every seed runs three checks: the default options,
  the default options with class support, and one rotating variant (by
  substring, the `--variants` selection also takes some plain variants).
  In all, 2.86M (square, d) pairs were checked.

  | runs | seeds |
  | --- | ---: |
  | default | 301 ran, 1 skipped by the oracle budget |
  | `--mode` 1, 2, 3 and 7 | 300 each |
  | `--mode 5` (hubs) | 36 with the oracle (4 over its budget); 300 with `--noracle` (1.63M pairs) |
  | `--mode 6` (saturated counters) | 300 with `--noracle` |
  | `--mode 4` (large sparse) | 300 with `--noracle` |
  | `--n` 4, 5, 7 and 8 | 150 each |
  | every variant in rotation (no `--variants`) | 298 (2 skipped) |
  | ARRANGE_DEBUG build | default 200; modes 1 and 7, `--n` 4 and 8: 100 each |
  | ARRANGE_DEBUG build, plain fuzz (with the poisoned words) | default and mode 7, 100 each |
  | build without GFNI / VPOPCNTDQ | default 149, mode 7 60 |
  | matrices build (class support a no-op there) | 100 |
  | portable build | 59 |
  | release build, plain fuzz | default 200, mode 7 100 |

  * Mode 4 has 300 seeds without 71600248. That seed's family (n = 8,
    4,373 vectors, 129-192 labels) kept the `--noracle` reference, the
    plain search with the default options, busy for over 10 minutes
    before any d-first check, when it was stopped.
  * Mode 4 with the oracle was dropped: its brute-force oracle spends
    minutes per seed. With `--noracle`, the plain search's squares are the
    reference, so modes 4-6 there test d-first and the class support
    against the plain search (which the plain fuzz tests against the
    oracle).
* **(c) msearch --diag-first on the full d lists, on vs off.** 39 sums with
  N = 451-3,626, 15 of them with a pair (632, 849, 836 and 12 more):
  * 125 semi-magic squares, 15 (square, SP traversal) pairs;
  * both runs pass the d-first vs plain-traversal differential (every pair
    of the plain squares found once, `set_count == sp_count`);
  * the multisets of (hash, d, partner, magic) of on and off are
    identical;
  * on/off nodes: 0.347 overall, 0.12-0.74 per sum;
  * the final binary repeats it (the same records and nodes).
* **(d) bench.** Quick / full / prod: 1,770,779 / 14,958,507 / 50,375,738
  nodes, every instance's nodes and hash equal to be625d8's, before and
  after the plain-root fix. prod took 13.76 / 13.94 s (13.43 s on the
  final binary) against be625d8's 13.86 / 14.14 s.
* **ctest.** `ctest -R fast_` passes all 49 tests. The new ones are
  `fast_msearch_dfirst_849_no_class_support`, a `--no-class-support` case
  in `fast_msearch_bad_args`, and `fast_class_support_plain_root`. Every
  changed file compiles with `-Wall -Werror`. arrange.c also compiles with
  `-Wextra -std=c17 -pedantic-errors`, as do msearch.c and
  class_support_test.c; for arrange.c this holds with ARRANGE_DEBUG,
  CLASS_PROF, CARRY_MAX_W=0, and without GFNI / VPOPCNTDQ or AVX-512BW.

**Paired measurements.** One binary, off / on / prototype alternated, min
of 2 (a1900 min of 5), `nice -n 10`. Each run is `msearch --diag-first
--diag-first-min-n 0 --d-stride k --d-offset o` with offset 13 (and 517 /
611 for the second offsets of c3648 / a2400, summed). CPU is the
process's (`cpu` of the dsum record). The load from the other workflow
was 1-3.

| P / S | N | labels | d | nodes off | nodes on | nodes on/off | off CPU-s | on CPU-s | prototype CPU-s | on / off | prototype / off | class share | class-free |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 13 7 4 3 1 1 / 1900 | 4,111 | 137 | 41 | 816,805 | 446,074 | 0.546 | 0.29 | 0.29 | 0.29 | 0.991 | 0.998 | | |
| 13 7 4 3 1 1 / 1950 | 5,896 | 153 | 30 | 1,449,997 | 838,119 | 0.578 | 0.59 | 0.50 | 0.51 | 0.851 | 0.869 | 30.4% | 0.59 |
| 12 6 3 2 1 1 / 988 | 6,671 | 159 | 34 | 1,969,071 | 993,003 | 0.504 | 0.86 | 0.70 | 0.71 | 0.815 | 0.833 | | |
| 13 7 4 3 1 1 / 2000 | 7,593 | 168 | 38 | 3,735,067 | 2,420,190 | 0.648 | 1.49 | 1.33 | 1.36 | 0.893 | 0.917 | 23.4% | 0.68 |
| 13 7 4 3 1 1 / 2100 | 11,306 | 192 | 29 | 5,387,694 | 3,254,658 | 0.604 | 2.46 | 1.99 | 2.10 | 0.810 | 0.854 | | |
| 12 6 3 2 1 1 / 1200 | 11,698 | 199 | 30 | 3,325,580 | 1,472,383 | 0.443 | 2.32 | 1.61 | 1.71 | 0.694 | 0.738 | 35.1% | 0.45 |
| 13 7 4 3 1 1 / 2200 | 15,199 | 211 | 31 | 11,426,894 | 7,205,201 | 0.631 | 7.41 | 5.65 | 6.02 | 0.763 | 0.812 | 14.7% | 0.65 |
| 11 6 4 3 2 1 / 2174 | 16,424 | 220 | 33 | 12,296,613 | 7,278,368 | 0.592 | 8.63 | 6.31 | 6.76 | 0.732 | 0.783 | 17.0% | 0.61 |
| 14 7 4 4 1 0 0 1 / 3648 | 20,538 | 252 | 42 | 10,656,143 | 4,968,486 | 0.466 | 9.55 | 6.30 | 6.76 | 0.659 | 0.708 | 24.6% | 0.50 |
| 9 6 4 3 1 1 1 1 / 2700 | 20,896 | 259 | 27 | 5,929,874 | 2,281,679 | 0.385 | 5.99 | 3.58 | 3.92 | 0.597 | 0.654 | 32.4% | 0.40 |
| 13 7 4 3 1 1 / 2400 | 22,993 | 245 | 39 | 25,476,926 | 14,930,420 | 0.586 | 20.79 | 14.51 | 15.41 | 0.698 | 0.741 | 14.4% | 0.60 |
| 13 7 4 3 1 1 / 2500 | 26,585 | 260 | 14 | 11,333,652 | 6,908,479 | 0.610 | 11.96 | 8.47 | 8.93 | 0.708 | 0.746 | | |
| 13 7 4 3 1 1 / 2650 | 31,743 | 279 | 11 | 10,221,345 | 5,634,085 | 0.551 | 12.64 | 8.31 | 8.66 | 0.657 | 0.685 | 11.6% | 0.58 |

Notes on the table:

* "class share" is CLASS_SUP's share of the search's cycles (the
  CLASS_PROF build, same samples).
* "class-free" is (on/off) x (1 - share): the ratio if the class support
  cost nothing, a bound no implementation reaches.
* The strides were 100 (a1900), 200 (a1950, b988, a2000), 400 (a2100,
  b1200), 500 (a2200, c2174), 1000 (c3648), 800 (c2700), 1200 (a2400),
  2000 (a2500) and 3000 (a2650).
* No sample had a pair.

* **Against the brief's targets.**
  * Required <= 0.80, measured: a2200 0.763, c3648 0.659, a2400 0.698,
    c2700 0.597, a2650 0.657.
  * Mean at N >= 20k: 0.664 (target <= 0.75).
  * a2000: 0.893 (target <= 0.95). a1950: 0.851 (target <= 1.00), so the
    rule needs no gate there.
  * Faster than the prototype on every sample: on/prototype 0.91-0.98,
    with 0.99 at a1900 (min of 5: 0.289 against 0.290 s; the first min of
    2 was 0.293 against 0.283, which is noise at 0.3 s).
  * At 15-23k: 0.60-0.76, mean 0.69, within 5% of the brief's
    zero-overhead estimate of about 0.73 everywhere. The worst is a2200
    at 0.763.
* **The gain grows with N.**
  * On 13 7 4 3 1 1 (ladder A, 8 sums, 4.1-31.7k): on/off = 0.971
    (N/4000)^-0.181 +- 0.017 (resid sd 0.033). Over all 13 sums: 0.947
    (N/4000)^-0.194 +- 0.032.
  * The node ratio has no trend (0.39-0.65, b = -0.03 +- 0.07).
  * What grows is the time each pruned (1,1) node would have cost. Its
    refutation costs ~N^1.17 per (1,1) node, while a call of the class
    support costs 6.7k cycles at 5.9k and 14.8k at 31.7k (~N^0.47).
  * Per node, on costs 1.16-1.82x off, falling with N on ladder A (1.82x
    at 4.1k, 1.19x at 31.7k).
* **Prototype.** prototype/off is 0.69-1.00 and on/prototype 0.91-0.99:
  the overhead removed here is 2-3% of the prototype's d-first time at
  5.9-7.6k and 4-9% at N >= 11k.

**d-first against the plain search, and the d-first law.**

The new d-first times are the integrated binary's per-sum d-first CPU
(section "Measurements on the integrated binary") times this on/off.

d-first / plain becomes:

| N | old | new |
| --- | ---: | ---: |
| 4.1k | 1.25 | 1.24 |
| 6.7k | 0.645 | 0.53 |
| 7.6k | 0.57 | 0.51 |
| 11.7k | 0.55 | 0.38 |
| 15-16k | 0.42-0.47 | 0.32-0.34 |
| 20.5k (252 labels) | 0.43 | 0.28 |
| 23k | 0.29 | 0.20 |
| 21-32k, > 256 labels | 0.32-0.43 | 0.19-0.31 |

Fits of t = a (N / 4000)^b CPU-s on the same 11 sums (4.1-31.7k):

| search | a | b | resid sd | t(16k) | t(32k) |
| --- | ---: | ---: | ---: | ---: | ---: |
| d-first, integrated (before) | 25.6 | 3.46 +- 0.13 | 0.27 | 3.1k | 34.3k |
| d-first with class support | 24.3 | 3.27 +- 0.17 | 0.34 | 2.2k | 21.6k |
| best of plain and the new d-first | 21.9 | 3.33 +- 0.17 | 0.34 | 2.2k | 22.3k |

* The residual sd grows because the gain is not uniform: 0.60 at c2700
  against 0.76 at a2200.
* (Written against be625d8, where v2 did not launch d-first units yet.
  On main, v2 launches them and fits a d-first time law and a d-first /
  plain ratio per engine, so the integration bumps msearch to engine 4
  and seeds engine 4's d-first law and ratio from these measurements:
  "Integration of the round-2 d-first changes". The plain search is
  unchanged.)
* CPU per expected (square, SP traversal) pair, best mode:
  * ~3.7-7.8k at 11.7-16.4k (was 5.4-10.6k);
  * 12.4k at 23k (was 17.8k);
  * 26.9k at 26.6k (was 38k);
  * 52.6k at 31.7k (was 80k).
* **Crossover.** Along 13 7 4 3 1 1, d-first / plain is 1.24 at 4.1k and
  0.66 at 5.9k (0.78 x 0.851), so it crosses 1 at ~4.6k (was ~4.9k).
  Pooled over the 11 sums, ln r = -0.062 - 0.773 ln(N/4000), N0 ~ 3.7k
  (3.9k before on the same sums).
* `--diag-first-min-n` stays 5000. Below it, the class support gains
  little or loses (next item), and the plain search's semi-magic squares
  feed the models.

**Negative results and limits.**

* **Small N.** On the gate (c) sums (full d lists, N = 0.45-3.6k), class
  support costs more than it saves. The d-first CPU is 1.36x off's in
  total (62.9 against 46.4 s; 1.0-1.68x per sum at N <= 3.0k) at
  0.12-0.74x the nodes. It is 0.90-0.97x at 3.1-3.6k, 0.99x at 4.1k and
  0.81-0.89x at 5.9-11.3k. A (1,1) node's subtree is cheap there, and the
  call's fixed costs (bounds, unions, keep masks, compaction) are not. The
  default `--diag-first-min-n` 5000 keeps d-first, and so the class
  support, above that. Gating the rule on |V_d| would be the fix if
  d-first is ever used lower.
* **The (1,2)/(2,1) layer** (the optional sub-step), on a2200 against the
  (1,1)-only time, min of 2:
  * The full fixpoint also at those children halves the nodes again
    (3,856,355 against 7,205,201), but takes 1.02x the time (5.75 against
    5.63 s). It makes 1.36M calls, 88% of them dead after 2.5 passes, and
    takes 30% of the search's cycles.
  * Capped at 1 / 2 / 3 passes there: 1.07x (5.43M nodes), 1.02x (3.98M),
    1.02x (3.86M).
  * The bar was <= 0.97x on a2200 and a2400, so it was not adopted. a2400
    was not measured once a2200 failed.
* Rejected organizations of the passes, all with the same nodes (see the
  table above):
  * Gauss-Seidel block order (fewer, costlier passes);
  * lazy per-group branches with union rescans;
  * a generic class loop (+9%);
  * starting on the other axis (+4%);
  * a histogram for the bounds;
  * a class-empty early death check (no gain: 18.0k against 17.5k cycles
    per call in the Gauss-Seidel build).
* Not done:
  * carrying the class bounds through the filters (<= 0.9%);
  * per-class unions out of the support filter (<= 0.6%);
  * the prototype's other levels (1 and 3-6: the scalar fixpoints and the
    fuller matchings), left out by the brief;
  * the matrix path (the V_d with more than 512 labels, or without
    AVX-512BW, get no class support).

## The star cover of the d loop (October 2026, branch c2/star)

The two diagonals of an n x n square of even order have no cell in common.
A magic square's n^2 numbers are distinct, so its two diagonals (both SP
traversals) have no number in common either. So for any number x, at most
one of them contains x, and the d-first loop finds the square from the
other one.

* The loop may skip every "star d", the d of the sum that contain a
  chosen x*. Every magic square is still found, once instead of twice when
  x* is on one of its diagonals.
* The square is still flagged: the partner test runs against the whole
  diagonal set, star d included.
* Only (square, SP traversal) pairs are lost. A systematic sample of the
  star d (every K-th, by rank in the unreduced list) estimates them.
* Not for odd n: the two diagonals share their center cell, which may be
  x*.

**Implementation** (src/c/dfirst.c, msearch):

* `dfirst_star_choose` picks x* from the whole unreduced list, so every
  `--d-range` unit and `--d-stride` sample of a sum agrees.
  * For every d, one pass over the postings of d's numbers gives |V_d| and
    the class sizes (R = the rarest class). No V_d is built.
  * x* = argmax over x of the sum, over the d containing x, of
    |V_d|^5.5 R^1.4, the per-d cost law of
    research/wip/complexity-profile.md. Ties go to the smaller x.
  * The CPU time is 0.036 s at N = 5.9k and 0.55 s at 31.7k: 0.03% and
    0.001% of the sum's d-first CPU. A d-range unit pays it again, which
    is still negligible. The fallback to the most frequent number was not
    needed.
* `dfirst_set_star(x, K, only)` sets the filter. With K >= 1 only the star
  d of rank % K == 0 are searched; with K = -1 none are.
* msearch options:
  * `--dfirst-star K` defaults to 4 for even n and to off for odd n.
    `--dfirst-star 0` gives be625d8's output. A nonzero K with odd n is
    refused.
  * `--dfirst-star-only` searches only the star d, for measurements. The
    sum is then not complete.
* The dsum record gets `star_x`, `star_k`, `nd_star`, `nd_star_skipped`,
  `nd_star_searched`, `pairs_star`, `nodes_star`, `cpu_star`,
  `pred_star_share`, `star_freq_rank`, `pred_share_top_freq` and
  `star_time`.
* The estimates (`est_pairs`, `est_nodes`, `est_time`) weigh a searched
  star d K x the d-stride weight. (On integ/round2 only `est_pairs` does:
  `est_nodes` and `est_time` are those of the run as made, and the
  K-weighted ones are `est_nodes_nostar` and `est_time_nostar`;
  "Integration of the round-2 d-first changes", "Fixes after the
  verification".)
  * With K = -1 they cover the d without x* only, and the schedulers do not
    count them.
  * `se_*` come from two strata (star d and the others), each a simple
    random sample of its d in the range. (Approximate: the star d are
    searched at a fixed rank % K, a systematic sample, which overweights
    the star stratum of `est_pairs` by K ceil(n/K) / n - 1 <= (K - 1) / n
    for n star d: 1.6x at n = 5, under 0.5% at the n of 600-1,150 of the
    sums measured.)
* `complete` still means that every magic square of the sum was found.
* Scheduler (both models): a complete sum's pairs count through
  `est_pairs`, once per sum (`scheduler.dsum_est_pairs`: totals
  `dfirst_pairs`, `dfirst_est_pairs`, `dfirst_est_sums`,
  `dfirst_star_sums`; SUMMARY_VERSION 4 in the prototype, 6 and then 7
  on integ/round2). A magic square is counted once,
  whether it came once or twice.

**Correctness.**

* (a) 39 small sums with N = 452-3,626 (cx/dfirst's gate-(b) list: 15 SP
  pairs, 2 of them on a star d) were each run with K = 0, 1, 4 and -1.
  * The K = 4 and K = -1 runs' "dsquare" records are byte-identical to the
    K = 0 run's records, filtered to the d without x* plus (K = 4) the star
    d of rank % 4 == 0. Star membership and ranks come from the K = 1 run's
    d log.
  * K = 1 equals K = 0: the same pairs, nodes and estimates.
  * The counts (nd + nd_star_skipped = N, nd_star_searched = ceil(nd_star
    / 4)) and est_pairs = pairs + 3 x pairs_star hold.
  * The default equals K = 4.
* (b) Magic exactness.
  * `fuzz_arrange --dfirst` now checks, per check:
    * x* against brute force (|V_d| and R recounted; the same weights
      summed in the same order; ties to the smaller number);
    * the filter with K = -1, a random K in 1-5, `only`, and a d-stride 3
      split, against the brute-force pairs filtered the same way, with the
      partner flags and the counts;
    * with x forced onto one diagonal of the planted magic pair (and not
      the other) and K = -1, the planted square must come exactly once,
      from the other diagonal, flagged.
  * 2,500 fresh seeds gave 0 fails. They covered: default modes 600;
    `--n 6` 300 (402 planted 6x6 magic pairs); `--n 4` 200; `--n 8` 100;
    mode 1 300; mode 7 200; the matrices build 300 plus mode 7 100;
    portable 200; no-GFNI 200.
  * Totals: 24,062 star runs, 737,723 kept pairs checked (198,642
    brute-force pairs, 111,404 flagged), and 4,062 planted magic squares
    with x* on one diagonal, all found and flagged.
  * Mode 5 (hubs) was stopped after 34 CPU-min inside its brute-force
    oracle. A rerun with `--budget 5000000` passed: 25 seeds ran (5 skipped
    by the oracle budget), with 48 planted squares.
  * Eight mutants of dfirst.c were each caught (14-120 fails in 60 seeds):
    the last number of d not tested, the rank offset, the exponent, ties
    to the larger number, the class count, `only` inverted, the skip
    counts, and the partner test made to ignore star d.
  * On real sums (a driver in the scratchpad, c2/star/real/star_real.c),
    222 magic pairs were planted in the box configuration: 3 per semi-magic
    square of 7 sums, N = 1.5-3k, the pairs added to the unreduced list.
    Two checks:
    * x* chosen on the augmented list (on a planted diagonal 66 times), K
      = -1, every planted square found and flagged;
    * x forced onto each planted diagonal: 444 checks, each found once from
      the other diagonal, flagged.
  * No result file holds a known magic square (none is known), so the
    "known magic squares" half of the gate is vacuous.
* (c) `--dfirst-star 0` is byte-identical to be625d8 on the 39 sums, except
  for the timing fields. Gates (a) and (c) were rerun on the final binary,
  with K = 4 as the default. arrange.c is untouched: bench quick / full /
  prod give 1,770,779 / 14,958,507 / 50,375,738 nodes, all ok.
* (d) `test_scheduler.py` `test_dfirst_star` runs synthetic star records
  through both schedulers.
  * A complete K = 4 sum: est_pairs, not pairs or K x pairs_star again.
  * A duplicate unit of it (K = 0), counted once.
  * A K = -1 sum: covered, its pairs found, no estimate.
  * A star-only part: not covered.
  * Magic squares that came once and twice: counted once each.
  * Fits unchanged, and incremental reads equal to a full read.
  * On msearch output: three `--d-range` units of a sum choose the same x*,
    and their counts, pairs and est_pairs add up to the whole sum's. K = 1
    equals K = 0, the default is K = 4, and odd n has no star cover.
* `ctest -R fast_` passes (48 tests). The new tests are
  `fast_msearch_dfirst_star` and five new refusals in
  `fast_msearch_bad_args`. dfirst.c and msearch.c compile with `-Wall
  -Wextra -Werror -std=c17 -pedantic-errors`.
* The independent x* check (scratchpad c2/star/gates/xstar.py: numpy from
  bin/enumerate's
  unreduced and reduced lists) agrees exactly at S = 849 and on a1950,
  a2000, b1200, c3648 and a2400: x*, its share, its star d and its
  frequency rank.

**Measured share.** s is the star d's share of the d loop's CPU.

* Method: per-d CPU from the d log, min of 2 alternating replicates, one
  heavy process at a time, ~10 CPU-min.
* The star d: all of them at 7.6k and 11.7k (`--dfirst-star 1
  --dfirst-star-only`), and 40-41 at 20-23k (K = 23 or 25).
* The other d: a d-stride sample of ~150 d.
* s = A / (A + B), with A the star d's CPU and B the other d's CPU
  estimated from the stride sample. The SE is by the delta method; A / T
  (T from the stride sample alone) agrees within 0.001.

| sum | N | x* (freq rank) | star d (share of d) | star d / mean d cost | predicted s | measured s | error | 1/(1 - 0.75 s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 13 7 4 3 1 1 / 2000 | 7,593 | 280 (1) | 627 (8.3%) | 1.14 | 0.0958 | 0.094 +- 0.005 | -2.0% | 1.076 |
| 12 6 3 2 1 1 / 1200 | 11,698 | 24 (4) | 707 (6.0%) | 1.20 | 0.0748 | 0.073 +- 0.004 | -2.0% | 1.058 |
| 14 7 4 4 1 0 0 1 / 3648 | 20,538 | 840 (8) | 1,018 (5.0%) | 1.36 | 0.0686 | 0.067 +- 0.005 | -2.1% | 1.053 |
| 13 7 4 3 1 1 / 2400 | 22,993 | 720 (31) | 919 (4.0%) | 1.51 | 0.0632 | 0.061 +- 0.005 | -3.9% | 1.048 |

* Every success criterion is met:
  * s >= 0.06 at 7-12k and >= 0.05 at 20-23k;
  * the prediction is within 2-4% (bar 30%), slightly high everywhere;
  * the choice of x* takes <= 0.03% of the CPU.
* **K = 4 is now the d-first default.** It searches a quarter of the star
  d, so the CPU per sum, and with it the CPU per expected magic square,
  falls by 1/(1 - 0.75 s): **1.076x at 7.6k, 1.058x at 11.7k, 1.053x at
  20.5k and 1.048x at 23k**.
  * K = -1 would give 1/(1 - s): 1.065-1.10x, but no estimate of the star
    d's pairs.
* The predicted share on all 12 sums of the round-2 profile is 0.104 at
  5.9k, 0.093-0.096 at 6.7-7.6k, 0.075-0.085 at 11-15k, 0.066-0.082 at
  16-21k and 0.062-0.068 at 23-32k.
  * It goes as ~0.110 (N/4000)^-0.28. With the measured 2-4% bias, K = 4
    saves 1.08x at 5.9k down to 1.048x at 31.7k.
  * So the gain is a constant factor, not a change of the N exponent. It
    shrinks slowly (+0.015 on the exponent of the time per sum), because
    x*'s share of the d falls (~N^-0.46) faster than its d's cost ratio
    rises (1.14 -> 1.51).
* The chosen x* is the most frequent number on 5 of the 12 sums, and of
  rank 4-31 on the others.
  * Choosing by cost gains up to 24% of share over the most frequent
    number (c3648: 0.069 against 0.056; a2650: 0.062 against 0.054). It
    gains nothing at a2400 (0.0632 against 0.0629, the most frequent
    number being 420).
* The pairs lose a little precision. A star d's pairs count K x, so the
  variance of est_pairs grows by about (K - 1) x the star d's share of the
  pairs (~1.2x at K = 4). Magic squares are unaffected.

**Limits (negative results).**

* The star is about the most an exact skip rule of this kind can take.
  * A skipped family of d is safe if no two of its members can be the two
    diagonals of a magic square. Without knowing the squares, that means
    pairwise intersecting families.
  * The other intersecting families, such as the d holding 2 of 3 fixed
    numbers or Hilton-Milner families, are much smaller than a star. 2 of
    3 numbers holds ~3 f^2 of the d, against f ~ 0.04-0.09 for a star.
    (By estimate, not measured.)
  * So the cover tops out at the largest star share, 6-10%.
* Two stars cannot be combined: a magic square with x1 on one diagonal and
  x2 on the other would be lost.

**For the integration** with the scheduler that launches d-first units
(e0ace54, after be625d8):

* Its d-first time law takes a part's share of the sum from
  `[d_lo, d_lo + nd)`. With the star cover nd counts only the d searched,
  so use `d_hi`.
* Its loop CPU is then the star-reduced cost. Either learn the law "as
  run" (every unit at K = 4; the old K = 0 records are 5-8% high), or
  scale K = 0 records by (1 - 0.75 pred_star_share), which those records
  lack (~0.07 at 20k).
* `dfirst_args` should pass `--dfirst-star 4` explicitly, so that a unit's
  cost does not depend on the binary's default.

CPU used: ~1.5 CPU-hours in all.

* Fuzz: ~55 min, of which mode 5 took 34.
* Gates: ~15 min.
* Measurements: ~10 min.
* x* checks: ~5 min.

## Pair rules on top of the class support (October 2026, branch c2/classhall: no-go at stage A, not built)

Question: the class support (`class_support = 1`, the level 2 of the
c2/vdclass prototype; "Class support in the V_d searches", branch
c2/classsup, code f46b203 and measurements a30a7fc) checks, at the (1,1) children of a V_d search, that
each candidate's cells and the classes it has to meet allow a matching:
every cell needs a class (the "bad" test) and every class needs a cell
(the "meet" test). The prototype's level 6 adds Hall's condition on pairs:
* cell side: two cells of the candidate that only the same class c covers
  (both in XO[c], the cells that class c alone covers outside the placed
  vectors);
* class side: two classes whose only cell of the candidate is the same
  cell.

On the same d, level 6 has 22-25% fewer nodes than level 2 (a2000
2,420,190 -> 1,876,913; a2200 7,205,201 -> 5,537,591; a2400 9,367,374 ->
7,060,610; 27% at a2650, 44% at c3648). It kills 2.6-3.3x more (1,1)
children. But the eager vectorized version (`CLASS_VEC=1 CLASS_LVL=6`)
is 1.24-1.58x slower.

The plan was a lazy version (`class_support = 2`): after the level-2
fixpoint, a masked test of 8 entries at a time against the OR of the XO[c]
and against the single-cell classes, with the exact test only on the lanes
that hit. Stage A was to decide by counting first:
* GO if the entries that need the exact test, (i) ∪ (ii) below, are at
  most 25% of the entries that pass level 2;
* and if the gated test is estimated at <= 5% of the d-first time.

**Answer: no-go on both counts.**
* The gate does not gate: (i) ∪ (ii), the entries that need the exact
  test, are 65-75% of the entries that pass level 2 (66-68% in the counts
  below, in all three sums), against the GO bar of 25%.
* The pair work is far over budget: to reach <= 0.90 at N >= 15k the
  whole pair stage may cost 8-15% of the level-2 time (below), while the
  eager pair work is 45-73% of it, i.e. 3-10x the budget. The lazy design
  saves at most half of that, and its estimate (~23% of the d-first time,
  against <= 5% wanted) is optimistic: it prices the exact test per entry
  that needs it, but the passes run on blocks of 8 entries, and with two
  thirds of the lanes hit nearly every block pays for all 8.
* The zero-cost ceiling is not the reason. With the pairs at no cost the
  d-first time would be 0.75-0.85 of level 2's (below), which is below
  the success bar (<= 0.90 at N >= 15k): the rules are worth having in
  principle; their cost rules them out.
* Stage B was not built.
* The instrumentation is in `research/wip/patches/c2-classhall/`: three
  patches on be625d8 (c2/profile's DPROF instrumentation, the c2/vdclass
  prototype, and the stage-A counts; `-DCLASS_SUP` only, default builds
  unchanged).

**Stage A counts.**
* Method: the scalar CLASS_FIX at level 6 (`CLASS_SUP=6 CLASS_LVL=6
  CLASS_VEC=0`), the (1,1) children, d-sampled at offset 13.
* Each pass tests every entry. The counts are per entry test, among the
  tests that pass level 2.
* The classes are the available classes other than the entry's. The cells
  are the entry's cells outside the placed vectors, other than its d cell.
* (i) a cell that exactly one class covers, i.e. a cell in some XO[c];
* (ii) a class with exactly one cell of the entry;
* (i2), (ii2): the actual violations, two such cells of one class or two
  such classes on one cell;
* (iii): the entries that the pair rules drop.

| sum (stride) | N | d | tests passing level 2 | (i) | (ii) | **(i) ∪ (ii)** | (i2) | (ii2) | (iii) dropped | calls with drops |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 13 7 4 3 1 1 / 2000 (200) | 7,593 | 38 | 126.7M | 51.0% | 47.2% | **66.4%** | 1.95% | 1.83% | 3.03% | 99.9% |
| 13 7 4 3 1 1 / 2200 (500) | 15,199 | 31 | 401.3M | 50.1% | 48.1% | **66.6%** | 1.92% | 1.85% | 3.05% | 99.97% |
| 13 7 4 3 1 1 / 2400 (1200) | 22,993 | 20 | 602.7M | 50.3% | 47.5% | **66.3%** | 1.87% | 1.77% | 2.96% | 99.99% |

* **Why so many.** At a (1,1) child, a candidate has k = 4 or 5 cells
  outside the placed vector, other than its d cell. It has exactly k
  classes to meet: k = 5 when its class is the placed vector's (37-54% of
  the tests), 4 otherwise.
  * So the matching is perfect and tight. A cell or class with a single
    partner is the norm, and an actual violation is rare: (iii) is 3%.
  * A cell covered by no class, or a class with no cell, is level 2's
    business.
* **The lazy order changes nothing.** `CLASS_STAGED=1` runs level 2 to
  its fixpoint, then the pair passes; nodes are the same. At a2000:
  * 82,895 of the 94,781 calls reach the pair stage;
  * that stage takes 7.1 passes per call (pair drops in 99.9% of the
    calls, which re-open level 2 on the other axis);
  * it makes 100.7M entry tests, against 134.1M in the level-2 stage;
  * (i) ∪ (ii) is 68.3% of the 94.4M tests that pass level 2.
* **Cheaper gates.**
  * (sc2), two single-class cells of any classes, is a one-mask popcount
    gate for the cell side, and it is 11.7% (a2000). But the cell side
    alone (level 3) gives -13.8% nodes at a2000 (2,085,773), against
    -22.4% for both sides.
  * The class side has no one-mask superset. "At most one cell per
    64-bit word" hits 80-87%.
* **Cost estimate** (a2200, by the go/no-go rule):
  * The level-2 pass of c2/classsup (its CLASS_PROF build, the same d)
    costs 9,234 cycles per call, 14.7% of the search. That is ~342 groups of 8 per
    call, ~3.4 cycles per entry slot.
  * The lazy pair stage needs ~1,080 gated tests per call (a2200's
    unstaged count x the a2000 staged / unstaged ratio 0.77).
  * At ~4x the per-entry cost of level 2, that is ~14.6k cycles per call,
    1.6x the whole level-2 pass: ~23% of the d-first time against <= 5%.

**The ceiling: the pair rules at zero cost.**
* Method: the prototype's vectorized level 2 against its eager level 6
  (`CLASS_VEC=1`), with TSC cycles around the class checks and the whole
  search (`clscyc`). Min of 2 alternating runs, the same d. The samples
  of c2/classsup: c3648 (14 7 4 4 1 0 0 1 / 3648) at offset 517, the
  others (a = 13 7 4 3 1 1) at offset 13.
* "rest" is the search without the class checks.
* The bound, (rest at level 6 + class checks at level 2) / (search at
  level 2), is the d-first time if the pairs cost nothing and level 2 cost
  no more than now.

| sum (stride) | N | nodes l2 -> l6 | search Gcyc l2 / eager l6 | class checks l2 / l6 | rest l2 -> l6 | bound | eager l6 / l2 |
| --- | ---: | --- | --- | --- | --- | ---: | ---: |
| a2000 (200) | 7,593 | 2,420,190 -> 1,876,913 (0.78) | 2.78 / 4.39 | 0.71 / 2.75 | 2.06 -> 1.64 (0.80) | 0.85 | 1.58 |
| a2200 (500) | 15,199 | 7,205,201 -> 5,537,591 (0.77) | 12.46 / 17.06 | 2.35 / 9.16 | 10.11 -> 7.91 (0.78) | 0.82 | 1.37 |
| c3648 (1000) | 20,538 | 2,027,116 -> 1,127,691 (0.56) | 5.69 / 7.97 | 1.74 / 5.44 | 3.94 -> 2.54 (0.64) | 0.75 | 1.40 |
| a2400 (1200) | 22,993 | 9,367,374 -> 7,060,610 (0.75) | 19.09 / 25.92 | 3.15 / 13.64 | 15.95 -> 12.28 (0.77) | 0.81 | 1.36 |
| a2650 (3000) | 31,743 | 5,634,085 -> 4,088,228 (0.73) | 17.95 / 22.18 | 2.70 / 10.73 | 15.25 -> 11.45 (0.75) | 0.79 | 1.24 |

* **The budget.** Hitting <= 0.90 needs the whole pair stage to cost <=
  8% (a2200), 9% (a2400), 11% (a2650) or 15% (c3648) of the level-2
  time. Hitting <= 1.02 at a2000 needs <= 17%.
* **The lazy design costs more than that.** The prototype's CLASS_VEC
  already runs the pairs only once level 2 is quiet ("the pairs last"),
  so the lazy design saves only the lanes the gate lets through:
  * the exact test still runs on 0.68 of the pair stage's lanes;
  * the gate test and the stage's extra level-2 passes run on all of
    them (at a2200, 2.58M passes against 1.35M at level 2);
  * an implementation like c2/classsup's might be ~1.34x cheaper than the
    prototype's (its level 2 takes 9,234 cycles per call against the
    prototype's 12,400 at a2200).
  * That is at best ~0.5 of the eager pair work. The eager pair work (the
    l6 - l2 class checks above) is 45-73% of the level-2 time, so the
    lazy version would still cost ~23-37%.
  * That gives ~1.0-1.1 at N >= 15k and ~1.2 at a2000.
* **Cell side alone.** Level 3 with the (sc2) gate would cost little, but
  at a2000 its ceiling is ~0.90: nodes x0.86 and rest ~x0.87 by the same
  elasticity, with level 2's class cost. So it is ~0.9 at best, and not
  tried.

**Takeaway.**
* The Hall pairs of the class matching kill real work: 22-44% of the
  nodes, 20-36% of the rest of the search.
* But they are a dense test. Two thirds of the candidates have a forced
  cell or class, and only 3% violate.
* So the test cannot be gated cheaply, and an exact per-entry test (a
  per-class popcount, or the per-label class masks with a table-driven
  Hall check) costs more than the pruning saves.
* This suggests that any further gain from the class structure must
  come from fewer candidates entering the (1,1) children, not from a
  stronger test there.

## Integration of the round-2 d-first changes (October 2026, branch integ/round2)

c2/classsup (f46b203, a30a7fc), c2/star (3e7519a) and c2/classhall (edfd553,
docs only) were cherry-picked onto main (189b49f: scheduler v2 with d-first
units), with the fixes their verifiers asked for. Measurements in this
section: one 4-core shared machine, `nice -n 10`, one heavy process at a
time, alternating arms of one binary, min of 2 rounds (scratch:
`integ2/gate/` and `integ2/e4/`).

### The gate of the class support

The class support was slower on small sums (gate (c) above: 1.36x the
d-first CPU at N <= 3k), and v2 may run d-first from N' ~3.4k on (and from
`--dfirst-min-n` 2000 with `--dfirst on`), below msearch's own 5000.

* **Measurement.** Per-d thread CPU from `--d-log`, class support on and
  off, the same d (`--d-stride`, offset 7, `--dfirst-star 0`), min of 2
  rounds per d, on 20 sums at N = 0.7-7.6k (gate (c)'s list, 13 7 4 3 1 1
  / 1900, 1950, 2000, 12 6 3 2 1 0 1 / 939, 951, 12 6 3 2 1 1 / 950).
* **What decides it is the V_d's label count, not |V_d|.** On / off by
  the V_d's labels (all 20 sums pooled):

  | V_d labels | d | on / off |
  | --- | ---: | ---: |
  | <= 116 | 13,958 | 1.465 |
  | 117-124 | 4,776 | 1.33-1.34 |
  | 125-128 | 1,773 | 1.118 |
  | 129-132 | 1,068 | 1.049 |
  | 133-136 | 2,856 | 0.997 |
  | 137-140 | 544 | 0.967 |
  | 141-160 | 4,929 | 0.85-0.93 |
  | > 160 | 3,369 | 0.884 |

  By |V_d| alone the two-word V_d (<= 128 labels) of 600-1,000 vectors
  were 1.36-1.43x, the wider ones of the same size 0.94-0.97x. A V_d of two
  words is cheap per node, and the class support's fixed costs per call
  (bounds, unions, keep masks, compaction) are not; the third word holds
  only the class labels at first.
* **The gate:** class support only in the V_d with at least 137 labels
  (`dfirst_set_class_min_labels`, msearch `--class-support-min-labels`,
  default `DFIRST_CLASS_MIN_LABELS` = 137; 0: every V_d). The label count
  is one pass over V_d's numbers (`vd_labels`, which the top-root choice
  already computes where the sum has few labels). Per sum, from the same
  per-d data, gated at L labels (on / off of the whole d loop):

  | gate | worst sum | s4031 (12 6 3 2 1 0 1 / 939, N 4.0k) | s4111 (13 7 4 3 1 1 / 1900) | N >= 5.7k |
  | --- | ---: | ---: | ---: | ---: |
  | none | 1.594 | 1.037 | 0.981 | 0.866-0.892 |
  | 129 (three words) | 1.035 | 1.035 | 0.981 | same |
  | 133 | 1.010 | 1.010 | 0.978 | same |
  | **137** | **1.000** | 1.000 | 0.992 | same |

  So no sum of the 20 is slower with the gate, the sums at N <= 3k are
  as without the class support, and nothing is lost above ~5k (all their
  V_d have more than 140 labels). The dsum record's `nd_class` counts the
  d whose search ran it (the gate passed, the top-label root on the
  carried path, a V_d of at least 2n vectors: `search_stats_t.class_used`;
  at first it counted every d that passed the gate, also under
  `--d-plain-root`, where the class support cannot run).
* Exactness: the class support is exact in any V_d, so any on/off pattern
  over the d gives the same pairs. `fuzz_arrange --dfirst` now runs the
  class support in every V_d on two seeds of three and, on the third, gates
  it at a random label count of the instance, which splits its V_d.
  `class_gate_test.sh` (ctest `fast_msearch_class_gate`) checks 16 5 4 2 /
  849 (102 labels: no class support, the nodes of `--no-class-support`)
  and d 0-100 of 13 7 4 3 1 1 / 1900 (V_d of 129-137 labels: some d with
  it; the gated nodes between off and on; the same pairs in all runs).

### Engine 4: the d-first time law and the ratio

msearch's `ENGINE_VERSION` and scheduler v2's `ENGINE` are now 4: the
d-first search changed (the gated class support and the star cover K = 4,
both on by default), the plain search did not (bench nodes and hashes
unchanged). Each engine's laws start from the previous engine's posterior
plus a shift (`amodel.time_prior`, and now `scheduler.current_ratio_level`
for the ratio: the engines' pairs are chained, each engine's level the
posterior of its own pairs under a prior at the previous engine's level
plus the shift).

**Paired measurement** (`integ2/e4/`): one binary, engine 3's settings
(`--no-class-support --dfirst-star 0`) against engine 4's (defaults; gate
at 129 labels in this run, which differs from 137 only on s4031 and s4111,
by ~1-3%), the same d (`--d-stride`, offset 7; whole d loops below 4k),
min of 2 alternating rounds, d-first CPU per sum (dsum `time`: the index
and the d loop), pairs equal in all arms:

| sum | N | labels | d searched | engine 3 s | engine 4 s | 4 / 3 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 16 5 4 2 / 849 | 1,515 | 102 | 1,515 | 0.55 | 0.51 | 0.933 |
| 10 6 4 2 1 1 / 855 | 2,444 | 110 | 2,444 | 3.42 | 3.07 | 0.899 |
| 12 6 3 2 1 0 1 / 900 | 2,994 | 121 | 2,994 | 6.07 | 5.19 | 0.854 |
| 10 6 3 2 1 / 1360 | 3,066 | 194 | 3,066 | 1.96 | 1.66 | 0.846 |
| 10 6 3 2 1 / 900 | 3,626 | 154 | 3,626 | 6.11 | 5.09 | 0.833 |
| 12 6 3 2 1 0 1 / 939 | 4,032 | 133 | 4,032 | 25.65 | 24.64 | 0.961 |
| 13 7 4 3 1 1 / 1900 | 4,111 | 137 | 2,055 | 15.40 | 14.07 | 0.914 |
| 12 6 3 2 1 1 / 950 | 5,704 | 149 | 1,426 | 22.37 | 18.49 | 0.826 |
| 13 7 4 3 1 1 / 1950 | 5,896 | 153 | 1,179 | 22.35 | 19.00 | 0.850 |
| 13 7 4 3 1 1 / 2000 | 7,593 | 168 | 633 | 22.01 | 17.79 | 0.808 |
| 12 6 3 2 1 1 / 1200 | 11,698 | 199 | 293 | 23.91 | 15.27 | 0.639 |
| 13 7 4 3 1 1 / 2200 | 15,199 | 211 | 102 | 24.33 | 17.50 | 0.719 |
| 14 7 4 4 1 0 0 1 / 3648 | 20,538 | 252 | 103 | 25.31 | 15.83 | 0.625 |
| 13 7 4 3 1 1 / 2400 | 22,993 | 245 | 58 | 30.37 | 21.46 | 0.707 |
| 13 7 4 3 1 1 / 2650 | 31,743 | 279 | 32 | 31.85 | 20.18 | 0.634 |

* Below 3k the star cover alone (the class support is gated off there)
  gives 0.85-0.93; with the gated class support alone (third arm, no star)
  the sums up to 5.9k were 0.90-1.04x engine 3 (0.96-1.04 where no V_d
  passes the gate: the run-to-run noise), and the star cover on top
  0.84-0.95x.
* **Fit**, N >= 3k (12 sums): ln(engine 4 / engine 3) = -0.138 (+- 0.029)
  - 0.162 (+- 0.027) ln(N/4000), residual sd 0.073 (N >= 4k: -0.097 -
  0.191 x; >= 5.5k: -0.138 - 0.164 x). This was the first shipped shift.
  Its verification found two flaws: the two P = 10 6 3 2 1 sums (S 900
  and 1360) lie outside v2's model grid (S 368-729 for this P), where the
  class support starts at a smaller N (154-194 labels at N 3.1-3.6k), and
  the sampled runs' `time` holds the whole x* choice and index but only
  1/d_stride of the loop (1.5-2.6% too high at 20-32k). Refit on N' over
  the 10 sums inside the grid: -0.099 - 0.194 x.
* **The c2 verifiers' numbers** (the class support's 0.924
  (N/4000)^-0.180, 11 sums at 5.9-31.7k, times the star cover's 1 - 0.75 s
  with s = 0.107 (N/4000)^-0.28) linearise over 4-32k to -0.161 - 0.162
  ln(N/4000).
* **The perf verifier's paired runs** (`integ2/perf/`): the 189b49f binary
  (B, engine 3) against the integrated one (E), identical d samples
  (`--d-stride k --d-offset 3`), whole-sum CPU = d_stride x the d loop +
  the fixed costs (not `est_time`), min of 2 alternating rounds, 13 sums;
  `N3` = the integrated binary with `--no-class-support --dfirst-star 0`,
  `CG` = with `--dfirst-star 0` (the class support alone):

  | sum | N | N' | E / B | N3 / B | CG / B |
  | --- | ---: | ---: | ---: | ---: | ---: |
  | 13 7 4 3 1 1 / 1880 | 3,485 | 3,444 | 0.900 | 0.990 | 0.987 |
  | 12 6 3 2 1 0 1 / 960 | 4,414 | 4,467 | 0.924 | 0.999 | 0.981 |
  | 10 6 4 2 1 1 / 930 | 4,685 | 4,352 | 0.910 | 1.002 | 0.996 |
  | 13 7 4 3 1 1 / 1950 | 5,896 | | 0.841 | | 0.891 |
  | 12 6 3 2 1 1 / 988 | 6,671 | | 0.813 | | 0.890 |
  | 13 7 4 3 1 1 / 2000 | 7,593 | | 0.780 | | 0.852 |
  | 12 6 3 2 1 1 / 1200 | 11,698 | | 0.676 | | 0.717 |
  | 13 7 4 3 1 1 / 2200 | 15,199 | | 0.701 | | 0.771 |
  | 11 6 4 3 2 1 / 2174 | 16,424 | | 0.730 | | 0.791 |
  | 14 7 4 4 1 0 0 1 / 3648 | 20,538 | | 0.589 | | 0.665 |
  | 9 6 4 3 1 1 1 1 / 2700 | 20,896 | | 0.558 | | 0.609 |
  | 13 7 4 3 1 1 / 2400 | 22,993 | | 0.662 | | 0.733 |
  | 13 7 4 3 1 1 / 2650 | 31,743 | | 0.652 | | 0.652 |

  N3 matches B within 0.99-1.03 with identical nodes on all 13 sums. Fit,
  x = ln(N'/4000): ln(E / B) = -0.105 (+- 0.032) - 0.210 (+- 0.028) x,
  residual sd 0.069; the class support alone -0.022 - 0.215 x. **Pooled**
  with the integrator's 10 in-grid sums (23 pairs): **-0.103 (+- 0.023) -
  0.204 (+- 0.020) x, residual sd 0.067.**
* **Shipped** (after the verification): `amodel.DFIRST_E4` = (-0.103,
  -0.204) (sds in `DFIRST_E4_SD`), the shift of both the d-first time law
  (`DFIRST_ENGINE_SHIFT[4]`: intercept 3.215 -> 3.112, slope 3.576 ->
  3.372, sd 0.25 kept) and the ratio (`DFIRST_RATIO_ENGINE_SHIFT[4]`: a0
  -0.028 -> -0.131, a1 -0.566 -> -0.770, sd 0.17 kept), since the plain
  search is unchanged (`ENGINE_TIME_SHIFT` has no entry for 4). With the
  7% calibration stream `--dfirst auto` switches to d-first at N' ~3.68k
  (engine 3: 4.29k; the first shift gave 3.49k). Against the verifier's
  13 sums: their own d-first law 3.181 (0.15) + 3.331 (0.13) x, rms 0.29
  around the shipped one (sum t / sum pred 1.00); ratio (E / plain)
  -0.113 (0.072) - 0.805 (0.062) x, rms 0.15 around the shipped one.
  Predicted engine 4 / 3: 0.927 at N' 3.5k, 0.831 at 6k, 0.599 at 30k
  (the first shift: 0.890, 0.816, 0.629; measured at 4.4-4.5k: 0.91-0.92).
* The d-first law's slope (engine 3's, shifted) fits the ladder (mostly 13
  7 4 3 1 1) but ran about 2x high for the calibration search's pool P at
  N' >= 12k (7 sums; joint fit of the 27 sums 3.37 + 3.23 x before the
  shift). The online refit learns only the intercept, so above 12k the
  planner charges pool P about 2x; that is 2% of E at 1-10 CPU-years and
  the "measured" forecast truth corrects it in the forecast. Not refitted.
* The scorer carries its engine (`AnalyticScorer.engine`), so that the
  slope matches the level. `current_ratio_level` chains the engines' pairs
  (above); at first it restarted engine 4 from its shipped prior once it
  had one pair, which dropped what engine 3's pairs had taught (a pair at
  the handed-over level moved it by 0.03-0.08); a test now checks that
  such a pair leaves the level unchanged.
* The calibration search in the target region measured engine 3's ratio
  level at about -0.08 (14 d-first sums, research/calibration-target.md
  3.6). It is not folded into the prior: the scheduler learns the level
  online from the streams, as before.

### The star cover with d-range units

c2/star was written before v2 launched d-first units. Its verifier's
required fixes, on main's units (research/scheduler-v2.md, "The star
cover"): dchunk records carry star_x and star_k; `Summary._dcover` keeps one
coverage per (x*, K) and never merges chunks of different x* or K (two x*
could each skip one diagonal of the same magic square); star-only chunks
(and star chunks without star_x) count towards no coverage; v2 passes
`--dfirst-star 4` to every d-first unit and, to a unit that continues a
sum, the K and x* (`--dfirst-star-x`, new) of the parts it continues, and
records both; a covered sum's pairs and their estimate come from its
group's chunks, each weighted by the share of its range new to the group,
so that duplicates count once and disjoint chunks add up to msearch's own
est_pairs; the d-first law takes a part's span as its d searched plus those
the star filter skipped, and learns from engine 4 only with v2's K.

The share weighting is exact for disjoint chunks and duplicates, which is
all the planner makes (it plans only the gaps of a group). A chunk that
partly overlaps its group (a hand-made unit) adds its pairs by the share of
its range that is new, as if they were spread evenly along d; on S = 589,
P = 12 6 3 2 1 with a duplicate unit 7 d longer the sum's est_pairs came
out 4.056 against msearch's 4 (dchunk records carry no per-d pairs).

### Fixes after the verification (integ/round2)

Three verifiers (correctness, performance, adversarial) passed the
integration with issues; one major, the rest minor. Each was checked first:

* **S-traversal prior (major, holds).** The first integration shipped
  band-marginal ratios (N' >= 3k +0.04 -> -0.023, a +0.036 6-12k effect)
  and called them the write-up's refit; the refit (calibration-target.md
  section 5, glm.out) has N' >= 3k +0.032 and 6-12k -0.018 and puts the
  3-6k deficit into k <= 5, ratio <= 1.1 and x < 0.1. The S refit alone
  would leave the shipped priors a mix (it lowers E; the P refit, not
  shipped, raises it: E at 10 CPU-years 0.251 against 0.283 with the whole
  refit), so now all three class-factor GLMs (squares, S, P) ship at the
  refit's posterior, every effect (means and sds; the squares' 6-12k
  effect is the same -0.078).
* **Engine-4 shift (holds).** Set to the perf verifier's pooled fit,
  (-0.103, -0.204), sds (0.023, 0.020) (above, "Engine 4"); the switch
  moves from N' 3.49k to 3.68k.
* **The ratio handoff (holds).** The ratio level is now chained across
  engines like the time laws (above).
* **`est_time` of a star run (holds; three verifiers).** It weighed a
  searched star d K x, i.e. estimated the CPU without the star cover (1.04-
  1.31x the run's own whole-sum CPU on d-stride samples). Now `est_time`
  and `est_nodes` are those of the run as made (each searched d weighs
  d_stride), with `se_time` over the star d actually searched (every K-th),
  and `est_time_nostar` / `est_nodes_nostar` keep the K-weighted figures.
  `est_pairs` is unchanged. Nothing in the scheduler used them;
  calibration-target/analyze.py's `t_dfirst_whole` would have.
* **calibrate.py and the star cover (holds; two verifiers).** Its d-first
  section weighed a dsquare by d_stride only, so on engine-4 output it
  undercounted the star d's pairs and the magic squares with x* on one
  diagonal. It now reads star_x / star_k / star_only from the file's dsum
  (or dchunk) records and weighs a dsquare whose dvec holds x* d_stride x
  K; K = -1 and star-only records count in the pairs, not in the
  estimates (`records_not_estimated`). New test `test_cli_star`; on
  msearch's own output for S = 589, P = 12 6 3 2 1 the weight is 4, as its
  est_pairs.
* **`nd_class` (holds; two verifiers).** It counted every d that passed the
  gate, also where the class support cannot run (`--d-plain-root`, the
  matrices path, V_d under 2n vectors). Now `search_stats_t.class_used`
  says whether it ran; `class_gate_test.sh` checks `--d-plain-root` (0).
* **Mixed engines in one ratio pair (holds).** A sum's d-first CPU for its
  ratio pair now takes only the parts of the engine of its first part
  (SUMMARY_VERSION 7).
* **A vacuous e2e assertion (holds).** `"3:dfirst" not in s.time` became
  vacuous with the engine bump; it now checks that no `:dfirst` key exists.
* **Partial overlaps in `_dcover` (documented, not changed).** The share
  weighting is exact for disjoint chunks and duplicates, which is all the
  planner makes; for a hand-made partly overlapping unit it is
  approximate (above). An exact count would need per-d pairs in the
  dchunk records.
* **Systematic star sample (documented).** The star d are searched at a
  fixed rank % K; the stratum's weight K overweights it by at most (K -
  1) / n for n star d (< 0.5% at the sums measured), and `se_*` are
  approximate. Left as is: a weight of n / ceil(n / K) per unit would
  make the units' estimates no longer add up to the whole sum's.
* **The d-first slope above 12k (documented).** Not refitted (above).
* **Measured truth on a learned state (documented)** in `--truth`'s help
  and `forecast_truth`: the online refit learns only the d-first level and
  one plain offset at N' >= 6k, so a learned state's "laws" stay above the
  measured CPU at N' >= 12k.
* **A_SQ 6 and PHI_SUM (documented)** next to `A_SQ`: the measured 0.41 is
  the spread beyond Poisson, so it also holds the per-sum overdispersion
  that PHI_SUM models (A ~8-9 without it); 6 is the refit's choice.
* **Docs**: the stale "0.12 / 0.25 / 0.39" forecast and "engine 3" text in
  README and scheduler-v2.md, the version number in the star section, the
  blank line before its heading, the gate default in a fuzz_arrange
  comment, the `--truth anchored` pointer of `forecast` (now `measured`),
  and the claim that engine 4's forecast equals the write-up's.
* **Plain-path timing (re-measured).** One verifier saw full at 1.016
  (contended); see the gates below.

### Gates of the integrated build (integ/round2)

After the verification fixes (the first integration's gates, on 6bc3124,
had the same outcome: 51/51, fuzz 0 fails on 3,900 seeds, bench nodes ok):

* Build: `cmake -S . -B build-r2 -D CMAKE_BUILD_TYPE=Release`; dfirst.c,
  msearch.c and arrange.c also compile with `-Wall -Wextra -Werror
  -std=c17 -pedantic-errors`.
* `ctest -R fast_`: 51 of 51 pass (`fast_msearch_class_gate` now also
  checks `nd_class` 0 under `--d-plain-root`).
* `fuzz_arrange`, fresh seeds (9600000-9670500), 0 fails: `--dfirst` 600,
  `--dfirst --n 6` 500, `--dfirst --mode 7` 500, `--variants cls` 500,
  plain mixed modes 600, `--mode 7` 500, and `--dfirst` 500 each on the
  matrices, portable and no-GFNI builds (3 seeds skipped by the harness).
  The d-first runs check every seed without the class support, with it in
  every V_d or gated at a random label count, and the star cover (K = -1,
  a random K, star-only, stride splits, planted magic pairs with x* on one
  diagonal): 326k (square, d) pairs, 44k star runs, 7,260 planted magic
  pairs with x* on a diagonal found and flagged.
* bench quick / full / prod: 1,770,779 / 14,958,507 / 50,375,738 nodes,
  hashes ok. Against 189b49f's build on an idle machine (alternating, 6
  rounds): min 0.992 / 0.982 / 0.998, median 1.045 / 1.000 / 1.002 of its
  time; the plain search is unchanged (a verifier's full at 1.016 was
  contended; the new `class_used` field of `search_stats_t` costs nothing
  measurable).
* `python3 scripts/test_scheduler.py` (48 s) and `test_calibrate.py`: all
  ok (new: `test_cli_star`, the chained ratio level, the refit priors).
* The paired engine 3 / engine 4 runs found the same pairs in every arm.
