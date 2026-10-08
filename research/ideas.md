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
