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
