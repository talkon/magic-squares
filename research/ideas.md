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
