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
