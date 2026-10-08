# Complexity in N: ideas, picks, prototype results (work in progress)

Raw results of the complexity-in-N workflow (profile in complexity-profile.md). Patches: patches/cx-*.

## Ideas

### Make the diagonal-first (traversal) decomposition the large-N magic search: its advantage over the plain search grows like N^-0.5

**description**

This decomposes the search by a vector d that will be a diagonal. For every d, run the existing semi-magic search on V_d = {v : |v & d| = 1}, then check the diagonals. This is already in the repo as src/c/dsearch.c.

What I added is a d-sampled driver (run every k-th d): research/decomp/dsample.c on branch cx/decomp, commit fd93d09, worktree /home/user/magic-squares/.claude/worktrees/cx-decomp. With stride 1 it reproduces dsearch exactly (622,221 nodes on 10 6 3 1 0 1, S = 391; at S = 836 of 12 6 3 2 1 0 1 it gives dsearch's 1 SP pair).

I measured it on the large-N sums, CPU-s per sum. Outputs are in /tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/cx/decomp/runs/ (ds_*.txt; same-binary plain baselines in pl_*.log).
- The ratio d-first / plain is ~0.47 at N = 20k (fit) and 0.43 at 23k (vs pooled plain).
- It finds every magic square: its two diagonals are SP traversals, so each is found twice.
- It does NOT output semi-magic squares in general. Only the ~1e-3 of them that have a vector (SP) traversal come out: 0 of them in all my large-N samples.

I also tried the next level of the same decomposition, fixing both diagonals (research/decomp/dsample2.c, random disjoint pairs). It is a negative result and should not be pursued:
- At N = 23k there are 2.2e8 disjoint pairs, with |V_{d1,d2}| ~ 489 and ~573 nodes per pair.
- Nearly all of that is the root loop over r1, which costs O(|V|^2) per pair.
- Estimated total: 1.27e11 nodes and ~1.0e5 CPU-s, about 2x the plain search, and still ~N^4.

**complexity_argument**

Measured d-first / plain CPU per sum. Plain is same-binary r1-sampled where marked, otherwise the profile's sampled totals with a ~5% overhead correction. d-first SE is 4-5% from 69-152 sampled d.

| sum | N | d-first CPU-s | plain CPU-s | ratio |
|---|---:|---:|---:|---:|
| 13 7 4 3 1 1, S = 2000 | 7.6k | 394 | 511 (same binary) | 0.77 |
| 12 6 3 2 1 1, S = 1080 | 8.9k | 462 | 656 | 0.74 |
| 13 7 4 3 1 1, S = 2100 | 11.3k | 1,505 | 2,484 | 0.64 |
| 13 7 4 3 1 1, S = 2200 | 15.2k | 5,236 | 10,300 (same binary, SE 15%; pooled ~10,800) | 0.48 |
| 14 7 4 4 1 0 0 1, S = 3648 | 20.5k | 8,086 | 16,830 | 0.50 |
| 13 7 4 3 1 1, S = 2400 | 23.0k | 17,817 | 38,900 (same binary, SE 18%; pooled with the two earlier samples ~41,500) | 0.43 |

- Node ratios agree: 0.47 at 23k (2.77e10 vs 5.95e10) and 0.48 at 15k.
- log-log slope of the ratio: -0.52.
- On ladder A, d-first time goes like N^3.49 and plain like N^4.03 (7.6k-23k).
- The V_d searches have the same bitset width as plain (236-245 of 245 labels at 23k), so this is not a width artifact.

Why it scales better, grounded in the profile:
- Profile section 5: an r1's cost grows like U^3.0-3.3 with its universe and is steeply super-linear. At fixed list sizes, larger-N sums have more work per square.
- d-first trades N roots for universes that shrink relative to N, because P(|v & d| = 1) falls as the label count grows: |V_d|/N = 0.217, 0.199, 0.196, 0.180, 0.156, 0.156 from 7.6k to 23k, ~N^-0.3.
- That is the measured mechanism; the N^-0.5 rate is the empirical fit.
- Per unit of magic expectation, cost then grows like ~N^(1.8-0.5) = N^1.3 instead of the profile's N^1.8. It is not output-sensitive.
- Per-d cost is nearly uniform: the largest single d is 0.76 s against a 0.35 s mean at 15k. So it samples and parallelizes trivially.

**expected_gain_at_N20k**

About 2.1x per sum at N = 20k (ratio 0.47 fitted, 0.48-0.50 measured at 15-20k), 2.3x at 23k, ~2.7x projected at 32k and ~3.4x at 50k if the N^-0.52 trend holds. It multiplies with any speedup of the inner arrange search (V_d runs the same code). That holds unless such a speedup depends on list sizes that V_d does not reach.

**cost_to_prototype**

The measurement is done: ~2.5 CPU-min of d-sampled runs plus ~2.5 CPU-min of plain baselines. Productionizing is an msearch mode:
- For each d, build V_d. It is O(N x 36) scalar now; a label -> vectors index would make it O(|V_d|).
- Call search_vectors on V_d. Per-d setup measured ~0.6 ms, 0.06 s per 100 d.
- Keep dsearch's diagonal check.
- About a day of work, plus a sampled plain run per sum if the semi-magic statistics are still wanted.

**risks**

- **Output changes.** Only magic squares (and the rare squares with an SP traversal) are found. The forecast and existence models need semi-magic counts and traversal rates, so plain r1-sampled runs must continue alongside, at a few percent of the CPU.
- Each magic square is found twice. Pruning on 'a second diagonal with larger index exists' is weak at the (2,2) layer: I estimate ~3 candidates there at 23k.
- **Uncertain denominators.** The ratio's denominators carry plain-sampling SE of 12-18% at 15-23k. Only ladder A has more than two points, and two other P agree (0.74 at 8.9k, 0.50 at 20.5k).
- Above 256 labels V_d keeps nearly all labels, so it also hits the matrix path (see idea 2).
- The trend beyond 23k is extrapolated.
- This is the collaborator's method. What is new here is that its advantage keeps growing past N = 6k (0.53-0.66 there on 14 7 5 3), and the negative two-diagonal result.

### Per-label-block bitset width (local universe of each r1), including the carried path for r1 below label 256 in sums with more than 256 labels

**description**

The r1 subproblem only involves labels <= x, where x is r1's largest label: vectors are sorted by label and every later vector has labels <= x.

The proposal: run SEARCH_ROOT per block of r1 with the bitsets re-packed to W_local = ceil((x+1)/64) words. Use the carried lists whenever x < 256, and the N x N matrices only for the top blocks.

I also tested whether relabelling deeper helps. LOCAL=2 build on cx/decomp counts the distinct labels in the lists, weighted by subtree time, at the (1,1) and (1,2)/(2,1) nodes. It does not:
- The local label set is essentially all labels <= x.
- Time-weighted mean: 210 of 245 labels at N = 23k, 180 of 211 at 15k.
- At most 192 labels for only 14% / 72% of the time, and the A nodes look the same.
- So per-r1 is the right granularity.

**complexity_argument**

The profile attributes N^0.3 of the N^4.1-4.3 growth to bitset width (x1.38 per extra word at identical nodes), plus x2.7 for the matrix path above 256 labels.

All numbers below come from the profile's own per-r1 logs (label and TSC cycles per r1), so no new runs were needed.

Time ratio from per-r1 width at x1.38 per word:

| sums | time ratio |
|---|---|
| a1900, a1950, a2000, a2100 | 0.73, 0.87, 0.96, 1.00 |
| a2200, a2300, a2400 | 0.78, 0.90, 0.96 |
| b920, b988, b1080, b1200, b1320, b1480, b1600 | 0.78, 0.94, 0.99, 0.73, 0.83, 0.92, 0.96 |
| c2146, c2174, c3648 | 0.98, 0.83, 0.97 |

The gain is largest just after L crosses a multiple of 64.

Above 256 labels:
- The r1 whose label is among the top 23 labels hold only 11-30% of the time on 7 sums with L = 211-252; the top 32 hold 39-64%.
- So for L ~ 279 (S = 2650, N = 31.7k), 70-89% of the time would run carried.
- Estimated time x (0.11-0.30 + (0.70-0.89)/2.7) = 0.37-0.63.

In growth terms:
- The time-weighted local label count trails L by ~31-35 labels.
- L grows ~N^0.34 on ladder A, so the 256-label wall moves from N ~ 25k to ~36k there.
- The width steps of section 8 are smoothed per block instead of hitting the whole sum at once.

**expected_gain_at_N20k**

1.0-1.37x on carried sums (mean ~1.12; e.g. 1.28x at a2200, 1.11x at a2300, 1.04x at a2400). 1.6-2.7x on sums with 257-~290 labels, which most sums with N >= 25k of 13 7 4 3 1 1 are. The gain fades above ~300 labels. It applies unchanged inside the V_d searches of idea 1.

**cost_to_prototype**

Small:
- A loop over label blocks in search_vectors / SEARCH_ROOT.
- s->bits re-packed per block, O(N W).
- Carried arrays allocated at the block's W.
- No new search code.
- Validation: ctest fast_, bench quick/full/prod, plus a > 256-label fuzz case.
- About half a day.

**risks**

- The x1.38/word factor was measured on whole sums at identical nodes; per block it may differ.
- The adaptive cross mode (xs_tick) would restart per block, so node counts may differ slightly; the squares will not.
- It is a constant factor, not an exponent change: it only delays the width steps.
- It overlaps with lifting the 256-label limit directly (carried W = 5-8). If that is done, this still saves 1.0-1.37x.

### r1-level cached relation for the long axis: build each (2,2) child's exactly-once list from its own entries instead of scanning the parent's long list

**description**

Shared sub-results across all (1,1) and A nodes of one r1.
- Build once per r1, lazily per row, the bitmask over cols1 positions of the cols meeting each row r exactly once (Lnk(r1, r)), plus cols1 x cols1 disjointness.
- A (2,2) child created by r2 under the A node (r1; c1, c2) gets its o-list as Lnk(r2) & Dc(c1) & Dc(c2).
- It then gathers only its ~157 kept entries (vpcompressb index extraction plus vpgatherqq per word).
- Today it scans the 1,300-1,400-entry long list and keeps 11% (filter o).
- After that, cross support, the support filter and the counts run unchanged on the carried lists. This is profile item (a), done via an r1-level cache instead of a per-node index.

Negative result that bounds this direction (prototype -DLOCAL11 on cx/decomp, fd93d09):
- I replaced the whole subtree below each (1,1) node by a search in node-local position space: exactly-once and disjoint bit matrices per (1,1) node, children by ANDs, forward checking and MRV via label -> position masks.
- It is correct: quick.txt and full.txt hashes all ok.
- But it is 3.3x slower at N = 7.6k (S = 2000, 38 r1) and 4.8x slower on full.txt. Forward checking alone lets 8x more nodes through (47.8M vs 6.0M; dead layers (3,3) 10.6M, (2,3) 5.9M), at ~300 cycles per node.
- Setup is 0.9M cycles per (1,1) node, 19% of the baseline's per-(1,1) cost there.
- Rows-first or cols-first branching orders do not help (24.7-29.7M nodes on a half sample).
- The label-union support test, which is what kills the (2,3) layer, has no cheap position-space form.
- So only this hybrid survives: masks to create lists, carried label tests after.

**complexity_argument**

Profile:
- Filter o is 21% of the cycles at N = 23k, and grows with the parent list (N^0.88; 2,441 cycles per (2,2) child in the rdtsc build).
- With the cache, the per-child cost follows the child's own list (N^0.53) plus three mask ANDs of |cols1|/64 words (29 words average at 23k).
- So the k12 x c22 term's filter-o part goes from ~N^1.35 to ~N^1.06 per A node.

Microbenchmark (research/decomp/mbench.c; W = 4, random 6-sets, keep 13%; TSC cycles per call):

| entries in the parent list | filter now | mask + gather |
|---:|---:|---:|
| 330 | 451 | 204 |
| 1,038 | 1,273 | 427 |
| 1,416 | 1,836 | 603 |
| 2,500 | 3,605 | 1,449 |

Setup per r1:
- |rows1| x |cols1| = 17M pair tests at 23k on average, at ~0.5 cycles per pair (measured in the prototype's builder): ~8.5M cycles.
- That is ~0.25% of the r1 subtree (~3.6e9 cycles per r1 at 23k).
- Plus 3.3M pairs for cols1 disjointness.

Saving: ~1,230 cycles x 29.5 children = 36k of ~290k plain cycles per A node at 23k, i.e. 12%. It grows with N, since the long list grows faster than the child's.

**expected_gain_at_N20k**

About 1.10-1.12x at 20k, ~1.13x at 23k, ~1.15-1.2x projected at 32k. More if the fixpoint of profile item (b) removes the survivors' share, which makes (2,2) creation a larger fraction. Multiplies with ideas 1 and 2.

**cost_to_prototype**

Medium (2-3 days):
- Carry a cols1-position index alongside the carried lists, or a mask of the A node's col list.
- A TRY_CHILD variant for the (2,2) children.
- Lazily built Lnk rows (2.2 MB per r1 at 23k as masks; L2-resident only for the active part).
- First check: replay the profile's (2,2) children through the gather path to confirm the microbenchmark in place.

**risks**

- The mask o-list is a superset of today's support-filtered list (support removals at the A node are not reflected), so cross support works on a few more entries.
- Gathers depend on cols1 staying in L2: 58 KB per word-array at 23k, 3,000-entry cols1 for the heavy r1.
- ideas.md found label -> position masks losing at N ~ 2-3k (LABEL_MASKS). The microbenchmark says this variant is already 2x at 330 entries, but in-place results there were much smaller than microbenchmarks before.
- It is a constant factor, not output-sensitivity.

The bigger lesson from the prototype: below the (1,2)/(2,1) layer, every local test that works on unions (forward check, cross, cell or pair support) weakens as lists grow. Position-space decompositions make nodes cheap but not fewer. None of the three ideas makes the work below an A node N-independent; idea 1 is the only one with a measured exponent change.

### Count-only pre-test: refute the five-deep children before materializing their lists (exact)

**description**

In TRY_CHILD, a child with 5 or more vectors placed (nearly always a (2,3)/(3,2) child of a searched (2,2) node) first runs FILTER_COUNT. That is FILTER_CARRY without the output: the exactly-once test on the other axis, then the disjoint + support test on the same axis, computing only the kept counts and a masked-OR union, with no compress and no store. It then applies the same count and forward checks that follow the real filters. If they fail, it returns; otherwise the normal FILTER_CARRY path runs. It kills exactly the children those checks would kill after the real filters, so the nodes and squares are unchanged: identical node counts on quick/full/prod and on every sampled run, fuzz_arrange 150 instances (1,221 squares, 0 fails) and --mode 6 pass. Prototype: -DPRETEST (env PRETEST_MIN, default 5) on branch cx/prune, commit 3dc9b05, worktree /home/user/magic-squares/.claude/worktrees/cx-prune; FILTER_COUNT plus ~25 lines in TRY_CHILD, compiled out by default. Next step (creation angle): a label->position index of each (2,2) survivor's lists, built once and shared by its k22 = 15-18 children, so the o-count and union become a few mask operations instead of a scan.

**complexity_argument**

The profile's fastest-growing term is the dead subtrees under the (2,2) survivors, k12*p22*k22*c23 (~N^2.2 per (1,2) node). Measured clean (sampled plain build, TSC, ladder A), it is 31% / 43% / 47% / 48% of all time at N = 5.9k / 11.3k / 15.2k / 23k: 12k to 136k cycles per searched (1,2) node, ~N^1.8. At 12 6 3 2 1 1, S=1080 it is 27%. That is more than the rdtsc profile's 10-32%. Of those (2,3)/(3,2) children, 96-99% die at exactly the checks after their two filters (profile: fc_o 18-29%, cnt_b 55-77%, fc_b 4-12%). Per 8 entries, FILTER_CARRY spends about 9 test operations, then per list 4 vpcompressq (2 uops each, port 5), 4 stores and 4 ORs. The pre-pass keeps the test and replaces the rest with 4 masked ORs. Derived from the 23k timing, the survivors' subtree cost falls by about 1/3 (136k to ~92k cycles per (1,2) node). Because the saving is a fixed fraction of a term whose share of the time grows with N, the gain grows with N: x0.935 at 5.9k, ~x0.91 at 15k, ~x0.845 at 23k (systematic r1 samples), and x0.92 / x0.81 / x0.79 on the heaviest r1 (plain builds, first 25M nodes). It does not change that term's exponent (still k22 x list per survivor). The per-(1,2)-node cost from 5.9k to 23k goes from ~N^1.46 to ~N^1.40.

**expected_gain_at_N20k**

13 7 4 3 1 1: CPU x0.89-0.93 at S=2200 (N=15.2k) and x0.83-0.86 at S=2400 (N=23k), sampled over all r1. 11 6 4 3 2 1 at S=2174 (N=16.4k): x0.85-0.89. Low-survival P: 14 7 4 4 1 0 0 1 at S=3648 (N=20.5k) x0.97-0.98 sampled (x0.855 on its heaviest r1), and 12 6 3 2 1 1 at S=1080 x0.98-1.01. Plain builds on the first 25M nodes: x0.81 (S=2200), x0.79 (S=2400), x0.92 (S=1950). Neutral at production sizes with the same nodes: full.txt 3.32-3.33 s vs 3.28-3.45 s, prod.txt 13.4-13.8 s vs 13.3-13.9 s. In short, 10-20% at N ~ 20k on the P/x that hold most squares, growing with N.

**cost_to_prototype**

Done: ~70 lines, measured (in-situ timings above, min of 2 alternating runs per variant). Productionizing on the carried path is about half a day (CMake/ctest variant, README). The matrix path (>256 labels, most sums with N > 25k) needs its own count-only version of filter_list + KEEP_COUNT, about 1 day. Alternatively it comes for free if the carried path's 256-label limit is lifted.

**risks**

The gain depends on the (2,2) survival p22 (high for long lists and x <= 0.5, low at high x or on some P). The shared machine adds +-3% noise (2 alternating runs each). The 1-4% of children that pass the pre-test pay one extra pass pair. This is a constant factor on the fastest-growing term, not output-sensitivity: the (2,2)-children creation (52-62% of the time after it) is untouched. An index-based (2,3) creation from another angle would subsume it. Under -DSCALE_PROF the early return skips some instrumentation counters, but the plain-build check confirms the gain.

### Iterate the (2,2) cascade to its fixpoint (or one more round) on its survivors, gated by N

**description**

After cross support and the support filter at a (2,2) child, a survivor reruns CROSS (both passes) and SUPPORT until no list changes, or for one extra round only. This is sound: every rule is a necessary condition, and no sampled productive node was killed. Prototype: -DSURV_FIX in cx/prune (env SURV_FIX_ROUNDS, SURV_FIX_MIN; FORCE_AFTER for the cascade order), combinable with -DPRETEST.

**complexity_argument**

It is the one pruning whose kill rate does not decay with N. On searched (2,2) nodes the cross+support fixpoint kills 95 / 91 / 86 / 90 / 92% at N = 4.1 / 5.9 / 11.3 / 15.2 / 23k, removing 92 / 83 / 77 / 82 / 82% of their subtree cycles. At the (2,2) children it cuts survival from 35-40% to 3.5-5% at N >= 11k. But each round rebuilds the 8 cell unions over ~150-entry lists and costs about 2/3 of the subtree it removes (estimated ~5.8k vs 8.8k cycles per survivor at 15.2k), and both scale with the same list lengths. So the net stays a constant rather than growing: nodes -55 to -75%, but CPU x0.93-0.95 at N 9-23k, x0.98 at 5.9k, x1.00-1.05 at 2.5k. One extra round matches the full fixpoint (x0.94 vs x0.92-0.96 at 15k) with half the node reduction. With idea 1 in place it adds 0-4% (pre-test + fixpoint vs pre-test alone): x0.83 vs 0.845 at 23k, 0.885 vs 0.91 at 15k, 0.935 vs 0.976 on 14 7 4 4 1 0 0 1, 0.96 vs 0.995 at S=1080, but 0.955 vs 0.935 at 5.9k.

**expected_gain_at_N20k**

Alone 5-8%. On top of idea 1, 0-4% (most where survival is low but the survivors are fat: 14 7 4 4 1 0 0 1, 12 6 3 2 1 1 at S=1080).

**cost_to_prototype**

Done (~15 lines). Needs a gate against the small-N loss. A list-length gate (sum of lists >= 150) was worse than always-on at 15k (x0.958 vs 0.93), so gate on N (> ~8k) or on a sampled survival estimate like the existing cross_after switch.

**risks**

Gains near the noise level (+-3%). The cost per round grows with list length. The obvious large-N reordering, the support filter first ('after' mode forced), is worse: x1.04-1.05 at 15k and x1.13 at 23k (x0.93-1.00 with the fixpoint), so the cascade must stay cross-first.

### Exactly-once (cell/pair) support from (1,1)-level relation matrices, at the (2,2) children and cost-adaptively at shallower layers: measured, does not scale (not recommended for pruning)

**description**

This is the angle's natural candidate, tested before anyone builds it. Once per (1,1) node, build the exactly-once relation between its row and col lists (bit matrices over list positions; every descendant's filter o then becomes an AND). Use it for 'cell support' (each candidate meets some other-axis candidate through every unmatched cell exactly once) and 'pair support', which are the exactly-once refinements of cross support that the profile said union tests cannot see. Apply them at the (2,2) children, and, when lists are long, at the (1,2)/(2,1) nodes. Also measured: magic-only diagonal pruning (a magic square needs two disjoint vectors meeting every placed vector exactly once). Analysed: nogoods across siblings, dominance between r1 subproblems, matching relaxations. Method: plain-C reference implementations iterated to their fixpoints, run on copies of the lists of sampled nodes (-DPRUNE_DIAG in cx/prune), N = 4.1-23k, 900-7,000 nodes per point.

**complexity_argument**

The build itself scales well: one orientation costs 8.7% / 6.0% / 3.9% / 3.9% of the (1,1) subtree at N = 4.1k / 5.9k / 15.2k / 23k (0.46-0.67 cycles per pair; 1.85k x 1.78k lists at 23k). The rules weaken with N as cross support does, because the expected number of exactly-once partners per (candidate, cell) also grows (m*q = 14*0.21 = 2.9 at 4.1k, 35*0.11 = 3.9 at 23k; P(none) ~ e^{-mq}). (2,2) children alive after one round + support: cross 15 / 25 / 38 / 40 / 37% vs cell 10 / 19 / 33 / 35 / 33% at N = 4.1 / 5.9 / 11.3 / 15.2 / 23k. At the fixpoint: cross 0.75 / 2.1 / 5.1 / 3.8 / 3.5%, cell 0.22 / 0.96 / 2.7 / 2.2 / 2.3%, all rules 0.06 / 0.40 / 1.2 / 1.2 / 1.5% (x25 from 4.1k to 23k). So no local rule makes the (2,2) layer output-sensitive. At searched (1,2)/(2,1) nodes, the cross/cell/pair fixpoints kill 3.2 / 5.6 / 1.2% of nodes at 4.1k but 0.6 / 0.8 / 0.0% at 15-23k, removing <= 0.6% of the subtree cycles at 4-6k and 0.000 at N >= 11k. They also leave the MRV count (the number of (2,2) children k12) at 0.98-1.00x. These nodes are locally consistent, more so as N grows, so no shallow test can cut k12, the factor that sets the cost per square. Magic-only: (1,2) nodes have 68 / 86 / 111 possible diagonals at 5.9k / 15.2k / 23k, and 99.9% of their children keep a disjoint pair (19-28 candidates each), so diagonals bite only ~7-8 vectors deep. Nogoods: a dead (2,2) node {r1,c1,c2,r2} cannot recur in its (1,1) subtree, because every other c2 shares c2's r1 cell and every other r2 shares the branching cell, and its death involves all four vectors, so there is nothing to reuse. Matching/Hall relaxations: the per-class domains are ~320 at (1,1) and ~8-10 (16 classes) at (2,2), so violations need near-empty classes (0.7% of (2,2) children at prod per ideas.md, fewer with longer lists).

**expected_gain_at_N20k**

For pruning: <= 2-3% at the (2,2) layer, and less once idea 1 is in. The cell/all-rule fixpoints remove another 6-10% of the survivors' subtree cycles over cross, for 1-2k-cycle passes. Zero at the (1,2)/(2,1) layer; ~0.1% for diagonal pruning. Not recommended as pruning. The relation matrices may still pay as a creation structure (filter o as an AND, 3.9% build cost per orientation at 15-23k), which is the other angle.

**cost_to_prototype**

The measurement is done (~350 lines of diagnostics, compiled out by default). A real implementation must carry list positions through every compress site (~6) plus the two matrix orientations, about 2-3 days, for a measured ceiling of a few percent.

**risks**

The kill rates come from fixpoints on list copies: sound (0 productive nodes killed), but only a few productive nodes fell in the samples. One-round kill rates depend on the pass order. Samples are 1/128 to 1/2048 of nodes. The negative conclusion covers union-based, exactly-once, pairwise, diagonal and matching rules at up to 4 placed vectors, not higher-order 4x4 tests; the opt2/residual solver results suggest those cost as much as the subtrees they replace.

### 1. Diagonal-first (V_d) search as the large-N production mode (recommended; measured 1.9-3x at N = 15-32k)

**description**

For each SP vector d, run the existing semi-magic search on V_d = { v : |v & d| = 1 }. This finds every (square, SP traversal) pair exactly once, so it finds every magic square (twice, once per diagonal) and nothing else of use. Above a crossover of N ~ 6k, use it instead of the plain search plus diagonal check. The code is essentially there: bin/dsearch, and my driver src/c/dsample.c on branch cx/magic (commit fbe1805, worktree .claude/worktrees/cx-magic). Building V_d with bitsets and the per-d setup is 0.2-0.9% of the time.

Measured d-first vs plain CPU per sum, from d-sampled runs:
- 13 7 4 3 1 1: 1.26x at N = 4.1k, 0.75 at 7.6k, 0.60 at 11.3k, 0.47 at 15.2k (4,944 vs 10,444 s), 0.40 at 23k (17,170 vs 42,433 s), 0.32 at 31.7k (119,600 vs 375,000 s; both on the matrix path, 279 labels).
- Other P: 12 6 3 2 1 1 at N = 6.7k: 0.75; at 11.7k: 0.66. 11 6 4 3 2 1 at 16.4k: 0.39. 14 7 4 4 1 0 0 1 at 20.5k: 0.53.
- Per-d cost is homogeneous: the top 10% of d hold 17-20% of the time; d-first SE 3-5% from 58-206 d.
- Plain at S = 2200 was re-measured on the same 76 r1 as cx/profile (10,444 vs 10,217 s).

Follow-ups:
- The dominant phase of d-first is the (1,1) -> (1,2)/(2,1) child creation: 14 children per node, each scanning the parent lists of ~370 + 360 entries and keeping 19% / 33%. The profile's idea (a), a label -> position index of the parent's lists, applies there directly.
- Keep a plain r1-sample (5-10% of CPU) per large-N sum for the calibration stream.

Data: scratchpad/cx/magic/runs (summary.tsv, d*.log, dsc*.json, dscp2200.json).

**complexity_argument**

Ladder A, 4.1k -> 23k: d-first time ~N^3.45 against plain ~N^4.1. The ratio falls like N^-0.67 (N^-0.52 between 15k and 32k).

Why, from the SCALE_PROF layer counts of the V_d searches (cx-profile instrumentation, scratchpad/cx/magic/runs/dsc*.json):
- |V_d| = 0.26N at 4k and 0.14N at 32k, so the lists stay short. At 15k the (1,2)/(2,1) lists are 119 against 589 in the plain search, with k12 = 3.8 against 24.6.
- With short lists the filters stay strong: **no (2,2) child survives in V_d at any N measured**. At 15k they die at the filters: fc o 35%, count b 21%, fc b 31%, support 12%. The plain search keeps 37% of them, and its dead subtrees below those survivors are 29-34% of its time. That term, the profile's fastest-growing one (~N^2.2 per (1,2) node), is gone.
- The price is redundancy. Each partial structure is re-explored under every d compatible with it: d-first has 12.5x the plain search's (1,2) nodes at 15k (1.44e9 vs 1.15e8) and 325x the (1,1) nodes.
- d-first (1,2) nodes grow ~N^3.1, its (2,2) children made ~N^3.45 (4.6e9 at 15k, against plain 2.84e9 (2,2) plus 1.6e10 (2,3)/(3,2) children), and cycles per node 600 -> 1,100 (the W step at 192 labels).

rdtsc profile at 15k:
- (1,2)/(2,1) creation 54%;
- (2,2) creation 19% (480 cycles per child);
- (1,1) creation 9%;
- (1,0) 2%.

Not output-sensitive: the output is ~0.4 (square, d) pairs per sum at 15k, so all the work is proving emptiness. The cost per expected magic square still grows ~N^1.1-1.2 instead of ~N^1.8.

**expected_gain_at_N20k**

Measured d-first/plain at N ~ 20k: ~0.42 interpolated on 13 7 4 3 1 1 (0.47 at 15k, 0.40 at 23k), 0.39 on 11 6 4 3 2 1 (16k), 0.53 on 14 7 4 4 1 0 0 1 (20.5k). So 1.9-2.5x lower cost per expected magic square at N = 20k, ~3x at 31.7k (0.32), and ~1.3x at 7-12k. Below N ~ 6k it loses (1.26x at 4.1k, 4-14x at N <= 1k). Applied to the existence frontier, whose cells near E = 1 have N ~ 21-22k, this cuts the CPU for one expected magic square by roughly 2-2.4x (central ~7,300 -> ~3,000-3,700 CPU-yr with the a = 3.92 time model). The cut is larger where the matrix path applies. Plain-side sampling SE is 22-34% at N >= 20k, so the ratios at 23k and 31.7k are +-30%.

**cost_to_prototype**

Prototype done: src/c/dsample.c (cx/magic, fbe1805), with output validated against dsearch: identical nodes on 13 6 3 2 at S = 517, and the 1 SP pair of 13 5 3 2 0 1 at S = 632 found. Total measurement CPU ~5 min. Production is ~1 day: an msearch mode (--diag-first, or automatic above N ~ 6k) that loops over d (in d-chunks as scheduling units, which are uniform in cost), outputs (square, d) records with the partner check (square_diag_stats), a fuzz/ctest case (plant a magic square: it must be found twice), and a refit of the scheduler's time model for this mode.

**risks**

- The gain is a constant factor that grows slowly with N (N^-0.67). It is not the output-sensitive cost per square the profile asks for.
- Each magic square is found twice. No cheap way to skip half of the d exists: the d left out would have to form a pairwise-intersecting family, at most ~3% of N.
- Its edge comes from having no dead (2,2) survivors. If the plain search gets the fixpoint or other (2,2) fixes from the other angles, plain gains 1.5-2x at 15-23k and the ratio shrinks. Re-measure then.
  - Conversely, d-first's main cost, the child creation over long parent lists, is the profile's idea (a), so it benefits from that work too.
- It no longer emits the ~1.7e-4^-1 x more semi-magic squares that calibrate the squares-per-sum and traversal models. Keep a plain r1-sample.
- The plain baselines at N >= 20k carry 22-34% sampling error. At 31.7k both runs are on the matrix path. The ratio after lifting the 256-label limit is unmeasured.
- Possible extension, not measured: restricting to magic squares whose largest label lies on a diagonal finds 1/3 of them at ~0.22 of the cost (u^3.45 averaged over the universes). That is ~1.5x per found square, but it gives up 2/3 of each sum's yield.

### 2. Both diagonals first (V_{d1,d2}), or second-diagonal / box-partner pruning inside V_d (measured; not recommended)

**description**

Fix two disjoint SP vectors (d1, d2) as the diagonals and search V_{d1,d2} = { v : |v & d1| = |v & d2| = 1 }. Each vector then has a type (a, b) in d1 x d2, and a magic square is 3 boxes (2 rows x 2 cols with d1 on the diagonal and d2 on the anti-diagonal of each box), i.e. pi^-1 pi' a fixed-point-free involution.

Measured with src/c/dpair.c (cx/magic): d1 sampled, 40-60 random disjoint d2 > d1 each, the generic search on V_{d1,d2}, 13 7 4 3 1 1.

| N | disjoint pairs | V_{d1,d2} | search-only CPU-s per sum | x plain | x d-first |
|---:|---:|---:|---:|---:|---:|
| 4.1k | 6.1e6 | 243 | 245 | 6.8 | 5.4 |
| 7.6k | 2.27e7 | 324 | 2,282 | 4.3 | 5.8 |
| 15.2k | 9.4-9.6e7 | 425-438 | 16,500-17,000 | 1.6 | 3.4 |
| 23k | 2.29e8 | 496 | 55,000 | 1.3 | 3.2 |

Building V_{d1,d2} and the per-pair setup double these totals (45,600 s at 15k).

The variant that keeps one diagonal at the root and prunes by the partner d2 inside V_d was also estimated. The candidate partners at a V_d (2,2) node are vectors disjoint from d that meet the 4 placed vectors once each: ~0.84 N x 0.16^4, about 8 at 15k and 11 at 23k. The set empties only at (3,3), and V_d searches never get past (2,2). ideas.md's 7% from diagonal pruning near S_min is consistent with this.

**complexity_argument**

- Pairs grow ~N^2.07 (0.83 N^2/2) and |V_{d1,d2}| ~N^0.45.
- Generic nodes per pair ~|V|^1.05 (275 -> 514 from 4.1k to 15.2k): the work per pair is dominated by its root loop.
- Total search time ~N^3.15, against d-first N^3.45 and plain N^4.1. So the ratio to d-first falls only like N^-0.3, from 5.4x to 3.2x between 4.1k and 23k, and would reach parity near N ~ 10^6.
- The partner structure can only cut deeper: an involution is 15 of 265 derangements (5.7%), but it constrains only structures with >= 3 vectors placed. A box-structured search therefore cannot remove the per-pair root work, which is O(|V_{d1,d2}|) to O(|V|^2) per pair, i.e. ~N^2.5-2.9 overall with a large constant.
- The redundancy grows with depth: each partial structure is compatible with ~N x 0.16^l diagonals (l vectors placed), and with N^2 x 0.16^(2l) diagonal pairs. Fixing more diagonal structure at the root multiplies the partial structures faster than it shortens the lists.

**expected_gain_at_N20k**

Negative: search-only it is ~1.4x plain and ~3.3x d-first at N = 20k (interpolated from 1.6x / 3.4x at 15.2k and 1.3x / 3.2x at 23k). Total with per-pair setup ~2.8x plain. A tailored implementation (precomputed V_d bitsets, box structure) would need to be >3x cheaper per pair than the generic one just to tie d-first at N <= 32k.

**cost_to_prototype**

Measurement prototype done (dpair.c, ~1 CPU-s per estimate). A tailored box-structured d-pair search would be ~2-3 days of work. I do not recommend it, given the redundancy argument and the measured slopes.

**risks**

- The generic V_{d1,d2} search is not the best possible d-pair search. A partner-structured one could be a few times cheaper, but most of the per-pair cost is root-level work that the partner structure does not touch.
- The estimates use 1,000-2,000 sampled pairs per sum (random d2), so they are noisy but unbiased.
- The ideas.md note ("same redundancy with a worse constant", 15x at N ~ 900) holds qualitatively. The gap narrows with N (6.8x -> 1.3x vs plain) but not fast enough to matter for N <= 32k.

### 3. Select (P, S) or squares by the mod-3 partner-pair congruences (measured; refuted), and other selection levers

**description**

heuristic.md shows that a partner pair can be magic only if u.(S, a) = 0 (mod 3) for every u orthogonal to the merged-row contrasts. When t = (S, v_p(P)) = 0 (mod 3), i.e. P a perfect cube and S = 0 (mod 3), nothing is forbidden. The lattice local CLT then predicts the full 3^d boost for every (square, tau): ~9x (5 primes) to ~27x (6 primes) over generic P, which would make cube P (e.g. 12 6 3 3 3, 6 6 3 3 3 3, S_min 1.4-2.6k) prime targets.

Tested at the level of residues, which is exact and model-free. For each (square, tau), count the sigma with Z1 + Z2 = 2t (mod 3) on the exponent coordinates (Y) and on all coordinates (XY). Script: scratchpad/cx/magic/py/mod3dist.py.

| squares | Y hit rate, as a multiple of 3^-k | XY hit rate, as a multiple of 3^-(k+1) |
|---|---|---|
| 400 generic ladder squares, 5 primes | 0.75-1.0 | 0.72-1.0 |
| 400 generic ladder squares, 6 primes | 0.6-0.8 | 0.67-0.85 |
| 24 cube-P squares (6 6 3 3 3, 9 3 3 3 3) | 0.46-1.2 | 0.46 at S = 0 mod 3 (vs 9x if boosted) |

On the cube squares at S = 0 (mod 3), the target residue is reached in 1% of the (square, tau), not in all of them. Their 12-cell product events (scratchpad/cx/magic/py/cube12.py): 24 events in groups of 8, i.e. 3 independent ones, against ~3.5 expected at the generic rate and 9-15 under the boost.

So the congruences are zero-sum. They move the probability between (square, tau) but raise it for no class of (P, S). Cube P are also poor in squares: 19 and 5 squares in 120 s and 89 s at N <= 2.7k.

I also checked selection within a sum by r1 universe (cheapest r1 first). From cx/profile's per-r1 logs, the cost per square is flat across U-quartiles: 0.93-1.16 relative in a1900, a1950, b920 and b988. There is no lever there either.

**complexity_argument**

- A necessary condition for magic is Z1 + Z2 = 2t (mod 3).
- Averaged over (square, tau, sigma), this holds with probability (0.6-1.0) x 3^-(k+1): the uniform value, not the boosted one.
- So the probability per square cannot gain on average from the congruences, and cube P show no boost. The conditional share given a hit (hit x 3^dim ~ 2.5) is offset exactly by how rarely the target residue is hit at all (1-12% of (square, tau)).
- The heuristic's 'allowed' test (lattice membership) is much weaker than an actual hit: it calls 100% of cube pairs allowed, but only 1-6% hit.
- Side result for the calibration work: the heuristic's net factor of 1.37 on E[magic] from the congruences is not supported at this level (~0.7-1.0). That lowers its magic estimates by ~1.4-2x.

**expected_gain_at_N20k**

None. Cube-P targeting gives no boost: measured 0.46-1.2x the uniform residue rate against a predicted 9-27x. The congruences cannot prune the search either: the merged rows need every column of the square, so the condition can only be checked on finished squares, where the diagonal check is already cheap. Within-sum selection by universe: none (cost per square flat across U).

**cost_to_prototype**

Done: ~3.5 CPU-min of msearch for the cube squares and ~1 CPU-min of Python. To firm up the cube depletion (0.46x, from 7-24 squares), ~10 CPU-min more of 6 6 3 3 3 / 9 6 3 3 3 would do.

**risks**

- The cube sample is small: 24 squares, of which 7 have S = 0 (mod 3). But the refutation is strong: the target residue is hit in 1 of 105 (square, tau) at S = 0 (mod 3), where the boost hypothesis needs most of them.
- The residue test assumes the exact-hit probability given the residue matches is the same for cube and generic squares (the same Gaussian part). Nothing suggests otherwise.
- I did not touch the heuristic itself. Its 1.37 congruence factor should be re-examined by its authors with this residue-level count.

## Picks

{
 "picks": [
  {
   "name": "Pick 1: diagonal-first (V_d) as the large-N production mode (msearch --diag-first), with a plain r1-sampled calibration stream",
   "spec": "GOAL: lower the cost per expected magic square at N >= ~6k by 2-3x. Merges idea 1 (cx/decomp) and idea 7 (cx/magic).\n\nBRANCH: cx/dfirst from fdb77fc, worktree /home/user/magic-squares/.claude/worktrees/cx-dfirst. Port src/c/dsample.c (cx/magic fbe1805) and the r1-sampling hooks SAMPLE_STRIDE/SAMPLE_OFFSET/SAMPLE_LOG (cx/profile 53252f1).\n\nIMPLEMENTATION\n1. msearch --diag-first [--diag-first-min-n N0, default 6000; auto mode switches on when reduced N >= N0].\n   - Loop d over the UNREDUCED SP vector list (a diagonal need not survive reduction, as in dsearch).\n   - Build V_d = {v in reduced list : |v & d| = 1} through a label -> reduced-vector posting index: bump a per-vector hit counter along d's 6 posting lists, keep the vectors with count == 1, then clear only the touched counters. Cost is O(sum of posting lengths), not O(36 N).\n   - Keep V_d as a subsequence of the reduced order, so label order and root symmetry breaking are unchanged.\n   - Call search_vectors(V_d). For each square, run square_diag_stats and emit a (square, d) JSONL record with sp_count, best_score and a magic flag.\n   - A magic square appears twice, once per diagonal: dedupe downstream by canonical square hash. Do not try to prune the second copy.\n2. Scheduling and sampling.\n   - Process d in contiguous chunks (~256 d) as scheduling and checkpoint units: per-d cost is homogeneous (top 10% of d hold 17-20% of the time).\n   - Keep dsample's --d-stride k / --d-offset o for unbiased per-sum estimates, and log per d: |V_d|, labels, nodes, squares, setup and search seconds.\n3. Calibration stream: --calib-r1-stride k runs the plain r1-sampled search at about 5-10% of the sum's d-first CPU. It emits semi-magic counts and traversal stats tagged with the stride, so the forecast and existence models keep their semi-magic inputs.\n   - Also report per sum the count of (square, SP-traversal) pairs from d-first. That is an exact count of the quantity magic depends on, a better calibration target than semi-magic counts x p_S p_P.\n4. Explicitly do NOT build:\n   - d-pair / V_{d1,d2} (measured 1.3-1.6x plain, 3.2-3.4x d-first);\n   - second-diagonal / partner pruning inside V_d (the candidate set empties only at (3,3), and V_d never gets past (2,2));\n   - the 'largest label on a diagonal' restriction. It keeps 1/3 of the yield at 0.22 of the cost. With the frontier's E(C) log-slope ~0.37, E_new(C) = E_old(4.5C)/3 ~ 0.58 E_old(C); it would need a slope >= 0.73 to break even.\n\nCORRECTNESS (all must pass)\n(a) At d-stride 1 it reproduces dsearch exactly:\n   - 622,221 nodes on 10 6 3 1 0 1, S = 391; identical nodes on 13 6 3 2, S = 517;\n   - the SP pairs of 13 5 3 2 0 1 / 632, 16 5 4 2 / 849 and 12 6 3 2 1 0 1 / 836.\n(b) Differential test against plain on >= 30 small sums (N 500-3000, e.g. those of bench quick/full): the multiset of (canonical square, d) from --diag-first must equal {(sq, t) : sq from the plain search, t a vector of the unreduced list that is a traversal of sq}.\n(c) New fuzz_arrange mode on abstract label sets:\n   - random instances plus random candidate diagonals; brute-force all (square, traversal-in-set) pairs and compare;\n   - planted squares with two traversals realizable as diagonal + antidiagonal must be reported twice and flagged;\n   - >= 150 instances, 0 fails.\n(d) ctest -R fast_ and the bench quick/full/prod hashes unchanged for the plain path.\n\nMEASUREMENT (same binary, same session, alternating; <= ~40 CPU-min)\n- d-first: d-sampled, relSE <= 5% (70-150 d).\n- plain: r1-sampled with STRATIFIED r1 (strata = r1-index quartiles, rate ~ the stratum's per-r1 cost, e.g. 4x denser in the first quartile, which holds ~90% of the time), relSE <= 10%. Today's plain denominators carry 15-34% SE.\n- Sums:\n  - 13 7 4 3 1 1 at S = 1900 / 2000 / 2200 / 2400 (N 4.1k / 7.6k / 15.2k / 23k);\n  - 12 6 3 2 1 1 at S = 988 (6.7k; full counts known);\n  - 11 6 4 3 2 1 at S = 2174 (16.4k);\n  - 14 7 4 4 1 0 0 1 at S = 3648 (20.5k).\n- Measure once on fdb77fc-equivalent search code. Re-measure on the integrated binary once picks 2 and 3 land, since their plain gains shrink the ratio somewhat.\n\nSUCCESS\n- d-first / plain CPU <= 0.55 on every sum with N >= 15k (>= 1.8x lower cost per expected magic square), and <= 0.50 geometric mean over the four N >= 15k sums.\n- Crossover N0 located (ratio = 1, expected ~5-6k) and used as the auto default.\n- If the ratio against the improved plain (picks 2+3) exceeds 0.65 at N ~ 20k, stop and re-evaluate before productionizing."
  },
  {
   "name": "Pick 2: carried bitsets above 256 labels (W = 5..8) plus per-r1 bitset width (the local universe of each first row)",
   "spec": "GOAL: remove the matrix-path penalty (x2.7 time, x2.8 nodes measured at N = 23k) above 256 labels, and the x1.38-per-word width cost, by giving each r1 block only the words it needs. Grafts idea 2 and profile item (c). It applies unchanged inside the V_d searches of pick 1, which also run on the matrix path at S = 2650 (279 labels).\n\nBRANCH: cx/wide from fdb77fc, worktree /home/user/magic-squares/.claude/worktrees/cx-wide.\n\nIMPLEMENTATION\n1. Carried path up to 512 labels.\n   - CARRY_MAX_W 8; placedw[2][SQ_MAX_N][CARRY_MAX_W]; raise COUNT_BYTES_MAX_W to 8 (check the cnt[] allocation of CNT_WORDS*W words per depth and axis).\n   - Instantiate arrange_core.h for W = 5, 6, 7; W = 8 exists.\n   - Width selection in search_vectors: Lw <= 320 -> 5, <= 384 -> 6, <= 448 -> 7, <= 512 -> 8.\n   - Audit every W-specific piece. The rotation packing of the exactly-once test in FILTER_CARRY/FILTER_COUNT (rotl until disjoint, 6 bits total) is generic. Also check COUNT_CARRY's bit transpose, carried CROSS, KEEP_CARRY/SUPPORT, PAD_LIST, and the label reconstruction from placedw (arrange.c ~l.441).\n   - The matrix path stays for > 512 labels and for CARRY_MAX_W=0 builds.\n2. Per-r1 width.\n   - Vectors are sorted by largest label, descending. SEARCH_ROOT calls TRY_CHILD(s, 0, ROW, bits(r1), r1, r1+1), so r1's subproblem uses only vectors with index > r1, whose labels are all <= x(r1).\n   - Partition the root range into maximal runs [i0, i1) with equal w(r1) = max(2, ceil((x(r1)+1)/64)).\n   - For each run: repack s->bits for the suffix [i0, N) at width w (O((N - i0) w)); fill vw[0] from that suffix; run the width-w root loop over r1 in [i0, i1). Allocate the carried arrays once, at the max width, and reuse them.\n   - Apply the r1-sampling selection to the GLOBAL r1 index, so the same r1 set is sampled as before.\n   - Carry the adaptive cross state (xs_tick / xs_n / cross_after) across runs, so its decisions track the single-width run.\n3. --min-words still forces a minimum width (testing).\n\nCORRECTNESS\n- With per-r1 width off: nodes and hashes identical on bench quick/full/prod.\n- With it on: squares and hashes identical; node differences only from cross adaptivity (report them).\n- Forced widths: quick.txt and full.txt with --min-words 5, 6, 7, 8 must give node counts AND hashes identical to the W = 4 carried run (the carried search is width-agnostic).\n- A CARRY_MAX_W=0 build gives identical squares.\n- fuzz_arrange: add instances with 257-500 labels plus a forced-width variant; >= 150 instances, 0 fails; ctest -R fast_ passes.\n\nMEASUREMENT (paired: identical r1 samples, old fdb77fc binary vs new, alternating, min of 2; <= ~45 CPU-min)\n- Above 256 labels:\n  - 9 6 4 3 1 1 1 1 at S = 2700 (N 20.9k, L = 259);\n  - 13 7 4 3 1 1 at S = 2650 (N 31.7k, L = 279). Sample r1 from the cheaper strata, plus 2-3 heavy r1 from r1 < 5,164, choosing ones that finish in < 60 s on the old binary. Compare per-r1 times on identical r1.\n  - Optionally a 13 7 4 3 1 1 sum at S ~ 2500 (N ~ 26-28k, L just over 256).\n- Carried sums: 13 7 4 3 1 1 at S = 2200 (predicted 1.28x) and 2400 (1.04x); 12 6 3 2 1 1 at S = 1200 (1.37x) and 1080 (~1.0x).\n- Inside V_d: d-sampled d-first at S = 2650, old vs new.\n\nSUCCESS\n- >= 2.0x lower CPU on the > 256-label sums at N 20-32k, plain and V_d alike (expected 2-2.7x).\n- >= 1.15x on S = 2200 and b1200.\n- No sum slower by > 2% (paired), and prod.txt within +-2%."
  },
  {
   "name": "Pick 3: count-only pre-test (PRETEST, exact) in production, depth-tuned per mode (plain depth >= 5, V_d depth 3-4), plus an optional gated extra cascade round",
   "spec": "GOAL: an exact constant-factor cut of the fastest-growing plain term, the dead subtrees under the (2,2) survivors. The term is 31% -> 48% of the time from N = 5.9k to 23k, and its share grows with N.\n\nBRANCH: cx/pretest from fdb77fc. Port FILTER_COUNT and the PRETEST block of TRY_CHILD from cx/prune 3dc9b05, without the PRUNE_DIAG / SCALE_PROF diagnostics.\n\nIMPLEMENTATION\n1. PRETEST on by default on the carried path, with a per-search threshold opts.pretest_min (vectors placed; plain default 5).\n   - A child at or above the threshold runs FILTER_COUNT on both axes: counts and masked-OR unions, no compress or store.\n   - It then applies the same count and forward checks that follow the real filters, and returns early if they fail.\n   - With pick 2 the matrix path only serves > 512 labels: no matrix-path FILTER_COUNT.\n2. Tune for V_d (pick 1) before choosing its value.\n   - In V_d, 87% of the (2,2) children die at exactly the checks PRETEST replicates (at 15k: fc o 35%, count b 21%, fc b 31%).\n   - The (1,1) -> (1,2)/(2,1) creation is 54% of d-first time.\n   - First record the death stages of the V_d (1,2)/(2,1) children with CHILD_PROF/SCALE_PROF (one d-sampled run, ~1 CPU-min). Then time pretest_min in {3, 4, 5, off} on d-sampled a2200 and a2400, and use the best value in --diag-first.\n3. Optional: SURV_FIX with one extra cross + support round at the (2,2) survivors (cx/prune), plain only.\n   - Gate it by sampled (2,2) survival over the last 256 (2,2) children (> ~20%), in the style of the cross_after switch. A list-length gate was worse than always-on at 15k.\n   - Keep the cascade cross-first: support first was x1.04-1.13 slower.\n   - Keep it only if it adds >= 3% on top of PRETEST on >= 2 of the N >= 15k sums; otherwise drop it.\n\nCORRECTNESS\n- PRETEST: node counts AND squares identical on bench quick/full/prod, on every paired sampled run (same r1, same d), on fuzz_arrange (>= 150 instances) and on pick 1's d-first fuzz; ctest -R fast_ passes.\n- SURV_FIX: squares identical; node counts lower.\n\nMEASUREMENT (paired, alternating, min of 2, identical r1/d samples; <= ~30 CPU-min)\n- plain: 13 7 4 3 1 1 at S = 2200 / 2400; 11 6 4 3 2 1 at S = 2174; 14 7 4 4 1 0 0 1 at S = 3648; 12 6 3 2 1 1 at S = 1080.\n- d-first: a2200 and a2400.\n- bench full/prod.\n\nSUCCESS\n- plain: >= 1.12x at a2400 (measured x0.845 there), >= 1.08x at a2200, geometric mean >= 1.08x over the four N >= 15k sums.\n- No regression > 1% on prod/full or on any sum.\n- V_d: adopt a depth-3/4 threshold only if it gives >= 1.08x on both d-first sums; otherwise leave V_d at the default with no regression."
  }
 ],
 "rationale": "I checked the claims against the repo before scoring:\n- Both d-first studies: cx/magic runs/summary.tsv and the cx/decomp numbers.\n- The carried-path limit in arrange.c: CARRY_MAX_W 4 is enforced only via placedw[..][4]. W = 5-7 are not instantiated, so 257-512 labels go to W = 8 on the matrix path.\n- The PRETEST code in cx/prune: exact, with the threshold set by the environment.\n- The SEARCH_ROOT suffix structure: TRY_CHILD(..., r1, r1+1) restricts each r1 to later vectors, all with labels <= x(r1). This is what makes per-r1 width sound.\n- The existence frontier: E = 1 cells sit at N ~ 20-25k, right at the 256-label wall. Examples: 9 6 4 3 1 1 1 1 at 20.9k with L = 259; 12 9 6 2 1 1 at 17k with L = 256; 13 7 4 3 1 1 with L > 256 from N ~ 25k.\n\nSCORES. Each line gives: gain at N 10-30k; confidence; prototype cost; risk to exactness; score out of 10.\n\n- Ideas 1 and 7, d-first (merged; pick 1):\n  - Gain: 2.1-2.5x per expected magic square at 20k; 3.1x at 31.7k, where both runs were on the matrix path.\n  - Confidence: high. Two independent agents agree: 0.47 vs 0.47 at 15k, 0.43 vs 0.40 at 23k. The mechanism is measured: V_d lists stay short, and no (2,2) child survives in V_d. The plain denominators carry 15-34% SE.\n  - Cost: ~1 day.\n  - Exactness: exact for magic squares, each found twice; the semi-magic output is lost, hence the calibration stream.\n  - Score: 9.\n  - It is the only idea with a measured exponent change: the ratio falls like N^-0.5 to N^-0.67. It is not output-sensitive: the cost per expected magic square still grows ~N^1.1-1.3.\n- Idea 2 plus profile item (c), wide carried bitsets and per-r1 width (pick 2):\n  - Gain: 1.6-2.7x on sums above 256 labels, which covers a large part of the frontier at N 20-32k. 1.0-1.37x (mean ~1.12) on carried sums. It multiplies with d-first, whose V_d searches hit the same wall.\n  - Confidence: medium-high. The x2.69 forced-matrix penalty and the x1.38 per word were measured; the per-r1 time shares come from the profile's logs.\n  - Cost: ~1 day.\n  - Exactness risk: none. It is the same search; forcing width must reproduce node counts exactly.\n  - Score: 8.\n- Idea 4, PRETEST (pick 3):\n  - Gain: x0.83-0.86 at 23k and x0.89-0.93 at 15k on high-survival P; about 1.0 on low-survival P.\n  - Confidence: high (measured in place, identical nodes). Cost: done; half a day to productionize. Exactness risk: none.\n  - Score: 6.5.\n  - It is picked because it is exact and cheap, and it serves the plain path (N < N0, the calibration stream, the honest d-first baseline). It may also help V_d at depth 3-4: 87% of V_d's (2,2) children die at exactly the checks it replicates. That is untested, so pick 3 tunes it there first.\n- Idea 5, SURV_FIX:\n  - Gain: 5-8% alone, 0-4% on top of PRETEST, near the noise. Needs a gate. Score: 3.5.\n  - Grafted into pick 3 as optional, with a keep-or-drop rule.\n- Idea 3, r1-level cached exactly-once relation for child creation:\n  - Gain: estimated 1.10-1.13x, plain only. V_d lists are short, with little reuse across siblings, and LABEL_MASKS history says gains measured in place are smaller than in microbenchmarks there.\n  - Confidence: medium-low (microbenchmark only). Cost: 2-3 days. Exactness: exact if order is preserved.\n  - Score: 4.5. Runner-up for a later round: replay the profile's (2,2) children through the gather path in place before building it.\n- Idea 6, cell/pair/diagonal/matching support: measured negative; do not build. Score: 1.\n  - Its value is diagnostic. No local rule cuts k12 or makes the (2,2) layer output-sensitive, which favours decomposition (pick 1) over more pruning.\n- Idea 8, both diagonals first (V_{d1,d2}), and partner pruning in V_d: measured negative, at 1.3-1.6x plain and 3.2-3.4x d-first. Score: 1.\n- Idea 9, mod-3 selection: refuted (0.46-1.2x against a predicted 9-27x); no search lever. Score: 1 for search.\n  - Forward the side result to the existence/heuristic owners: the heuristic's 1.37 congruence factor is not supported (~0.7-1.0), which lowers its magic estimates by ~1.4-2x.\n\nCOMBINED EXPECTATION\n- N 15-23k (mostly carried): pick 1 gives ~2.1-2.5x. Pick 2 adds 1.0-1.3x and pick 3 1.0-1.1x inside V_d, so about 2.3-3x lower cost per expected magic square.\n- N 25-32k, above 256 labels: about 2.7-3x from d-first times ~2x from pick 2, roughly 5-6x.\n- The calibration plain stream gains ~1.1-1.3x from picks 2 and 3 on carried sums, ~2x above 256 labels.\n\nNone of these is output-sensitive. The cost per square still rises with N, but more slowly, which moves the E = 1 frontier's ~3,000-14,000 CPU-years down by roughly 2-4x.\n\nMeasurement discipline for all three picks:\n- paired runs on identical r1 or d samples, alternating, min of 2;\n- stratified r1 sampling for plain, to get under 10% SE;\n- one heavy process at a time, each pick within ~45 CPU-min."
}

## Prototype results and verification

```
{
 "impl": {
  "name": "pick-1-diagonal-first-v-: msearch --diag-first (V_d per SP vector d, posting index, top-label root)",
  "branch": "cx/dfirst (worktree /home/user/magic-squares/.claude/worktrees/cx-dfirst)",
  "commits": "c66450a \"msearch --diag-first: d-first search with a posting index and a top-label root\" (on fdb77fc). New: src/c/dfirst.c, src/c/dfirst.h. Changed: msearch.c (--diag-first, --diag-first-min-n default 3000, --d-stride/--d-offset/--d-range/--d-chunk/--d-log, --calib-r1-stride, --d-plain-root, --r1-stride/--r1-strata/--r1-log/--sample-seed, SAMPLE_* env; dsquare/dchunk/dsum/csquare/csum records), arrange.h/arrange.c/arrange_core.h (r1 sampling incl. stratified in search_opts_t with estimates and SEs; top_numbers/top_root_only; default root loop unchanged), bench.c (SAMPLE_* env), fuzz_arrange.c (--dfirst mode), CMakeLists.txt (fast_fuzz_dfirst, fast_fuzz_dfirst_matrices, fast_msearch_dfirst_849), src/c/README.md, research/ideas.md (new section \"Diagonal-first search in msearch\"). Scratch: /tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/cx/dfirst/ (meas/ with run.sh, an.py, all JSONL and per-r1/per-d logs; gate_a/, gate_b*/, gate_b.py, fz/).",
  "results": "Beyond the spec, I added a new root for the V_d searches. Every square in V_d contains all 6 numbers of d, so d's numbers get the top labels (the rarest class on top) and only the first rows through that number are roots, about 1/6 of V_d. On the same d this gives 0.51-0.57x the nodes and 0.64-0.70x the time at every N from 4.1k to 23k, with identical pairs. It is the default; --d-plain-root gives the spec's plain root.\n\nMeasurement on fdb77fc search code, thread CPU-s per sum, all estimated. Plain search: stratified r1 (quartile strides k:2k:6k:40k), relSE 2-9%. d-first: d-sampled, 120-304 d per sum, relSE 3-6%. Two replicates of each, run alternately, same binary, same session. Columns: plain / d-first / ratio / spec's plain-root ratio.\n- 13 7 4 3 1 1 S=1900 (N 4.1k): 36 / 31 / 0.86 / 1.28\n- S=1950 (5.9k): 187 / 115 / 0.61 / 0.91\n- 12 6 3 2 1 1 S=988 (6.7k): 274 / 156 / 0.57 / 0.84\n- S=2000 (7.6k): 509 / 282 / 0.56 / 0.82\n- S=2100 (11.3k): 2,978 / 1,019 / 0.34 / 0.51\n- S=2200 (15.2k): 9,395 / 3,242 / 0.35 / 0.50\n- 11 6 4 3 2 1 S=2174 (16.4k): 11,598 / 4,262 / 0.37 / 0.54\n- 14 7 4 4 1 0 0 1 S=3648 (20.5k): 18,031 / 4,825 / 0.27 / 0.41\n- 13 7 4 3 1 1 S=2400 (23.0k): 45,614 / 12,672 / 0.28 / 0.45\n\nAt N >= 15k the geometric mean is 0.31 (plain root 0.47). Every N >= 15k sum is <= 0.37, which beats the 0.55 / 0.50 success bar. That is 2.7-3.7x lower cost per expected magic square at 15-23k, and about 1.8x at 6-8k. The plain search's cost per semi-magic square at 15-23k is 4.2-33 CPU-s; d-first spends the equivalent of 1.4-8.8.\n\nAt large x d-first still wins at large N: 0.46 at x=0.53 (N 11.7k) and 0.37 at x=0.74 (N 15k).\n\nCrossover: the ratio falls like N^-0.65 and reaches 1 at N ~ 3,000 for the new root; the plain root crosses at ~5,500, as expected. The auto default is 3000.\n\nNot kept:\n- Forward checking on d's numbers: 1% fewer nodes, time within noise, so the core was left unchanged.\n- Reducing each V_d before its search: removes 0-2% of V_d and is slower.\n\nThe V_d posting index plus search setup are 0.2-1.2% of d-first time. Plain-path regression check on bench/prod against fdb77fc's bench: identical nodes, 13.70-13.87 s vs 13.86-14.04 s.\n\nTotal CPU about 32 min.",
  "correctness_ok": true,
  "recommended": true,
  "risks": "- d-first does not enumerate semi-magic squares. Its sums write dsquare/dchunk/dsum records (and csquare/csum for the calibration stream), and scheduler.py does not read these yet: it will neither mark those sums covered nor fit its models to them. Integration is still to do.\n- With default --diag-first the cutover is at N >= 3000 vectors after reduction, so the semi-magic squares above it come only from the --calib-r1-stride sample. Raise --diag-first-min-n if more full semi-magic output is wanted near 3-5k, where the gain is 0-15%.\n- Every measured sum has <= 256 labels (carried path). The matrix path (> 256 labels, most sums with N >= 25k) was not measured.\n- The spec's re-measure on the integrated binary with picks 2 and 3 was not possible here. Those picks speed up both arms, but V_d lists are about 6x shorter, so the ratio may move; the 0.65-at-20k stop rule still needs checking then.\n- The machine was shared during the measurements (load 4-8 on 4 cores). Both arms were affected alike and alternated, but absolute times carry contention.\n- d-first relSE was 5-6% on four sums, a little over the 5% target. With the new root the per-d CV is 0.54-0.61, so about 150 d are needed for 5%.\n- At small N with large x (e.g. 10 6 3 2 1 / 1360: N 3.1k, x ~1.1) d-first is 2.6x slower, but such sums essentially never have squares.\n- The checkpoint-and-resume path (--d-range from the last dchunk) was tested by split-and-merge, not by killing a real run."
 },
 "verdict": {
  "name": "pick-1-diagonal-first-v-: msearch --diag-first (V_d per SP vector d, posting index, top-label root)",
  "correctness_ok": true,
  "gain_confirmed": true,
  "recommend_merge": true,
  "problems": "VERIFIED (copies of c66450a and fdb77fc built in scratchpad/cx/review-dfirst/{new,base}; the branch itself was not touched):\n- Exactness of the top-label root. Every square of V_d contains all 6 numbers of d. The rarest of them gets label L-1, and vectors are sorted by descending label tuple, so the 2 vectors of any square that hold L-1 sort before every other vector. The square's lowest-index vector is therefore always below root_limit, so restricting the roots this way is exact. V_d keeps exactly the vectors with |v&d|=1 (hit counter, bitmap). The diagonal set is always the whole unreduced list (dfirst_new(red,..,all,start_i,raw)), so set_count and partner stay correct under --d-range and --d-stride.\n- Memory and portability. inv[128] is enough for SQ_MAX_N=8 (at most 105 involutions). Every search_opts_t is created through search_opts_default. Clean under -Wall -Wextra (gcc and clang). ASan+UBSan: 200 dfirst fuzz seeds plus msearch d-first, d-range, d-stride, calib and strata runs all clean.\n- Gates:\n  - ctest -R fast_: 33/33 pass (fast_scheduler needed the bin/ symlinks in my copy).\n  - fuzz_arrange --dfirst, fresh seeds: 1000 native, 1000 matrices, 1000 cascadelake, 500 cascadelake matrices and 400 portable, all 0 fails. Plain fuzz: 1000 native + 1000 cascadelake + 500 cascadelake-matrices, 0 fails.\n  - bench quick/full/prod: hashes and nodes identical to fdb77fc. prod 14.08 s vs 14.17 s.\n- Gate (a): --d-plain-root reproduces dsearch's nodes on all 32 sums I ran (622,221 at 391 and 1,723,937 at 517, among them).\n- Gate (b): on my 32 sums, both roots match a Python brute-force oracle (traversals of the plain squares with sum S and product P).\n- Stronger real-data test (scripts/dtest.c and real_oracle.py). Diagonal set = random traversals of the plain squares plus planted diagonal/anti-diagonal pairs, run against the full plain outputs at N 4.1k, 6.0k, 6.7k and 7.9k. All 1,266 pairs and 296 partner flags are exact: no misses, no duplicates.\n- Splits: d-range splits, d-stride 3 splits and dchunk records all add up to the full run.\n\nGAIN, re-measured (same binary, plain = stratified r1 with my own seeds, plus a third replicate with a different 8-strata design; d-first = d-sampled with 80 d, alternating). Ratio is d-first / plain CPU:\n- 13 7 4 3 1 1 S=2200 (N 15.2k): 0.307 +- 0.022. Claimed 0.35.\n- 14 7 4 4 1 0 0 1 S=3648 (N 20.5k): 0.346 +- 0.025. Claimed 0.27.\n- 13 7 4 3 1 1 S=2400 (N 23k): 0.255 +- 0.021. Claimed 0.28.\n- NEW, matrix path, which the prototype did not measure: 9 6 4 3 1 1 1 1 S=2700 (N 20.9k, 259 labels; 38 of 70 V_d also have > 256 labels): 0.257 +- 0.038.\n- Geometric mean over the four is 0.29. Every N >= 15k sum is <= 0.35, well inside the <= 0.55 / 0.50 bar.\n- Full, unsampled runs at N 4.1k (13 7 4 2 0 1 / 1166): plain 34.2 s, d-first 27.4 s, ratio 0.80, which fits the N0 ~ 3000 crossover. The sampled estimates were within 3-4% of these full runs.\n- V_d build plus setup is <= 0.3% of the time. The top 10% of d hold 20% of the time.\n\nPROBLEMS / CAVEATS:\n1. Per-sum ratios move with machine load. My 3648 d-first estimate is 19% above the prototype's (2.6 sigma); 2200 moved the other way. So individual ratios carry +-0.05-0.08 of contention noise, and the qualitative result (about 3x at 15-23k) is the robust part.\n2. Gate (b) on real sums is nearly vacuous. SP pairs are extremely rare: 3 pairs in my 32 sums, and 0 in the full runs at 988, 1166, 1302 and 3231. Exactness rests on the fuzz and my planted-traversal real-data test, both of which pass.\n3. Integration hazard, beyond the noted \"scheduler ignores dsum\": msearch still writes a \"done\" record covering d-first sums. Those sums have no \"sum\" or \"square\" records, and the scheduler's Results.covered treats a done-covered sum with no sum record as searched and empty. If --diag-first were ever passed in scheduler units, the yield and time models would silently mis-fit. Also, with --d-range units, a sum below --diag-first-min-n is fully plain-searched by every unit, duplicating work and records. Fix both before production.\n4. The default --diag-first-min-n 3000 sits at the crossover. Gain at 3-5k is 0-20%, and those sums lose their semi-magic output for the models. Consider 5000-6000 (the spec's value) or a mandatory --calib-r1-stride.\n5. Minor:\n   - Strata with fewer than 2 samples (Q4 almost always) add nothing to the SE, so SE is slightly understated. Negligible: Q4 holds < 1% of the time.\n   - bench.c never fcloses the SAMPLE_LOG file.\n   - gen_inv's inv[128] would overflow silently if SQ_MAX_N were ever raised to >= 9; add a static assert.\n6. Not done (not possible here): re-measuring with picks 2 and 3 integrated.\n\nCPU used: about 35 min."
 }
}
```

```
{
 "impl": {
  "name": "pick-2-carried-bitsets-a",
  "branch": "cx/wide (worktree /home/user/magic-squares/.claude/worktrees/cx-wide)",
  "commits": "272cbc2 arrange: carried bitsets up to 512 labels, each r1 at the words it needs (on fdb77fc)",
  "results": "Paired vs fdb77fc (identical r1 sets via the same -DR1_SAMPLE patch on both; per-r1 thread CPU, min of 2 alternating runs). >256 labels (base on the matrix path): 13 7 4 3 1 1 S=2650 (N 31.7k, L 279; 4 heavy r1 <5164 + 7 at stride 4096) 64.66 -> 22.80 s = 2.84x (5-word r1 2.25x, 4-word r1 3.49x, per r1 2.2-4.2x), nodes 3.33x fewer; S=2500 (N 26.6k, L 260, 6 r1) 3.13x; 9 6 4 3 1 1 1 1 S=2700 (N 20.9k, L 259, stride 800) 53.28 -> 18.58 s = 2.87x, nodes 5.1x fewer; V_d searches at S=2650 (dsample, 16 d at stride 2000, |V_d|~4.5k, 274-279 labels) 54.78 -> 25.18 s = 2.18x (1.93-2.44 per d). Carried sums (per-r1 width only, identical nodes): 13 7 4 3 1 1 S=2200 (N 15.2k) 1.225x; 12 6 3 2 1 1 S=1200 (N 11.7k) 1.31x; S=2400 (N 23k) 1.001x; b1080 (N 8.9k) 1.005x; b988 (N 6.7k) 1.09x; 14 7 4 4 1 0 0 1 S=3231 (N 7.9k) 1.06x; 12 9 6 2 1 1 S=3500 (N 17.1k, L 256) 1.02x; prod.txt (N 1.5-3k) base/new 1.004 total (within +-3% per instance). Same-width r1 within noise (one r1 per process, 5 interleaved rounds: 0.996). Cost per square (whole-sum estimates; squares identical so ratio = time ratio): S=2200 4.3 -> 3.5 s/sq, b1200 5.7 -> 4.4, S=2700 ~27 -> 9.3, S=2400 13.4 -> 13.4, S=2650 ~70 -> ~25 (square count 4600+-3400); local slope along 13 7 4 3 1 1 above 2e4 drops from 5.2 to roughly 3-3.4 (rough, 6-15 r1 samples). Width cost at identical nodes (full.txt, forced 2..8 words): 1.00,1.45,1.98,2.43,2.93,3.43,3.98 (time ~ proportional to words). Cross unions chunked by 4 words: 1.6% (min)/7% (median) less on 5-word r1, unchanged at <=4 words. Nodes: quick 1,770,779 / full 14,958,507 / prod 50,375,738 identical with per-r1 width on/off and --min-words 4-8; matrices build node-identical to fdb77fc. Squares and hashes identical on every sampled r1 with squares (S=2700, 2500, 2200, 2400, 1200, 3231).",
  "correctness_ok": true,
  "recommended": true,
  "risks": "Measurements on a shared 4-core machine at load 3-10 (run-to-run spread 5-12%); S=2650 heavy r1 were hand-picked to finish <60 s on base (biased toward cheaper heavy r1; per-r1 ratios 2.2-4.2x). Whole-sum base estimate at S=2650 is scaled from the existence study's older code (x0.9). Memory: carried lists are allocated at full N per depth (53 MB at N=31.7k W=5; ~270 MB at N=100k W=8), still below the matrices. Three more instantiations (compile time, code size). CTest: 8-word gather and plain-C cross tests moved to bench_matrices since the native build now carries 8 words; new tests fast_bench_quick_wide5/7, wide6_fixed, matrices_wide8, fast_fuzz_arrange_wide, fast_fuzz_arrange_matrices_wide. >512 labels still on the matrix path. This removes an implementation penalty, not the algorithmic N^1.5-per-(1,2)-node growth: cost per square still rises with N. Narrowing the width inside an r1 and relabelling a node's labels into fewer words were judged not worth it (recorded in research/ideas.md); -DR1_SAMPLE sampling hooks are compiled out by default."
 },
 "verdict": {
  "name": "pick-2-carried-bitsets-a",
  "correctness_ok": true,
  "gain_confirmed": true,
  "problems": "No blocking problems. The branch was left untouched: I worked on git-archive copies of fdb77fc and 272cbc2 in scratchpad/cx/review-wide. I used about 40 CPU-minutes, one heavy process at a time.\n\nCODE REVIEW (exactness, memory safety, portability)\n- The run partition in search_root covers [0, N) exactly, in ascending r1 order. Each r1 is searched exactly once, with the same min_v = r1+1 on both axes, so symmetry breaking is unchanged and nothing can be duplicated.\n- The narrower width is valid. Vectors are sorted by their descending label tuples and bits/lab are built after the sort, so for every v > r1, lab[v][0] <= lab[r1][0] < 64*w(r1).\n- Entries of vw[0] below i0 are never read. PAD_LIST at count and the cap of count+64 still hold.\n- Words above w in cells, uni, placedw and cnt can hold stale data, but they are never read at width w. report_bits reads only W words.\n- The adaptive cross state (xs_*, cross_after) lives in s and is not reset between runs, so nodes match a single-width run exactly.\n- FILTER_CARRY's rotation packing is generic in W: at most a*b <= 16 bad rotations when a+b <= 8.\n- CROSS_CHUNK covers words 0..3 and 4..W-1, with 64-byte-aligned stores into Ub[CROSS_MAXY*W].\n- Byte counters are 8 words per word of labels, inside the CNT_WORDS*W_ allocation.\n- The width override only applies when (Lw+63)/64 <= CARRY_MAX_W, and every width search_root_w can be asked for is instantiated. Non-AVX512BW builds stay on the matrix path and never call r1_words.\n- Spec deviations, both equivalent or better: COUNT_BYTES is forced whenever CARRY is set instead of raising COUNT_BYTES_MAX_W, so the matrix W=8 path keeps its bit-sliced counters and CARRY_MAX_W=0 builds stay node-identical to fdb77fc. And bits are read with stride bw instead of being repacked.\n\nGATES (all run by me)\n- ctest -R fast_: 36 of 36 pass. fast_scheduler first failed only because new/bin/enumerate did not exist yet (an environment step, not a code issue); it passes once the binaries are linked into bin/.\n- fuzz_arrange on fresh seeds, all variants, 0 fails everywhere:\n  - native: 1000 mixed (71000+) and 1000 mode 7 (72000+);\n  - -march=cascadelake: 1000 + 1000;\n  - -DARRANGE_DEBUG (padding asserts): 300 + 300;\n  - ASan+UBSan: 200 mode 7, clean;\n  - clang-18: 300 + 300.\n- Mode 7 really exercises every width: with an R1_SAMPLE fuzz build, 300 seeds found squares under r1 at every width 2-8 (4793, 2086, 1097, 533, 259, 73 and 9 squares).\n- bench quick/full: nodes and hashes identical per instance for base, new, --no-r1-width, --min-words 3/5/6/7/8 and --no-r1-width --min-words 6.\n- bench prod: identical nodes, all ok.\n- cascadelake and clang builds: identical nodes. Matrix builds (CARRY_MAX_W=0): base and new node-identical on quick/full, plain, with --min-words 8 and with --no-cross.\n- At S=2700, the same 6-r1 list gives identical per-r1 nodes, squares and hash for native, --no-r1-width (W5), --min-words 8 and cascadelake.\n\nMEASUREMENTS (my own r1 samples, chosen to differ from the prototype's; base fdb77fc vs new, both built with the same R1_SAMPLE patch, alternating, min of 2 per r1)\n| Sum | Sample | Base -> new (CPU-s) | Gain | Notes |\n|---|---|---|---|---|\n| 13 7 4 3 1 1, S=2650 (N 31.7k, L 279) | 13 r1: the 4 heavy ones 964, 2164, 2764, 4564 (964 is the heaviest known, 64 s on base) plus 9 light ones from 6000 to 30000 | 207.4 -> 70.1 | 2.96x | 5-word r1 2.54x, 4-word 3.53x; nodes 2.6x fewer; 3 squares, same hash |\n| same sum, whole-sum estimate | strata weighted 1291 / 2953 | 374k -> 119k | 3.16x | base estimate matches the existence study's 375k CPU-s |\n| 9 6 4 3 1 1 1 1, S=2700 (N 20.9k) | stride 800, offset 684 (disjoint from the prototype's) | 75.7 -> 21.6 | 3.50x | 2 squares, same hash |\n| 13 7 4 3 1 1, S=2500 (N 26.6k) | 5 r1 | 55.6 -> 13.7 | 4.06x | no 5-word r1 in this sample |\n| V_d search (dsample), S=2650 | 13 d, stride 2500, offset 1234 | 48.87 -> 22.40 | 2.18x | |\n| 13 7 4 3 1 1, S=2200 | | | 1.226x | |\n| 12 6 3 2 1 1, S=1200 | | | 1.293x | |\n| 13 7 4 3 1 1, S=2400 | | | 1.076x | |\n| 12 6 3 2 1 1, S=1080 | | | 1.006x | |\n| 12 6 3 2 1 1, S=988 | two runs | | 1.032x / 1.079x | |\n| bench/prod.txt | 3 alternations | | 1.004x | |\n\n- On the carried sums (L <= 256) nodes, squares and hashes were identical.\n- All success criteria are met.\n- One apparent slowdown was noise: in the first S=988 run the 3-word r1 looked 6% slower, but a focused rerun of 12 heavy 3-word r1 (5 interleaved rounds) and the repeated pair both gave 1.00.\n- The cascadelake-built code shows the same gain on this CPU: 3.2x at S=2700, same squares.\n\nCAVEATS (minor)\n- arrange.o text grows from 163 KB to 300 KB and compile time from 5.5 to 8.7 s (three more instantiations). No i-cache penalty showed up at equal width.\n- The V_d gain (2.18x, with only 1.24x fewer nodes) is just above the 2.0x target.\n- Per-r1 width gives nothing where the heavy r1 need all the words (S=2400).\n- The matrix-path base is memory-bound and noisy under load (base run-to-run spread up to 20% at S=2650), so ratios use the min of 2 runs.\n- Real Skylake-X / Cascade Lake hardware is untested; only the -march=cascadelake build was run, on this newer Xeon.\n- Per-r1 node counts in sampled runs depend on which r1 are in the sample, because the adaptive cross state carries from one r1 to the next. That predates this change, and the comparisons above used identical r1 sets.\n- This removes an implementation penalty, not the algorithmic growth: cost per square still rises with N (roughly 81 -> 26 s per square at S=2650).",
  "recommend_merge": true
 }
}
```

```
{
 "impl": {
  "name": "pick-3-count-only-pre-te (PRETEST: count-only filters for the deep children, apply-from-masks for survivors)",
  "branch": "cx/pretest (worktree /home/user/magic-squares/.claude/worktrees/cx-pretest, from fdb77fc)",
  "commits": "24b60b7 arrange: pretest of the deep children (count-only filters, opts.pretest_min). Changes: src/c/arrange_core.h (FILTER_COUNT, FILTER_APPLY, noinline PT_APPLY, pretest block in TRY_CHILD), arrange.c/.h (opts.pretest_min default 5, 0 = off; ptmask buffers; CHILD_PROF pretest phase), bench.c and msearch.c (--pretest-min K), fuzz_arrange.c (pt_* variants plus a node-equality check against pretest off), src/c/README.md, research/ideas.md (new section \"The pretest of the deep children\", which also records the negative results). Scratch: /tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/cx/pretest/ holds pair.sh/an.py (r1-paired), dpair.sh/dan.py (d-paired), mkdf.sh (cx/dfirst + patch), runs/, df/, and survfix_src/ (the dropped SURV_FIX code).",
  "results": "Change from the cx/prune port: survivors no longer run the two filters again. FILTER_COUNT stores one byte of keep mask per 8 entries, and the lists of the survivors are compressed from those masks. This was x1.037 faster than the re-filter version at a2400 in the same session.\n\nPLAIN, final code vs fdb77fc. Paired, alternating, per-r1 min of 2 runs, same r1 sample, thread CPU of the sampled r1:\n- 13 7 4 3 1 1/1850 (N 2.5k, every r1): 2.91 -> 2.88 s, x1.009\n- 13 7 4 3 1 1/2000 (N 7.6k, stride 40): 13.29 -> 12.50, x1.063\n- 12 6 3 2 1 1/1080 (N 8.9k, stride 80): 9.10 -> 8.60, x1.058\n- a2200 (N 15.2k, stride 800): 14.37 -> 12.98, x1.107\n- 11 6 4 3 2 1/2174 (N 16.4k): 14.49 -> 12.93, x1.120\n- 14 7 4 4 1 0 0 1/3648 (N 20.5k, only 12 r1): 9.19 -> 8.96, x1.026 (the earlier re-filter build measured x1.060)\n- a2400 (N 23k): 32.14 -> 27.46, x1.170\n\nGeometric mean over the four N >= 15k sums: x1.105. Squares and nodes are identical, so cost per square falls by the same factor: 4.0 -> 3.6 CPU-s at a2200, 10.8 -> 9.6 at 2174, 15.5 -> 13.2 at 2400. This is a constant factor that grows slowly with N, not a lower exponent.\n\nBench, min of 4 alternating runs, same nodes and hashes:\n- prod: 13.85 -> 13.49 s (-2.6%)\n- full: 3.359 -> 3.333 s (-0.8%)\n- quick (--repeat 10, 5 rounds): 0.2719 -> 0.2753 s (+1.3%). Its two smallest instances are +3-6%. The same binary with the pretest off gave 0.2733.\n- --pretest-min 4 in the plain search: x0.969 at a2200, so plain stays at 5.\n\nV_D (cx/dfirst + this patch, d-sampled, per-d min of 2 runs, against the cx/dfirst binary):\n- Death stages at a2200 (60 d, CHILD_PROF): (1,2) fc o 6.4%, fc b 0.4%, support 69.1%, live 24.2%. (2,1) fc o 11.9%, fc b 1.3%, support 49.1%. (2,2) 88% at the checks the pretest repeats. Creating the (1,2)/(2,1) children is 61% of the cycles.\n- Thresholds at a2200: off 1.020, 1 0.952, 3 0.987, 4 1.053, 5 1.015.\n- Thresholds at a2400: off 0.980, 3 0.941, 4 1.062, 5 1.020.\n- 4 is best but below the 1.08 bar on both sums, so it is not adopted. V_d keeps the default 5, with no regression.\n\nSURV_FIX, gated at more than 1/5 survival over the last ~256 (2,2) children, on top of the pretest:\n- a2200 x1.008, c2174 x1.021, d3648 x1.019, a2400 x1.039 (always-on x1.050).\n- Nodes fall 23-49%, but only 1 of 4 sums gains 3% or more, so it was dropped.\n\nSUCCESS CRITERIA: a2400 x1.170 (needed 1.12) met. a2200 x1.107 (needed 1.08) met. Geometric mean x1.105 (needed 1.08) met. No sum regressed and prod/full got faster. The only slowdown measured is quick.txt +1.3%, which is not on the regression list. CPU used: about 1.3 CPU-hours.",
  "correctness_ok": true,
  "recommended": true,
  "risks": "1. A sampled run in the plain table has 12-19 r1 at N >= 15k, on a shared machine (load 2-7). Session-to-session noise is about 5%: a2400 gave x1.10 in one session and x1.19 in another with the re-filter build, and 3648 gave x1.026 vs x1.060.\n2. Small N regresses slightly: quick.txt is +1.3%, and its two 12-50 ms instances are +3-6%. A list-length gate might remove this; it was not tried.\n3. V_d was tested on a scratch build of cx/dfirst with this patch applied by hand. Three hunks needed manual fixing: sstate_t, the opts default, and the normalization in search_vectors. Integration must keep `s.pretest_min = opts->pretest_min > 0 && opts->forward_check ? opts->pretest_min : INT_MAX` in search_vectors. Without it, s.pretest_min stays 0 from the memset and the pretest runs at every depth. That is still correct, but about 4% slower in V_d. My own scratch script hit exactly this, and a few intermediate d-first comparisons were invalid until it was fixed; the reported numbers come after the fix.\n4. The pretest exists only on the carried path. cx/wide extends that path to 512 labels, and the pretest should apply there unchanged because FILTER_COUNT is generic in W, but this was not tested with W = 5..8.\n5. TRY_CHILD now has a `goto filtered` into the block after the b filter. It skips the count/fc checks the pretest already made. It is legal C because W is a constant and nothing there is a variable-length array.\n6. Whether to give V_d threshold 4 (x1.05-1.06) is left to the integrator. The bar was 1.08, so it was not adopted."
 },
 "verdict": {
  "name": "pick-3-count-only-pre-te (PRETEST: count-only filters for the deep children, apply-from-masks for survivors) -- cx/pretest 24b60b7",
  "correctness_ok": true,
  "gain_confirmed": true,
  "problems": "VERIFIED (my own copies via git archive in /tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/cx/verify-pretest/; the branch is untouched, HEAD 24b60b7, clean).\n\nEXACTNESS (diff review): FILTER_COUNT runs the same tests as FILTER_CARRY (same rotation and once/excl setup, same popcount or (c, c-1) test). Its union (a masked or of x) equals FILTER_CARRY's or of the maskz-compressed lanes. The pretest checks are the same: count o, then fc o, then pex = ~(ncell_o|uo) with sup ? pex : NULL, count b, fc b. It only runs when forward_check is on (normalized to INT_MAX otherwise), so the same children die. The cross_after sampling counter (xs_tick) only advances for children that reach cross, so its decisions are unchanged. After `goto filtered` only dead (= 0), nvalid, uni and vw[d+1] are read, never the skipped ko, kb or excl. vw[d+1] gets byte-identical contents (maskz compress, 8-lane stores at k). ptmask is consumed by PT_APPLY before the recursion, so one buffer per search is safe. Memory: ptmask is cap/8+1 bytes with cap >= count+64, and at most ceil(cnt/8) <= cap/8 groups are written. The FILTER_APPLY stores reach k+8 <= cap. All of it is freed, and free(NULL) is fine on the matrix path. Under ASan+UBSan (pt_* variants, 300 seeds): clean.\n\nGATES: ctest -R fast_: 30/30 pass. fast_scheduler first failed only because my archive copy lacked the bin/ symlinks that build.sh creates; it passes once they exist. fuzz_arrange on 1000 fresh seeds (731000-731999): 0 fails native and 0 fails with -march=cascadelake (no vpopcntq in the binary, so the non-VPOPCNTDQ branch was tested). Mutation check, 200 seeds each: fc o check dropped (338 fails), count b off by one (1148), wrong keep mask for the b list (first seed) and support excl dropped (2, by the node check). The fuzz plus the node-equality check is sensitive. cx/dfirst + this patch, rebuilt by hand (same 3 rejected arrange.c hunks + arrange.h, as reported): d-first fuzz 600 seeds with the pretest at depths 1-5, 0 fails; plain fuzz 300 seeds, 0 fails. bench prod/full/quick: identical nodes and all instances ok.\n\nGAIN (r1-sampled, alternating base/new, per-r1 min of 2 runs, FRESH offsets different from the prototype's; nodes and squares identical at every r1):\n- a2200 (stride 800, off 501): x1.100, bootstrap90 [1.05, 1.14]. Claimed 1.107.\n- a2400 (stride 1600, off 977): x1.205 [1.09, 1.34]. Second sample (off 1450): x1.158 [1.05, 1.26]. Pooled over 28 r1: x1.182. Claimed 1.170.\n- c2174 (off 450): x1.188 [1.12, 1.24]. Claimed 1.120.\n- d3648 (stride 1200, off 555): x1.068 [1.02, 1.12]. Claimed 1.026.\n- Geometric mean over the four N >= 15k sums: x1.13 (claimed 1.105).\n- Mid and small N: b1080 x1.053 (claimed 1.058); a1850 (every r1) x0.995 (claimed 1.009).\n- V_d (dfirst+patch vs dfirst, a2200, 60 d, fresh seed 7, default 5): x0.997, identical nodes. Neutral, no regression.\nAll plain success bars are met, but the margin is thin at single-sample level: each a2400 sample's bootstrap lower bound is below 1.12, a2200's is below 1.08.\n\nPROBLEMS / RISKS:\n1. PORTABILITY (real, trivial fix): the `filtered:` label is directly followed by a declaration (`const int cross = ...`). That is a C23 feature. clang 18 warns by default (-Wc23-extensions, 3 warnings), gcc -std=gnu17 -pedantic-errors gives 3 errors where base gives 0, and GCC <= 10 and clang <= 15 reject it outright. The cluster's compiler is unknown. Fix: `filtered:;` (checked: 0 pedantic errors, 0 clang warnings).\n2. Small-N cost is neutral to slightly negative, not the claimed gains. Sum of per-instance mins: quick +0.6% and a1850 -0.5%. full showed +2.0% over 9 rounds, but in isolation its biggest instance (14.7.5.3, the one at +3%) was even (1.892 vs 1.918 min), and prod -0.3%. Everything is within the machine's ~2-5% noise. The README's \"0-3% on prod/full\" is optimistic; it should say neutral within noise.\n3. Diagnostic nit: under -DPROFILE the pretest's early return skips prof_created/prof_kids/prof_pruned, so the deep-depth created/pruned counts would be undercounted. prof_print has no callers, so this has no practical effect.\n4. Integration with cx/dfirst needs hand-merging of the arrange.c hunks (limits.h, default, s.pretest_min normalization) and arrange.h. The dfirst fuzz loop also needs the `if (V->pretest)` line, otherwise the pt_* variants are not exercised in --dfirst mode.\n5. Not tested with cx/wide W=5..8 (as stated by the prototype).\n6. The gain is a constant factor (~1.1-1.2x at N 15-23k, growing slowly with N), not a lower exponent.\n\nCPU used by this verification: ~30 CPU-min (about 10 of it a hung mutant).",
  "recommend_merge": true
 }
}
```
