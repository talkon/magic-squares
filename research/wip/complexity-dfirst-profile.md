# Profile of the diagonal-first (V_d) search at large N

Part 1 of round 2 of the complexity study (c2/profile), 2026-10-09.

* Code: be625d8 (engine 3) plus instrumentation that is compiled out by default.
  * Branch `c2/profile`, commit cd9faac, worktree `.claude/worktrees/c2-profile`.
  * The default build is unchanged: bench quick has 1,770,779 nodes, and `fuzz_arrange --dfirst` passes.
* Instrumentation builds: `cmake -DMAGIC_DEFS=...`, in `cmake-build-dprof` and `cmake-build-cprof`.
  * `-DDPROF` gives per-layer counts of the children made and searched, and the TSC cycles spent
    creating them (without their subtrees) and selecting them. It also records list lengths, the
    vectors placed in searched nodes, and how the support filter kills. One line per d goes to
    `DPROF_LOG`.
  * `-DDPROF -DCHILD_PROF` adds the round-1 phase timers and death stages.
* Checks: both instrumented builds visited the same nodes as the release build on every sampled
  d of all 12 sums. `fuzz_arrange --dfirst` (60 seeds) and the plain fuzz pass with DPROF.
* Overhead: DPROF adds ~9% CPU (estimated CPU per sum 1.01-1.09x the ideas.md values).
  CHILD_PROF adds ~80%, and its phase shares over-weight the short phases.
* Machine: 4-core Xeon at 2.1 GHz TSC, with AVX-512 (VPOPCNTDQ, VBMI, GFNI). Shared, load ~2.
  Everything ran under `nice -n 10`, one heavy process at a time.
* Sampling: d-sampled with `--d-stride k --d-offset 13 % k`, 106-152 d per sum, the same d in
  every build. Per-d CV is 0.55-0.66, so the means have a relSE of ~5-6%.
* CPU used: ~45 CPU-min in total, ~25 of it on the two ladder passes.

## Bottom line

**1. A V_d search is a three-layer tree, and nothing below it survives.**
* Roots: the r1 through T, the rarest of d's numbers in V_d. R ~ 0.075 |V_d| per d.
* Then the (1,1) nodes, and then the (1,2)/(2,1) children. About 30-40% of those survive.
* Each survivor makes 2.4-3.2 (2,2) children, and **none survived on any sampled d of any sum**.
  * That is 0 of ~10^7 (2,2) children made across the samples.
  * The (1,3)/(3,1) nodes searched are ~0.007 per d.
* About 10^-5 of V_d's vectors are ever placed at depth >= 4.
* Every searched node of every V_d is dead: the search proves emptiness at depth 3.
  * The expected output is ~2-4x10^-5 (square, d) pairs per d.

**2. Where the cycles go** (DPROF, exclusive cycles; 12 sums, N = 5.9k-31.7k):
* **Creating the (1,2)/(2,1) children: 64-75% of all cycles**, plus 1-3% for their selection.
* Creating the (2,2) children: 5-16%.
* Creating the (1,1) children: 8-13%.
* The roots: 0.5-1.4%.
* Setup: 0.1-1.6%, falling with N. This covers V_d from the posting index, the search's
  relabelling and bitsets, and the depth-0 list fills.
* Unattributed: 3-8%.

CHILD_PROF splits the (1,2)/(2,1) creation into phases:
* the exactly-once filter o: 14-20%;
* the disjointness filter b: 14-26%;
* the support filter: 20-23%;
* count and selection: 2-6%.

**3. Time per sum is N^3.44 on ladder A** (13 7 4 3 1 1, 5.9-31.7k). It factors as follows:

| factor | N exponent | value at 5.9k / 15.2k / 31.7k |
|---|---:|---|
| number of d | 1 | N |
| roots per d (R) | 0.79 | 92 / 203 / 355 |
| (1,1) children per root (k1) | 0.49 | 17.6 / 30 / 40 |
| cycles per searched (1,1) node | 1.17 | 28k / 85k / 189k |
| of which: (1,2)/(2,1) children per (1,1) node (k11) | 0.39 | 17.5 / 26 / 34 |
| of which: cycles per (1,2)/(2,1) child | 0.83 | 1.0k / 2.1k / 4.0k |
| = cycles per d | 2.44 | 0.021 s / 0.25 s / 1.28 s |

* All 12 sums together give cost per d ~N^2.35.
* The **output**, from the expected (square, SP traversal) pairs per sum (amodel E_semi x SQ12 x
  720 p_SP), grows only N^0.71: 0.16 / 0.64 / 0.50 at 5.9k / 15.2k / 31.7k.
  * It is N^1.08 up to 23k and falls above.
  * **Pairs per d are flat to falling (N^-0.3).**
  * So the cost per pair grows ~N^2.7, which matches ideas.md's N^2.8-3.2 per expected magic
    square.

**4. The terms that grow faster than the output.**
* **(a) The number of refuted (d, r1, c1) triples, i.e. searched (1,1) nodes: N^2.28 per sum
  against output N^0.7.**
  * Their productivity falls like N^-1.6.
  * This is redundancy across d. Measured on the vector lists, a random pair (r1, c1) meeting
    once lies in the top-root tree of 26 / 39 / 52 different d at N = 7.6k / 15.2k / 23k
    (~N^0.6).
  * A random (r1; c1, c2) lies in 5.4 / 7.0 / 8.4 trees.
  * d-first searches 67x / 123x / 146x the plain search's (1,1) nodes at 5.9k / 15.2k / 23k.
  * It searches 5-7.7x the plain search's (1,2)/(2,1) nodes.
* **(b) The cost of refuting one (1,1) node: N^1.17.** Its dominant part is the creation of k11
  third vectors (N^0.39) at a cost per child of N^0.83:
  * **Each child's two filters scan the parent's whole (1,1) lists** (435 -> 1,164 entries,
    N^0.59) and keep 24-28% of them.
    * The child's own lists grow N^0.35-0.6: 58+63 -> 104+165.
    * The 26-34 siblings each rescan the same parent lists.
    * Filters o + b are 28-43% of all cycles (CHILD_PROF).
  * **The bitset width** costs x1.25 per 64-bit word inside V_d, at identical nodes (paired,
    S = 2000: 3 -> 4 words, S = 2200: 4 -> 5).
    * V_d keeps 99% of the sum's labels, so it steps to 4 words at ~13k and to 5 words at ~25k.
    * That is N^0.27 of the ladder's exponent. Without the steps, cost per d would be ~N^2.17.
  * **Most children die in the support filter** (47-85%, after filters o and b). Each call runs
    4.0-5.5 passes and costs ~460-1,260 cycles.
    * 79-82% of these kills are a dead cell, the rest a count.
    * In 63% of the dead-cell kills, a cell of the new vector itself (c2 or r2) has lost all
      support.
    * d's numbers are among the dead cells in only 4-7% (1-2% alone).
    * So the class structure of V_d (one row and one col per number of d) is not what kills.
* **(c) The (2,2) children of the (1,2)/(2,1) survivors.**
  * There are k12 ~ 2.4-3.2 per survivor (N^0.18). All die, 88-94% at the forward check or
    count right after their two filters, i.e. at the pretest's checks (pretest_min is 5; these
    have 4 vectors).
  * They take 5-16% of the cycles, ~N^1.27 per (1,1) node. They are smaller than (a) and (b).
* Not super-linear: setup, the root layer and selection (each <= 1.6%, falling).

**5. Per d: the cost is extremely uneven in |V_d|, but so is the output, so ordering d is no
lever.**
* Within and across sums, ln(cycles per d) = 5.65 ln|V_d| + 1.43 ln R - 2.68 ln N + c (R^2 0.97,
  1,701 d).
* The top 10% of d hold 20-22% of the time.
* The transversals of 145 plain-search squares at S = 1950 sit on the expensive d:
  * mean |V_t| is 1,464 against 1,378 over the d, and mean R is 134 against 97;
  * the cheapest half of d by predicted cost holds 23% of the cost and 18% (S-traversals) to 25%
    (all transversals) of that output proxy.
* So yield per cost is about homogeneous across d. The proxy assumes P(t is SP) does not depend on
  |V_t|.

**6. V_d cannot be shrunk statically.**
* 75-91% of V_d's vectors are placed in some searched node, and 66-89% in a depth-3 node.
* The expected share of V_d in a square is ~10^-7.
* So a static reduction of V_d has nothing to remove (ideas.md's reduce on V_d: 0.05%). The waste
  is in how often each vector is combined, not in which vectors exist.

## 1. Per sum

d-sampled with the DPROF build. "est CPU" is mean per-d CPU x N, with ~9% instrumentation
overhead. "E pairs" is the model's expected pairs per sum.

| sum | N | labels | d | mean V_d | V_d / N | labels of V_d | roots R / d | R / V_d | nodes / d | CPU / d (s) | est CPU (s) | ideas.md d-first CPU | E pairs | top 10% of d | CV |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 7 4 3 1 1 / 1950 | 5,896 | 153 | 148 | 1,385 | 0.235 | 152 | 95 | 0.069 | 5.06e4 | 0.0216 | 127 | 116 | 0.16 | 0.21 | 0.62 |
| 12 6 3 2 1 1 / 988 | 6,671 | 159 | 148 | 1,452 | 0.218 | 158 | 114 | 0.078 | 5.65e4 | 0.0258 | 172 | 149-155 | 0.14 | 0.21 | 0.58 |
| 13 7 4 3 1 1 / 2000 | 7,593 | 168 | 152 | 1,660 | 0.219 | 166 | 121 | 0.073 | 8.98e4 | 0.0411 | 312 | 251 | 0.26 | 0.22 | 0.60 |
| 13 7 4 3 1 1 / 2100 | 11,306 | 192 | 151 | 2,230 | 0.197 | 190 | 161 | 0.072 | 2.01e5 | 0.100 | 1,133 | 1,019 | 0.42 | 0.21 | 0.55 |
| 12 6 3 2 1 1 / 1200 | 11,698 | 199 | 150 | 2,151 | 0.184 | 197 | 168 | 0.078 | 1.21e5 | 0.0864 | 1,011 | 870 | 0.08 | 0.22 | 0.59 |
| 13 7 4 3 1 1 / 2200 | 15,199 | 211 | 152 | 2,713 | 0.179 | 210 | 206 | 0.076 | 3.43e5 | 0.246 | 3,739 | 3,432 | 0.64 | 0.21 | 0.56 |
| 11 6 4 3 2 1 / 2174 | 16,424 | 220 | 150 | 2,818 | 0.172 | 218 | 222 | 0.079 | 3.61e5 | 0.269 | 4,414 | 4,056 | 0.46 | 0.22 | 0.64 |
| 14 7 4 4 1 0 0 1 / 3648 | 20,538 | 252 | 129 | 3,214 | 0.157 | 250 | 242 | 0.075 | 2.84e5 | 0.263 | 5,396 | 5,407 | 0.11 | 0.20 | 0.58 |
| 9 6 4 3 1 1 1 1 / 2700 | 20,896 | 259 | 131 | 3,277 | 0.157 | 256 | 237 | 0.072 | 2.21e5 | 0.232 | 4,841 | 4,457 | 0.05 | 0.22 | 0.60 |
| 13 7 4 3 1 1 / 2400 | 22,993 | 245 | 128 | 3,569 | 0.155 | 243 | 288 | 0.081 | 6.48e5 | 0.559 | 12,859 | 11,785 | 0.66 | 0.21 | 0.62 |
| 13 7 4 3 1 1 / 2500 | 26,585 | 260 | 107 | 3,981 | 0.150 | 258 | 307 | 0.077 | 8.42e5 | 0.923 | 24,542 | 23,366 | 0.61 | 0.22 | 0.66 |
| 13 7 4 3 1 1 / 2650 | 31,743 | 279 | 106 | 4,445 | 0.140 | 276 | 359 | 0.081 | 9.80e5 | 1.28 | 40,474 | 39,917 | 0.50 | 0.20 | 0.57 |

* No pair occurred in the samples, consistent with 2-4x10^-5 expected per d.
* Every sum's V_d have 99% of the sum's labels. The top-label root makes the roots 7-8% of V_d
  (the rarest of 6 classes).
* Setup:
  * the posting index takes 1-5 ms per sum;
  * V_d from the postings takes 0.02-0.18% of the cycles;
  * the search setup of a V_d (labels, sort, bitsets) takes 0.08-1.4%;
  * the depth-0 list fills take 0.01-0.06%.
  * Together that is 1.6% at 5.9k and 0.09% at 31.7k.

## 2. The tree per d

Per d:
* "made" counts children created; "live" counts those that survived and were searched.
* k1, k11 and k12 are the children per searched node at (1,0), (1,1) and (1,2)/(2,1).
* p12 is the survival of the (1,2)/(2,1) children: rows-placed (2,1) and cols-placed (1,2)
  separately.

| sum | R live | k1 | (1,1) live | k11 | (1,2)+(2,1) made | p12 (1,2) / (2,1) | k12 | (2,2) made | (2,2) live | (1,3)+(3,1) made |
|---|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|
| a1950 | 92 | 17.6 | 1,610 | 17.5 | 28.1k | 0.25 / 0.32 | 2.36 | 19.6k | 0 | 1.2k |
| b988 | 110 | 19.3 | 2,120 | 17.9 | 37.9k | 0.12 / 0.20 | 2.29 | 15.1k | 0 | 1.3k |
| a2000 | 118 | 20.8 | 2,440 | 19.4 | 47.4k | 0.25 / 0.35 | 2.54 | 37.7k | 0 | 2.2k |
| a2100 | 158 | 25.8 | 4,040 | 23.6 | 95.4k | 0.29 / 0.40 | 2.81 | 96.3k | 0 | 4.7k |
| b1200 | 165 | 24.8 | 4,050 | 21.3 | 86.2k | 0.06 / 0.17 | 2.36 | 27.9k | 0 | 2.6k |
| a2200 | 203 | 30.0 | 6,060 | 26.2 | 159k | 0.28 / 0.41 | 2.94 | 171k | 0 | 7.9k |
| c2174 | 219 | 31.0 | 6,730 | 26.7 | 180k | 0.21 / 0.35 | 2.98 | 164k | 0 | 9.8k |
| c3648 | 238 | 30.5 | 7,230 | 26.0 | 188k | 0.07 / 0.20 | 2.66 | 81.7k | 0 | 6.7k |
| c2700 | 232 | 28.5 | 6,560 | 25.4 | 167k | 0.06 / 0.13 | 2.54 | 44.1k | 0 | 3.1k |
| a2400 | 284 | 35.4 | 10,000 | 30.2 | 302k | 0.24 / 0.38 | 3.14 | 321k | 0 | 14.9k |
| a2500 | 304 | 38.6 | 11,700 | 32.2 | 376k | 0.26 / 0.40 | 3.21 | 435k | 0 | 18.9k |
| a2650 | 355 | 40.3 | 14,200 | 33.6 | 478k | 0.19 / 0.35 | 3.18 | 464k | 0 | 23.4k |

* MRV at a root does not branch on T. Pairs of class-T vectors meeting only at T are ~93 per root
  at 15k, against k1 = 30.
* So the (1,1) nodes are pairs (r1, c_y) through r1's most constrained cell y.

**Against the plain search** (round-1 layer counts, whose nodes are unchanged by engine 3), per
sum:

| sum | N | d-first (1,1) searched | plain (1,1) | ratio | d-first (1,2)+(2,1) searched | plain | ratio | E pairs | d-first (1,1) nodes per E pair |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| a1950 | 5,896 | 9.5e6 | 1.41e5 | 67 | 4.9e7 | 9.7e6 | 5.0 | 0.16 | 5.9e7 |
| a2000 | 7,593 | 1.85e7 | 2.22e5 | 84 | 1.13e8 | 1.86e7 | 6.0 | 0.26 | 7.3e7 |
| a2100 | 11,306 | 4.6e7 | 4.46e5 | 102 | 3.9e8 | 5.21e7 | 7.5 | 0.42 | 1.1e8 |
| a2200 | 15,199 | 9.2e7 | 7.50e5 | 123 | 8.8e8 | 1.15e8 | 7.7 | 0.64 | 1.4e8 |
| a2400 | 22,993 | 2.3e8 | 1.57e6 | 146 | 2.35e9 | 3.53e8 | 6.7 | 0.66 | 3.5e8 |

**Redundancy across d** (`redund2.py` on the unreduced lists of a2000 / a2200 / a2400; 3,000
random structures each; T(d) is the rarest number of d in V_d):

| structure | 7.6k | 15.2k | 23k | slope |
|---|---:|---:|---:|---:|
| d with r1, c1 in V_d (plain root) | 472 | 724 | 911 | N^0.6 |
| ... and r1 a root (T(d) in r1: the top root) | 25.7 | 38.8 | 52.0 | N^0.64 |
| d with r1, c1, c2 in V_d (plain root) | 106 | 137 | 154 | N^0.34 |
| ... top root | 5.4 | 7.0 | 8.4 | N^0.40 |
| share of (r1; c1, c2) with >= 1 such d | 0.89 | 0.92 | 0.92 | |

## 3. Cycles

Exclusive DPROF cycles as a share of the search's cycles. "(1,2)/(2,1)" is the creation of those
children; "sel" is MRV plus selection at the searched nodes of the layer. The per-child cycles are
the creation cost per child made.

| sum | roots + sel | (1,1) create | (1,1) sel | (1,2)/(2,1) create | (1,2)/(2,1) sel | (2,2) create | deeper | unattr. | V_d + setup (of all) | cycles per (1,1) child | per (1,2)/(2,1) child | per (2,2) child | cycles per searched (1,1) node |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| a1950 | 1.3% | 8.8% | 1.2% | 64.9% | 3.4% | 11.8% | 0.5% | 8.0% | 1.6% | 2,420 | 1,020 | 268 | 27.6k |
| b988 | 1.4% | 10.2% | 1.3% | 70.0% | 2.3% | 7.4% | 0.5% | 6.9% | 1.5% | 2,550 | 984 | 261 | 25.1k |
| a2000 | 1.1% | 9.0% | 1.0% | 64.7% | 3.3% | 12.9% | 0.5% | 7.5% | 1.0% | 3,120 | 1,160 | 292 | 34.9k |
| a2100 | 0.8% | 8.7% | 0.9% | 63.9% | 3.3% | 15.3% | 0.5% | 6.6% | 0.6% | 4,470 | 1,400 | 333 | 51.7k |
| b1200 | 1.2% | 12.8% | 1.1% | 72.3% | 1.4% | 6.2% | 0.4% | 4.6% | 0.6% | 5,670 | 1,510 | 397 | 44.5k |
| a2200 | 0.7% | 8.9% | 0.7% | 65.8% | 2.7% | 15.8% | 0.5% | 5.1% | 0.3% | 7,520 | 2,140 | 476 | 85.0k |
| c2174 | 0.7% | 9.4% | 0.7% | 67.6% | 2.3% | 13.9% | 0.6% | 4.8% | 0.3% | 7,770 | 2,120 | 477 | 83.5k |
| c3648 | 0.9% | 11.9% | 0.9% | 73.3% | 1.4% | 7.6% | 0.4% | 3.7% | 0.3% | 8,990 | 2,140 | 511 | 76.0k |
| c2700 | 1.1% | 13.4% | 1.0% | 75.4% | 1.0% | 4.7% | 0.2% | 3.3% | 0.4% | 9,780 | 2,190 | 514 | 73.9k |
| a2400 | 0.6% | 8.7% | 0.6% | 68.1% | 2.3% | 14.9% | 0.4% | 4.4% | 0.2% | 10,200 | 2,650 | 543 | 117k |
| a2500 | 0.5% | 8.1% | 0.6% | 68.4% | 2.1% | 16.2% | 0.4% | 3.6% | 0.1% | 13,400 | 3,530 | 723 | 166k |
| a2650 | 0.6% | 8.6% | 0.6% | 71.5% | 2.0% | 13.1% | 0.4% | 3.3% | 0.1% | 16,100 | 4,010 | 756 | 189k |

**Phases** (CHILD_PROF, share of its rdtsc total, inflated short phases):

| sum | (1,0)+(1,1) children | (1,2)/(2,1): filter o | filter b | support | count + kids | all | (2,2): filter o | filter b | cross + support | all | cycles per (1,2)/(2,1) child: filter o, filter b, support |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| a1950 | 9.4% | 13.7% | 14.4% | 22.7% | 5.6% | 56.4% | 5.6% | 2.9% | 0.7% | 9.2% | 280, 295, 464 |
| b988 | 10.9% | 16.2% | 17.2% | 22.8% | 4.0% | 60.2% | 3.5% | 1.9% | 0.3% | 5.7% | 331, 350, 464 |
| a2000 | 9.2% | 13.9% | 14.9% | 22.3% | 5.5% | 56.6% | 6.3% | 3.3% | 0.8% | 10.4% | 348, 374, 558 |
| a2100 | 8.5% | 13.8% | 16.0% | 21.7% | 5.8% | 57.3% | 7.6% | 4.0% | 0.8% | 12.4% | 395, 457, 622 |
| b1200 | 13.1% | 17.7% | 22.3% | 21.4% | 3.0% | 64.4% | 3.4% | 1.8% | 0.3% | 5.5% | 520, 653, 628 |
| a2200 | 8.9% | 14.7% | 19.2% | 21.2% | 5.3% | 60.4% | 8.3% | 4.6% | 0.8% | 13.7% | 594, 779, 860 |
| c2174 | 9.5% | 15.5% | 20.3% | 22.1% | 4.6% | 62.5% | 7.4% | 4.0% | 0.6% | 12.0% | 576, 750, 819 |
| c3648 | 12.0% | 18.1% | 24.8% | 22.3% | 2.8% | 68.0% | 4.3% | 2.1% | 0.2% | 6.6% | 588, 808, 725 |
| c2700 | 14.2% | 20.1% | 26.0% | 21.5% | 1.9% | 69.5% | 2.7% | 1.3% | 0.1% | 4.1% | 642, 834, 690 |
| a2400 | 8.8% | 15.8% | 21.5% | 21.8% | 4.4% | 63.5% | 8.3% | 4.3% | 0.7% | 13.3% | 681, 924, 940 |
| a2500 | 8.5% | 16.1% | 24.2% | 20.1% | 4.4% | 64.8% | 9.2% | 4.7% | 0.8% | 14.7% | 915, 1375, 1139 |
| a2650 | 9.0% | 17.4% | 26.1% | 20.8% | 3.9% | 68.2% | 7.6% | 3.8% | 0.5% | 11.9% | 1048, 1576, 1260 |

The (1,1) children's filter o scans rows1 (2,049 entries at 15k) and keeps 25%. It is 4.5% of the
cycles there.

## 4. Lists

Means per layer:
* rows1 and cols1 are a root's lists.
* "(1,1) lists" are a searched (1,1) node's lists.
* "in" is the parent's lists scanned per (1,2)/(2,1) child (o, b).
* "after o / b" is the child's lists after its two filters.
* (2,2) gives the same for the (2,2) children.

| sum | rows1 | cols1 | (1,1) lists | (1,2)/(2,1) in o + b | after o / b | searched (1,2)/(2,1) lists | (2,2) in o / b | after o / b |
|---|---:|---:|---|---:|---|---:|---|---|
| a1950 | 966 | 275 | 247 / 187 | 435 | 58 / 63 | 63 | 78 / 55 | 20 / 6.1 |
| a2000 | 1,190 | 312 | 298 / 213 | 509 | 64 / 75 | 73 | 94 / 61 | 22 / 6.8 |
| a2100 | 1,650 | 377 | 395 / 264 | 652 | 76 / 99 | 92 | 121 / 72 | 26 / 7.6 |
| a2200 | 2,040 | 431 | 479 / 305 | 775 | 84 / 118 | 109 | 148 / 81 | 28 / 8.4 |
| a2400 | 2,720 | 510 | 622 / 367 | 976 | 95 / 144 | 133 | 189 / 93 | 31 / 8.3 |
| a2500 | 3,070 | 551 | 694 / 399 | 1,075 | 102 / 163 | 150 | 216 / 101 | 34 / 9.6 |
| a2650 | 3,430 | 589 | 757 / 427 | 1,164 | 104 / 165 | 156 | 225 / 103 | 33 / 8.2 |

## 5. Exponents

Least squares of ln q on ln N:
* ladder A: all 7 points (5.9-31.7k);
* ladder A's 5 points with <= 256 labels (5.9-23k);
* ladder A at 15-31.7k;
* all 12 sums.

| quantity | A: 7 | A: <= 256 labels | A: 15-31.7k | all 12 |
|---|---:|---:|---:|---:|
| CPU per sum | 3.44 | 3.43 | 3.27 | 3.35 |
| cycles per d | 2.44 | 2.43 | 2.27 | 2.35 |
| V_d | 0.69 | 0.70 | 0.67 | 0.70 |
| labels of V_d | 0.35 | 0.35 | 0.37 | 0.37 |
| roots R | 0.79 | 0.82 | 0.75 | 0.76 |
| k1 | 0.49 | 0.51 | 0.41 | 0.47 |
| searched (1,1) per d | 1.28 | 1.33 | 1.16 | 1.23 |
| searched (1,1) per sum | 2.28 | 2.33 | 2.16 | 2.23 |
| k11 | 0.39 | 0.41 | 0.35 | 0.38 |
| p12 (survival) | 0.05 | 0.12 | -0.19 | 0.04 |
| searched (1,2)/(2,1) per sum | 2.72 | 2.86 | 2.31 | 2.64 |
| k12 | 0.18 | 0.21 | 0.12 | 0.18 |
| cycles per (1,1) child | 1.14 | 1.09 | 1.04 | 1.13 |
| cycles per (1,2)/(2,1) child | 0.83 | 0.73 | 0.87 | 0.80 |
| cycles per (2,2) child | 0.65 | 0.55 | 0.66 | 0.63 |
| parent lists scanned per (1,2)/(2,1) child | 0.59 | 0.60 | 0.56 | 0.59 |
| its lists after filter o / b | 0.35 / 0.59 | 0.37 / 0.61 | 0.30 / 0.48 | 0.36 / 0.58 |
| cycles per searched (1,1) node | 1.17 | 1.10 | 1.12 | 1.14 |
| of which (1,2)/(2,1) creation (k11 x c12) | 1.22 | 1.14 | 1.22 | 1.18 |
| of which (2,2) creation | 1.27 | 1.29 | 0.94 | 1.22 |
| **E pairs per sum (output)** | **0.71** | 1.08 | -0.29 | 0.51 |
| E pairs per d | -0.29 | 0.08 | -1.29 | -0.49 |
| rows1 / cols1 | 0.75 / 0.45 | | | |
| (1,1) lists rows / cols | 0.67 / 0.49 | | | |

* The bitset width costs x1.25 per word in V_d, at identical nodes.
  * Measured with DPROF_MINW, min of 2 alternating runs, the same d.
  * S = 2000 at 3 -> 4 words: 7.08 -> 8.83 s on 13.65M nodes.
  * S = 2200 at 4 -> 5 words: 13.06 -> 16.34 s on 16.55M nodes.
  * Ladder A's V_d run at 3 words to 11k, 4 words at 15-23k and 5 words at 26-32k. That is
    1.25^2 over the ladder, N^0.27 of its exponent.
* The model's output saturates above 23k on ladder A: E_semi rises ~N^1.9 but 720 p_SP falls ~N^-1.2.
  The exponent of the output over 15-31.7k is therefore negative. On the other P (fixed x ~
  0.25), pairs per sum are 0.05-0.46 at 16-21k.

## 6. How the children die

Share of the children made (CHILD_PROF):

| sum | layer | made / d | fc o | count b | fc b | cross | support | live |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| a1950 | (1,2) | 10.6k | 7.9% | 0.0% | 0.8% | | 66.6% | 24.8% |
| a1950 | (2,1) | 17.6k | 12.1% | 0.0% | 1.9% | | 53.7% | 32.3% |
| a1950 | (2,2) | 19.6k | 36.3% | 24.0% | 31.2% | 6.6% | 1.9% | 0 |
| a2200 | (1,2) | 49.6k | 6.1% | 0.0% | 0.3% | | 65.9% | 27.6% |
| a2200 | (2,1) | 109k | 11.2% | 0.0% | 1.1% | | 47.0% | 40.7% |
| a2200 | (2,2) | 171k | 28.6% | 22.3% | 35.9% | 1.8% | 11.3% | 0 |
| c3648 | (1,2) | 53.8k | 7.2% | 0.0% | 1.0% | | 84.9% | 6.9% |
| c3648 | (2,1) | 134k | 14.0% | 0.0% | 2.9% | | 63.0% | 20.1% |
| c3648 | (2,2) | 81.7k | 32.2% | 30.4% | 31.1% | 2.3% | 4.1% | 0 |

(All 12 sums are in `cprof_tables.txt`. The pattern is the same everywhere: (1,2) children die
65-85% at support and (2,1) 47-65%, more where p12 is low; (2,2) children die 86-94% at fc o +
count b + fc b.)

**What the support filter kills on** (DPROF, d-sampled at stride 150-450):

| sum | layer | passes per call | killed | by count | by a dead cell | dead cell in the new vector | dead cell is one of d's numbers | only d's numbers |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| a2000 | (1,2) | 4.46 | 75% | 19% | 81% | 61% | 5.0% | 1.5% |
| a2000 | (2,1) | 4.17 | 62% | 17% | 83% | 61% | 7.4% | 1.9% |
| a2200 | (1,2) | 5.22 | 70% | 21% | 79% | 64% | 4.2% | 1.1% |
| a2200 | (2,1) | 4.64 | 53% | 18% | 82% | 63% | 6.8% | 1.7% |
| a2400 | (1,2) | 5.48 | 70% | 22% | 78% | 63% | 4.0% | 1.1% |
| a2400 | (2,1) | 4.70 | 50% | 18% | 82% | 63% | 6.7% | 1.6% |
| c3648 | (1,2) | 4.00 | 91% | 25% | 75% | 64% | 4.0% | 1.0% |
| c3648 | (2,1) | 4.12 | 72% | 20% | 80% | 63% | 6.7% | 1.3% |

* The support filter iterates 4-5.5 passes over the ~90-165-entry lists before the child dies.
* The usual reason is that a cell of the new vector (c2 or r2) has no supported partner left.
  * So the child is dead because of how the new vector meets the (1,1) node's lists, not because
    of d's class structure.
  * Forward checking on d's numbers was already measured at 1% (ideas.md).

## 7. Per d, and output per d

* **Cost law.** Per d, pooled over the 12 sums (1,701 d):
  * ln cycles = 5.65 ln|V_d| + 1.43 ln R - 2.68 ln N + c, with resid sd 0.29 and R^2 0.967.
  * |V_d| alone gives slope 3.49 (R^2 0.66). Within each sum the |V_d|, R slopes are 5.1-6.1 and
    1.25-1.5.
  * Cost per d is homogeneous in the scheduling sense: the top 10% of d hold 20-22%, CV 0.55-0.66.
* **Output per d.** The model's pairs per d (E pairs / N) are 2.7e-5 at 5.9k, 4.2e-5 at 15.2k,
  2.9e-5 at 23k and 1.6e-5 at 31.7k (~N^-0.3 along A).
  * Within a sum, the transversals of 145 squares of an r1-sampled plain run at a1950 were placed
    by their predicted V_t cost (`trav.py`):

  | cheapest share of d by predicted cost | share of the cost | share of all transversals (104k) | share of S-traversals (100) |
  |---:|---:|---:|---:|
  | 25% | 6% | 12% | 9% |
  | 50% | 23% | 25% | 18% |
  | 75% | 50% | 39% | 33% |
  | 90% | 75% | 52% | 43% |

  * The yield follows the cost.
  * Skipping or ordering d by predicted cost does not lower the cost per pair, under the stated
    assumption that P(SP) is independent of |V_t|.
* **Vector usage** (DPROF_VEC):

  | | a1950 | a2200 | a2400 | a2650 | c2700 |
  |---|---:|---:|---:|---:|---:|
  | share of V_d placed in a searched node | 0.76 | 0.90 | 0.91 | 0.91 | 0.75 |
  | share placed in a depth-3 node | 0.74 | 0.89 | 0.89 | 0.89 | 0.66 |
  | share placed at depth >= 4 | 0 | 1e-5 | 0 | 0 | 0 |

  * In a square: none in the samples. The expected share is 12 x (pairs per d) / |V_d| ~ 2x10^-7.

## 8. Synthesis: the d-first time law

Per sum, with N_d ~ N diagonals:

    T ~ N x R x k1 x [ c11 + k11 x ( c12 + p12 x k12 x c22 ) ]

| term | growth | what it is |
|---|---|---|
| N x R x k1 = searched (1,1) nodes, i.e. refuted (d, r1, c1) triples | N^2.28 | each (r1, c1) is re-refuted under ~N^0.6 diagonals (26 -> 52 from 7.6k to 23k); output per triple falls N^-1.6 |
| k11 | N^0.39 | candidates through the MRV cell of a (1,1) node |
| c12, cycles per (1,2)/(2,1) child | N^0.83 | parent lists scanned (N^0.59) + bitset width (N^0.27) |
| p12 k12 c22, the (2,2) layer | N^0.8 per child, 5-16% of the time | all dead, 88-94% at the pretest's checks |
| c11 | N^1.14 | rows1 scan (N^0.75) + width, 8-13% of the time |
| output: pairs per sum | N^0.71 (N^1.1 to 23k, falling above) | E_semi x 720 p_SP |

* **The excess over the output is N^2.7**: N^1.57 from the number of (1,1) triples (the
  re-exploration of the same (r1, c1) under many d) and N^1.17 from the refutation cost of each.
  * The latter is N^0.39 more children, N^0.59 longer parent lists per child and N^0.27 bitset
    width.
* Unlike the plain search, d-first has no layer whose node count tracks the output.
  * Its lists stay short (~100-165 after filtering at 31.7k), so the (2,2) filters keep killing
    everything.
  * But the redundancy across d means that even the (1,1) layer grows faster than the output.
* Ranked by cycles at N 15-32k, the targets are:
  * **(i) the (1,2)/(2,1) creation, 66-75%.** Within it, filters o + b (28-43%) are the part
    that scales with the parent's list rather than the child's. The 26-34 siblings scan the same
    780-1,160 entries. The support cascade (20-23%) is the other part: 4-5.5 passes, killing on a
    cell of the new vector.
  * **(ii) the (2,2) children, 13-16%** on high-survival P, with the pretest's checks (fc o,
    count b, fc b) killing 88-94%.
  * **(iii) the (1,1) children, 8-13%.**
  * **(iv) the width steps**, x1.25 per word, 4 -> 5 words above 256 labels.
* None of these changes the N^2.28 of term (a). Only a decomposition that shares the refutation
  of an (r1, c1) or (r1; c1, c2) structure across its ~7-52 compatible d would. Note that
  91-92% of (r1; c1, c2) triples lie in at least one top-root V_d tree.

## Files

Scratch directory: `/tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/c2/profile/`
* `sums.txt` lists the 12 sums with their d strides; `run.sh build tags...` runs them.
* Per-sum runs in `runs/`:
  * `dprof.<tag>.dp` has the per-d lines;
  * `.dlog` is msearch's d log;
  * `.jsonl` holds the dsum records;
  * `cprof.<tag>.err` has the CHILD_PROF tables;
  * `sup.<tag>.err` has the support-kill counts;
  * `plain.a1950.jsonl` holds the 145 sampled squares.
* Analysis:
  * `an.py` and `rep.py` produce `rep_dprof.txt`;
  * `expo.py` produces `expo.txt`;
  * `cprof.py` produces `cprof_tables.txt`;
  * `cmp.py` produces `cmp.txt`;
  * `epairs.py` produces `epairs.out` (the model's pairs);
  * `redund2.py` computes the multiplicity across d;
  * `trav.py` produces `trav_a1950.txt`.
* Per-d line: `D i |V_d| labels roots nodes pairs vd_cyc setup_cyc rootfill_cyc search_cyc searches
  cpu C <class sizes ascending> U <vused k=1..17> L<9r+c>:made,live,mk,sel,kids,in_o,in_b,out_o,
  n_o,out_b,n_b,lst_row,lst_col ...`
* Env: `DPROF_LOG=file`, `DPROF_VEC=1` (the vector marks; slower), `DPROF_MINW=k` (forced width).
