# How the arrangement search's work grows with N

Part 1 of the "cost per square at large N" study (cx/profile), 2026-10-08.

* Code: fdb77fc (`claude/magic-search-speedups`) plus instrumentation that is compiled out by
  default, on branch `cx/profile` (commit d5daef6, worktree `.claude/worktrees/cx-profile`).
  With the default flags the build is unchanged: same nodes on `quick.txt` (1,770,779) and
  `full.txt` (14,958,507).
* Machine: the shared 4-core Xeon with AVX-512 (VPOPCNTDQ, VBMI, GFNI). Native build
  (`-O3 -march=native -flto`), run one process at a time.
* Search CPU used: ~1,250 s, all of it r1-sampled.

## Bottom line

**1. Time per sum grows like N^4.1-4.3 on the carried-bitset path (up to 256 labels).**
* Three ladders from N = 2k to 23k give the same answer:
  * 13 7 4 3 1 1 at S = 1850-2400;
  * 12 6 3 2 1 1 at S = 840-1200;
  * x = ln(S/S_min) ~ 0.25 across 7 P.
* Nodes grow like N^3.6-4.0.
* N^0.3 of the growth is an implementation step, not the algorithm. Each extra 64-bit word of
  label bitset costs x1.38 time at identical nodes, and the width steps up at 128 and 192
  labels.
* **The existence study's "local slope 5.2 above 2x10^4" is the matrix path.**
  * The S = 2650 sum had 279 labels, so it ran on the matrix path: no cross support and an
    8-slice MRV.
  * On the same r1 of S = 2400 (N = 23k), forcing that path costs x2.7 time and x2.8 nodes.
  * On the carried path the slope stays ~4 up to N = 23k. N^4 extrapolates the 2400 sum to
    ~150k CPU-s at N = 31.7k; times 2.7 that is ~400k, against the 375k measured.

**2. The search is output-proportional down to the (1,2)/(2,1) layer. All of the super-linear
excess is below it.**
* Searched (1,2)/(2,1) nodes (r1, c1 and one more row or col) grow like **N^2.7** on all three
  ladders (2.67-2.70).
* A roughly constant fraction of them, **1-3x10^-5 at every N**, has a square below it:
  2.1e-5 at N = 4.1k, 3.0e-5 at 7.6k, 2.4e-5 at 15k, 1.3e-5 at 23k (3-35 squares per sample).
* So squares are ~2x10^-5 x (searched (1,2)/(2,1) nodes). Everything down to and including
  the creation of those nodes is 3-7% of the time.
* **The cycles spent per searched (1,2)/(2,1) node grow x12 from N = 4.1k to 23k**: 27.6k to
  340k TSC cycles in the profiling build, i.e. **~N^1.5**. That is the growth of the cost per
  square.

**3. That per-node cost has two parts.**

* **(a) Creating the (2,2) children: 48-64% of the time** (82% at x = 0.74).
  * Each (1,2)/(2,1) node has k12 ~ N^0.53 of them (12 to 30).
  * Each costs ~N^0.82 cycles (1,440 to 5,940).
  * The exactly-once filter scans the parent's long list, which grows like N^0.88 (330 to 1,420
    entries). It keeps 20% of it at N = 4k and 11% at N = 23k. Filter o alone goes from 373 to
    2,441 cycles per child.
  * Cross support then works on the filtered lists, which grow like N^0.53 (65 to 157 entries):
    713 to 2,023 cycles.
  * The number of (2,2) children itself grows like N^3.2. That is close to the output: 0.4-3.6
    million (2,2) children per square.
* **(b) The dead subtrees under the surviving (2,2) children.**
  * The share of (2,2) children that survive cross support and the support filter rises from
    3.6% at N = 2.5k to 34-39% at N >= 11k (13 7 4 3 1 1).
  * At x ~ 0.25 and N = 15-20k it is 11-37% depending on P. It is near 0 at x >= 0.7.
  * Each survivor creates k22 ~ N^0.6 children (4.6 to 18), and nearly all of them die at once.
  * The (2,3)/(3,2) children go from 13% to 85% of all nodes, and from 10% (N = 4.1k) to 32% of
    the cycles.
  * The dead work below the survivors grows x41 per (1,2)/(2,1) node between N = 4.1k and 23k
    (~N^2.2), the fastest-growing term.

**4. Why the (2,2) filters weaken with N.**
* Cross support and the support filter test a candidate against unions of other-axis
  candidates.
  * At N = 23k there are ~140-160 candidates per axis after the filters. Every candidate passes
    through one of the 4 unmatched cells of each placed vector, so there are ~35-40 through each
    cell.
  * With unions that large, a candidate's free labels almost always meet every one of them.
    This is an interpretation; the coverage of the unions was not measured.
* Cross support kills 50% of the (2,2) children at N = 2.5k and 10% at N = 23k.
* The survivors are consistent under one round of the filters, but not globally:
  * Cross support and the support filter, iterated to their common fixpoint, kill **86-95% of the
    searched (2,2) nodes at every N**. This is sound: none of the productive ones died.
  * Only 1x10^-5 to 1x10^-6 of the searched (2,2) nodes have a square below them, a fraction
    that falls like ~N^-1.3.
* Their children die by count (`cnt_b`: 47% to 77% of them):
  * Placing a third col (or row) cuts the other axis from ~60-130 candidates (N = 4-23k) to the
    ~12-15 that meet it exactly once.
  * Their union then supports only 0.3-1.5 of the remaining same-axis candidates on average
    (3 needed).
  * Same-axis disjointness is not what binds: 42-65 same-axis candidates are disjoint from the
    new vector.

**5. Per r1.**
* **Work per r1 is ~kids1^0.85 x cols1^3.0-3.3**, with R^2 = 0.98-1.00 in every sum.
  * cols1 is the number of vectors meeting r1 exactly once; kids1 is the number of choices of
    c1.
  * Squares per r1 go like cols1^2.1 kids1^0.5 (pooled Poisson fit). So the cost per square of
    an r1 grows like ~cols1^1 kids1^0.4.
  * Pooled across sums at fixed list sizes, larger-N sums have fewer squares per r1 (N^-1.4)
    and more work (N^0.4).
* **The work is not concentrated in a few r1.**
  * The top 1% of r1 by cost hold 5-10% of the time, the top 10% ~40% and the top half 93-97%.
  * The first half by index (the largest universes) holds ~92-94% of the time.
  * The r1 with a square are 2% (N = 4k) to 10-20% (N = 11-23k) of all r1, and they hold 5-50%
    of the time.
* Within a sum, (2,2) survival goes from 0 at the r1 with short lists to 70-85% at the r1 with
  the longest lists. At a fixed list size it falls with N.

**6. Implication for an output-sensitive search.** The work below a (1,2)/(2,1) node has to stop
growing with N. In decreasing order of the time at stake at N = 15-23k:
* **Create the (2,2) children in time proportional to their own lists, not the parent's.**
  * At N = 23k, 89% of the parent's entries scanned per child are discarded.
  * An example would be a label -> position index of the parent's lists, shared by its 12-30
    children.
  * The round-2 label -> position masks (opt2/residual) were tried at N ~ 2-3k, with lists of
    ~200 and 8-12 children. They lost there to the folded exactly-once test. At N = 23k the
    parent lists are 7x longer and there are 30 children per node.
  * Filter o is 16-21% of the time. With filter b, cross support and the support filter, the
    (2,2) creation is 48-64%.
* **Kill the dead survivors before their 15-18 children.**
  * The fixpoint catches ~90% of them. research/ideas.md found it neutral at N ~ 2-3k.
  * But the survivors (count, select and children) take 29-34% of the time at N = 15-23k,
    against 11% at N = 4k and ~3% on bench/prod.txt, where ideas.md measured it.
  * The remaining ~10% need a test that sees the "meets exactly once" structure between the
    two axes, which unions cannot.
* **Change the decomposition below (1,2)/(2,1).**
  * That is the last layer whose productive fraction holds with N.
  * A completion method whose cost per (1,2)/(2,1) node does not depend on N would make the
    cost per square constant.
* **Lift the 256-label limit of the carried path**: a factor 2.7 at N = 23k. Most sums with
  N >= 25k have more than 256 labels.

## Method

* **r1 sampling.**
  * The root loop runs r1 = off, off + k, off + 2k, ... (env `SAMPLE_STRIDE`, `SAMPLE_OFFSET`),
    as in the existence study, but on fdb77fc.
  * k x (sampled totals) is unbiased for the sum.
  * Every diagnostic run uses a subset of the same r1 as the ladder run.
  * Check: the two halves of a stride-2 run add up to the full search (same squares). The nodes
    differ by 3 in 960k, because the adaptive cross mode (`xs_tick`) depends on the order of
    the children.
* **Builds** (`src/c/arrange.c`, `arrange_core.h` on cx/profile):
  * `-DSCALE_PROF`: counters only, ~8% overhead. These give the totals, layers and lists.
  * `-DSCALE_PROF -DCHILD_PROF`: rdtsc per phase, ~80% overhead. These give the phase shares and
    the cycles per call. The shares include rdtsc overhead, which inflates the short phases.
  * `-DSCALE_FIX`: the fixpoint test on a copy of the lists of each searched (2,2) node, and
    the disjoint-only count at the (2,3)/(3,2) children.
* **Sampling error.** Standard errors of the time and node totals (simple-random-sampling
  formula, conservative for systematic samples) are 3-5% at N <= 8.6k, 9-18% at 11-19k, and
  22-34% at 20-23k (29-76 r1).
  * Ratios within a sample (layer shares, list sizes, cost per node) are much tighter.
  * Square counts at N >= 11k rest on 1-24 squares per sample. Where better counts exist they
    are used: two full runs, and the existence study's samples of S = 2200 and 2400.
* **Notation.**
  * (r,c) is a node or child with r rows and c cols placed.
  * "made" counts children created; "searched" counts nodes that survived their creation.
  * k1, k11, k12 and k22 are children per searched node at (1,0), (1,1), (1,2)/(2,1) and (2,2).
  * o is the axis that must meet the new vector exactly once; b is the axis of the new vector
    (disjoint from it).

## 1. Per sum: time, nodes, squares

P = 13 7 4 3 1 1 (S_min 1718, x = 0.07-0.33). At S = 2200 and 2400 the squares come from the
pooled samples (2,533 from 38 squares in 228 r1, and 2,745 from 8 in 67):

| P | S | N | labels | W | r1 sampled | nodes | CPU-s | squares (sample) | squares est. | CPU-s / square | nodes / square | cycles / node |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 7 4 3 1 1 | 1850 | 2,523 | 115 | 2 | 2,523 | 1.05e7 | 3.05 | 10 | 10 | 0.305 | 1.05e6 | 609 |
| 13 7 4 3 1 1 | 1900 | 4,111 | 137 | 3 | 4,111 | 8.86e7 | 35.9 | 77 | 77 | 0.466 | 1.15e6 | 850 |
| 13 7 4 3 1 1 | 1950 | 5,896 | 153 | 3 | 1,474 | 5.05e8 | 182 | 85 | 340 | 0.536 | 1.49e6 | 757 |
| 13 7 4 3 1 1 | 2000 | 7,593 | 168 | 3 | 759 | 1.38e9 | 528 | 54 | 540 | 0.977 | 2.56e6 | 801 |
| 13 7 4 3 1 1 | 2100 | 11,306 | 192 | 3 | 226 | 6.76e9 | 2,484 | 24 | 1,201 | 2.07 | 5.63e6 | 774 |
| 13 7 4 3 1 1 | 2200 | 15,199 | 211 | 4 | 76 | 1.92e10 | 10,217 | 10 | 2,533 | 4.03 | 7.58e6 | 1,119 |
| 13 7 4 3 1 1 | 2300 | 18,957 | 228 | 4 | 47 | 4.49e10 | 24,178 | 7 | 2,823 | 8.56 | 1.59e7 | 1,137 |
| 13 7 4 3 1 1 | 2400 | 22,992 | 245 | 4 | 29 | 7.03e10 | 42,433 | 5 | 2,745 | 15.5 | 2.56e7 | 1,269 |

P = 12 6 3 2 1 1 (S_min 706, x = 0.17-0.82). S = 988 and 1080 use the full runs' counts.
* Beyond S = 1200 (x > 0.55) this P's squares vanish, which the sampled counts show.
* The nodes fall while the cycles per node rise: the (2,2) creation regime of section 4.

| P | S | N | labels | W | r1 sampled | nodes | CPU-s | squares (sample) | squares est. | CPU-s / square | nodes / square | cycles / node |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 12 6 3 2 1 1 | 840 | 2,812 | 114 | 2 | 2,812 | 1.88e7 | 5.54 | 32 | 32 | 0.173 | 5.88e5 | 632 |
| 12 6 3 2 1 1 | 920 | 4,803 | 141 | 3 | 4,803 | 1.47e8 | 74.9 | 108 | 108 | 0.694 | 1.36e6 | 1,075 |
| 12 6 3 2 1 1 | 988 | 6,671 | 159 | 3 | 1,111 | 5.98e8 | 273 | 52 | 245 | 1.11 | 2.44e6 | 964 |
| 12 6 3 2 1 1 | 1080 | 8,891 | 179 | 3 | 444 | 1.24e9 | 656 | 16 | 275 | 2.39 | 4.50e6 | 1,180 |
| 12 6 3 2 1 1 | 1200 | 11,697 | 199 | 4 | 183 | 3.17e9 | 2,207 | 6 | 384 | 5.75 | 8.25e6 | 1,483 |
| 12 6 3 2 1 1 | 1320 | 13,565 | 216 | 4 | 106 | 3.04e9 | 3,034 | 1 | 128 | 23.7 | 2.38e7 | 2,103 |
| 12 6 3 2 1 1 | 1480 | 15,074 | 236 | 4 | 78 | 2.13e9 | 3,106 | 0 | 0 | n/a | n/a | 3,085 |
| 12 6 3 2 1 1 | 1600 | 15,813 | 250 | 4 | 66 | 1.76e9 | 2,806 | 0 | 0 | n/a | n/a | 3,375 |

x ~ 0.25 across P. b920 and a2200 are reused here. Squares at N >= 16k rest on 1-4 per sample:

| P | S | N | labels | W | r1 sampled | nodes | CPU-s | squares (sample) | squares est. | CPU-s / square | nodes / square | cycles / node |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 7 3 2 1 | 796 | 2,049 | 106 | 2 | 2,049 | 5.32e6 | 1.49 | 4 | 4 | 0.371 | 1.33e6 | 591 |
| 12 8 4 2 1 | 1115 | 3,164 | 124 | 2 | 3,164 | 2.84e7 | 7.81 | 29 | 29 | 0.269 | 9.81e5 | 577 |
| 12 6 3 2 1 1 | 920 | 4,803 | 141 | 3 | 4,803 | 1.47e8 | 74.9 | 108 | 108 | 0.694 | 1.36e6 | 1,075 |
| 13 7 4 3 2 | 2146 | 8,618 | 173 | 3 | 539 | 1.82e9 | 787 | 34 | 544 | 1.45 | 3.35e6 | 909 |
| 13 7 4 3 1 1 | 2200 | 15,199 | 211 | 4 | 76 | 1.92e10 | 10,217 | 10 | 2,533 | 4.03 | 7.58e6 | 1,119 |
| 11 6 4 3 2 1 | 2174 | 16,424 | 220 | 4 | 55 | 2.19e10 | 12,950 | 4 | 1,194 | 10.8 | 1.83e7 | 1,243 |
| 14 7 4 4 1 0 0 1 | 3648 | 20,538 | 252 | 4 | 34 | 1.72e10 | 16,830 | 1 | 604 | 27.9 | 2.84e7 | 2,061 |

Fitted exponents d ln q / d ln N (least squares on log-log), over all points and over the 4
largest-N points. Ladder B is fitted on its 5 points with x <= 0.53.
* The square counts at the top of each ladder rest on few sampled squares, and on 13 7 4 3 1 1
  they level off as x grows. So the "top 4" exponents of squares and of CPU per square are
  rough.

| quantity | A: all 8 points | A: top 4 | B (x <= 0.53): all 5 | B: top 4 | C: all 7 | C: top 4 |
|---|---:|---:|---:|---:|---:|---:|
| nodes | 4.02 | 3.37 | 3.60 | 3.37 | 3.81 | 2.94 |
| CPU time | 4.31 | 4.03 | 4.12 | 3.73 | 4.25 | 3.75 |
| squares | 2.51 | 1.17 | 1.75 | 1.34 | 2.40 | 0.51 |
| CPU per square | 1.80 | 2.86 | 2.37 | 2.39 | 1.85 | 3.25 |
| searched (1,2)+(2,1) nodes | 2.68 | 2.69 | 2.67 | 2.51 | 2.70 | 2.48 |
| k12 ((2,2) children per (1,2)/(2,1) node) | 0.53 | 0.40 | 0.51 | 0.44 | 0.53 | 0.38 |
| (2,2) children | 3.21 | 3.09 | 3.18 | 2.95 | 3.22 | 2.86 |
| (2,2) children per square | 0.70 | 1.93 | 1.42 | 1.61 | 0.82 | 2.35 |
| (2,2) survival | 0.89 | -0.06 | 0.36 | 0.25 | 0.66 | -0.42 |
| searched (2,2) nodes | 4.10 | 3.03 | 3.54 | 3.20 | 3.88 | 2.44 |
| k22 | 0.61 | 0.40 | 0.63 | 0.59 | 0.61 | 0.51 |
| (2,3)+(3,2) children | 4.71 | 3.42 | 4.17 | 3.79 | 4.49 | 2.94 |
| parent o-list of a (2,2) child | 0.88 | 0.79 | 0.90 | 0.86 | 0.88 | 0.82 |
| its o-list after the filter | 0.55 | 0.45 | 0.53 | 0.50 | 0.53 | 0.45 |
| cycles per node | 0.29 | 0.66 | 0.53 | 0.38 | 0.44 | 0.82 |
| labels | 0.34 | 0.34 | 0.39 | 0.39 | 0.36 | 0.42 |

## 2. Where the nodes are: the tree layer by layer

Searched nodes, children per node, (2,2) survival, and the dead (2,3)/(3,2) layer:
* "share of nodes" is (2,3)+(3,2) children as a share of all nodes.
* "(2,2) live" is the share of (2,2) children that survive their filters and are searched.

P = 13 7 4 3 1 1:

| N | searched (1,0) | k1 | (1,1) | k11 | (1,2)+(2,1) | k12 | (2,2) made | (2,2) live | (2,2) searched | k22 | (2,3)+(3,2) made | share of nodes | (2,2) made / square |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,523 | 2,260 | 13 | 29,380 | 31.4 | 9.24e5 | 8.86 | 8.18e6 | 0.0359 | 2.93e5 | 4.61 | 1.35e6 | 0.129 | 8.18e5 |
| 4,111 | 3,813 | 18.4 | 70,232 | 48.7 | 3.42e6 | 12.3 | 4.19e7 | 0.149 | 6.23e6 | 6.93 | 4.32e7 | 0.487 | 5.44e5 |
| 5,896 | 5,528 | 25.5 | 1.41e5 | 68.9 | 9.70e6 | 15.4 | 1.50e8 | 0.248 | 3.72e7 | 9.3 | 3.46e8 | 0.684 | 4.40e5 |
| 7,593 | 7,253 | 30.6 | 2.22e5 | 83.8 | 1.86e7 | 17.6 | 3.28e8 | 0.295 | 9.68e7 | 10.7 | 1.04e9 | 0.749 | 6.08e5 |
| 11,306 | 10,856 | 41 | 4.46e5 | 117 | 5.21e7 | 21.5 | 1.12e9 | 0.363 | 4.06e8 | 13.8 | 5.59e9 | 0.827 | 9.31e5 |
| 15,199 | 15,199 | 49.3 | 7.50e5 | 154 | 1.15e8 | 24.6 | 2.84e9 | 0.372 | 1.06e9 | 15.3 | 1.62e10 | 0.846 | 1.12e6 |
| 18,957 | 18,554 | 56.9 | 1.06e6 | 195 | 2.06e8 | 27.6 | 5.69e9 | 0.39 | 2.22e9 | 17.6 | 3.90e10 | 0.869 | 2.01e6 |
| 22,992 | 22,992 | 68.1 | 1.57e6 | 225 | 3.53e8 | 28.2 | 9.97e9 | 0.336 | 3.35e9 | 17.9 | 5.99e10 | 0.853 | 3.63e6 |

P = 12 6 3 2 1 1:

| N | searched (1,0) | k1 | (1,1) | k11 | (1,2)+(2,1) | k12 | (2,2) made | (2,2) live | (2,2) searched | k22 | (2,3)+(3,2) made | share of nodes | (2,2) made / square |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,812 | 2,562 | 14.1 | 36,045 | 35 | 1.26e6 | 9.7 | 1.22e7 | 0.0825 | 1.01e6 | 5.2 | 5.25e6 | 0.279 | 3.82e5 |
| 4,803 | 4,478 | 22.6 | 1.01e5 | 59.2 | 6.00e6 | 13.4 | 8.04e7 | 0.0992 | 7.98e6 | 7.6 | 6.06e7 | 0.412 | 7.44e5 |
| 6,671 | 6,263 | 29.2 | 1.83e5 | 79.2 | 1.45e7 | 16.3 | 2.36e8 | 0.155 | 3.67e7 | 9.45 | 3.47e8 | 0.58 | 9.65e5 |
| 8,891 | 8,470 | 35.8 | 3.03e5 | 97.3 | 2.95e7 | 18.1 | 5.34e8 | 0.115 | 6.13e7 | 11 | 6.74e8 | 0.544 | 1.94e6 |
| 11,697 | 11,313 | 41.1 | 4.65e5 | 120 | 5.59e7 | 20 | 1.12e9 | 0.139 | 1.55e8 | 12.9 | 1.99e9 | 0.629 | 2.91e6 |
| 13,565 | 13,309 | 45.2 | 6.01e5 | 125 | 7.53e7 | 20.3 | 1.53e9 | 0.0706 | 1.08e8 | 13.3 | 1.44e9 | 0.472 | 1.20e7 |
| 15,074 | 13,721 | 47.4 | 6.51e5 | 132 | 8.61e7 | 19.9 | 1.71e9 | 0.0143 | 2.45e7 | 13.4 | 3.29e8 | 0.155 | n/a |
| 15,813 | 15,573 | 41.2 | 6.42e5 | 133 | 8.56e7 | 18.8 | 1.61e9 | 3.41e-3 | 5.48e6 | 12.9 | 7.06e7 | 0.04 | n/a |

x ~ 0.25 across P:

| N | searched (1,0) | k1 | (1,1) | k11 | (1,2)+(2,1) | k12 | (2,2) made | (2,2) live | (2,2) searched | k22 | (2,3)+(3,2) made | share of nodes | (2,2) made / square |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,049 | 1,818 | 11.4 | 20,786 | 25 | 5.20e5 | 7.95 | 4.13e6 | 0.037 | 1.53e5 | 4.13 | 6.31e5 | 0.119 | 1.03e6 |
| 3,164 | 2,865 | 14.8 | 42,493 | 36.3 | 1.54e6 | 10.5 | 1.62e7 | 0.115 | 1.87e6 | 5.69 | 1.06e7 | 0.374 | 5.59e5 |
| 4,803 | 4,478 | 22.6 | 1.01e5 | 59.2 | 6.00e6 | 13.4 | 8.04e7 | 0.0992 | 7.98e6 | 7.6 | 6.06e7 | 0.412 | 7.44e5 |
| 8,618 | 8,234 | 35 | 2.88e5 | 98 | 2.82e7 | 18.9 | 5.34e8 | 0.218 | 1.17e8 | 10.8 | 1.26e9 | 0.691 | 9.82e5 |
| 15,199 | 15,199 | 49.3 | 7.50e5 | 154 | 1.15e8 | 24.6 | 2.84e9 | 0.372 | 1.06e9 | 15.3 | 1.62e10 | 0.846 | 1.12e6 |
| 16,424 | 15,827 | 63.3 | 1.00e6 | 160 | 1.61e8 | 25.1 | 4.03e9 | 0.283 | 1.14e9 | 15.5 | 1.77e10 | 0.809 | 3.37e6 |
| 20,538 | 19,934 | 61.2 | 1.22e6 | 190 | 2.31e8 | 25.7 | 5.94e9 | 0.113 | 6.71e8 | 16.4 | 1.10e10 | 0.64 | 9.84e6 |

* **The upper tree's branching factors are universal.** On all three ladders:
  * k1 ~ N^0.74-0.77, k11 ~ N^0.86-0.89 and k12 ~ N^0.51-0.53;
  * (1,2)/(2,1) nodes ~ N^2.7 and (2,2) children ~ N^3.2.
* **The layer below is not universal.** The survival p22 depends on P and on x, and k22 grows
  like N^0.6.
* At N >= 11k on 13 7 4 3 1 1, 83-87% of all nodes are dead (2,3)/(3,2) children.
* Everything down to the creation of the (1,2)/(2,1) nodes is 3-7% of the time.

## 3. The productive fraction per layer

The fraction of searched nodes with at least one square in their subtree (`-DSCALE_FIX` runs,
13 7 4 3 1 1):

| 13 7 4 3 1 1 | N | (1,0) = r1 | (1,1) | (1,2)+(2,1) | (2,2) | (2,3)+(3,2) |
|---|---:|---:|---:|---:|---:|---:|
| S = 1900 | 4,111 | 1.8e-02 | 1.0e-03 | 2.1e-05 | 1.1e-05 | 2.1e-02 |
| S = 1950 | 5,896 | 4.7e-02 | 2.0e-03 | 2.8e-05 | 7.6e-06 | 2.2e-02 |
| S = 2000 | 7,593 | 6.8e-02 | 2.5e-03 | 3.0e-05 | 5.9e-06 | 2.2e-02 |
| S = 2100 | 11,306 | 9.3e-02 | 2.5e-03 | 2.1e-05 | 2.6e-06 | 1.7e-02 |
| S = 2200 | 15,199 | 1.6e-01 | 3.7e-03 | 2.4e-05 | 2.8e-06 | 2.8e-02 |
| S = 2200 (2nd sample) | 15,199 | 1.6e-01 | 2.8e-03 | 1.8e-05 | 1.8e-06 | 1.9e-02 |
| S = 2400 | 22,992 | 2.0e-01 | 3.1e-03 | 1.3e-05 | 1.1e-06 | 2.3e-02 |

* The fraction at (1,2)+(2,1) stays at 1-3x10^-5 from N = 4k to 23k.
* Between that layer and the searched (2,2) nodes it drops 2x at N = 4k and 12x at N = 23k.
  The leaf-side layer, (2,3)+(3,2), stays at ~2%.
* So the loss of selectivity with N is in the step from (1,2)/(2,1) to (2,2): more children
  (N^0.53), more of them surviving, and no change in how many lead to squares.

## 4. Cost per (1,2)/(2,1) node, and phase shares

Cycles per searched (1,2)/(2,1) node, by what they buy:
* These are TSC cycles of the rdtsc build (`-DCHILD_PROF`), whose overhead inflates the short
  phases.
* "below the survivors" is the creation of all children with five or more vectors placed.
* "above, amortized" is the (1,0) and (1,1) layers, the creation of the (1,2)/(2,1) nodes and
  their MRV and selection.

| sum (rdtsc build) | N | cycles per searched (1,2)/(2,1) node | creating its (2,2) children | k12 | cycles per (2,2) child | (2,2) survivors: count + select | below the survivors | above, amortized | unattributed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 7 4 3 1 1, S = 1900 | 4,111 | 27,554 | 17,632 | 12.3 | 1,439 | 503 | 2,632 | 1,962 | 4,824 |
| S = 2000 | 7,593 | 72,803 | 40,999 | 17.6 | 2,327 | 1,639 | 14,727 | 3,263 | 12,175 |
| S = 2100 | 11,306 | 126,609 | 61,152 | 21.6 | 2,836 | 3,139 | 36,175 | 5,203 | 20,940 |
| S = 2200 | 15,199 | 205,903 | 111,362 | 24.3 | 4,582 | 4,162 | 55,202 | 8,769 | 26,409 |
| S = 2400 | 22,992 | 339,846 | 175,548 | 29.5 | 5,944 | 5,608 | 108,939 | 11,430 | 38,322 |
| 12 6 3 2 1 1, S = 1480 | 15,074 | 78,875 | 64,282 | 19.9 | 3,231 | 149 | 1,370 | 7,413 | 5,661 |

Share of all cycles by phase:

| sum (rdtsc build) | N | (1,0), (1,1) and MRV at (1,2)/(2,1) | creating (1,2)/(2,1) | (2,2) children: filter o | filter b | cross | support + after | (2,2) survivors: count + select | creating (2,3)/(3,2) | unattributed | cycles per (2,2) child | cycles per (2,3)/(3,2) child |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 13 7 4 3 1 1, S = 1900 | 4,111 | 1.4% | 5.7% | 16.6% | 8.6% | 27.6% | 11.3% | 1.8% | 9.6% | 17.5% | 1,439 | 208 |
| S = 2000 | 7,593 | 0.5% | 4.0% | 17.4% | 6.4% | 22.0% | 10.6% | 2.3% | 20.2% | 16.7% | 2,327 | 270 |
| S = 2100 | 11,306 | 0.4% | 3.7% | 15.6% | 5.0% | 18.1% | 9.5% | 2.5% | 28.6% | 16.5% | 2,836 | 321 |
| S = 2200 | 15,199 | 0.3% | 4.0% | 17.9% | 6.9% | 19.5% | 9.7% | 2.0% | 26.8% | 12.8% | 4,582 | 428 |
| S = 2400 | 22,992 | 0.2% | 3.2% | 21.2% | 5.5% | 16.1% | 8.8% | 1.7% | 32.1% | 11.3% | 5,944 | 572 |
| 12 6 3 2 1 1, S = 1480 | 15,074 | 0.6% | 8.8% | 35.8% | 11.0% | 22.9% | 11.8% | 0.2% | 1.7% | 7.2% | 3,231 | 385 |

Cycles per call at the (2,2) children:

| N | filter o | filter b | cross (both passes) | support | after (cross + support again) | count (survivors) |
|---:|---:|---:|---:|---:|---:|---:|
| 4,111 | 373 | 193 | 713 | 339 | 878 | 124 |
| 7,593 | 718 | 264 | 1,010 | 469 | 1,313 | 147 |
| 11,306 | 916 | 297 | 1,184 | 551 | 1,572 | 183 |
| 15,199 | 1,520 | 583 | 1,812 | 792 | 2,386 | 225 |
| 22,992 | 2,441 | 635 | 2,023 | 930 | 2,678 | 252 |
| 15,074 (12 6 3 2 1 1) | 1,420 | 437 | 1,490 | 599 | 1,446 | 232 |

* The (2,2) children's creation stays at 48-64% of the time. Within it, the exactly-once filter
  over the parent's list grows from 17% to 21%, while cross support shrinks relatively.
* The (2,3)/(3,2) creation grows from 10% to 32%.
* At x = 0.74 (12 6 3 2 1 1, S = 1480) almost nothing survives, and creating the (2,2)
  children is 82% of the time: 36% for filter o alone.
* The cost per (2,2) child grows like N^0.82 and the cost per (2,3)/(3,2) child like N^0.59.
  Both follow the lengths of the lists they scan (section 5).
* The cycles per node of the whole search drift only N^0.3. The node mix shifts toward the cheap
  (2,3)/(3,2) children, which offsets the longer lists; the drift that remains is the width
  steps of section 8.

## 5. List lengths

Means at each layer:
* rows1 and cols1 are the depth-1 lists after r1 is placed; (1,1) gives that node's lists.
* "(1,2)/(2,1) short / long": the axis with two vectors of the other axis placed, and the other
  axis.
* "(2,2) child: parent o / b": the lists that the (2,2) children's filters scan, then their
  sizes after filters o and b.
* "(2,3)": the same for the children of the searched (2,2) nodes.

P = 13 7 4 3 1 1:

| N | rows1 | cols1 | (1,1) rows | (1,1) cols | (1,2)/(2,1) short | long | (2,2) child: parent o | parent b | after filter o | after filter b | searched (2,2) lists | (2,3) parent lists | (2,3) after o | (2,3) after b |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,523 | 870 | 369 | 284 | 266 | 84.8 | 185 | 203 | 93 | 46.2 | 38.9 | 34.9 | 40.7 | 10.2 | 1.54 |
| 4,111 | 1,477 | 532 | 434 | 403 | 112 | 299 | 327 | 123 | 64.6 | 57.6 | 49.6 | 57.4 | 11.9 | 1.56 |
| 5,896 | 2,198 | 700 | 581 | 548 | 136 | 422 | 462 | 148 | 81 | 73.8 | 64.2 | 73.5 | 13.1 | 1.45 |
| 7,593 | 2,918 | 841 | 723 | 677 | 155 | 531 | 581 | 169 | 93.3 | 86.2 | 73.4 | 84.1 | 13.7 | 1.36 |
| 11,306 | 4,529 | 1,123 | 999 | 944 | 187 | 748 | 816 | 204 | 115 | 107 | 91.7 | 104 | 14.5 | 1.15 |
| 15,199 | 6,102 | 1,350 | 1,280 | 1,194 | 215 | 964 | 1,038 | 232 | 132 | 124 | 102 | 116 | 14.8 | 0.968 |
| 18,957 | 7,913 | 1,648 | 1,496 | 1,414 | 238 | 1,187 | 1,275 | 255 | 148 | 138 | 116 | 132 | 15.3 | 0.886 |
| 22,992 | 9,479 | 1,828 | 1,613 | 1,574 | 241 | 1,319 | 1,416 | 258 | 157 | 138 | 117 | 132 | 14.6 | 0.66 |

P = 12 6 3 2 1 1:

| N | rows1 | cols1 | (1,1) rows | (1,1) cols | (1,2)/(2,1) short | long | (2,2) child: parent o | parent b | after filter o | after filter b | searched (2,2) lists | (2,3) parent lists | (2,3) after o | (2,3) after b |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,812 | 967 | 397 | 309 | 289 | 90.4 | 205 | 224 | 98.8 | 50.3 | 43.8 | 38.1 | 44.5 | 10.7 | 1.67 |
| 4,803 | 1,789 | 583 | 490 | 463 | 117 | 342 | 373 | 128 | 68 | 59.3 | 51.8 | 59.5 | 11.4 | 1.12 |
| 6,671 | 2,555 | 737 | 641 | 604 | 139 | 462 | 504 | 151 | 82.3 | 73.5 | 62.9 | 71.7 | 12.1 | 1.02 |
| 8,891 | 3,517 | 896 | 801 | 753 | 155 | 581 | 632 | 169 | 92.1 | 81.1 | 71.8 | 81.6 | 12.2 | 0.775 |
| 11,697 | 4,718 | 1,095 | 982 | 938 | 175 | 740 | 808 | 189 | 107 | 93.5 | 83.6 | 94.8 | 12.8 | 0.682 |
| 13,565 | 5,510 | 1,189 | 1,046 | 987 | 178 | 778 | 863 | 195 | 108 | 91.8 | 86.8 | 97.3 | 12.3 | 0.503 |
| 15,074 | 6,310 | 1,278 | 1,172 | 1,141 | 181 | 851 | 930 | 196 | 108 | 86 | 87.3 | 96.4 | 11.5 | 0.341 |
| 15,813 | 6,545 | 1,277 | 1,164 | 1,122 | 181 | 857 | 945 | 196 | 107 | 81 | 84.8 | 92.8 | 10.9 | 0.266 |

x ~ 0.25:

| N | rows1 | cols1 | (1,1) rows | (1,1) cols | (1,2)/(2,1) short | long | (2,2) child: parent o | parent b | after filter o | after filter b | searched (2,2) lists | (2,3) parent lists | (2,3) after o | (2,3) after b |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,049 | 675 | 316 | 232 | 220 | 75.2 | 146 | 162 | 82.8 | 40.4 | 33.9 | 31.3 | 36.5 | 9.92 | 1.65 |
| 3,164 | 1,097 | 447 | 354 | 327 | 101 | 234 | 256 | 112 | 55.9 | 49.6 | 42.6 | 49.8 | 11.7 | 1.81 |
| 4,803 | 1,789 | 583 | 490 | 463 | 117 | 342 | 373 | 128 | 68 | 59.3 | 51.8 | 59.5 | 11.4 | 1.12 |
| 8,618 | 3,366 | 884 | 793 | 740 | 159 | 572 | 618 | 172 | 93.9 | 85.8 | 70.7 | 80.5 | 12.7 | 0.963 |
| 15,199 | 6,102 | 1,350 | 1,280 | 1,194 | 215 | 964 | 1,038 | 232 | 132 | 124 | 102 | 116 | 14.8 | 0.968 |
| 16,424 | 6,508 | 1,362 | 1,339 | 1,224 | 215 | 980 | 1,045 | 229 | 129 | 120 | 101 | 116 | 14.2 | 0.798 |
| 20,538 | 8,768 | 1,629 | 1,466 | 1,365 | 223 | 1,151 | 1,256 | 241 | 136 | 119 | 106 | 119 | 13.1 | 0.422 |

How the lists scale:

| list | exponent |
|---|---|
| rows1 | ~N^1.1 |
| cols1 | ~N^0.7 |
| (1,1) lists | ~N^0.8 |
| long list at (1,2)/(2,1) (scanned by every (2,2) child) | ~N^0.85-0.9 |
| short list at (1,2)/(2,1) | ~N^0.45, saturating at ~180-240 |
| (2,2) child's lists after filters | ~N^0.53 |
| searched (2,2) node's lists | ~N^0.5 |

* After the third row or col, the other axis collapses to 10-15 entries. The same axis then
  has fewer than 2 entries left; 3 are needed.

## 6. How the (2,2) children die, and why the survivors multiply

Share of the (2,2) children dying at each stage, and of the (2,3)/(3,2) children:
* `cnt` means too few candidates; `fc` means an unmatched cell lost all of its candidates.
* `after` is the "cross support after the support filter" mode on its sampled children.

P = 13 7 4 3 1 1:

| N | (2,2) cnt_o | (2,2) fc_o | (2,2) cnt_b | (2,2) fc_b | (2,2) cross | (2,2) support | (2,2) after | (2,2) live | (2,3) fc_o | (2,3) cnt_b | (2,3) fc_b | (2,3) support | (2,3) live |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,523 | 0.0% | 0.9% | 0.4% | 3.3% | 50.4% | 37.1% | 4.4% | 3.6% | 39.12% | 46.80% | 9.90% | 2.94% | 0.02% |
| 4,111 | 0.0% | 0.2% | 0.1% | 1.0% | 36.6% | 43.0% | 4.1% | 14.9% | 29.28% | 54.85% | 11.73% | 3.70% | 0.01% |
| 5,896 | 0.0% | 0.1% | 0.1% | 0.6% | 24.6% | 45.8% | 3.9% | 24.8% | 23.57% | 61.07% | 11.61% | 3.53% | 0.00% |
| 7,593 | 0.0% | 0.1% | 0.1% | 0.4% | 20.8% | 45.6% | 3.7% | 29.5% | 21.56% | 63.92% | 11.07% | 3.28% | 0.00% |
| 11,306 | 0.0% | 0.0% | 0.0% | 0.2% | 14.8% | 45.2% | 3.4% | 36.3% | 18.27% | 69.83% | 9.32% | 2.48% | 0.00% |
| 15,199 | 0.0% | 0.0% | 0.0% | 0.1% | 12.6% | 46.7% | 3.3% | 37.2% | 17.77% | 72.87% | 7.48% | 1.80% | 0.00% |
| 18,957 | 0.0% | 0.0% | 0.0% | 0.1% | 9.5% | 48.1% | 3.3% | 39.0% | 16.24% | 75.38% | 6.76% | 1.55% | 0.00% |
| 22,992 | 0.0% | 0.0% | 0.0% | 0.1% | 9.8% | 52.5% | 3.9% | 33.6% | 17.93% | 77.00% | 4.23% | 0.77% | 0.00% |

P = 12 6 3 2 1 1:

| N | (2,2) cnt_o | (2,2) fc_o | (2,2) cnt_b | (2,2) fc_b | (2,2) cross | (2,2) support | (2,2) after | (2,2) live | (2,3) fc_o | (2,3) cnt_b | (2,3) fc_b | (2,3) support | (2,3) live |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,812 | 0.0% | 0.6% | 0.3% | 2.1% | 48.1% | 36.6% | 4.1% | 8.3% | 35.27% | 48.42% | 11.70% | 3.72% | 0.02% |
| 4,803 | 0.0% | 0.2% | 0.2% | 1.2% | 34.7% | 49.5% | 4.4% | 9.9% | 30.89% | 59.35% | 7.52% | 1.78% | 0.00% |
| 6,671 | 0.0% | 0.1% | 0.1% | 0.6% | 28.3% | 51.2% | 4.1% | 15.5% | 26.95% | 64.09% | 7.07% | 1.61% | 0.00% |
| 8,891 | 0.0% | 0.1% | 0.1% | 0.5% | 28.9% | 54.8% | 4.1% | 11.5% | 26.50% | 67.69% | 4.69% | 0.88% | 0.00% |
| 11,697 | 0.0% | 0.0% | 0.1% | 0.4% | 23.9% | 57.5% | 4.2% | 13.9% | 24.28% | 70.87% | 3.97% | 0.70% | 0.00% |
| 13,565 | 0.0% | 0.1% | 0.1% | 0.6% | 26.0% | 61.9% | 4.2% | 7.1% | 26.60% | 70.63% | 2.25% | 0.33% | 0.00% |
| 15,074 | 0.0% | 0.1% | 0.2% | 0.9% | 29.0% | 64.3% | 4.1% | 1.4% | 30.87% | 67.78% | 1.00% | 0.11% | 0.00% |
| 15,813 | 0.0% | 0.1% | 0.3% | 1.5% | 25.1% | 69.0% | 3.7% | 0.3% | 34.90% | 64.17% | 0.55% | 0.06% | 0.00% |

x ~ 0.25:

| N | (2,2) cnt_o | (2,2) fc_o | (2,2) cnt_b | (2,2) fc_b | (2,2) cross | (2,2) support | (2,2) after | (2,2) live | (2,3) fc_o | (2,3) cnt_b | (2,3) fc_b | (2,3) support | (2,3) live |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2,049 | 0.0% | 1.3% | 0.8% | 4.6% | 52.1% | 33.4% | 4.1% | 3.7% | 39.95% | 44.83% | 10.56% | 3.31% | 0.03% |
| 3,164 | 0.0% | 0.5% | 0.2% | 1.8% | 42.4% | 39.4% | 4.1% | 11.5% | 31.89% | 49.46% | 13.38% | 4.67% | 0.01% |
| 4,803 | 0.0% | 0.2% | 0.2% | 1.2% | 34.7% | 49.5% | 4.4% | 9.9% | 30.89% | 59.35% | 7.52% | 1.78% | 0.00% |
| 8,618 | 0.0% | 0.0% | 0.0% | 0.3% | 21.3% | 52.4% | 4.1% | 21.8% | 24.38% | 67.18% | 6.73% | 1.51% | 0.00% |
| 15,199 | 0.0% | 0.0% | 0.0% | 0.1% | 12.6% | 46.7% | 3.3% | 37.2% | 17.77% | 72.87% | 7.48% | 1.80% | 0.00% |
| 16,424 | 0.0% | 0.0% | 0.0% | 0.1% | 15.1% | 52.8% | 3.7% | 28.3% | 19.19% | 74.00% | 5.56% | 1.15% | 0.00% |
| 20,538 | 0.0% | 0.0% | 0.0% | 0.4% | 19.1% | 65.0% | 4.1% | 11.3% | 23.82% | 74.06% | 1.76% | 0.24% | 0.00% |

* **Cross support is what loses its grip.** It kills 50% of the (2,2) children at N = 2.5k, 21%
  at 7.6k and 10% at 23k.
* **The support filter kills 37-53% on ladder A** (up to 69% at high x).
* The survivors of both, searched, have lists of ~100-120 per axis.

**Fixpoint and productivity of the searched (2,2) nodes** (`-DSCALE_FIX`):

| sum | N | searched (2,2) nodes (est. per sum) | killed by cross + support iterated to their fixpoint | searched (2,2) nodes with a square below | searched (1,2)+(2,1) nodes with a square below |
|---|---:|---:|---:|---:|---:|
| 13 7 4 3 1 1, S = 1900 | 4,111 | 6.23e+06 | 95.3% | 35 of 3,114,003 (1.1e-05) | 35 of 1,701,579 (2.1e-05) |
| S = 1950 | 5,896 | 3.59e+07 | 90.8% | 17 of 2,248,247 (7.6e-06) | 17 of 598,918 (2.8e-05) |
| S = 2000 | 7,593 | 9.55e+07 | 89.0% | 28 of 4,779,815 (5.9e-06) | 28 of 937,678 (3.0e-05) |
| S = 2100 | 11,306 | 4.3e+08 | 86.0% | 11 of 4,295,497 (2.6e-06) | 11 of 526,502 (2.1e-05) |
| S = 2200 | 15,199 | 1.01e+09 | 90.5% | 7 of 2,536,963 (2.8e-06) | 7 of 286,257 (2.4e-05) |
| S = 2200 (2nd sample) | 15,199 | 1.3e+09 | 90.6% | 3 of 1,622,014 (1.8e-06) | 3 of 167,549 (1.8e-05) |
| S = 2400 | 22,992 | 4.01e+09 | 91.5% | 3 of 2,613,666 (1.1e-06) | 3 of 237,637 (1.3e-05) |
| 12 6 3 2 1 1, S = 1480 | 15,074 | 2.66e+07 | 99.7% | 0 of 68,886 (0.0e+00) | 0 of 260,814 (0.0e+00) |

* At every N, ~90% of the searched (2,2) nodes are dead under cross support and the support
  filter iterated to their common fixpoint. This is sound: no productive node is killed.
* It tests only the unions, so it is cheap relative to the subtree it would save:
  * a survivor's subtree is ~k22 x 570 + 500 ~ 11k cycles at N = 23k;
  * one round of cross support (both passes) is ~2,000 cycles, and the fixpoint needs one or
    more extra rounds plus support passes.
* research/ideas.md ("The support / cross cascade") measured this as neutral at N ~ 2-3k. There,
  everything below the searched (2,2) nodes was ~3% of the time. Here, the survivors take
  29-34% at N = 15-23k.

**What kills the (2,3)/(3,2) children** (`-DSCALE_FIX`):

| sum | N | (2,3)/(3,2) children reaching the b filter | same-axis candidates disjoint from the new vector | kept with the support test (3 needed) |
|---|---:|---:|---:|---:|
| 13 7 4 3 1 1, S = 1950 | 5,896 | 15,933,085 | 42.2 (short of 3: 0.00%) | 1.45 |
| 13 7 4 3 1 1, S = 2200 | 15,199 | 19,101,439 | 65.1 (short of 3: 0.00%) | 0.88 |
| 12 6 3 2 1 1, S = 1480 | 15,074 | 638,837 | 57.9 (short of 3: 0.00%) | 0.34 |

* The binding constraint is not pairwise disjointness within an axis.
* It is the exactly-once relation across the axes.
  * After a third col is placed, only ~12-15 row candidates meet it exactly once (and the same
    holds for a third row).
  * The remaining col candidates must lie within the placed rows and the union of those few
    rows. On average fewer than 1.5 do.
* The searched (2,2) nodes have long lists on both axes, but the bipartite "meets exactly once"
  relation between them is sparse.
  * Each new col meets ~11-13% of the row candidates exactly once.
  * The chance that a 4 x 4 completion (rows pairwise disjoint, cols pairwise disjoint, all 16
    pairs exactly once) exists is ~1e-6, as measured.
  * The union-based filters cannot see this.

## 7. Per-r1 subproblems

By position decile of r1. The vectors are sorted by label, so early r1 have the largest
universes U.
* n11, n12, c22 and s22 are per r1: searched (1,1) nodes, searched (1,2)+(2,1) nodes, (2,2)
  children made, and (2,2) nodes searched.
* "cyc/c22" is the total cycles of the r1 per (2,2) child made.

P = 13 7 4 3 1 1, S = 1900 (N = 4.1k, all r1):

```
== a1900 N 4111: per-r1 by position decile (r1 index / N)
 decile  U/N   alive  rows1  cols1  kids1  n11/r1   n12/r1   c22/r1   s22/r1  cyc/r1(M) time%  sq   sq/r1   cyc/c22
   0    0.95   0.87    2782    979   10.2      8.8      729    11650     5294     37.06  20.3   17  0.0413    3181
   1    0.85   0.96    2488    880   18.5     17.8     1276    19965     6194     50.42  27.5   25  0.0608    2526
   2    0.75   0.96    2227    747   24.9     24.0     1485    21024     2656     38.99  21.3   13  0.0316    1854
   3    0.65   0.97    1912    660   26.5     25.8     1458    18708      861     26.82  14.6   10  0.0243    1434
   4    0.55   0.97    1619    554   25.0     24.4     1204    13237      128     15.33   8.4    6  0.0146    1158
   5    0.45   0.97    1298    472   23.7     22.9      978     9412       14      9.12   5.0    6  0.0146     969
   6    0.35   0.96     986    384   19.8     19.0      665     5226        1      3.99   2.2    0  0.0000     763
   7    0.25   0.95     681    288   15.9     15.1      378     2220        0      1.18   0.6    0  0.0000     532
   8    0.15   0.94     385    186   10.5      9.9      135      486        0      0.22   0.1    0  0.0000     460
   9    0.05   0.72     135     79    4.7      3.2        5        8        0      0.01   0.0    0  0.0000    1403
  ln cyc ~ rows1: coef 2.77 sd 0.83 R2 0.90
  ln cyc ~ cols1: coef 3.38 sd 0.93 R2 0.88
  ln cyc ~ U: coef 3.06 sd 0.83 R2 0.90
  ln cyc ~ rows1+cols1+kids1: coef 0.89 1.96 0.82 sd 0.32 R2 0.99
```

P = 13 7 4 3 1 1, S = 2200 (N = 15.2k, 76 r1):

```
== a2200 N 15199: per-r1 by position decile (r1 index / N)
 decile  U/N   alive  rows1  cols1  kids1  n11/r1   n12/r1   c22/r1   s22/r1  cyc/r1(M) time%  sq   sq/r1   cyc/c22
   0    0.95   1.00   11689   2485   35.1     35.1     8315   245242   170743   3321.73  24.7    1  0.1250   13545
   1    0.85   1.00   10530   2179   73.1     73.1    15627   454050   305213   5232.36  34.1    4  0.5714   11524
   2    0.75   1.00    9253   1945   64.4     64.4    11304   322271   133142   2455.16  18.3    2  0.2500    7618
   3    0.65   1.00    8128   1607   67.1     67.1    12635   331573    75034   1722.70  11.2    1  0.1429    5196
   4    0.55   1.00    6797   1452   55.1     55.1     8476   192985    19174    729.28   5.4    0  0.0000    3779
   5    0.45   1.00    5375   1275   72.9     72.9     9614   193478     5642    575.63   4.3    1  0.1250    2975
   6    0.35   1.00    4120   1048   51.4     51.4     5675    92253      227    213.86   1.4    0  0.0000    2318
   7    0.25   1.00    2891    795   39.2     39.2     3216    41985        3     71.20   0.5    1  0.1250    1696
   8    0.15   1.00    1706    500   27.4     27.4     1485    12844        0     11.24   0.1    0  0.0000     875
   9    0.05   1.00     543    202   10.5     10.2      250     1286        0      0.78   0.0    0  0.0000     607
  ln cyc ~ rows1: coef 2.99 sd 0.92 R2 0.93
  ln cyc ~ cols1: coef 3.82 sd 1.00 R2 0.91
  ln cyc ~ U: coef 3.24 sd 0.90 R2 0.93
  ln cyc ~ rows1+cols1+kids1: coef 0.34 2.86 0.83 sd 0.41 R2 0.99
```

* The first half of the r1 (U > N/2) hold ~92-94% of the time and nearly all the squares.
* The (2,2) survival s22/c22 in the top deciles is 45% at N = 4.1k and 70% at 15k. That is why
  the cost per (2,2) child there is 3-13k cycles, against 0.5-2k for the late r1.

Per-r1 regressions (all alive r1 of each sum; ln of TSC cycles):

| sum | N | top 1% / 10% / 50% of r1: share of time | r1 with a square (sampled) | their share of time | ln cyc ~ cols1, kids1 | R^2 |
|---|---:|---|---|---:|---|---:|
| a1850 | 2,523 | 0.06 / 0.40 / 0.96 | 9 / 2,523 | 0.01 | cols1^2.82 kids1^0.91 | 0.98 |
| a1900 | 4,111 | 0.05 / 0.39 / 0.96 | 76 / 4,111 | 0.05 | cols1^3.02 kids1^0.87 | 0.98 |
| a1950 | 5,896 | 0.06 / 0.41 / 0.96 | 78 / 1,474 | 0.14 | cols1^3.15 kids1^0.87 | 0.99 |
| a2000 | 7,593 | 0.06 / 0.41 / 0.96 | 50 / 759 | 0.16 | cols1^3.19 kids1^0.83 | 0.98 |
| a2100 | 11,306 | 0.05 / 0.42 / 0.96 | 22 / 226 | 0.25 | cols1^3.32 kids1^0.86 | 0.99 |
| a2200 | 15,199 | 0.07 / 0.44 / 0.96 | 9 / 76 | 0.32 | cols1^3.28 kids1^0.84 | 0.99 |
| a2300 | 18,957 | 0.10 / 0.38 / 0.96 | 5 / 47 | 0.18 | cols1^3.57 kids1^0.91 | 0.99 |
| a2400 | 22,992 | 0.22 / 0.50 / 0.95 | 5 / 29 | 0.52 | cols1^3.33 kids1^0.77 | 1.00 |
| b920 | 4,803 | 0.06 / 0.39 / 0.96 | 101 / 4,803 | 0.06 | cols1^3.05 kids1^0.89 | 0.99 |
| b1200 | 11,697 | 0.07 / 0.40 / 0.95 | 5 / 183 | 0.07 | cols1^3.12 kids1^0.82 | 0.99 |
| b1600 | 15,813 | 0.08 / 0.38 / 0.95 | 0 / 66 | 0.00 | cols1^2.85 kids1^0.91 | 1.00 |

Pooled fits:

| fit | ladder A (13 7 4 3 1 1) | ladder B (12 6 3 2 1 1, x <= 0.53) |
|---|---|---|
| ln cyc per r1 | 3.13 ln U (R^2 0.73) | 3.11 ln U (R^2 0.71) |
| ln cyc per r1, with N | 3.03 ln cols1 + 0.88 ln kids1 + 0.42 ln N (sd 0.37, R^2 0.98) | 3.03 ln cols1 + 0.88 ln kids1 + 0.33 ln N (sd 0.36, R^2 0.99) |
| (2,2) children per r1 | ~U^2.42 (N/U)^-0.08: a function of U only | ~U^2.38 (N/U)^-0.17 |
| squares per r1 (Poisson) | ~U^1.81 +- 0.09; ~cols1^2.08 kids1^0.51; with N: cols1^3.29 kids1^0.79 N^-1.42 +- 0.25 | ~U^1.81 +- 0.13; ~cols1^2.20 kids1^0.56; with N: cols1^4.13 kids1^0.96 N^-2.68 +- 0.36 |

* Within a sum, an r1 with universe U costs ~U^3.0-3.3 and yields ~U^1.8-2.0 squares. So its cost
  per square grows ~U^1.1-1.5.
* Summed over the r1, with U uniform in 0..N, this gives time ~N^4.0-4.3, close to what is
  observed.

**(2,2) survival as a function of the r1's col list** (`s22 / c22` pooled over the r1 in each
bin of cols1, with the number of r1 in parentheses). It rises steeply with the list size within
a sum, and the threshold moves up with N:

```
cols1        a1900        a1950        a2000        a2100        a2200        a2300        a2400
 400-700     2.2% (1321)  0.3% ( 355)  0.1% ( 150)  0.0% (  31)  0.0% (   8)  0.0% (   4)  0.0% (   3)
 700-1000   25.5% ( 974) 10.2% ( 341)  3.9% ( 157)  0.3% (  36)  0.0% (   9)  0.0% (   6)  0.0% (   2)
1000-1400   58.3% ( 169) 52.0% ( 317) 35.1% ( 180)  9.6% (  47)  2.5% (  15)  0.3% (   6)  0.0% (   5)
1400-2000                85.0% (  25) 74.1% (  93) 53.7% (  60) 27.8% (  20) 11.2% (  12)  3.3% (   4)
2000-3000                                          86.6% (  15) 70.7% (  16) 60.1% (  13) 29.4% (   8)
3000-5000                                                                    85.8% (   2) 70.0% (   5)
```

## 8. Implementation steps on top of the algorithm

| comparison (same r1 sample) | nodes | time |
|---|---:|---:|
| 13 7 4 3 1 1, S = 2100 (192 labels): W = 3 -> forced W = 4 (`--min-words 4`) | same (71.1M) | x1.38 (745 -> 1,037 cycles per node) |
| 12 6 3 2 1 1, S = 840 (114 labels): W = 2 -> forced W = 3 | same (9.4M) | x1.38 (645 -> 875) |
| 13 7 4 3 1 1, S = 2400 (245 labels, N = 23k): carried bitsets -> matrix path (`--min-words 8`, as for > 256 labels: N x N matrices, no cross support, 8-slice MRV) | x2.79 | x2.69 (10 r1: 30.3 s -> 81.4 s; x1.5-3.4 per r1) |

* Ladder A crosses 128 labels at N ~ 3k and 192 at N ~ 13k. Those two steps (x1.9) are the
  N^0.3 drift of the cycles per node.
* The matrix path costs x2.7 at N = 23k (x1.5-3.4 per r1). Whether the factor grows with N
  was not measured.
* **Many of the sums the existence study cares about have more than 256 labels**: N >= 25k for
  13 7 4 3 1 1.

## 9. Synthesis: what grows faster than N^3, and why

On 13 7 4 3 1 1, S = 1850-2400, the time is roughly
`T ~ A x [k12 x c22 + k12 x p22 x k22 x c23] x w(W)`, with:

| quantity | measured growth |
|---|---|
| A, searched (1,2)/(2,1) nodes | N^2.7, about proportional to the squares (productive fraction ~2e-5 at every N) |
| k12, (2,2) children per A node | N^0.53 |
| c22, cycles per (2,2) child | N^0.82; it scans a parent list ~N^0.88 and cross-tests lists ~N^0.53 |
| p22, survival | 0.04 -> 0.37 (P- and x-dependent): ~N^1.5 from 2.5k to 11k, then flat |
| k22 | N^0.6 |
| c23 | N^0.6 |
| w(W) | x1.38 per word: N^0.3 on this ladder |

* **Super-linear relative to the output:**
  * per-child filter costs that scale with the parent's list instead of the child's (N^0.8);
  * a (2,2) layer whose union-based filters stop discriminating as the lists grow, so that
    dead survivors and their children grow like N^4.1-4.7 overall (N^3.0-3.4 over the top 4
    sums);
  * more (2,2) children per square (N^0.7-1.4 over the three ladders);
  * bitset width.
* **Not super-linear:** the upper tree (r1, (1,1), (1,2)/(2,1)), which is 3-7% of the time. A
  tracks the output.
* **Cost per square:**
  * On the carried path it rises from 0.2-0.5 CPU-s at N = 2-5k to 4-16 CPU-s at N = 15-23k
    (ladder A; 11-28 on the other P at N = 16-20k, where the counts are poor).
  * That is ~N^1.8 on ladder A, and N^1.85 on the fixed-x ladder.
  * The matrix path above 256 labels multiplies it by another ~2.7.

## Files

Scratch: `/tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/cx/profile/`
* `runs/`: per sum, `<tag>.log` (per-r1 lines), `<tag>.json` (layer counters) and `.inst`.
  * `a*`: ladder A.
  * `b*`: ladder B.
  * `c*`: fixed x.
  * `p*`: rdtsc phase builds.
  * `f*` and `d*`: fixpoint and disjointness diagnostics.
  * `m2400`: matrix path.
  * `w*` and `v*`: bitset width.
* `an.py`: loading, layer tables, per-r1 tables and regressions, fixpoint, phases and per-node
  costs.
* `report.py`: the section 1-6 tables. `report_out*.txt` and `tables2.txt` are its output.
* `run1.sh`-`run5.sh`: the exact runs (strides, offsets).
* `bin/`: the builds.
  * `bench_sc` and `bench_sc2`: `-DSCALE_PROF`.
  * `bench_scp2`: with `-DCHILD_PROF`.
  * `bench_fx2`: with `-DSCALE_FIX`.
  * `bench_plain2`: default flags.
* Per-r1 log line: `r1 label U alive rows1 cols1 kids1 squares nodes tsc_cycles cpu_ns n11
  n12+n21 c22 s22`.
