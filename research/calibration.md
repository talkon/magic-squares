# Calibration ladder: the rare events of the squares found, predicted vs observed (October 2026)

Every semi-magic square found so far was checked for the rarer events on the
way to a magic square:

* traversals (possible diagonals) with the magic sum S, the magic product P,
  or both (SP);
* pairs of disjoint diagonals of each type, from S+S up to SP+SP. An SP+SP
  pair is a magic square;
* **sub-events**: only some of a diagonal's coordinates (its sum and the
  exponents of the primes of P) pinned at their targets, on one diagonal or
  on both diagonals of a pair. They are 10-10^5 times commoner than the rungs
  and test the predictors below the magic rung.

For each square the exact counts were compared with what three predictors
say the square "should count for", rung by rung up the ladder. The aim is to
test the rungs we can observe before trusting the one we cannot (SP+SP).

* Tool: `scripts/calibrate.py`. Per-square heuristic:
  `scripts/square_heuristic.py`. Test: `scripts/test_calibrate.py`.
* The data (msearch output of October 2026) are not in the repo. The tables
  below come from one run over all of it (see "Reproduce").
* This is the second version, revised after two reviews; the changes are
  listed at the end.

## Short answer

**Data.** 7,021 distinct 6x6 semi-magic squares (deduplicated by hash) on
113 values of P. The counts recomputed from the grids match every record's
s_count, p_count, sp_count and best_score.

* By source: 4,920 from the 36 hand-picked seed P, searched from S_min; 901
  from the 40 units the scheduler chose on fresh P (sched40); 1,200 from the
  large-N, many-divisor and existence runs. By P: 5,594 squares are on the
  36 seed P.
* **Two thirds of the squares are not new evidence for the SP rung.** 4,741
  squares (62 P) lie in sums the first search (758,949 squares, per-P totals
  only) had already covered, on P chosen from that search: 34 of the 36 seed
  P were searched by it, and 4,527 of the 5,594 seed-P squares lie in its
  range. At least one is literally the same square (7aa6fdb4ec9759c3, line
  32 of its top-100 list). The other 2,280 squares (55 P) are new.

**Do the predictions line up rung by rung?** (Full table below.)

| rung | observed | regression E, O/E [90%] | own-count null E, O/E | heuristic E, O/E [90%] |
|---|---:|---|---|---|
| S traversals | 8,670 | 7,530: 1.15 [1.13, 1.17] | (= obs) | 8,704: 1.00 [0.98, 1.01] |
| P traversals | 2,106 | 1,978: 1.06 [1.03, 1.10] | (= obs) | 1,874: 1.12 [1.08, 1.16] |
| S+S pairs | 147 | 90.7: 1.62 [1.41, 1.86] | 123: 1.19 [1.04, 1.37] | 124: 1.18 [1.03, 1.35] |
| S+P pairs | 48 | 42.9: 1.12 [0.87, 1.42] | 51.6: 0.93 [0.72, 1.18] | 46.3: 1.04 [0.80, 1.32] |
| P+P pairs | 8 | 7.48: 1.07 [0.53, 1.93] | 8.70: 0.92 [0.46, 1.66] | 8.38: 0.95 [0.48, 1.72] (in-sample) |
| SP traversals | 14 | 6.17: 2.27 [1.37, 3.55] | (= obs) | 8.76: 1.60 [0.97, 2.50] |
| ... inside the first search's range | 11 | 4.63: 2.37 [1.33, 3.93] | | 6.75: 1.63 [0.91, 2.70] |
| ... new squares | 3 | 1.54: 1.95 [0.53, 5.05] | | 2.01: 1.49 [0.41, 3.86] |
| squares with best pair >= SP+S | 2 | 0.186 | 0.534 (p = 0.10) | 0.299 |
| SP+SP (magic) | 0 | 7.2e-5 | 0 | 2.5e-4 (corrected: 1.8e-4) |

* **The heuristic** (from each square's 36 entries; its Edgeworth
  truncation was chosen on the P+P rung and one other event, so P+P is
  in-sample for it):
  * S is right (1.00). P is 12% low, a 5-sigma deficit (p = 8e-8) that
    varies between bins (0.76-1.35).
  * S+P and P+P line up. S+S is 1.18 [1.03, 1.35], an excess the own-count
    null shows too and nothing explains.
  * SP is 1.60 [0.97, 2.50] (p = 0.06). The rungs above it (SP+0 = 15 x 14
    - 2 pairs, SP+S, SP+P, the top rung) are the same 14 SP traversals and
    are not separate confirmations. Given the SP count the top rung is 2 vs
    0.48-0.53 squares (p = 0.09-0.10).
  * On the new squares SP is 3 vs 2.01 (p = 0.33): no information.
* **The regression** (the scheduler's model) is low on every rung, and more
  so higher up: S 1.15, S+S 1.62, SP 2.27, top rung 2 vs 0.19. It is
  in-sample on the seed squares.

**What the sub-events add** (new; 10-10^5 times more events than the rungs):

* **The pair factor is observable below the magic rung.** With d exponents
  pinned on both diagonals (12 million down to 65 pairs for d = 1..5), the
  data's pair factor given the counts stays at 1.00 while the heuristic's
  grows to 1.38. The ratio falls: 0.99, 0.95, 0.92, 0.86 [0.82, 0.89], 0.72
  [0.57, 0.87].
* **Its congruence part is right, its correlation part is not.** Split by
  the exact mod-3 class of each pair type: forbidden pair types have exactly
  0 pairs at every d, and boosted ones show the predicted ~3x (2.5-3.1).
  The shortfall is the correlation between the two diagonals: in pair types
  without extra congruences the data show 1.00, 0.97, 0.95, 0.88 where the
  heuristic has 1.02-1.07. Partner diagonals are slightly anti-clustered.
* **The default Edgeworth variant is supported out of sample**: its
  single-diagonal deficit (up to 12%) cancels that, and two-diagonal totals
  come out at 0.95-1.02 up to 10 pinned coordinates. The "mixed" variant is
  1.25-1.6x too high there, the others much worse.
* **The heuristic's sum-product coupling falls short, but only on the seed
  P inside the first search's range.** Pinning the sum with j exponents,
  over the j exponents alone: 1.00, 1.01, 1.03, 1.12, 1.35 for j = 0..4, then
  1.60 at SP. On re-found seed-P squares it reaches 1.57 [1.18, 1.99] at
  j = 4. On the new squares it is 0.95 [0.78, 1.12] at j = 3 (149 events, 50
  times their SP count) and 0.94 at j = 4; on sched40 0.99 and 0.90.

**Is there a trend with rarity?**

* Against the heuristic, no significant trend on the ladder. The SP excess
  is real in the sub-events but specific to the seed P near S_min, inside
  the first search's range: P selected on the first search, on partly the
  same squares.
* Against the regression, the climb is the traversal-class rate errors
  multiplied together, mainly rho near S_min, which the first search showed
  already.
* Pairs add no clumping beyond the class rates except S+S (about 1.2). The
  top rung pooled with the first search (new squares only) is 0.83-0.96,
  consistent with 1.

**P(magic | square).** E[magic] summed over the 7,021 squares, and over the
901 sched40 squares (the ones most like the forecast's):

| estimate | all: E[magic] | x regression | sched40: x regression |
| --- | ---: | ---: | ---: |
| regression (the scheduler's model) | 7.2e-5 | 1 | 1 |
| heuristic as it stands | 2.5e-4 | 3.45 | 2.07 |
| **heuristic, pair factor corrected by the sub-events** | **1.8e-4 (1.6-2.3e-4)** | **2.5 (2.2-3.2)** | **1.63 (1.5-2.1)** |
| ... x the pooled SP factor^2 (2.6 [1.0, 5.9]) | 4.7e-4 [1.8e-4, 1.1e-3] | 6.6 [2.6, 15] | 4.3 [1.7, 9.6] |

* Zero magic squares observed is consistent with every row.
* The pair factor at magic goes from 1.74 to about 1.26 (1.1-1.6).

**What it means for `research/forecast.md`.** Compared like with like on
sched40 (scheduler units on fresh P), the heuristic's P(magic | square) is
**about 0.75x the forecast's central assumption (0.6-0.9x)**.

* The 10 CPU-year central value of 0.03-0.045 becomes about 0.02-0.035.
* If the pooled SP excess applied to fresh P, it would be 1.5-1.8x
  [0.6, 4]. The new squares' sub-events argue against that, and the SP rung
  on fresh P (2 vs 1.29) cannot tell.

**What to measure next.** SP traversals and the sub-events on fresh P and
(P, S) outside the first search's range, near S_min (forecast.md's item 1).
The sum + 3-4 exponent sub-events give 10-50 times the SP events, so about
50,000 new squares would pin the coupling to a few percent and the SP rate
to 15%. Rerun `scripts/calibrate.py` on the new units.

## Data, events, predictors

**Events.** These are exactly those of `square_diag_stats` (src/c/square.c).

* Traversals: the 720 permutations. Each is S (sum S, including SP), P
  (product P, including SP) or SP (both). The exclusive classes S_only and
  P_only are also used.
* Pairs: the 5,400 unordered partner pairs {sigma, sigma o tau}, where tau
  is one of the 15 fixed-point-free involutions.
* A pair's type uses the exclusive classes, and its score is the sum of the
  traversal scores (S = 2, P = 3, SP = 7): S+S 4, S+P 5, P+P 6, SP+0 7,
  SP+S 9, SP+P 10, SP+SP 14 (magic).
* "Both sums S" (S+S + SP+S + SP+SP) and "both products P" are the
  inclusive versions.
* Each square's best_score gives the "squares with best pair >= ..." rows.
* **Sub-events.** A traversal's coordinates are its sum X and the exponents
  of the k primes of P. A sub-event pins d of them at their targets: the k
  cyclic windows of d consecutive primes, with or without X
  (`square_heuristic.sub_sets`), on one diagonal or on both diagonals of a
  partner pair. The full sets are the rungs.

**Splits of the data.**

* **First search: searched / new.** The first search covered every sum
  from S_min to max_S for 1,138 values of P (stats_short.txt; max_S = 99999
  means all sums). A square is "searched" if its (P, S) lies in that range.
  Such squares were within its reach and are partly the same squares; only
  its top-100 list survives to check identity.
* **P source.** The first source with a square of that P, so "seed P"
  means the 36 hand-picked P, whichever run found the square. The new
  squares on seed P come from only 6 of them (2 the first search never
  searched, 4 above its max_S).

**Predictors.**

* **regression**: the scheduler's model (`DEFAULT_MODEL`, which is
  `forecast/state/model_6.json`).
  * p_S and p_P per square come from its (P, S, S_min).
  * q_SP = rho p_S p_P with rho = 2.148.
  * Traversal classes are i.i.d. over the 720 traversals, so an X+Y pair
    count is 5400 (2) q_X q_Y. Magic is 5400 q_SP^2, as in the scheduler.
  * It is in-sample on the seed squares: `scheduler.fit_model` fitted p_S,
    p_P and rho with each seed sum as a row and with the first search's
    per-P totals, which contain most of those squares again. They are
    diluted only by the first search's squares far from S_min.
* **null** (own-count null): the square's own numbers of S_only, P_only and
  SP traversals, placed on a random subset of the 720. It reproduces the
  traversal rungs by construction, so it tests only the pair rungs "given
  the counts". For a sub-event it is 15/719 C(n1, 2) for the square's n1
  hitting traversals.
* **heuristic**: `square_heuristic.predict(grid, S, P)` (see "The
  heuristic"). It uses only the square's entries and never looks at which
  traversals hit S or P. No parameter is fitted to counts, but **its
  Edgeworth truncation was chosen** because it fits the P+P rung and the
  12-cell product event, so those two are in-sample for it.
  * Its variants (`heur-mixed`, `heur-marg`, `heur-gauss`, `heur-nolat`),
    computed from the same components, measure its model uncertainty.
  * Its conditional mode (`--heuristic-cond`) predicts each rung from the
    realized one below it. It gives SP = 9.8 from the realized S and P
    counts.
* `calibrate.py` also reports two in-sample predictors, perP and eb. They
  are used only to split the regression's misses into parts.

**Intervals.**

* "[90%]" is the exact (Garwood) Poisson interval of the observed count,
  divided by E. E is treated as exact.
* Sub-event intervals are bootstrap intervals over squares (the windows of
  a square are nested and its pairs cluster).
* On the ladder, a bootstrap over squares agrees with the Poisson intervals
  except for SP+0, where 15 pairs share one SP traversal. A bootstrap over
  the 113 P gives nearly the same intervals for S, P, SP and S+S, so
  clustering by P does not matter (review).
* p-values are one-sided Poisson tails.

## The ladder (all 7,021 squares)

| rung | observed | regression E | O/E [90%] | null E | O/E [90%] | heuristic E | O/E [90%] |
|---|---:|---:|---|---:|---|---:|---|
| S traversals | 8,670 | 7,530 | 1.15 [1.13, 1.17] | (= obs) | - | 8,704 | 1.00 [0.98, 1.01] |
| P traversals | 2,106 | 1,978 | 1.06 [1.03, 1.10] | (= obs) | - | 1,874 | 1.12 [1.08, 1.16] |
| S+S pairs | 147 | 90.7 | 1.62 [1.41, 1.86] | 123 | 1.19 [1.04, 1.37] | 124 | 1.18 [1.03, 1.35] |
| S+P pairs | 48 | 42.9 | 1.12 [0.87, 1.42] | 51.6 | 0.93 [0.72, 1.18] | 46.3 | 1.04 [0.80, 1.32] |
| P+P pairs | 8 | 7.48 | 1.07 [0.53, 1.93] | 8.70 | 0.92 [0.46, 1.66] | 8.38 | 0.95 [0.48, 1.72] |
| pairs with both sums S | 148 | 90.8 | 1.63 [1.42, 1.87] | 124 | 1.20 [1.04, 1.37] | 125 | 1.19 [1.03, 1.36] |
| pairs with both products P | 9 | 7.52 | 1.20 [0.62, 2.09] | 8.80 | 1.02 [0.53, 1.78] | 8.45 | 1.07 [0.56, 1.86] |
| SP traversals | 14 | 6.17 | 2.27 [1.37, 3.55] | (= obs) | - | 8.76 | 1.60 [0.97, 2.50] |
| SP+0 pairs | 208 | 92.3 | 2.25 [2.00, 2.53] | 210 | 0.99 [0.88, 1.11] | 131 | 1.59 [1.41, 1.78] |
| SP+S pairs | 1 | 0.144 | 6.9 [0.36, 33] | 0.438 | 2.3 [0.12, 11] | 0.231 | 4.3 [0.22, 21] |
| SP+P pairs | 1 | 0.0447 | 22 [1.2, 106] | 0.104 | 9.6 [0.49, 46] | 0.0723 | 14 [0.71, 66] |
| squares with best pair >= SP+S | 2 | 0.186 | 10.7 [1.9, 34] | 0.534 | 3.7 [0.67, 11.8] | 0.299 | 6.7 [1.2, 21] |
| SP+SP pairs (magic) | 0 | 7.2e-5 | - | 0 | - | 2.5e-4 | - |

The "best pair at least" rows (observed; regression / null / heuristic /
heuristic-cond):

| at least | observed | regression | null | heuristic | heuristic-cond |
| --- | ---: | ---: | ---: | ---: | ---: |
| S+S | 201 | 141 | 188 | 183 | 182 |
| S+P | 68 | 55.3 | 72.5 | 62.4 | 67.5 |
| P+P | 22 | 13.6 | 22.6 | 17.1 | 19.9 |
| SP | 14 | 6.16 | 14 | 8.75 | 9.70 |
| SP+S | 2 | 0.186 | 0.534 | 0.299 | 0.783 |
| SP+P | 1 | 0.045 | 0.103 | 0.072 | 0.346 |

**SP and the top rung by split** (O/E with the 90% interval; p = P(X >= O)):

| split | squares | SP obs | regression | heuristic | p (heur) | top rung obs | reg / null / heur E |
|---|---:|---:|---|---|---:|---:|---|
| all | 7,021 | 14 | 2.27 [1.37, 3.55] | 1.60 [0.97, 2.50] | 0.062 | 2 | 0.19 / 0.53 / 0.30 |
| inside the first search's range | 4,741 | 11 | 2.37 [1.33, 3.93] | 1.63 [0.91, 2.70] | 0.082 | 0 | 0.14 / 0.39 / 0.23 |
| new | 2,280 | 3 | 1.95 [0.53, 5.05] | 1.49 [0.41, 3.86] | 0.33 | 2 | 0.043 / 0.14 / 0.065 |
| seed P (by P) | 5,594 | 12 | 2.52 [1.45, 4.08] | 1.71 [0.99, 2.77] | 0.054 | 1 | 0.14 / 0.43 / 0.24 |
| other P | 1,427 | 2 | 1.42 [0.25, 4.47] | 1.15 [0.20, 3.62] | 0.52 | 1 | 0.044 / 0.10 / 0.063 |
| inside range, seed P | 4,527 | 11 | 2.51 [1.41, 4.15] | 1.71 [0.96, 2.83] | 0.063 | 0 | 0.13 / 0.39 / 0.22 |
| new, seed P (6 P) | 1,067 | 1 | 2.64 [0.14, 12.5] | 1.71 [0.09, 8.1] | 0.44 | 1 | 0.008 / 0.041 / 0.016 |
| new, other P | 1,213 | 2 | 1.73 [0.31, 5.44] | 1.41 [0.25, 4.43] | 0.42 | 1 | 0.035 / 0.10 / 0.049 |
| sched40 | 901 | 2 | 1.83 [0.33, 5.76] | 1.55 [0.28, 4.88] | 0.37 | 1 | 0.034 / 0.10 / 0.046 |

Notes on the ladder:

* The regression's best-pair probabilities are exact for i.i.d. classes on
  the partner graph, up to the Monte Carlo error of its tables. The null
  uses the same exact graph tables. The heuristic uses Poisson clumping; on
  the exact graph at the heuristic's own rates (pairs independent) it would
  give 176 / 60.9 / 15.6 / 8.75 / 0.293 / 0.058, within 10% of its own
  values except SP+P (0.058 vs 0.072).
* Tail probabilities:
  * SP: 14 or more is p = 0.0046 under the regression, 0.062 under the
    heuristic and 0.12 under the conditional heuristic.
  * The rungs above SP are not independent of it: SP+0 is 15 x 14 - 2
    pairs, and SP+S, SP+P and the top rung each sit on one of the 14 SP
    traversals. Unconditionally the top rung is p = 0.015 (regression) and
    0.037 (heuristic); given the SP count it is 2 vs 0.48-0.53 (class
    factors of either predictor, or the null), p = 0.09-0.10.
  * The two top-rung squares are both new squares (null E 0.14, p = 0.009
    on that post-hoc split; 0 vs 0.39 inside the first search's range).
    Neither half should be over-read.
  * P rung under the heuristic: p = 8e-8. Its deficit is systematic but
    small.
  * S+S: p = 0.02-0.03 under both the null and the heuristic.

## The trend with rarity: per-class rate factors

A pair rung X+Y should be off by f_X f_Y if the only errors are in the
traversal-class rates f (O/E for S_only, P_only and SP) and the pairs are
otherwise right. The residual is what is left:

| rung | regression O/E | f_X f_Y | residual | heuristic O/E | f_X f_Y | residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| class rates (S_only, P_only, SP) | 1.15, 1.06, 2.27 | | | 1.00, 1.12, 1.60 | | |
| S+S | 1.62 | 1.32 | 1.23 | 1.18 | 0.99 | 1.19 |
| S+P | 1.12 | 1.22 | 0.92 | 1.04 | 1.12 | 0.93 |
| P+P | 1.07 | 1.13 | 0.95 | 0.95 | 1.26 | 0.76 |
| SP+0 | 2.25 | 2.27 | 0.99 | 1.59 | 1.60 | 0.99 |
| SP+S | 6.9 | 2.61 | 2.7 (1 event) | 4.3 | 1.59 | 2.7 (1 event) |
| SP+P | 22 | 2.41 | 9.3 (1 event) | 14 | 1.79 | 7.7 (1 event) |
| squares best >= SP+S | 2 vs 0.186 | | 2 vs 0.48, p = 0.085 | 2 vs 0.299 | | 2 vs 0.50, p = 0.090 |

* For the regression, the climb up the ladder is the class rates. The same
  split with each P's own rates (perP) as the middle step gives:
  * the P's rates: S+S 1.32, SP 2.27;
  * spread of counts across squares of a P: 0.98-1.03 on the rungs with
    many events, so traversal counts are barely over-dispersed;
  * pairs given counts: as the null column above.
* **The top rung with the first search.** Its totals contain 10 such
  squares (9 SP+S, 1 SP+P, from the report). Pooled at the P level with
  only the 2,280 squares here outside its range (the others are already in
  its totals): 12 observed against 12.34 (each P's own i.i.d. rates) + 0.14
  (null) = 12.5, O/E 0.96 [0.55, 1.56].
  * On these squares, per-P i.i.d. rates give only 0.87x the null's SP+S +
    SP+P, because SP squares are S-richer. Scaling the first search's
    expectation by that gives 14.4 and O/E 0.83 [0.48, 1.35].
  * So 0.83-0.96: no extra pair clumping at the top rung.
* The report's estimates of 4.58 and 1.46 miss a factor of 2: an unordered
  pair can have its SP diagonal on either side.
* The only pair-level excess on the ladder is S+S, at 1.19 [1.04, 1.37]
  against both the null and the heuristic. The heuristic's own correlation
  between the two sums adds only 1.01. The sub-events show the same excess
  with one exponent pinned too (1.25 [0.94, 1.58]); it stays unexplained.

## Sub-events: the predictors below the rungs

`calibrate.py` pins d of a diagonal's coordinates and counts exactly, per
square, the traversals (one diagonal) and the partner pairs (both
diagonals) that hit them, against the heuristic and the own-count null.
Intervals are bootstraps over squares. The full tables, also by P source,
first-search range and source, are in the tool's `ladder.md`.

**Both diagonals, d exponents pinned on each (all squares):**

| d (coordinates) | pairs | O/E default [90%] | mixed | marg | gauss | O/null (data's pair factor) | heuristic's kappa | O/(null kappa) [90%] |
|---|---:|---|---:|---:|---:|---:|---:|---|
| 1 (2) | 11,966,313 | 0.997 | 0.997 | 0.979 | 0.842 | 1.000 | 1.014 | 0.986 [0.986, 0.987] |
| 2 (4) | 836,388 | 0.996 [0.992, 1.000] | 0.957 | 0.860 | 0.630 | 1.000 | 1.052 | 0.951 [0.95, 0.95] |
| 3 (6) | 47,420 | 1.023 [1.011, 1.036] | 0.907 | 0.686 | 0.435 | 1.009 | 1.094 | 0.922 [0.91, 0.93] |
| 4 (8) | 2,056 | 1.020 [0.98, 1.06] | 0.804 [0.77, 0.84] | 0.474 | 0.264 | 1.005 | 1.174 | 0.856 [0.82, 0.89] |
| 5 (10) | 65 | 0.949 [0.74, 1.17] | 0.638 [0.50, 0.78] | 0.274 | 0.136 | 0.989 | 1.375 | 0.719 [0.57, 0.87] |

* **The pair factor can be observed below the magic rung.** Given the
  counts, the data's pair factor stays at 1.00 for up to 5 exponents (10
  coordinates), while the heuristic's kappa grows to 1.38.
  * As O/null = kappa^beta, beta = 0.02 [-0.05, 0.09]: the data show none
    of the heuristic's pair factor *in total*.
  * The moment residual (O2/E2)/(O1/E1)^2 agrees: 0.98, 0.95, 0.93, 0.86,
    0.76.
  * The trend is the same in every group: at d = 4 it is 0.88 on seed P,
    0.79 on sched40, 0.73 on verify-recompute, 0.85 on existence, 0.86
    inside the first search's range and 0.85 outside.
* **The default Edgeworth variant is right in total out of sample**: 0.95-
  1.02 up to 10 coordinates. Its one-diagonal rates are 1-12% low and its
  pair factor too high, and the two cancel. "mixed" over-predicts these
  events by 1.25x at 8 and 1.6x at 10 coordinates; "marg" and "gauss" much
  more.

**The same pairs by the exact congruence class of their pair type tau:**

| d | forbidden: pairs (null) | no extra congruence: pairs, O/null [90%], kappa | extra congruences: pairs, O/null [90%], kappa |
|---|---|---|---|
| 1 | 0 (69,038) | 11,889,081, 1.002, 1.016 | 77,232, 2.54 [2.51, 2.57], 2.46 |
| 2 | 0 (27,750) | 761,493, 0.972 [0.970, 0.975], 1.022 | 74,895, 2.97 [2.95, 2.99], 3.13 |
| 3 | 0 (4,986) | 36,037, 0.945 [0.935, 0.955], 1.032 | 11,383, 2.94 [2.89, 2.99], 3.11 |
| 4 | 0 (600) | 937, 0.883 [0.83, 0.93], 1.071 | 1,119, 2.92 [2.77, 3.08], 3.30 |
| 5 | 0 (45) | 2, -, 3.9 | 63, 3.06 [2.44, 3.70], 4.33 |

* "Forbidden" means the target is outside the coset of the pair's lattice;
  "extra congruences" means the pair lattice is finer than the two
  diagonals' lattices, and the heuristic multiplies those pairs by the
  index (3, 9, ...).
* **The congruence part is confirmed**: no forbidden pair type ever has a
  pair (exact), and the boosted ones show about 3x, as predicted, up to
  d = 4.
* **The correlation part is not.** Without extra congruences, kappa is
  only the Gaussian-Edgeworth correlation of the two diagonals (1.02-1.07).
  The data show the opposite sign and a growing effect: 0.97, 0.95, 0.88.
  Partner diagonals are slightly anti-clustered, and the heuristic is
  0.95, 0.92, 0.83 too high there. In the boosted class it is too high by
  about the same factor (0.95, 0.94, 0.88, then 0.71 at d = 5).

**One diagonal, the sum X with j exponents pinned** (the coupling is O/E of
X with the j exponents over O/E of the same j exponents alone; j = k is the
SP rung):

| j | exponents alone: O, O/E | X + j exponents: O, O/E | coupling, all [90%] | inside the first search's range, seed P | new squares | sched40 |
|---|---|---|---|---|---|---|
| 0 | - | 8,670, 0.996 | 0.996 [0.98, 1.01] | | | |
| 1 | 5,883,743, 1.007 | 10,362, 1.014 | 1.007 [0.99, 1.03] | | | |
| 2 | 1,367,724, 1.025 | 2,473, 1.053 | 1.03 [0.99, 1.07] | 1.04 [0.99, 1.10] | 0.98 [0.91, 1.07] | 0.96 [0.84, 1.09] |
| 3 | 298,401, 1.051 | 598, 1.18 | 1.12 [1.03, 1.22] | 1.20 [1.07, 1.33] | 0.95 [0.78, 1.12] | 0.99 [0.73, 1.29] |
| 4 | 59,689, 1.087 | 150, 1.46 | 1.35 [1.06, 1.64] | 1.57 [1.18, 1.99] | 0.94 [0.48, 1.44] | 0.90 [0.22, 1.81] |
| k (SP) | | 14, 1.60 | | 11 vs 6.44: 1.71 | 3 vs 2.01: 1.49 | 2 vs 1.29: 1.55 |

* The heuristic's coupling between the sum and the exponents falls short
  more and more as exponents are pinned. That is the SP excess, and it is
  not 14-event noise: 598 and 150 events at j = 3, 4.
* **But it is confined to the re-found squares on seed P**, and within
  them to S near S_min (S/S_min < 1.2: 1.22 [1.08, 1.38] and 1.66 [1.24,
  2.15]; at 1.2 S_min and above 1.05 and 1.20).
* On the 2,280 new squares (and on sched40) it is consistent with 1, with
  149 events at j = 3: 0.95 [0.78, 1.12] excludes 1.2. The new squares on
  seed P (6 P) give 0.98 [0.74, 1.22], the other P 0.92 [0.70, 1.16].
* The seed P were chosen as the best of the first search, near S_min, so a
  coupling excess there, on partly the same squares, is what selection
  produces.

## The SP rung

* **Against the regression** SP is 2.27x. The heuristic explains part of
  it: its own prediction is 1.42x the regression's (S rate 1.16x, P rate
  0.95x, and an effective rho of 2.80 per square against the model's
  2.15). Conditioning on the realized S and P counts (`heuristic-cond`)
  gives E = 9.8, so O/E = 1.43 [0.87, 2.24] and P(f <= 1) = 0.10.
* **rho depends on P and on S/S_min**, from first principles. The heuristic's
  rho is 4.8 below 1.1 S_min, then 3.4, 2.9, 2.4 and 2.0 above 1.5; 3.0 on
  the seed P and 2.2 on sched40.
* **The first search does not replicate the excess independently.**
  * Its rho is 2.20 over 519 SP traversals, mostly far from S_min, and
    2.87 on the 34 seed P it searched (82 SP traversals in 91,700
    squares). The seed P were picked from it, so that is the same
    selection, not a replication.
  * With its rho of 2.20, 7.6 SP traversals would be expected here
    (p = 0.025). That tests mostly a subset of its own data (11 of the 14
    events are inside its range) against its average.
  * The excess is also not independent of forecast.md's rho of 2.6, which
    comes from the first search's squares near S_min. Applying both would
    count it twice.
* **Posterior of the SP factor f** (O/E against the heuristic, Jeffreys
  gamma, 90%):
  * pooled: 1.62 [1.01, 2.43], f^2 = 2.6 [1.0, 5.9];
  * inside the first search's range: 1.65 [0.97, 2.60];
  * new squares: 1.58 [0.54, 3.51], f^2 = 2.5 [0.3, 12]; sched40 1.69
    [0.44, 4.3]; other P 1.25 [0.33, 3.2]: no information;
  * against the regression: 2.30 [1.44, 3.45], f^2 = 5.3 [2.1, 11.9].
* **Conclusion.** The SP excess over the heuristic is a real deficit in its
  sum-product coupling, seen clearly only on the re-found seed-P squares
  near S_min. On new squares the sub-events, with 10-50 times the SP events,
  show no such deficit.

## P(magic | square)

How the heuristic's E[magic] relates to the regression's: ratio = rates x
rho x pairs, where rates are the heuristic's S and P rates with the model's
rho; rho is the heuristic's S-P correlation; pairs is its SP+SP over
5400 q_SP^2, split into the mod-3 congruences and the correlation between
the diagonals.

| squares | ratio | rates | rho | pairs (congruences x correlation) | heuristic rho |
| --- | ---: | ---: | ---: | ---: | ---: |
| all | 3.45 | 1.22 | 1.63 | 1.74 (1.35 x 1.29) | 2.80 |
| S/S_min [1, 1.1) | 25 | 1.93 | 4.77 | 2.75 (1.87 x 1.47) | 4.79 |
| [1.1, 1.2) | 6.4 | 1.20 | 2.58 | 2.07 (1.49 x 1.39) | 3.44 |
| [1.2, 1.3) | 2.8 | 1.00 | 1.84 | 1.49 (1.23 x 1.21) | 2.92 |
| [1.3, 1.5) | 2.2 | 1.15 | 1.23 | 1.52 (1.26 x 1.21) | 2.36 |
| [1.5, 2) | 2.0 | 1.65 | 0.88 | 1.36 (1.17 x 1.16) | 2.01 |
| 4 primes | 1.8 | 1.23 | 1.18 | 1.21 (1.11 x 1.09) | 2.30 |
| 5 primes | 3.2 | 1.16 | 1.72 | 1.63 (1.32 x 1.24) | 2.80 |
| 6 primes | 8.8 | 1.40 | 2.42 | 2.60 (1.69 x 1.53) | 3.28 |
| seed P | 3.9 | 1.14 | 1.90 | 1.79 (1.37 x 1.31) | 2.99 |
| sched40 | 2.1 | 1.37 | 1.06 | 1.43 (1.22 x 1.17) | 2.17 |
| new squares | 3.1 | 1.35 | 1.24 | 1.84 (1.41 x 1.30) | 2.49 |

* **Rates.** The S and P rungs test this factor. The heuristic is right on
  S and 12% low on P.
* **rho.** The SP rung and the sum + exponent sub-events test it. They ask
  for more on the re-found seed-P squares, not on new ones.
* **Pairs.** The sub-events test it below the magic rung. They confirm the
  congruence part and not the correlation part.

**The pair factor corrected by the sub-events.** The data's pair factor
over the heuristic's, r(d), is measured for d = 1..5 exponents on both
diagonals. A magic pair pins the k exponents of P (plus the sum), so each
square gets r(k):

* r(4) = 0.856 and r(5) = 0.719 are measured. r(6) = 0.66 and r(7) = 0.61
  are extrapolated from r(5) with the fitted slope (-0.08 per coordinate).
* Weighted by the heuristic's E[magic] (4, 5 and 6 primes: 0.18, 0.48,
  0.34), this is 0.73.
* Alternatives give the range:
  * keep the confirmed congruences and set the correlation part to 1
    (pairs given the counts independent otherwise): 0.78;
  * also apply the data's own anti-clustering (0.88 at d = 4, extrapolated):
    0.64;
  * the sum coordinate is not corrected; its S+S residual of 1.19 would
    raise any of these by up to 1.19.
* The pair factor at magic thus goes from 1.74 to about 1.26 (1.1-1.6).

**Model uncertainty: the Edgeworth variants** (E[magic] over all squares;
each "x f^2" uses that variant's own SP factor):

| variant | E[magic] | x regression | SP O/E | two-diagonal sub-events at 8 / 10 coordinates | x f^2 |
| --- | ---: | ---: | ---: | --- | ---: |
| pair (default) | 2.48e-4 | 3.45 | 1.60 | 1.02 / 0.95 | 6.5e-4 (9.1x) |
| mixed | 3.78e-4 | 5.27 | 1.30 | 0.80 / 0.64 | 6.6e-4 (9.2x) |
| marg | 2.12e-3 | 29.5 | 1.10 | 0.47 / 0.27 | 2.6e-3 (36x) |
| gauss | 4.42e-3 | 62 | 0.78 | 0.26 / 0.14 | 2.8e-3 (38x) |
| nolat (no pair congruences) | 1.80e-4 | 2.51 | 1.60 | (congruences confirmed) | 4.7e-4 (6.6x) |

* The variants differ by up to 18x at magic. The sub-events select the
  default: "mixed" (1.5x higher at magic) is 1.25-1.6x too high on
  two-diagonal events with 8-10 coordinates, and "marg" and "gauss" are
  rejected outright.
* Calibrated on the SP rung, "pair" and "mixed" agree (6.5e-4 vs 6.6e-4),
  so the SP-calibrated figure is the more robust one *if* the SP excess is
  general. It is not on the new squares (above).

**Estimates of E[magic] over these 7,021 squares** (1.0e-8 per square under
the regression):

| estimate | E[magic] | x regression |
| --- | ---: | ---: |
| regression | 7.2e-5 | 1 |
| regression x observed SP factor^2 (independent pairs) | 3.8e-4 [1.5e-4, 8.5e-4] | 5.3 [2.1, 11.9] |
| heuristic, pairs independent (no congruences either) | 1.4e-4 | 2.0 |
| heuristic, congruences and the data's anti-clustering | 1.6e-4 | 2.2 |
| **heuristic, pair factor x r(k) (sub-events)** | **1.8e-4** | **2.5** |
| heuristic, congruences only | 1.9e-4 | 2.7 |
| the last two x the sum part 1.19 | 2.1-2.3e-4 | 3.0-3.2 |
| heuristic as it stands | 2.5e-4 | 3.45 |
| heuristic "mixed" variant (disfavoured) | 3.8e-4 | 5.3 |
| heuristic x r(k) x pooled SP factor^2 | 4.7e-4 [1.8e-4, 1.1e-3] | 6.6 [2.6, 15] |
| heuristic as it stands x pooled SP factor^2 | 6.5e-4 [2.5e-4, 1.5e-3] | 9.1 [3.5, 20] |
| forecast.md's assumption (rho 2.6, pair 1.1) | 1.15e-4 | 1.6 |

My reading: the per-square P(magic) of squares like these is **about 2.5x
the regression (2.2-3.2x)** once the pair factor is corrected. It is up to
about 7x only if the SP excess of the re-found seed-P squares is general,
which the new squares do not support.

## What it means for research/forecast.md

* The forecast's central rows take P(magic | square) = 1.61x the model:
  rho 2.6, giving (2.6 / 2.148)^2 = 1.46, times a pair factor of 1.1.
* Its "fresh P x0.6" factor contains the fresh-P traversal rate, 0.35 x
  1.26^2. So on fresh P its P(magic | square) is 2.56x the model. With
  sched40's observed P rate (1.06x the model) as well, it would be 2.87x.
* The squares the forecast schedules are like sched40: scheduler units on
  fresh P, all at 1.2-2 S_min. The all-squares mix (mostly re-found seed-P
  squares near S_min) does not describe them, so it is not used here.
* **Like with like on sched40** (per-square P(magic) relative to the
  model):

| heuristic on sched40 | x model | / forecast (2.56) | with the observed S and P rates on both sides (x1.25 / 2.87) |
| --- | ---: | ---: | ---: |
| pairs independent (no congruences) | 1.45 | 0.57 | 0.63 |
| **pair factor x r(k)** | **1.63** | **0.64** | **0.71** |
| congruences only | 1.78 | 0.69 | 0.77 |
| either, x the sum part 1.19 | 1.94-2.12 | 0.76-0.83 | 0.85-0.92 |
| as it stands | 2.07 | 0.81 | 0.90 |
| "mixed" (disfavoured by the sub-events) | 2.83 | 1.10 | 1.23 |
| r(k) x pooled SP factor^2 (2.6 [1.0, 5.9]) | 4.3 [1.7, 9.6] | 1.67 [0.65, 3.75] | 1.49 [0.58, 3.34] (instead of the rates) |

* "Observed rates" scales the heuristic by its sched40 S and P misses
  ((1.02 x 1.10)^2 = 1.25) and the forecast by the model's (2.87 instead of
  2.56). That is the like-with-like column.
* **Net: E(C) is about 0.75x forecast.md's central values (0.6-0.9x).**
  * The 10 CPU-year central value of 0.03-0.045 becomes about 0.02-0.035.
  * It replaces only the forecast's rho and pair factors; its other factors
    (squares on fresh P, large N, many-divisor P) are untouched.
  * The pair factor is the main change from the previous version (0.8x
    before): the sub-events cut the heuristic's 1.43 on sched40 to about
    1.13.
* **Upper scenario.** If the pooled SP excess applied to fresh P, it would
  be 1.5-1.7x [0.6, 3.8]. The fresh-P SP rung cannot tell (2 vs 1.29, O/E
  1.55 [0.28, 4.88]). The new squares' sum + exponent sub-events, with
  10-50 times more events, show no coupling deficit (0.95 [0.78, 1.12]),
  so this scenario is now disfavoured, not merely untested.
* **Targeting.** By the heuristic, the regression undervalues near-S_min
  sums (6-25x below 1.2 S_min), 6-prime P (9x) and many-divisor P (10-22x
  for tau >= 4,000).
  * These ratios include its pair factor (2.1-2.8 near S_min and for 6
    primes), which the sub-events cut by about 0.7, and its rho near S_min,
    which is untested on fresh P (sched40 has no squares below 1.2 S_min).
  * If the rest holds, a scheduler ranking by the heuristic's P(magic)
    could add more than the forecast's E. Using the heuristic for unseen
    (P, S) would need a regression on its per-square output, because it
    needs the squares.

## Bins

Traversal rungs: O/E with 90% intervals; bins with no SP event show the SP
expectation. Pair rungs: S+S under each predictor, the top rung (observed /
regression / null / heuristic E) and E[magic]. `calibrate.py` writes the
same tables for every event.

* The heuristic gets S right in the S/S_min, tau, primes and P-source bins
  (0.89-1.06, within noise); N 4,000-8,000 is 1.07 [0.99, 1.17].
* The P deficit is not the same everywhere: 1.02-1.16 by S/S_min, 1.30
  [1.08, 1.55] at N 4,000-8,000, 0.76 [0.50, 1.11] at N >= 8,000, 1.35
  [1.11, 1.63] on verify-region-pool and 1.00 on existence.
* The regression's S rate has the wrong shape in S/S_min: 1.40x low below
  1.1 S_min and 1.30x at 1.5-2, against 1.06-1.16 in between.
* The heuristic's SP excess is 1.56-1.69 in every S/S_min bin that has SP
  events. That is mostly seed P everywhere (12 of 14 events).

**S / S_min**

| S/S_min | squares | S: reg | S: heur | P: reg | P: heur | SP obs | SP: reg | SP: heur |
|---|---:|---|---|---|---|---:|---|---|
| [1, 1.1) | 276 | 1.40 [1.26, 1.54] | 1.01 [0.91, 1.11] | 0.88 [0.66, 1.15] | 1.02 [0.77, 1.33] | 0 | (E 0.094) | (E 0.27) |
| [1.1, 1.2) | 2,808 | 1.12 [1.09, 1.15] | 0.99 [0.96, 1.01] | 1.07 [1.00, 1.15] | 1.12 [1.05, 1.19] | 6 | 2.94 [1.28, 5.80] | 1.69 [0.74, 3.33] |
| [1.2, 1.3) | 1,159 | 1.06 [1.01, 1.10] | 0.99 [0.95, 1.04] | 1.06 [0.98, 1.15] | 1.14 [1.05, 1.23] | 3 | 2.25 [0.61, 5.80] | 1.68 [0.46, 4.33] |
| [1.3, 1.5) | 1,429 | 1.16 [1.12, 1.21] | 1.00 [0.96, 1.04] | 1.08 [1.00, 1.16] | 1.16 [1.08, 1.25] | 3 | 1.81 [0.49, 4.68] | 1.56 [0.42, 4.03] |
| [1.5, 2) | 1,308 | 1.30 [1.24, 1.35] | 1.02 [0.98, 1.07] | 1.06 [0.98, 1.15] | 1.09 [1.01, 1.18] | 2 | 1.96 [0.35, 6.16] | 1.68 [0.30, 5.29] |
| >= 2 | 41 | 1.52 [1.11, 2.03] | 0.90 [0.66, 1.20] | 1.12 [0.66, 1.78] | 0.93 [0.55, 1.48] | 0 | (E 0.018) | (E 0.029) |

| S/S_min | S+S obs | S+S: reg | S+S: null | S+S: heur | best >= SP+S: obs / reg / null / heur | E[magic] reg | heur | heur / reg |
|---|---:|---|---|---|---|---:|---:|---:|
| [1, 1.1) | 5 | 2.60 [1.02, 5.46] | 1.34 [0.53, 2.82] | 1.41 [0.55, 2.96] | 0 / 0.002 / 0 / 0.008 | 4.1e-7 | 1.0e-5 | 25 |
| [1.1, 1.2) | 68 | 1.62 [1.31, 1.98] | 1.28 [1.03, 1.56] | 1.23 [1.00, 1.51] | 1 / 0.06 / 0.27 / 0.12 | 1.8e-5 | 1.1e-4 | 6.4 |
| [1.2, 1.3) | 20 | 1.29 [0.85, 1.87] | 1.07 [0.71, 1.56] | 1.11 [0.73, 1.61] | 0 / 0.041 / 0.062 / 0.059 | 1.7e-5 | 4.7e-5 | 2.8 |
| [1.3, 1.5) | 27 | 1.30 [0.92, 1.79] | 0.94 [0.66, 1.29] | 0.91 [0.64, 1.26] | 0 / 0.055 / 0.12 / 0.072 | 2.5e-5 | 5.3e-5 | 2.2 |
| [1.5, 2) | 27 | 2.63 [1.86, 3.63] | 1.45 [1.03, 2.01] | 1.54 [1.08, 2.12] | 1 / 0.027 / 0.082 / 0.039 | 1.1e-5 | 2.2e-5 | 2.0 |
| >= 2 | 0 | 0 [0, 24] | 0 [0, 20] | 0 [0, 8.5] | 0 / 0.0003 / 0 / 0.0008 | 1.0e-7 | 4.3e-7 | 4.4 |

**tau(P)**

| tau(P) | squares | S: reg | S: heur | P: reg | P: heur | SP obs | SP: reg | SP: heur |
|---|---:|---|---|---|---|---:|---|---|
| [0, 2,000) | 1,731 | 1.27 [1.23, 1.31] | 1.00 [0.97, 1.04] | 1.08 [1.02, 1.15] | 1.11 [1.04, 1.18] | 4 | 1.98 [0.68, 4.52] | 1.66 [0.57, 3.79] |
| [2,000, 4,000) | 2,594 | 1.05 [1.01, 1.08] | 0.98 [0.95, 1.01] | 1.10 [1.04, 1.15] | 1.17 [1.11, 1.23] | 8 | 2.79 [1.39, 5.03] | 2.01 [1.00, 3.63] |
| [4,000, 8,000) | 2,393 | 1.17 [1.13, 1.20] | 1.01 [0.98, 1.04] | 0.98 [0.90, 1.08] | 1.02 [0.93, 1.11] | 2 | 1.64 [0.29, 5.17] | 0.89 [0.16, 2.81] |
| [8,000, 16,000) | 274 | 1.27 [1.09, 1.47] | 0.93 [0.80, 1.08] | 0.91 [0.69, 1.19] | 1.32 [0.99, 1.71] | 0 | (E 0.049) | (E 0.11) |
| >= 16,000 | 29 | 1.74 [1.22, 2.42] | 0.89 [0.62, 1.23] | 0 [0, 0.60] | 0 [0, 3.5] | 0 | (E 0.008) | (E 0.009) |

| tau(P) | S+S obs | S+S: reg | S+S: null | S+S: heur | best >= SP+S: obs / reg / null / heur | E[magic] reg | heur | heur / reg |
|---|---:|---|---|---|---|---:|---:|---:|
| [0, 2,000) | 44 | 1.85 [1.42, 2.38] | 1.12 [0.86, 1.44] | 1.14 [0.87, 1.46] | 1 / 0.065 / 0.20 / 0.092 | 2.9e-5 | 6.2e-5 | 2.2 |
| [2,000, 4,000) | 48 | 1.51 [1.17, 1.92] | 1.34 [1.04, 1.71] | 1.28 [0.99, 1.63] | 0 / 0.086 / 0.27 / 0.13 | 3.6e-5 | 1.1e-4 | 3.1 |
| [4,000, 8,000) | 55 | 1.59 [1.25, 1.99] | 1.17 [0.92, 1.46] | 1.17 [0.92, 1.46] | 1 / 0.035 / 0.062 / 0.077 | 7.1e-6 | 7.2e-5 | 10 |
| [8,000, 16,000) | 0 | 0 [0, 7.6] | 0 [0, 4.4] | 0 [0, 3.8] | 0 / 0.0005 / 0 / 0.0019 | 9.2e-8 | 2.0e-6 | 22 |
| >= 16,000 | 0 | 0 [0, 37] | 0 [0, 9.0] | 0 [0, 9.1] | 0 / 0.0001 / 0 / 0.0002 | 2.1e-8 | 0 | 0 |

**N** (the number of vectors of the sum)

| N | squares | S: reg | S: heur | P: reg | P: heur | SP obs | SP: reg | SP: heur |
|---|---:|---|---|---|---|---:|---|---|
| [0, 1,000) | 155 | 1.18 [1.06, 1.30] | 0.94 [0.85, 1.04] | 0.97 [0.79, 1.18] | 1.00 [0.81, 1.22] | 1 | 3.4 [0.17, 16] | 2.4 [0.12, 11] |
| [1,000, 2,000) | 4,351 | 1.17 [1.15, 1.20] | 1.00 [0.98, 1.02] | 1.08 [1.03, 1.12] | 1.12 [1.08, 1.17] | 12 | 2.62 [1.51, 4.25] | 1.88 [1.09, 3.05] |
| [2,000, 4,000) | 1,671 | 1.10 [1.06, 1.14] | 0.97 [0.94, 1.01] | 1.09 [1.00, 1.19] | 1.15 [1.05, 1.25] | 1 | 0.91 [0.05, 4.34] | 0.57 [0.03, 2.69] |
| [4,000, 8,000) | 544 | 1.17 [1.07, 1.27] | 1.07 [0.99, 1.17] | 1.01 [0.84, 1.21] | 1.30 [1.08, 1.55] | 0 | (E 0.15) | (E 0.15) |
| >= 8,000 | 300 | 1.15 [1.02, 1.29] | 1.05 [0.94, 1.18] | 0.64 [0.43, 0.93] | 0.76 [0.50, 1.11] | 0 | (E 0.058) | (E 0.055) |

| N | S+S obs | S+S: reg | S+S: null | S+S: heur | best >= SP+S: obs / reg / null / heur | E[magic] reg | heur | heur / reg |
|---|---:|---|---|---|---|---:|---:|---:|
| [0, 1,000) | 12 | 3.24 [1.87, 5.25] | 2.30 [1.33, 3.73] | 2.03 [1.17, 3.28] | 0 / 0.012 / 0.081 / 0.020 | 6.5e-6 | 1.8e-5 | 2.7 |
| [1,000, 2,000) | 95 | 1.58 [1.33, 1.88] | 1.11 [0.93, 1.32] | 1.13 [0.95, 1.34] | 2 / 0.14 / 0.43 / 0.22 | 5.6e-5 | 1.8e-4 | 3.3 |
| [2,000, 4,000) | 36 | 1.54 [1.14, 2.03] | 1.29 [0.96, 1.70] | 1.19 [0.89, 1.58] | 0 / 0.031 / 0.021 / 0.056 | 8.3e-6 | 4.5e-5 | 5.4 |
| [4,000, 8,000) | 2 | 0.86 [0.15, 2.69] | 0.65 [0.12, 2.04] | 0.73 [0.13, 2.30] | 0 / 0.0025 / 0 / 0.0027 | 4.8e-7 | 1.1e-6 | 2.2 |
| >= 8,000 | 2 | 1.62 [0.29, 5.09] | 1.50 [0.27, 4.72] | 1.29 [0.23, 4.06] | 0 / 0.0009 / 0 / 0.0010 | 1.2e-7 | 3.5e-7 | 3.0 |

* Sums with N >= 4,000 (844 squares) expect only about 0.2 SP traversals,
  so they cannot test the SP rate. 12 of the 14 SP traversals are at
  N = 1,000-2,000.

**Number of primes of P**

| primes | squares | S: reg | S: heur | P: reg | P: heur | SP obs | SP: reg | SP: heur |
|---|---:|---|---|---|---|---:|---|---|
| 4 | 1,081 | 1.16 [1.11, 1.23] | 1.01 [0.96, 1.06] | 1.09 [1.02, 1.16] | 1.13 [1.06, 1.20] | 3 | 1.93 [0.53, 4.98] | 1.63 [0.44, 4.21] |
| 5 | 2,830 | 1.11 [1.08, 1.14] | 0.99 [0.96, 1.02] | 1.07 [1.02, 1.13] | 1.13 [1.07, 1.20] | 8 | 2.62 [1.30, 4.73] | 1.91 [0.95, 3.44] |
| 6 | 3,081 | 1.19 [1.16, 1.22] | 1.00 [0.97, 1.02] | 1.02 [0.95, 1.11] | 1.10 [1.02, 1.18] | 3 | 1.93 [0.53, 4.99] | 1.11 [0.30, 2.86] |
| 8 | 29 | 1.74 [1.22, 2.42] | 0.89 [0.62, 1.23] | 0 [0, 0.60] | 0 [0, 3.5] | 0 | (E 0.008) | (E 0.009) |

| primes | S+S obs | S+S: reg | S+S: null | S+S: heur | best >= SP+S: obs / reg / null / heur | E[magic] reg | heur | heur / reg |
|---|---:|---|---|---|---|---:|---:|---:|
| 4 | 12 | 1.55 [0.90, 2.52] | 1.05 [0.60, 1.70] | 1.14 [0.66, 1.85] | 1 / 0.047 / 0.14 / 0.061 | 2.6e-5 | 4.5e-5 | 1.8 |
| 5 | 68 | 1.72 [1.39, 2.10] | 1.35 [1.09, 1.65] | 1.35 [1.09, 1.65] | 0 / 0.095 / 0.29 / 0.14 | 3.7e-5 | 1.2e-4 | 3.2 |
| 6 | 67 | 1.55 [1.25, 1.90] | 1.10 [0.89, 1.34] | 1.06 [0.86, 1.30] | 1 / 0.045 / 0.10 / 0.095 | 9.5e-6 | 8.4e-5 | 8.8 |
| 8 | 0 | 0 [0, 37] | 0 [0, 9.0] | 0 [0, 9.1] | 0 / 0.0001 / 0 / 0.0002 | 2.1e-8 | 0 | 0 |

**Source** (of the square) **and P source** (first source with that P)

| source | squares | S: reg | S: heur | P: reg | P: heur | SP obs | SP: reg | SP: heur |
|---|---:|---|---|---|---|---:|---|---|
| seed36 | 4,920 | 1.12 [1.10, 1.15] | 0.99 [0.97, 1.01] | 1.08 [1.03, 1.13] | 1.13 [1.08, 1.18] | 12 | 2.62 [1.51, 4.25] | 1.75 [1.01, 2.84] |
| sched40 | 901 | 1.26 [1.20, 1.32] | 1.02 [0.97, 1.07] | 1.06 [0.98, 1.15] | 1.10 [1.01, 1.19] | 2 | 1.83 [0.33, 5.76] | 1.55 [0.28, 4.88] |
| verify-recompute | 152 | 1.50 [1.26, 1.78] | 0.90 [0.76, 1.07] | 0.65 [0.42, 0.98] | 0.95 [0.60, 1.42] | 0 | (E 0.033) | (E 0.094) |
| verify-region-pool | 489 | 1.16 [1.06, 1.26] | 1.08 [0.98, 1.17] | 1.04 [0.85, 1.26] | 1.35 [1.11, 1.63] | 0 | (E 0.14) | (E 0.14) |
| verify-model-form | 36 | 1.24 [0.77, 1.91] | 1.09 [0.67, 1.67] | 1.20 [0.76, 1.83] | 1.37 [0.86, 2.08] | 0 | (E 0.014) | (E 0.014) |
| existence | 514 | 1.24 [1.16, 1.33] | 0.99 [0.93, 1.06] | 0.96 [0.81, 1.14] | 1.00 [0.84, 1.19] | 0 | (E 0.30) | (E 0.36) |
| other | 9 | 1.05 [0.61, 1.70] | 0.77 [0.44, 1.24] | 1.36 [0.47, 3.12] | 1.30 [0.45, 2.98] | 0 | (E 0.013) | (E 0.016) |
| *P of seed36* | 5,594 | 1.13 [1.10, 1.15] | 1.00 [0.98, 1.02] | 1.07 [1.03, 1.12] | 1.13 [1.08, 1.18] | 12 | 2.52 [1.45, 4.08] | 1.71 [0.99, 2.77] |
| *P of verify-recompute* | 284 | 1.34 [1.17, 1.54] | 0.91 [0.79, 1.05] | 0.82 [0.61, 1.07] | 1.27 [0.95, 1.67] | 0 | (E 0.054) | (E 0.12) |
| *P of existence* | 235 | 1.30 [1.19, 1.42] | 0.96 [0.88, 1.05] | 1.10 [0.90, 1.34] | 1.10 [0.90, 1.33] | 0 | (E 0.25) | (E 0.31) |

| source | S+S obs | S+S: reg | S+S: null | S+S: heur | best >= SP+S: obs / reg / null / heur | E[magic] reg | heur | heur / reg |
|---|---:|---|---|---|---|---:|---:|---:|
| seed36 | 111 | 1.55 [1.32, 1.82] | 1.21 [1.03, 1.41] | 1.17 [1.00, 1.37] | 1 / 0.14 / 0.43 / 0.23 | 5.2e-5 | 2.0e-4 | 3.9 |
| sched40 | 20 | 1.81 [1.20, 2.64] | 1.05 [0.70, 1.53] | 1.17 [0.78, 1.71] | 1 / 0.034 / 0.10 / 0.046 | 1.5e-5 | 3.2e-5 | 2.1 |
| verify-recompute | 0 | 0 [0, 10.5] | 0 [0, 3.6] | 0 [0, 3.5] | 0 / 0.0004 / 0 / 0.0018 | 7.5e-8 | 1.9e-6 | 26 |
| verify-region-pool | 2 | 0.88 [0.16, 2.77] | 0.68 [0.12, 2.13] | 0.76 [0.13, 2.39] | 0 / 0.0024 / 0 / 0.0025 | 4.5e-7 | 8.6e-7 | 1.9 |
| verify-model-form | 0 | 0 [0, 66] | 0 [0, 144] | 0 [0, 49] | 0 / 0.0002 / 0 / 0.0002 | 6.7e-8 | 8.7e-8 | 1.3 |
| existence | 14 | 2.66 [1.61, 4.17] | 1.69 [1.02, 2.64] | 1.58 [0.95, 2.46] | 0 / 0.0097 / 0 / 0.015 | 3.8e-6 | 1.1e-5 | 2.9 |
| other | 0 | 0 [0, 16] | 0 [0, 24] | 0 [0, 8.6] | 0 / 0.0006 / 0 / 0.0009 | 2.9e-7 | 3.4e-7 | 1.2 |
| *P of seed36* | 115 | 1.53 [1.31, 1.79] | 1.20 [1.02, 1.40] | 1.17 [0.99, 1.36] | 1 / 0.14 / 0.43 / 0.24 | 5.2e-5 | 2.0e-4 | 3.9 |

(The P of sched40 are exactly its squares; "P of other" holds 7 squares.)

**First-search range**

| squares | squares | S: reg | S: heur | P: reg | P: heur | SP obs | SP: reg | SP: heur |
|---|---:|---|---|---|---|---:|---|---|
| inside its range | 4,741 | 1.12 [1.10, 1.15] | 0.98 [0.96, 1.00] | 1.09 [1.04, 1.13] | 1.14 [1.09, 1.18] | 11 | 2.37 [1.33, 3.93] | 1.63 [0.91, 2.70] |
| new | 2,280 | 1.23 [1.19, 1.28] | 1.03 [1.00, 1.07] | 1.02 [0.95, 1.09] | 1.10 [1.02, 1.17] | 3 | 1.95 [0.53, 5.05] | 1.49 [0.41, 3.86] |
| new, seed P | 1,067 | 1.19 [1.13, 1.26] | 1.06 [1.01, 1.12] | 0.98 [0.85, 1.12] | 1.07 [0.93, 1.23] | 1 | 2.6 [0.14, 12.5] | 1.7 [0.09, 8.1] |
| new, other P | 1,213 | 1.27 [1.21, 1.33] | 1.01 [0.96, 1.06] | 1.03 [0.95, 1.11] | 1.11 [1.02, 1.20] | 2 | 1.73 [0.31, 5.44] | 1.41 [0.25, 4.43] |

| squares | S+S obs | S+S: reg | S+S: null | S+S: heur | best >= SP+S: obs / reg / null / heur | E[magic] reg | heur | heur / reg |
|---|---:|---|---|---|---|---:|---:|---:|
| inside its range | 112 | 1.59 [1.35, 1.86] | 1.22 [1.04, 1.43] | 1.18 [1.00, 1.38] | 0 / 0.14 / 0.39 / 0.23 | 5.5e-5 | 2.0e-4 | 3.6 |
| new | 35 | 1.73 [1.28, 2.30] | 1.10 [0.82, 1.46] | 1.18 [0.87, 1.57] | 2 / 0.043 / 0.14 / 0.065 | 1.7e-5 | 5.3e-5 | 3.1 |
| new, seed P | 15 | 1.77 [1.09, 2.73] | 1.34 [0.83, 2.06] | 1.36 [0.84, 2.10] | 1 / 0.008 / 0.041 / 0.016 | 1.6e-6 | 1.8e-5 | 11 |
| new, other P | 20 | 1.70 [1.13, 2.47] | 0.98 [0.65, 1.42] | 1.07 [0.71, 1.56] | 1 / 0.035 / 0.10 / 0.049 | 1.5e-5 | 3.5e-5 | 2.2 |

## What to measure next

1. **SP traversals and sum + exponent sub-events on fresh P near S_min,
   outside the first search's range** (forecast.md's item 1). Require P the
   first search did not pick from and (P, S) above its max_S.
   * The heuristic expects about 1.25e-3 SP traversals per square near
     S_min, so 50,000 squares give 60-100 SP events and pin the SP factor
     to +-15%.
   * The sum + 3 and + 4 exponent sub-events come 50 and 10 times as often
     as SP (on the 2,280 new squares: 149 and 31 events vs 3). They give the
     coupling to a few percent with the same squares, and show whether the
     re-found seed-P excess near S_min recurs anywhere else.
   * They also test the heuristic where it disagrees most with the
     regression (S/S_min < 1.2, 6 primes, tau >= 4,000). sched40 has no
     squares below 1.2 S_min.
   * Run `scripts/calibrate.py` on the new units; the sub-event tables come
     with it.
2. **The pair factor at d >= 5.** r(5) rests on 65 pairs and r(6), r(7) are
   extrapolated. The same 50,000 squares would give about 500 pairs at
   d = 5 and the first events at d = 6.
3. **The S+S excess** (1.19x, also with one exponent pinned). It needs no
   new runs, only a model of the correlation between the two sums beyond
   Hoeffding's covariance and its Edgeworth terms. A model of the
   anti-clustering of the exponents would serve the same purpose: the
   heuristic's Gaussian pair has the wrong sign there.

## The heuristic (scripts/square_heuristic.py)

For a traversal sigma, Z(sigma) = sum_i z[i, sigma(i)], with
z[i, j] = (a_ij, v_p1(a_ij), ..., v_pk(a_ij)). Sigma is S-type iff
Z_0 = S, and P-type iff the remaining coordinates equal v(P).

1. **Exact moments.** For a semi-magic square, E[Z] over a uniform random
   permutation is exactly t = (S, v(P)), so the targets are the centre of
   the distribution. Hoeffding's formula gives
   Cov Z = 1/(n-1) sum_ij d_ij d_ij^T, with d = z - t/n.
   * For a partner pair (sigma, sigma o tau), the cross-covariance is
     C_tau = 1/(n-1) sum_i sum_c d_ic d_tau(i)c^T. Averaged over the 15 tau
     it is exactly -Sigma/5.
   * Both were checked against enumeration (`test_calibrate.py`; the review
     to 6e-16 on 60 squares).
2. **Lattice local CLT at the mean.** Z lies on Z(id) + L, where L is
   generated by the 2x2 interaction contrasts.
   * P(Z_I = t_I) ~ [t_I in Z(id)_I + L_I] covol(L_I) / ((2 pi)^(r/2)
     sqrt(det Cov Z_I)).
   * There is no exponential factor, because t is the mean.
3. **Edgeworth factors.** Every coordinate gets its univariate first-order
   factor, and every pair of coordinates its bivariate factor, combined in a
   Kirkwood product and each clamped to [0.25, 4].
   * This truncation was chosen because it fits the P+P rung and the
     12-cell product event (904 observed vs 953), so those are in-sample.
   * The sub-events support it out of sample on two-diagonal totals, where
     the other variants over-predict increasingly.
   * On one diagonal it is 1-12% low, more so as more exponents are pinned;
     the P rung's 12% deficit is one case of this.
4. **Pairs.** The pair (Z1, Z2) has the exact covariance and the exact
   lattice of the pair. Classes come from inclusion-exclusion over the point
   events {X = S} and {Y = v(P)}.
   * The Gaussian density at the mean is higher for correlated pairs
     whatever the sign of the correlation. The data show partner diagonals
     slightly anti-clustered in their exponents, so this part of the pair
     factor is too high (see the sub-events).
5. **Congruences of partner pairs.** For every square and every partner
   involution tau, Z1 + Z2 satisfies at least k + 1 - 4 independent linear
   congruences mod 3, where k is the number of primes of P.
   * The constant of each congruence is u . t, while a magic pair needs
     u . 2t. So a pair of type tau can be magic only if u . (S, v(P)) = 0
     mod 3 for all such u.
   * Proof: the merged rows R_p = z_{i_p} + z_{i_p'} have equal column sums.
     So mod 3 every 2x2 contrast is ±(gamma_c - gamma_c') with
     gamma_c = R_1c - R_2c, and sum_c gamma_c = 0 bounds the rank by 4.
     For u orthogonal to that span, u . R is separable, and the row and
     column sums give u . (Z1 + Z2) = u . t.
   * This rules out magic for 90% of all (square, tau). The review checked
     the rank bound on 105,315 (square, tau) and the zeros on an integer
     lattice: 2,379 of the 7,021 squares (34%) provably cannot be permuted
     into a magic square, the SP+S square among them.
   * The surviving pairs are 3^d times likelier. The sub-events confirm
     both the zeros and the size of the boost.
6. **best_score** uses Poisson clumping.
7. **Conditional mode.**
   * S and P are the realized counts.
   * SP = n_S n_P / 720 x rho_square.
   * Pairs are the own-count null times the model's per-tau pair
     correlation, with the exclusive classes: an SP diagonal's S_only
     partners are (SP, S inclusive) minus (SP, SP). The first version took
     them inclusive, a bias of about q_SP / p_S (0.2%).
8. **Sub-events** (`components(..., sub_events=True)`): the same machinery
   for every window of coordinates, with exact counts and the split by
   congruence class.

It takes about 0.1 s per square, 0.25 s with the sub-events.

## Approximations, explicitly

1. **Regression and null.** Traversal classes are independent across the
   720 traversals of a square. For the null they are uniform given the
   counts.
   * The regression is in-sample on the seed squares, which it counts twice
     (as rows and inside the first search's totals). Its comparison on the
     re-found squares inside the first search's range is in-sample too.
2. **Best-pair probabilities, regression and null.** These use W(x, y) on the
   exact partner graph:
   * exact for y <= 2 and for (x, y) = (0, 3);
   * Monte Carlo elsewhere, with s.e. <= 0.0022;
   * a Poisson approximation outside the table, only in cells of
     multinomial mass <= 1e-11;
   * the regression's mixture over the class counts is truncated at
     a <= 40, b <= 25, c <= 6, with captured mass >= 1 - 1e-11.
3. **Heuristic.**
   * The local CLT uses only 6 terms.
   * The pairwise Edgeworth (Kirkwood) product is an ad-hoc truncation,
     chosen on P+P and the 12-cell product event. Its variants are listed as
     model uncertainty.
   * The size of the congruence boost assumes smoothness on the coarse
     lattice (confirmed by the sub-events up to d = 4).
   * Best scores use Poisson clumping.
   * It extrapolates two rungs, from SP to SP+SP.
4. **Pair-factor correction.** The sub-events have pairs for up to 5
   pinned exponents per diagonal; magic pins k = 4-6 exponents and the sum.
   * r(6) and r(7) are extrapolated (log-linear in d from r(5)).
   * The sum coordinate is not corrected: S+S shows the opposite residual
     (1.19), so it enters as a range.
   * The correction is applied to the heuristic's E[magic] square by square
     by its number of primes. It assumes the residual, measured on windows
     of consecutive primes over all squares, holds for the full set of
     primes of each square.
5. **Intervals** treat E as exact.
   * Poisson intervals ignore the clustering of pairs within a square. The
     bootstrap in the tool's output covers it; it matters only for SP+0.
   * The SP-factor posterior assumes one common factor for all squares; the
     sub-events show it is not common (seed P near S_min only).
6. **First search.**
   * Only per-P totals are available, so its checks are at the P level.
   * Its best-pair counts (9 SP+S, 1 SP+P) are from the report.
   * The regression is evaluated at each P's midpoint sum, as the scheduler
     does for legacy totals.
   * "Inside its range" means it covered that (P, S). The square itself can
     be checked against it only through its top-100 list.
   * The pooled top rung uses each P's i.i.d. rates for the first search,
     which run about 13% below the own-count null on these squares; both
     are given.
7. **Data.** Records come from the October 2026 msearch runs (not in the
   repo) and are deduplicated by hash (a grid-based canonical form gives the
   same classes).
   * The existence runs were still growing during the work. This analysis
     is frozen to 7,021 hashes (`--hashes`).
   * By the end there were 7,058 squares. The 37 new ones have no SP
     traversal and change none of the conclusions.
   * N comes from the msearch sum records. S_min comes from the scheduler's
     pinfo caches, with `bin/enumerate` as the fallback; both agree for all
     113 P.
8. **Forecast translation.** It rescales only the rho and pair factors of
   forecast.md's central rows, on squares like sched40, and assumes the
   forecast's future squares resemble them. sched40 has no squares below
   1.2 S_min, where the heuristic's rho is highest and least tested.

## Reproduce

```
python3 scripts/test_calibrate.py                     # ~15 s
python3 scripts/calibrate.py --jobs 4                 # data/sched/units -> data/calibrate/
python3 scripts/calibrate.py --jobs 4 --heuristic-cond --eb \
    --state data/sched seed=runs/seed36 fresh=data/sched/units   # labels name the sources
```

* Output, in `--out`, default `data/calibrate`:
  * `ladder.md`: every table, overall and binned by S/S_min, tau(P), N,
    number of primes, source, P source and first-search range; the
    sub-event tables (overall and by P source, first-search range and
    source); E[magic] for every predictor, variant and correction; the first
    search's P-level checks and the pooled top rung;
  * `results.json`: every number;
  * `per_square.jsonl`: per square, observed counts and every predictor.
* The heuristic, its variants and the sub-events are cached per square,
  keyed by the hash of `square_heuristic.py`'s code. `--no-sub-events` skips
  the sub-events (2.5x cheaper).
* For the 7,021 squares the heuristic with sub-events took about 30
  CPU-minutes; the rest takes under a minute.
* The run for this write-up used `--hashes` to restrict it to the 7,021
  squares that existed when the analysis was frozen.

## Changes in this version (after two reviews)

* **Sub-events** in `calibrate.py` and `square_heuristic.py`, with the
  split by congruence class. They show that the pair factor can be observed
  below the magic rung; that the heuristic's congruence part holds and its
  correlation part does not; that the default Edgeworth variant is right in
  total; and that the sum-product coupling deficit is confined to re-found
  seed-P squares. E[magic] now comes with the corrected pair factor:
  2.5e-4 becomes 1.8e-4 (1.6-2.3e-4).
* **First-search range and P source** splits in every SP and top-rung
  table: 11 of the 14 SP traversals are on re-found squares from P selected
  on the first search. "Seed P" is now by P (5,594 squares; outside them SP
  is 2 vs 1.74).
* **Model uncertainty**: the Edgeworth truncation was chosen on P+P and the
  12-cell event; "nothing fitted" is gone, and E[magic] is given for every
  variant with the evidence for preferring the default.
* **Forecast translation** on sched40 only, like with like. It is now
  0.75x (0.6-0.9x) rather than 0.8-1.6x, with the upper scenario disfavoured
  instead of 2-4x. The "not too optimistic" sentence is gone.
* **Double counting**: the SP-containing rungs are not separate
  confirmations; the top rung is given conditional on the SP count.
* **Pooled top rung** with the first search: only the squares outside its
  range, with the bias of its per-P expectation: 0.83-0.96.
* **Bins**: tau, N, primes, source, P-source and first-search tables in
  full, S/S_min >= 2 listed, bin claims qualified.
* **Conditional mode**: SP+S and SP+P use the exclusive classes.
* **Regression**: described as in-sample on the seed squares.
