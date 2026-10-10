#ifndef DFIRST_H
#define DFIRST_H

#include <stdint.h>
#include <stdio.h>

#include "arrange.h"

/*
 * Diagonal-first search ("d-first"): the squares that have a given vector d
 * as a traversal (one number in every row and every col: a possible
 * diagonal) are exactly the semi-magic squares built from
 *   V_d = { v : |v & d| = 1 },
 * so running the semi-magic search on V_d for every d of a "diagonal set"
 * finds every (square, traversal in the set) pair exactly once. With the
 * vectors of one (P, S) as the diagonal set (the unreduced list: a diagonal
 * need not survive the reduction), these are the (square, SP traversal)
 * pairs, and a magic square is a square with two of them that fit on the
 * two diagonals together, so it is found (at least) twice, once per
 * diagonal (see research/ideas.md, "Diagonal-first search"; at least
 * once with the star filter below).
 *
 * V_d is built from a number -> vector posting index over the search's
 * vectors: a hit counter per vector is bumped along the postings of d's
 * numbers, the vectors hit once are kept, and only the touched counters are
 * cleared (O(sum of the posting lengths) per d, not O(n^2 N)). V_d keeps the
 * order of the vector list (the search relabels and sorts it anyway).
 */

typedef struct dfirst_s dfirst_t;

/* one (square, d) pair */
typedef struct {
  const square_t *sq;
  const uint64_t *d; /* the diagonal: n numbers */
  size_t d_index;    /* its index in the diagonal list */
  int set_count;     /* traversals of sq in the diagonal set (>= 1, d) */
  int partner;       /* 1 if some other traversal t of sq in the diagonal set
                        goes on the other diagonal when d is the main one
                        (rows and cols permuted): with the vectors of one
                        (P, S) as the set, sq is then a magic square */
} dsquare_t;

/* called on every pair; return nonzero to stop */
typedef int (*dsquare_cb)(const dsquare_t *ds, void *ctx);

typedef struct {
  uint64_t nd;       /* diagonals searched */
  uint64_t vd_total; /* sum of |V_d| */
  uint64_t nodes;    /* search nodes, summed over d */
  uint64_t pairs;    /* (square, d) pairs found */
  uint64_t partners; /* of them with partner = 1 */
  double vd_seconds, setup_seconds, search_seconds; /* V_d, search setup,
                                                        search (wall) */
  double cpu_seconds; /* thread CPU time of the whole loop */
  double d_cpu, d_cpu2, d_pairs2; /* sums over d of the per-d CPU seconds,
                                     their squares and the squares of the
                                     pairs (for standard errors) */
  int truncated;      /* some V_d search hit opts->node_limit */
  int stopped;        /* the callback asked to stop */
  /* the star filter (dfirst_set_star; all 0 without it): the d of the
   * sample that it skipped (star d not of a searched rank, and with only
   * the d without the star number), and of the d searched (counted in nd,
   * nodes, pairs, d_cpu.. above too) the star d, with their nodes, pairs
   * and the sums over them of the per-d CPU, its squares and the squares
   * of the pairs */
  uint64_t nd_star_skipped, nd_other_skipped, nd_star;
  uint64_t nodes_star, pairs_star;
  double d_cpu_star, d_cpu2_star, d_pairs2_star;
  /* the d whose search ran the class support (opts->class_support, the
   * gate of dfirst_set_class_min_labels, and search_stats_t.class_used:
   * the top-label root on the carried path, V_d of at least 2n vectors) */
  uint64_t nd_class;
} dfirst_stats_t;

/*
 * The index over the search's vectors vecs[vstart, vstart + vcount) and the
 * diagonal set diags[dstart, dstart + dcount) (both with the same n; the
 * lists must outlive the index).
 */
dfirst_t *dfirst_new(const vec_list_t *vecs, size_t vstart, size_t vcount,
                     const vec_list_t *diags, size_t dstart, size_t dcount);
void dfirst_free(dfirst_t *df);

/* V_d of diagonal i (0-based within the diagonal set), appended to out in
 * the order of vecs; returns |V_d| */
size_t dfirst_vd(dfirst_t *df, size_t i, vec_list_t *out);

/*
 * The root of the V_d searches. Every square of V_d contains all n numbers
 * of d (its rows meet d once each and are disjoint), so by default (on = 1)
 * d's numbers get the largest labels, the one in the fewest vectors of V_d
 * the largest (search_opts_t top_numbers), and only the first rows through
 * it are searched (top_root_only): the roots are that class of V_d (~1/n
 * of it) instead of all of V_d. 0.6-0.7x the nodes and time of the plain
 * root (on = 0, as bin/dsearch). The (square, d) pairs are the same.
 * With on = 1, a V_d searched with the intersection matrices of 8 words
 * (257-512 labels in a build with -DCARRY_MAX_W below 8) gets the plain
 * root: the top root took 1.27x its CPU there (and 0.68x with the
 * matrices of up to 4 words; research/ideas.md).
 */
void dfirst_set_top_root(dfirst_t *df, int on);

/*
 * The gate of the class support (search_opts_t class_support, which the
 * caller's opts turn on): it runs only in the V_d with at least min_labels
 * labels (0: in every V_d). Default DFIRST_CLASS_MIN_LABELS = 137, from the
 * paired per-d CPU of 20 sums at N = 0.7-7.6k (on / off, by the V_d's
 * labels): 1.33-1.47x at <= 124 labels (two 64-bit words; every |V_d|
 * from 200 to 1,000), 1.12x at 125-128, 1.05x at 129-132, 1.00x at
 * 133-136, 0.97x at 137-140 and 0.85-0.93x above. |V_d| alone does not
 * separate them (two-word V_d of 600-1,000 vectors: 1.36-1.43x; wider
 * ones of the same size: 0.94-0.97x). Gated there, no sum of the 20 was
 * slower (worst 1.000; at 129: 1.035), and the sums at N >= 5.7k kept their
 * gain (research/ideas.md, "Integration of the round-2 d-first changes").
 * The gate costs one pass over V_d's numbers (vd_labels, which the
 * top-root choice computes anyway where the sum has few labels).
 */
#define DFIRST_CLASS_MIN_LABELS 137
void dfirst_set_class_min_labels(dfirst_t *df, uint32_t min_labels);

/* distinct numbers among the search's vectors */
uint32_t dfirst_num_labels(const dfirst_t *df);

/* is the n-set t (any order) in the diagonal set? */
int dfirst_in_set(const dfirst_t *df, const uint64_t *t);

/*
 * The star cover (research/ideas.md, "The star cover of the d loop"). The
 * two diagonals of a magic square of even order n have no cell in common,
 * so (its n^2 numbers being distinct) no number in common: whatever the
 * number x, at most one of them contains x, and the square is found from
 * the other one. So the d-first loop still finds every magic square (once
 * instead of twice when one diagonal contains x), and flags it (the
 * partner test runs against the whole diagonal set), without the "star
 * d", the diagonals that contain x. They lose only (square, d) pairs,
 * which a sample of them (every k-th) estimates. Not for odd n: the two
 * diagonals share the center cell.
 *
 * dfirst_star_choose picks x = x*, the number whose star d hold the most
 * predicted cost: x* = argmax_x sum over the d containing x of
 * |V_d|^5.5 R_d^1.4 (R_d: the rarest class of V_d, i.e. the fewest vectors
 * of V_d through one number of d, which the top-label root searches; the
 * per-d cost law of research/wip/complexity-profile.md), ties to the
 * smaller x, over the whole diagonal set (so every --d-range unit and
 * --d-stride sample of a sum agrees). One pass over the postings of every
 * d (no V_d is built).
 */
typedef struct {
  uint64_t x;          /* x* (0 if no number of the set is in vecs) */
  double share;        /* predicted cost share of its star d:
                          sum over them of the weights / sum over all d */
  uint64_t nstar;      /* d of the set that contain x */
  int freq_rank;       /* x's rank (1 = most) by the number of d of the
                          set containing it, ties counted ahead */
  uint64_t top_freq;   /* the number in the most d (ties: smaller) */
  double top_freq_share; /* its predicted cost share */
  double seconds;      /* thread CPU of the choice */
} dfirst_star_t;
dfirst_star_t dfirst_star_choose(dfirst_t *df);
/* the same stats for the given number x (x* forced, msearch
 * --dfirst-star-x: a scheduler continuing a sum whose parts were searched
 * with that x*); x = 0: dfirst_star_choose */
dfirst_star_t dfirst_star_choose_x(dfirst_t *df, uint64_t x);

/*
 * The star filter of dfirst_search, with x the star number: the star d
 * are ranked by their index in the diagonal set (0, 1, ..); k >= 1:
 * search only the star d of rank % k == 0 (k = 1: all of them); k = -1:
 * none; k = 0: no filter (only = 0). With only, the d without x are not
 * searched (a measurement of the star d's cost; k = 0 is then k = 1). The
 * d that the filter skips are counted in nd_star_skipped /
 * nd_other_skipped, not searched, not logged.
 */
void dfirst_set_star(dfirst_t *df, uint64_t x, int k, int only);
/* the star d with index in [lo, hi) of the diagonal set (0 without a star
 * filter) */
uint64_t dfirst_star_count(const dfirst_t *df, size_t lo, size_t hi);

/*
 * Search V_d for the diagonals i = lo + off, lo + off + stride, ... < hi
 * (stride 1, off 0: all of [lo, hi); stride > 1 with a random off: an
 * unbiased sample, stride x the totals estimates the full loop), except the
 * d that the star filter skips. If dlog is set, one line per d searched:
 * "i |V_d| labels nodes pairs vd_s setup_s search_s cpu_s", plus " star"
 * (1 for a star d, else 0) with a star filter.
 */
dfirst_stats_t dfirst_search(dfirst_t *df, size_t lo, size_t hi, size_t stride,
                             size_t off, const search_opts_t *opts,
                             dsquare_cb cb, void *ctx, FILE *dlog);

#endif // !DFIRST_H
