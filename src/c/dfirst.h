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
 * diagonal (see research/ideas.md, "Diagonal-first search").
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
  double d_cpu2, d_pairs2; /* sums over d of the squares of the per-d CPU
                              seconds and pairs (for standard errors) */
  int truncated;      /* some V_d search hit opts->node_limit */
  int stopped;        /* the callback asked to stop */
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
 */
void dfirst_set_top_root(dfirst_t *df, int on);

/* distinct numbers among the search's vectors */
uint32_t dfirst_num_labels(const dfirst_t *df);

/* is the n-set t (any order) in the diagonal set? */
int dfirst_in_set(const dfirst_t *df, const uint64_t *t);

/*
 * Search V_d for the diagonals i = lo + off, lo + off + stride, ... < hi
 * (stride 1, off 0: all of [lo, hi); stride > 1 with a random off: an
 * unbiased sample, stride x the totals estimates the full loop). If dlog is
 * set, one line per d: "i |V_d| labels nodes pairs vd_s setup_s search_s
 * cpu_s".
 */
dfirst_stats_t dfirst_search(dfirst_t *df, size_t lo, size_t hi, size_t stride,
                             size_t off, const search_opts_t *opts,
                             dsquare_cb cb, void *ctx, FILE *dlog);

#endif // !DFIRST_H
