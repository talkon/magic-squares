#ifndef ARRANGE_H
#define ARRANGE_H

#include <stdint.h>
#include <stdio.h>

#include "enumerate.h"

/*
 * Arrangement search: find every way to choose n "rows" and n "cols" among a
 * list of vectors (n-sets of numbers) such that rows are pairwise disjoint,
 * cols are pairwise disjoint, and every row meets every col in exactly one
 * number. Each such configuration is a semi-magic square (when the vectors all
 * have the same sum and product). Squares are found exactly once up to
 * transposition and row/col permutation.
 */

#define SQ_MAX_N 8

typedef struct {
  /* rows[i][j] = entry of the square; cols are implied */
  uint64_t rows[SQ_MAX_N][SQ_MAX_N];
  uint64_t cols[SQ_MAX_N][SQ_MAX_N];
  int n;
} square_t;

typedef struct {
  uint64_t nodes;      /* search nodes visited */
  uint64_t squares;    /* semi-magic squares found */
  double seconds;      /* wall time of the search (excluding setup) */
  double setup_seconds;/* relabelling + intersection tables */
  int num_labels;      /* distinct numbers among the vectors */
  int truncated;       /* 1 if the node limit was hit */
  int class_used;      /* the class support ran (opts.class_support with
                          the carried path, the top-label root and every
                          vector in one class) */
  /* r1 sampling (search_opts_t r1_*; all 0 otherwise): the number of first
   * rows searched, and estimates of the full search's totals (sum over the
   * sampled r1 of stride x value) with their standard errors. est_squares
   * is unbiased; est_nodes and est_seconds are unbiased up to the adaptive
   * cross-support state, which carries over from one r1 to the next (so a
   * sampled r1 can take +-1 node from its history: ~1e-5 relative; exact
   * with opts.cross = 0)
   * (stratified simple-random-sampling formula, conservative for the
   * systematic samples taken here); seconds are thread CPU time.
   * r1_strata_nose: the strata with fewer than 2 sampled r1 (out of more),
   * whose variance cannot be estimated and is missing from se_* (which
   * are then too small; with 0 sampled, the stratum is also missing from
   * est_*): use more r1 per stratum, or fewer strata, when it is not 0 */
  uint64_t r1_sampled;
  double est_nodes, est_squares, se_squares, est_seconds, se_seconds;
  int r1_strata_nose;
} search_stats_t;

/* called on every square found; return nonzero to stop the search */
typedef int (*square_cb)(const square_t *sq, void *ctx);

typedef struct {
  uint64_t node_limit; /* 0 = unlimited */
  int forward_check;   /* prune nodes with an uncoverable cell (default 1) */
  int mrv;             /* branch on the most constrained cell (default 1) */
  int min_words;       /* use label bitsets of at least this many 64-bit
                          words (for testing; default 0 = smallest that fits) */
  int gather;          /* when filtering candidate lists with the
                          intersection matrices, use gathers even when N is
                          small enough for permutes (for testing; default 0) */
  int support;         /* drop candidates with a cell that no candidate of the
                          other axis covers (default 1; needs forward_check) */
  int cross;           /* drop candidates that miss every candidate of the
                          other axis through some unmatched cell (default 1:
                          where it pays off, see arrange_core.h; 2: always;
                          needs support) */
  int r1_width;        /* with carried bitsets, search each first row r1
                          with only the words of labels <= its largest label
                          (default 1; 0: every r1 at the width of all the
                          labels, for testing). The same nodes and squares
                          either way. */
  int pretest_min;     /* the children with at least this many vectors
                          placed first run the count and forward checks of
                          their two filters on counts and unions only, and
                          are only filtered for real if they pass (the same
                          nodes and squares; candidate lists that carry their
                          bitsets only; 0 = never; default 5, i.e. from the
                          (2,3)/(3,2) children on) */
  /* r1 sampling (for measurements; default r1_stride 0: every r1, no
   * plan; stride 1: every r1 with the sampling statistics): search only the
   * first rows r1 = r1_offset, r1_offset + r1_stride, ... of the root list
   * (r1 is the square's lowest-index vector, so each r1 is an independent
   * subproblem and stride x the sampled totals is unbiased for the sum),
   * or with r1_nstrata > 0, split the root list into r1_nstrata equal
   * index ranges and sample range h with r1_sstride[h] from r1_soffset[h]
   * (stratified: the early r1, whose universes are the largest, hold most
   * of the time), or with r1_nlist > 0, exactly the r1 of r1_list (for
   * paired measurements; the estimates are then the sampled totals). The
   * root list is [0, N), or with top_root_only its prefix of roots; the r1
   * are indices into it, whatever width each is searched at (r1_width).
   * r1_log: one line per sampled r1 (see arrange.c). */
  uint32_t r1_stride, r1_offset;
  int r1_nstrata;
  uint32_t r1_sstride[8], r1_soffset[8];
  FILE *r1_log;
  const uint32_t *r1_list;
  uint32_t r1_nlist;
  /* for the diagonal-first search (dfirst.c): give these numbers the
   * largest labels, the first the largest (the others by the usual order,
   * ignoring them), and with top_root_only search only the first rows r1
   * through the largest label: exact when every square contains
   * top_numbers[0], as the squares of V_d all contain every number of d */
  const uint64_t *top_numbers;
  int n_top, top_root_only;
  /* class support (the V_d searches of dfirst.c; default 0 = off): when
   * top_numbers holds n numbers (n_top == n) and every vector has exactly
   * one of them, its class (as every vector of V_d has one of d's
   * numbers), the rows of a square have the n classes, a row of class j
   * meets the col of class j at that number and every other col at a
   * number outside top_numbers. So, with 1, the children with one row and
   * one col placed drop the candidates with a cell that neither a placed
   * vector of the other axis nor an other-axis candidate of the right
   * class (its own for its top number, another for the rest) covers, or
   * that meet no other-axis candidate of some class not yet placed there
   * outside the placed vectors, on both axes to the fixpoint (see
   * CLASS_SUP in arrange_core.h). The same squares, fewer nodes; only on
   * the carried path (a no-op on the matrices, or when some vector has
   * not exactly one top number). */
  int class_support;
} search_opts_t;

void search_opts_default(search_opts_t *o);

/* will a search of `labels` distinct numbers with opts carry its bitsets
 * in the candidate lists (1), or use the N x N intersection matrices (0)?
 * (see CARRY_MAX_W in arrange.c) */
int search_carries(uint32_t labels, const search_opts_t *opts);

/*
 * Search vectors l[start, start+count). All vectors must have l->n elements.
 */
search_stats_t search_vectors(const vec_list_t *l, size_t start, size_t count,
                              const search_opts_t *opts, square_cb cb,
                              void *ctx);

#endif // !ARRANGE_H
