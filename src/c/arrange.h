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
} search_stats_t;

/* called on every square found; return nonzero to stop the search */
typedef int (*square_cb)(const square_t *sq, void *ctx);

typedef struct {
  uint64_t node_limit; /* 0 = unlimited */
  int forward_check;   /* prune nodes with an uncoverable cell (default 1) */
  int mrv;             /* branch on the most constrained cell (default 1) */
  int min_words;       /* use label bitsets of at least this many 64-bit
                          words (for testing; default 0 = smallest that fits) */
  int gather;          /* filter candidate lists with gathers even when N is
                          small enough for permutes (for testing; default 0) */
} search_opts_t;

void search_opts_default(search_opts_t *o);

/*
 * Search vectors l[start, start+count). All vectors must have l->n elements.
 */
search_stats_t search_vectors(const vec_list_t *l, size_t start, size_t count,
                              const search_opts_t *opts, square_cb cb,
                              void *ctx);

#endif // !ARRANGE_H
