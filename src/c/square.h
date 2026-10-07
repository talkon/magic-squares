#ifndef SQUARE_H
#define SQUARE_H

#include <stdint.h>
#include <stdio.h>

#include "arrange.h"

/* the n x n grid of a semi-magic square: grid[i][j] = rows[i] & cols[j] */
typedef struct {
  uint64_t g[SQ_MAX_N][SQ_MAX_N];
  int n;
} grid_t;

/* diagonal scores, as in postprocess.py */
#define DIAG_NONE 0
#define DIAG_S 2
#define DIAG_P 3
#define DIAG_SP 7

typedef struct {
  uint64_t sum, product;
  /* number of traversals (permutations) with magic sum / product / both */
  int s_count, p_count, sp_count;
  /* best score of a pair of diagonals (14 = magic square) */
  int best_score;
  /* rows/cols reordered so the best pair is on the main diagonals */
  grid_t best;
} diag_stats_t;

/* build the grid from a square's rows and cols */
void square_to_grid(const square_t *sq, grid_t *grid);

/*
 * canonical form, independent of row/col order and transposition: sorted rows
 * (each sorted ascending) and sorted cols; the lexicographically smaller of
 * the two is put first. out must hold 2 * n * n values.
 */
void square_canonical(const square_t *sq, uint64_t *out);

/* order-independent hash of a square (hash of its canonical form) */
uint64_t square_hash(const square_t *sq);

/* traversal and diagonal statistics */
void square_diag_stats(const square_t *sq, diag_stats_t *ds);

/* print grid as JSON array of arrays */
void grid_print_json(FILE *fp, const grid_t *g);

#endif // !SQUARE_H
