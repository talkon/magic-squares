#ifndef TYPES_H
#define TYPES_H

#include <stdlib.h>

#include "bitset.h"

#define VEC_SIZE 6

/* number of slots reserved for valid_rows and valid_cols */
#define NUM_VALID_SLOTS 64

/* vectors = rows = cols = (s, p)-sets. starts with max, then descending */
typedef size_t vec_t[VEC_SIZE];

/* sum of vector */
static inline size_t vec_sum(vec_t vec) {
  size_t sum = 0;
  for (size_t i = 0; i < VEC_SIZE; i++) {
    sum += vec[i];
  }
  return sum;
}

/* all vectors with given product */
typedef struct {
  /* all the vectors */
  vec_t *vecs;
  size_t num_vecs;
  /* vectors, indexed by sum */
  struct {
    /* vectors with same sum, in ascending order of max element */
    vec_t *vecs;
    size_t length;
  } *sums;
  size_t num_sums;
} vecs_t;

/* global search state, does not change between search_aux calls */
typedef struct {
  /* label_to_elt[label] = original element */
  size_t *label_to_elt;
  size_t num_labels;
  /* vecs, with labels instead of elts */
  vec_t *vecs;
  size_t num_vecs;
  /* bitarrays of each vec */
  bitset_t *bitarrays;
  /* inters[i][j] = size of vecs[i] & vecs[j] */
  unsigned char **inters;
  /* arena-allocated valid_rows and valid_cols slots */
  size_t *row_idx_slots[NUM_VALID_SLOTS];
  size_t *col_idx_slots[NUM_VALID_SLOTS];
} global_t;

/* search state for rows/cols */
typedef struct row_table {
  /* vecs taken */
  size_t vecs[VEC_SIZE];
  /* number of vecs taken */
  size_t num_vecs;
  /* vecs that could go in */
  size_t *valid;
  size_t num_valid;
} row_table;

/* backtracking state, what changes between search_aux calls */
typedef struct search_table {
  /* table elements that aren't in both a row and a col */
  bitset_t unmatched;
  /* number of tables tested */
  size_t *num_searched;
} search_table;

#endif // !TYPES_H
