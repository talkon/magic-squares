#ifndef NEW_INTERS_H
#define NEW_INTERS_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "types.h"


// Not using bitsets here, don't want extra pointer jumping
typedef struct inter_row{
    int size;
    // number of 1-bit vector intersections in each bitarray that are rows and not cols
    int num_rows;
    uint64_t* inters_0;
    uint64_t* inters_1;
    // for xor computations
    uint64_t* scratch;
    int* vecs;
} inter_row;

typedef struct inter_grid{
    inter_row* rows;
    int* vec_to_idx;
} inter_grid;

uint64_t* inter_bitarray(inter_row row, int index, int num_inters){
    return ((num_inters == 0) ? row.inters_0 : row.inters_1) + index * (1 + (row.size - 1) / 64);
}

/* set ith bit to 1; no bounds checking */
static inline void inter_set(inter_row row, int vec_index, int main_index, int num_inters) {
    uint64_t* array = inter_bitarray(row, vec_index, num_inters);
  size_t shiftedi = main_index / 64;
  array[shiftedi] |= ((uint64_t)1) << (main_index % 64);
}

/* set ith bit to 0; no bounds checking */
static inline void inter_unset(inter_row row, int vec_index, int main_index, int num_inters) {
    uint64_t* array = inter_bitarray(row, vec_index, num_inters);
  size_t shiftedi = main_index / 64;
  array[shiftedi] &= ~(((uint64_t)1) << (main_index % 64));
}


#endif



