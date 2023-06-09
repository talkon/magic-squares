#ifndef ARRANGEMENT_H
#define ARRANGEMENT_H

#include <stdint.h>

#include "bitset.h"

#define max(a, b) (a > b ? a : b)

#define VEC_SIZE 6

#define ROW 0
#define COL 1
#define ROWCOL 2
#define NONE 3

typedef struct vec {
    int elts[VEC_SIZE];
} vec;

typedef struct vecgroupinfo {
    vec *start;
    int length;
} vecgroupinfo;

typedef struct vecgroup {
    vec *vecs;
    vecgroupinfo *infos;
    int numvecs;
    int numsums;
} vecgroup;

typedef struct relabeling {
    vec *vecs;
    int num_vecs;
    int *label_to_elt;
    bitset_t *bitarrays;
} relabeling;

// row state
typedef struct row_table {
    // current rows and cols taken
    int vecs[6];
    // number of rows and cols already taken
    int numvecs;
    // rows and columns that can feasibly go in
    int *valid;
    int num_valid;
} row_table;

// backtracking state
typedef struct search_table {
    row_table rows;
    row_table cols;
    // table elements that aren't in both a row and a col
    bitset_t unmatched;
    int num_unmatched;
    // number of tables tested
    int *num_searched;
    int *row_idx_slots[64];
    int *col_idx_slots[64];
} search_table;

typedef struct indexed_int {
    int elt;
    int freq;
} indexed_int;

int vec_sum(vec r);

vecgroup read_vecs(char *filename);

void search_sum(vecgroup group, int sum);

#endif // !ARRANGEMENT_H
