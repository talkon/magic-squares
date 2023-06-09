#include <stdlib.h>

#include "bitset.h"
#include "types.h"

typedef struct indexed_int {
  int elt;
  int freq;
} indexed_int;

int index_cmp_reverse(const void *a, const void *b) {
  int v1 = (*(indexed_int *)a).freq, v2 = (*(indexed_int *)b).freq;
  return (v1 < v2) - (v1 > v2);
}

int int_reverse(const void *a, const void *b) {
  return (*(int *)a < *(int *)b) - (*(int *)a > *(int *)b);
}

int vec_reverse(const void *a, const void *b) {
  int v1 = (*((vec_t *)a))[0], v2 = (*((vec_t *)b))[0];
  return (v1 < v2) - (v1 > v2);
}

global_t elt_relabeling(vecs_t vecs, size_t sum) {
  global_t g;
  vec_t *vec_arr = vecs.sums[sum].vecs;
  size_t num_vecs = vecs.sums[sum].length;
  indexed_int *elt_freqs = calloc(sum, sizeof(indexed_int));
  for (int i = 0; i < sum; i++) {
    elt_freqs[i].elt = i;
  }
  for (int i = 0; i < num_vecs; i++) {
    for (int j = 0; j < VEC_SIZE; j++) {
      elt_freqs[vec_arr[i][j]].freq++;
    }
  }
  qsort(elt_freqs, sum, sizeof(indexed_int), index_cmp_reverse);
  size_t *elt_to_label = malloc(sum * sizeof(size_t));
  size_t *label_to_elt = malloc(sum * sizeof(size_t));
  for (size_t i = 0; i < sum; i++) {
    elt_to_label[elt_freqs[i].elt] = i;
    label_to_elt[i] = elt_freqs[i].elt;
  }
  vec_t *new_vecs = malloc(num_vecs * sizeof(vec_t));
  for (size_t i = 0; i < num_vecs; i++) {
    for (size_t j = 0; j < VEC_SIZE; j++) {
      new_vecs[i][j] = elt_to_label[vec_arr[i][j]];
    }
    qsort(&(new_vecs[i][0]), VEC_SIZE, sizeof(size_t), int_reverse);
  }
  qsort(new_vecs, num_vecs, sizeof(vec_t), vec_reverse);
  int max_elt = new_vecs[0][0];
  g.bitarrays = malloc(sizeof(bitset_t) * num_vecs);
  for (size_t i = 0; i < num_vecs; i++) {
    g.bitarrays[i] = bitset_create(max_elt + 1);
    for (size_t j = 0; j < VEC_SIZE; j++) {
      bitset_set(g.bitarrays[i], new_vecs[i][j]);
    }
  }
  g.vecs = new_vecs;
  g.num_vecs = num_vecs;
  g.label_to_elt = label_to_elt;
  g.inters_0 = malloc(num_vecs * sizeof(bitset_t));
  g.inters_1 = malloc(num_vecs * sizeof(bitset_t));

  for (int i = 0; i < num_vecs; i++) {
    g.inters_0[i] = bitset_create(num_vecs+64);
    g.inters_1[i] = bitset_create(num_vecs+64);
    for (int j = 0; j < num_vecs; j++) {
      int num_inters = bitset_and_count(g.bitarrays[i], g.bitarrays[j]);
      if(num_inters == 0){
        bitset_set(g.inters_0[i], j);
      }
      if(num_inters == 1){
        bitset_set(g.inters_1[i], j);
      }
    }
  }
  free(elt_to_label);
  free(elt_freqs);
  return g;
}
