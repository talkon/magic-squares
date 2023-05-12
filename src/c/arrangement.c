#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "labelling.h"
#include "search.h"
#include "types.h"

double get_wall_time() {
  struct timeval time;
  if (gettimeofday(&time, NULL)) {
    //  Handle error
    return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

// IMPORTANT: the arrays pointed to by the returned vecgroup must be freed after
// use
vecs_t read_vecs(char *filename) {
  FILE *fp = fopen(filename, "r");
  vecs_t vecs;
  vecs.num_vecs = 0;
  char buf[1024];
  while (fgets(buf, 1024, fp) != NULL)
    vecs.num_vecs++;

  vecs.vecs = calloc(vecs.num_vecs, sizeof(vec_t));
  fclose(fp);
  fp = fopen(filename, "r");

  vecs.num_vecs = 0;
  while (fgets(buf, 1024, fp) != NULL) {
    size_t *vec = vecs.vecs[vecs.num_vecs];
    sscanf(buf, "%ld %ld %ld %ld %ld %ld %ld", &vec[0], &vec[1], &vec[2], &vec[3],
           &vec[4], &vec[5], &vec[6]);
    vecs.num_vecs++;
  }
  vecs.num_sums = vec_sum(vecs.vecs[vecs.num_vecs - 1]);
  size_t cursum = 0;

  // all infos are initialized to 0 length and null starts
  vecs.sums = calloc(vecs.num_sums + 1, sizeof(*vecs.sums));
  for (size_t i = 0; i < vecs.num_vecs; i++) {
    vec_t *vec = &vecs.vecs[i];
    size_t vecsum = vec_sum(*vec);
    if (vecsum != cursum)
      vecs.sums[vecsum].vecs = vecs.vecs + i;
    vecs.sums[vecsum].length++;
    cursum = vecsum;
  }

  fclose(fp);
  return vecs;
}

// Exhaustively search for magic squares with the given sum, 
// and returns the total number of search nodes
uint64_t search_sum(vecs_t vecs, int sum) {
  size_t total_vecs = vecs.sums[sum].length;
  if (total_vecs == 0)
    return 0;
  printf("sum %d nvecs %ld\n", sum, total_vecs);

  global_t g = elt_relabeling(vecs, sum);
  // printf("inters calculated\n");

  search_table table;
  row_table rows;
  row_table cols;

  rows.num_vecs = 0;
  cols.num_vecs = 0;
  rows.valid = malloc(total_vecs * sizeof(size_t));
  cols.valid = malloc(total_vecs * sizeof(size_t));
  for (size_t i = 0; i < total_vecs; i++) {
    rows.valid[i] = i;
    cols.valid[i] = i;
  }
  rows.num_valid = total_vecs;
  cols.num_valid = total_vecs;

  for (size_t i = 0; i < 3 * VEC_SIZE; i++) {
    g.row_idx_slots[i] = calloc(total_vecs, sizeof(size_t));
    g.col_idx_slots[i] = calloc(total_vecs, sizeof(size_t));
  }

  size_t num_searched = 0;
  table.num_searched = &num_searched;
  table.unmatched = bitset_create(64 * g.bitarrays[0].size);

  search_aux(g, rows, cols, table);
  printf("num searched: %ld\n", *table.num_searched);
  fflush(stdout);

  uint64_t ret = *table.num_searched;

  // freeing
  for (int i = 0; i < g.num_vecs; i++) {
    bitset_free(g.inters_0[i]);
    bitset_free(g.inters_1[i]);
  }
  free(g.inters_0);
  free(g.inters_1);
  free(g.vecs);
  free(g.label_to_elt);
  for (int i = 0; i < g.num_vecs; i++) {
    bitset_free(g.bitarrays[i]);
  }
  free(g.bitarrays);
  free(rows.valid);
  free(cols.valid);
  bitset_free(table.unmatched);
  for (int i = 0; i < 3 * VEC_SIZE; i++) {
    free(g.row_idx_slots[i]);
    free(g.col_idx_slots[i]);
  }

  return ret;
}

void search(char *filename, int sum, int min_sum, int max_sum, uint64_t count_cutoff) {
  vecs_t vecs = read_vecs(filename);
  printf("read\n");
  double before = get_wall_time();
  if (sum != -1) {
    search_sum(vecs, sum);
  } else {
    if (min_sum == -1)
      min_sum = vec_sum(vecs.vecs[0]);
    if (max_sum == -1)
      max_sum = vecs.num_sums;
    uint64_t total_count = 0;
    for (int i = min_sum; i < max_sum; i++) {
      total_count += search_sum(vecs, i);
      if (count_cutoff >= 0 && total_count > count_cutoff) {
        printf("terminated: total count %ld exceeds cutoff %ld\n", total_count, count_cutoff);
        break;
      }
    }
  }
  double difference = get_wall_time() - before;
  printf("completed in %.5f secs\n", (double)difference);
  free(vecs.vecs);
  free(vecs.sums);
}

int main(int argc, char *argv[]) {
  char *filename;
  int sum = -1;
  int min_sum = -1;
  int max_sum = -1;
  uint64_t count_cutoff = -1;
  int opt;

  while (1) {
    static struct option long_options[] = {
        {"file", required_argument, 0, 'f'},
        {"sum", required_argument, 0, 's'},
        {"min-sum", required_argument, 0, 'n'},
        {"max-sum", required_argument, 0, 'x'},
        {"count", required_argument, 0, 'c'},
        {0, 0, 0, 0}};

    opt = getopt_long(argc, argv, "f:s:n:x:", long_options, NULL);
    if (opt == -1)
      break;

    switch (opt) {
    case 'f':
      filename = optarg;
      break;
    case 's':
      sum = strtol(optarg, NULL, 10);
      break;
    case 'n':
      min_sum = strtol(optarg, NULL, 10);
      break;
    case 'x':
      max_sum = strtol(optarg, NULL, 10);
      break;
    case 'c':
      count_cutoff = (uint64_t) strtoll(optarg, NULL, 10);
    case '?':
      break;
    default:
      abort();
    }
  }

  search(filename, sum, min_sum, max_sum, count_cutoff);
}
