#include "arrangement.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdint.h>
#include "arrangement.h"
#include <time.h>

int sum(vec r) {
  int s = 0;
  for (int i = 0; i < VEC_SIZE; i++) {
    s += r.elts[i];
  }
  return s;
}

// IMPORTANT: the arrays pointed to by the returned vecgroup must be freed after
// use
vecgroup read_vecs(char *filename) {
  FILE *fp = fopen(filename, "r");
  int numvecs = 0;
  char buf[1024];
  while (fgets(buf, 1024, fp) != NULL) {
    numvecs++;
  }
  vec *vec_arr = malloc(sizeof(vec) * numvecs);
  fp = fopen(filename, "r");
  int veccount = 0;
  int maxsum;
  while (fgets(buf, 1024, fp) != NULL) {
    vec slot;
    sscanf(buf, "%d %d %d %d %d %d", &slot.elts[0], &slot.elts[1],
           &slot.elts[2], &slot.elts[3], &slot.elts[4], &slot.elts[5]);
    vec_arr[veccount] = slot;
    veccount++;
    maxsum = sum(slot);
  }
  int cursum = 0;
  // all infos are initialized to 0 length and null starts
  vecgroupinfo *infos = calloc(maxsum + 1, sizeof(vecgroupinfo));
  for (int i = 0; i < numvecs; i++) {
    vec curvec = vec_arr[i];
    if (sum(curvec) != cursum) {
      infos[sum(curvec)].start = vec_arr + i;
    }
    infos[sum(curvec)].length++;
    cursum = sum(curvec);
  }
  vecgroup r;
  r.vecs = vec_arr;
  r.infos = infos;
  r.numvecs = numvecs;
  r.numsums = maxsum + 1;
  return r;
}

int index_cmp_reverse(const void *a, const void *b) {
  int v1 = (*(indexed_int *)a).freq, v2 = (*(indexed_int *)b).freq;
  return (v1 < v2) - (v1 > v2);
}

int int_reverse(const void *a, const void *b) {
  return (*(int *)a < *(int *)b) - (*(int *)a > *(int *)b);
}

int vec_reverse(const void *a, const void *b) {
  int v1 = (*((vec *)a)).elts[0], v2 = (*((vec *)b)).elts[0];
  return (v1 < v2) - (v1 > v2);
}

// free relabeling pointers after use
relabeling elt_relabeling(vecgroup group, int sum) {
  vec *vec_arr = group.infos[sum].start;
  int total_vecs = group.infos[sum].length;
  relabeling r;
  indexed_int *elt_freqs = calloc(sum, sizeof(indexed_int));
  for (int i = 0; i < sum; i++) {
    elt_freqs[i].elt = i;
  }
  for (int i = 0; i < total_vecs; i++) {
    for (int j = 0; j < VEC_SIZE; j++) {
      elt_freqs[vec_arr[i].elts[j]].freq++;
    }
  }
  qsort(elt_freqs, sum, sizeof(indexed_int), index_cmp_reverse);
  int *elt_to_label = malloc(sum * sizeof(int));
  int *label_to_elt = malloc(sum * sizeof(int));
  for (int i = 0; i < sum; i++) {
    elt_to_label[elt_freqs[i].elt] = i;
    label_to_elt[i] = elt_freqs[i].elt;
  }
  vec *new_vecs = malloc(total_vecs * sizeof(vec));
  for (int i = 0; i < total_vecs; i++) {
    for (int j = 0; j < VEC_SIZE; j++) {
      new_vecs[i].elts[j] = elt_to_label[vec_arr[i].elts[j]];
    }
    qsort(&(new_vecs[i].elts[0]), VEC_SIZE, sizeof(int), int_reverse);
  }
  qsort(new_vecs, total_vecs, sizeof(vec), vec_reverse);
    int max_elt = new_vecs[0].elts[0];
    r.bitarrays = malloc(sizeof(uint64_t*) * total_vecs);
    printf("%d\n", max_elt);
    for(int i = 0; i < total_vecs; i++){
        r.bitarrays[i] = calloc((max_elt / 64) + 1, sizeof(uint64_t));
        for(int j = 0; j < VEC_SIZE; j++){
            *((r.bitarrays[i]) + (new_vecs[i].elts[j] / 64)) |= (uint64_t)1 << (new_vecs[i].elts[j] % 64);
        }
    }
  r.vecs = new_vecs;
  r.num_vecs = total_vecs;
  r.label_to_elt = label_to_elt;
    r.bitarray_size = (max_elt / 64) + 1;
  free(elt_to_label);
  free(elt_freqs);
  return r;
}

int num_inters(vec a, vec b) {
  int p1 = 0;
  int p2 = 0;
  int tot = 0;
  while (p1 < VEC_SIZE && p2 < VEC_SIZE) {
    if (a.elts[p1] > b.elts[p2]) {
      p1++;
    } else {
      if (a.elts[p1] == b.elts[p2]) {
        p1++;
        tot++;
      }
      p2++;
    }
  }
  return tot;
}

unsigned char **intersections(vec *vecs, int numvecs) {
  unsigned char **inters = malloc(numvecs * sizeof(char *));
  for (int i = 0; i < numvecs; i++) {
    inters[i] = malloc(numvecs * sizeof(char));
    for (int j = 0; j < numvecs; j++) {
      inters[i][j] = num_inters(vecs[i], vecs[j]);
    }
  }
  return inters;
}

int max_unmatched(search_table table) {
  int m = 0;
  for (int i = 0; i < table.num_unmatched; i++) {
    if (m < table.unmatched[i]) {
      m = table.unmatched[i];
    }
  }
  return m;
}

void record(search_table table, relabeling r){
    *table.num_searched += 1;
    if(table.numrows * table.numcols == VEC_SIZE * VEC_SIZE){
        printf("solution found\n");
        for(int i = 0; i < VEC_SIZE; i++){
            for(int j = 0; j < VEC_SIZE; j++){
                printf("%d ", r.label_to_elt[r.vecs[table.rows[i]].elts[j]]);
            }
            printf("\n");
        }
    }
}

int fill_valids(int *valid_rows, int last_row, int inter_val,
                unsigned char **inters, relabeling r, int minvec, int numvalids,
                int *oldrows) {
  int row_ptr = 0;
  int last_max = r.vecs[last_row].elts[0];
  for (int i = 0; i < numvalids; i++) {
    int test_row = oldrows[i];
    bool is_valid_row =
        (inters[last_row][test_row] == inter_val) && test_row > minvec;
    if (is_valid_row) {
      valid_rows[row_ptr] = test_row;
      row_ptr++;
    }
  }
  return row_ptr;
}

void search_aux(relabeling r, search_table table, unsigned char **inters,
                int last_appended, int minvec) {
  search_table prev_table = table;
  int max_row_elt = 0;
  if (table.num_valid_rows > 0) {
    max_row_elt = r.vecs[table.valid_rows[0]].elts[0];
  }
  int max_col_elt = 0;
  if (table.num_valid_cols > 0) {
    max_col_elt = r.vecs[table.valid_cols[0]].elts[0];
  }
  int max_elt = (max_row_elt > max_col_elt) ? max_row_elt : max_col_elt;
  /*if(max_unmatched(table) > max_elt){
      //printf("a");

      return;
  }*/
  if (!((table.numrows + table.num_valid_rows >= 6) &&
        (table.numcols + table.num_valid_cols >= 6))) {
    return;
  }
  record(table, r);
  int last_row = 0;
  if (table.numrows > 0) {
    last_row = table.rows[table.numrows - 1];
  }
  int last_col = 0;
  if (table.numcols > 0) {
    last_col = table.cols[table.numcols - 1];
  }
  int *valid_rows = malloc(sizeof(int) * table.num_valid_rows);
  int *valid_cols = malloc(sizeof(int) * table.num_valid_cols);

  // printf("w");
  // printf("%d ", minvec);
  //   printf("%d.", table.num_valid_rows);
  //     printf("%d ", table.num_valid_cols);
  if (last_appended == ROW || last_appended == ROWCOL) {
    /*int row_ptr = 0;
    int last_max = r.vecs[last_row].elts[0];
    for(int i = 0; i < table.num_valid_rows; i++){
        int test_row = table.valid_rows[i];
        bool is_valid_row = (inters[last_row][test_row] == 0) && test_row >
    minvec; if(is_valid_row){ valid_rows[row_ptr] = test_row; row_ptr++;
        }
    }*/

    table.num_valid_rows =
        fill_valids(valid_rows, last_row, 0, inters, r, minvec,
                    table.num_valid_rows, table.valid_rows); // row_ptr;

    /*int col_ptr = 0;
    for(int i = 0; i < table.num_valid_cols; i++){
        int test_col = table.valid_cols[i];
        bool is_valid_col = (inters[last_row][test_col] == 1) && test_col >
    minvec; if(is_valid_col){ valid_cols[col_ptr] = test_col; col_ptr++;
        }
    }*/

    table.num_valid_cols =
        fill_valids(valid_cols, last_row, 1, inters, r, minvec,
                    table.num_valid_cols, table.valid_cols); // col_ptr;
    table.valid_rows = valid_rows;
    table.valid_cols = valid_cols;
  }
  if (last_appended == COL || last_appended == ROWCOL) {
    /*int row_ptr = 0;
    int last_max = r.vecs[last_col].elts[0];
    for(int i = 0; i < table.num_valid_rows; i++){
        int test_row = table.valid_rows[i];
        bool is_valid_row = (inters[last_col][test_row] == 1) && test_row >
    minvec; if(is_valid_row){ valid_rows[row_ptr] = test_row; row_ptr++;
        }
    }*/
    table.num_valid_rows =
        fill_valids(valid_rows, last_col, 1, inters, r, minvec,
                    table.num_valid_rows, table.valid_rows); // row_ptr;

    /*int col_ptr = 0;
    for(int i = 0; i < table.num_valid_cols; i++){
        int test_col = table.valid_cols[i];
        bool is_valid_col = (inters[last_col][test_col] == 0) && test_col >
    minvec; if(is_valid_col){ valid_cols[col_ptr] = test_col; col_ptr++;
        }
    }*/
    table.num_valid_cols =
        fill_valids(valid_cols, last_col, 0, inters, r, minvec,
                    table.num_valid_cols, table.valid_cols); // col_ptr;
    table.valid_rows = valid_rows;
    table.valid_cols = valid_cols;
  }
  //    printf("%d.", table.num_valid_rows);
  //    printf("%d !", table.num_valid_cols);

    int max_unmatched = unmatched_max(table.unmatched, r.bitarray_size);


    for(int i = 0; i < table.num_valid_rows; i++){
        if(table.numrows == VEC_SIZE) break;
        int max_elt = r.vecs[table.valid_rows[i]].elts[0];
        if(max_elt < max_unmatched){
            break;
        }
        else{
            
            table.rows[table.numrows] = table.valid_rows[i];
            table.numrows++;
            
            int roww = table.valid_rows[i];
            for(int j = 0; j < r.bitarray_size; j++){
                table.unmatched[j] ^= r.bitarrays[roww][j];
            }
            search_aux(r, table, inters, ROW, table.valid_rows[i]);
            for(int j = 0; j < r.bitarray_size; j++){
                table.unmatched[j] ^= r.bitarrays[roww][j];
            }
            
            //table.unmatched = prev_unmatched;
            //table.num_unmatched = prev_num_unmatched;
            table.numrows--;
            //free(unmatched);
        }
    }

  for (int i = 0; i < table.num_valid_cols; i++) {
    if (table.numcols == VEC_SIZE)
      break;
    if (table.numrows == 0 && table.numcols == 0) {
      break;
    }
    int max_elt = r.vecs[table.valid_cols[i]].elts[0];
    if (max_elt < max_unmatched) {
      break;
    } else {

            table.cols[table.numcols] = table.valid_cols[i];
            table.numcols++;
         
            int coll = table.valid_cols[i];
            for(int j = 0; j < r.bitarray_size; j++){
                table.unmatched[j] ^= r.bitarrays[coll][j];
            }
            search_aux(r, table, inters, COL, table.valid_cols[i]);
            for(int j = 0; j < r.bitarray_size; j++){
                table.unmatched[j] ^= r.bitarrays[coll][j];
            }


            table.numcols--;


        }
    }
    free(valid_rows);
    free(valid_cols);
    return;
}

void search_sum(vecgroup group, int sum) {
  vec *start = group.infos[sum].start;
  int total_vecs = group.infos[sum].length;
  if (total_vecs == 0) {
    return;
  }
  relabeling r = elt_relabeling(group, sum);
  unsigned char **inters = intersections(r.vecs, total_vecs);
  printf("inters calculated\n");
  int *valid_row_indices = malloc(total_vecs * sizeof(int));
  int *valid_col_indices = malloc(total_vecs * sizeof(int));
  for (int i = 0; i < total_vecs; i++) {
    valid_row_indices[i] = i;
    valid_col_indices[i] = i;
  }
  search_table table;
  table.numrows = 0;
  table.numcols = 0;
  table.num_valid_rows = total_vecs;
  table.num_valid_cols = total_vecs;
  table.valid_rows = valid_row_indices;
  table.valid_cols = valid_col_indices;
  table.num_unmatched = 0;
  int num_searched = 0;
  table.num_searched = &num_searched;
    table.unmatched = calloc(r.bitarray_size, sizeof(uint64_t));
    
  search_aux(r, table, inters, NONE, 0);
  printf("num searched: %d\n", *table.num_searched);
  free(inters);
  free(r.vecs);
  free(r.label_to_elt);
    for(int i = 0; i < r.num_vecs; i++){
        free(r.bitarrays[i]);
    }
    free(r.bitarrays);
  free(valid_row_indices);
  free(valid_col_indices);
    free(table.unmatched);
}

void test_print_vecs(char *filename) {
  vecgroup group = read_vecs(filename);
  printf("read\n");
  clock_t before = clock();
  search_sum(group, 327);
  clock_t difference = clock() - before;
  printf("completed in %d msecs\n", difference);
  int msec = difference / CLOCKS_PER_SEC;
  free(group.vecs);
  free(group.infos);
}

int main(int argc, char *argv[]) {
  int opt;
  char filename[1024];
  while ((opt = getopt(argc, argv, "f:")) != -1) {
    test_print_vecs(optarg);
  }
}
