#include <getopt.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "arrangement.h"

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


int vec_sum(vec r){
    int s = 0;
    for(int i = 0; i < VEC_SIZE; i++){
        s += r.elts[i];
    }
    return s;
}

// IMPORTANT: the arrays pointed to by the returned vecgroup must be freed after use
vecgroup read_vecs(char *filename){
    FILE *fp = fopen(filename,"r");
    int numvecs = 0;
    char buf[1024];
    while(fgets(buf, 1024, fp) != NULL){
        numvecs++;
    }
    vec *vec_arr = malloc(sizeof(vec) * numvecs);
    fp = fopen(filename,"r");
    int veccount = 0;
    int maxsum;
    while(fgets(buf, 1024, fp) != NULL){
        vec slot;
        sscanf(buf, "%d %d %d %d %d %d", &slot.elts[0], &slot.elts[1], &slot.elts[2], &slot.elts[3], &slot.elts[4], &slot.elts[5]);
        vec_arr[veccount] = slot;
        veccount++;
        maxsum = vec_sum(slot);
    }
    int cursum = 0;
    // all infos are initialized to 0 length and null starts
    vecgroupinfo* infos = calloc(maxsum + 1, sizeof(vecgroupinfo));
    for(int i = 0; i < numvecs; i++){
        vec curvec = vec_arr[i];
        if(vec_sum(curvec) != cursum){
            infos[vec_sum(curvec)].start = vec_arr + i;
        }
        infos[vec_sum(curvec)].length++;
        cursum = vec_sum(curvec);
    }
    vecgroup r;
    r.vecs = vec_arr;
    r.infos = infos;
    r.numvecs = numvecs;
    r.numsums = maxsum + 1;
    return r;

}

int index_cmp_reverse(const void *a, const void *b){
    int v1 = (*(indexed_int*)a).freq, v2 = (*(indexed_int*)b).freq;
    return (v1 < v2) - (v1 > v2);
}

int int_reverse(const void *a, const void *b){
    return (*(int*)a < *(int*)b) - (*(int*)a > *(int*)b);
}

int vec_reverse(const void *a, const void *b){
    int v1 = (*((vec*)a)).elts[0], v2 = (*((vec*)b)).elts[0];
    return (v1 < v2) - (v1 > v2);
}

// free relabeling pointers after use
relabeling elt_relabeling(vecgroup group, int sum){
    vec *vec_arr = group.infos[sum].start;
    int total_vecs = group.infos[sum].length;
    relabeling r;
    indexed_int *elt_freqs = calloc(sum, sizeof(indexed_int));
    for(int i = 0; i < sum; i++){
        elt_freqs[i].elt = i;
    }
    for(int i = 0; i < total_vecs; i++){
        for(int j = 0; j < VEC_SIZE; j++){
            elt_freqs[vec_arr[i].elts[j]].freq++;
        }
    }
    qsort(elt_freqs, sum, sizeof(indexed_int), index_cmp_reverse);
    int *elt_to_label = malloc(sum * sizeof(int));
    int *label_to_elt = malloc(sum * sizeof(int));
    for(int i = 0; i < sum; i++){
        elt_to_label[elt_freqs[i].elt] = i;
        label_to_elt[i] = elt_freqs[i].elt;
    }
    vec *new_vecs = malloc(total_vecs * sizeof(vec));
    for(int i = 0; i < total_vecs; i++){
        for(int j = 0; j < VEC_SIZE; j++){
            new_vecs[i].elts[j] = elt_to_label[vec_arr[i].elts[j]];
        }
        qsort(&(new_vecs[i].elts[0]), VEC_SIZE, sizeof(int), int_reverse);
    }
    qsort(new_vecs, total_vecs, sizeof(vec), vec_reverse);
    int max_elt = new_vecs[0].elts[0];
    r.bitarrays = malloc(sizeof(bitset_t) * total_vecs);
    for(int i = 0; i < total_vecs; i++){
        r.bitarrays[i] = bitset_create(max_elt);
        for(int j = 0; j < VEC_SIZE; j++){
            bitset_set(r.bitarrays[i], new_vecs[i].elts[j]);
        }
    }
    r.vecs = new_vecs;
    r.num_vecs = total_vecs;
    r.label_to_elt = label_to_elt;
    free(elt_to_label);
    free(elt_freqs);
    return r;
}

unsigned char** intersections(relabeling r, int numvecs){
    unsigned char** inters = malloc(numvecs * sizeof(char*));
    for(int i = 0; i < numvecs; i++){
        inters[i] = malloc(numvecs * sizeof(char));
        for(int j = 0; j < numvecs; j++){
            inters[i][j] = bitset_and_count(r.bitarrays[i], r.bitarrays[j]);
        }
    }
    return inters;
}

bool any_match(vec v, int target){
  for(int i = 0; i < VEC_SIZE; i++){
    if(v.elts[i] == target){
      return true;
    }
  }
  return false;
}

int unmatched_max(uint64_t *unmatched, int num_64s){
    for(int i = num_64s - 1; i >= 0; i--){
        if(unmatched[i] != 0){
            return 64 * i + 63 - __builtin_clzll(unmatched[i]);
        }
    }
    return -1;
}

void record(search_table table, relabeling r){
    *table.num_searched += 1;
    if(table.rows.numvecs * table.cols.numvecs == VEC_SIZE * VEC_SIZE){
      //return;
        printf("solution found\n");
        for(int i = 0; i < VEC_SIZE; i++){
            printf("%d    ", table.rows.vecs[i]);
            for(int j = 0; j < VEC_SIZE; j++){
                printf("%d ", r.label_to_elt[r.vecs[table.rows.vecs[i]].elts[j]]);
            }
            printf("\n");
        }
        printf("\n");
        for(int i = 0; i < VEC_SIZE; i++){
            printf("%d    ", table.cols.vecs[i]);
            for(int j = 0; j < VEC_SIZE; j++){
                printf("%d ", r.label_to_elt[r.vecs[table.cols.vecs[i]].elts[j]]);
            }
            printf("\n");
        }
    }

}

int fill_valids(int *rows_valid, int last_row, int inter_val, unsigned char** inters, relabeling r, int numvalids, int *oldrows, int minvec){
    int row_ptr = 0;
    int last_max = r.vecs[last_row].elts[0];
    for(int i = 0; i < numvalids; i++){
        int test_row = oldrows[i];
        bool is_valid_row = (inters[last_row][test_row] == inter_val) && test_row >= minvec;
        int is_valid = (int) is_valid_row;
        rows_valid[row_ptr] = test_row * is_valid + rows_valid[row_ptr] * (1 - is_valid);
        row_ptr += is_valid;
    }
    return row_ptr;
}


static inline search_table search_prep(relabeling r, search_table new_table, unsigned char** inters, int to_add, int new_minvec, int last_vec) {
    if(new_table.rows.numvecs > 1 || new_table.cols.numvecs > 1){
      new_minvec = 0;
    }
    int *rows_valid = new_table.row_idx_slots[new_table.rows.numvecs+new_table.cols.numvecs+1];//malloc(sizeof(int) * table.rows.num_valid);
    int *cols_valid = new_table.col_idx_slots[new_table.rows.numvecs+new_table.cols.numvecs+1];//malloc(sizeof(int) * table.cols.num_valid);
    new_table.rows.num_valid = fill_valids(rows_valid, last_vec, to_add, inters, r, new_table.rows.num_valid, new_table.rows.valid, new_minvec);//row_ptr;
    new_table.cols.num_valid = fill_valids(cols_valid, last_vec, 1-to_add, inters, r, new_table.cols.num_valid, new_table.cols.valid, new_minvec);//col_ptr;
    new_table.rows.valid = rows_valid;
    new_table.cols.valid = cols_valid;
    return new_table;
}


void search_aux(relabeling r, search_table table, unsigned char** inters, int last_appended, int minvec){
    record(table, r);

    if(!((table.rows.numvecs + table.rows.num_valid >= 6) && (table.cols.numvecs + table.cols.num_valid >= 6))){
        return;
    }

    int max_unmatched = bitset_maximum(table.unmatched);

    for (int to_add = ROW; to_add <= COL; to_add++) {
        if (to_add == COL && table.rows.numvecs == 0 && table.cols.numvecs == 0) break;
        row_table active = to_add == ROW ? table.rows : table.cols;
        for (int i = 0; i < active.num_valid; i++) {
            if (active.numvecs == VEC_SIZE) break;
            int max_elt = r.vecs[active.valid[i]].elts[0];
            if(max_elt < max_unmatched){
                break;
            }
            if(!any_match(r.vecs[active.valid[i]], max_unmatched) && (to_add == COL || table.rows.numvecs > 0)){
              continue;
            }
            search_table new_table = table;
            row_table* new_active = to_add == ROW ? &new_table.rows : &new_table.cols;
            new_active->vecs[active.numvecs] = active.valid[i];
            new_active->numvecs++;
            int roww = active.valid[i];
            bitset_inplace_xor(new_table.unmatched, r.bitarrays[roww]);
            new_table = search_prep(r, new_table, inters, to_add, active.valid[i], active.valid[i]);
            search_aux(r, new_table, inters, to_add, active.valid[i]);
            bitset_inplace_xor(new_table.unmatched, r.bitarrays[roww]);
        }
    }

    return;
}

void search_sum(vecgroup group, int sum){
    vec *start = group.infos[sum].start;
    int total_vecs = group.infos[sum].length;
    if(total_vecs == 0){
        return;
    }
    printf("sum %d nvecs %d\n", sum, total_vecs);
    relabeling r = elt_relabeling(group, sum);
    unsigned char** inters = intersections(r, total_vecs);
    //printf("inters calculated\n");
    
    int *valid_row_indices = malloc(total_vecs * sizeof(int));
    int *valid_col_indices = malloc(total_vecs * sizeof(int));
    for(int i = 0; i < total_vecs; i++){
        valid_row_indices[i] = i;
        valid_col_indices[i] = i;
    }
    search_table table;
    for(int i = 0; i < 3*VEC_SIZE; i++){
        table.row_idx_slots[i] = malloc(total_vecs * sizeof(int));
        table.col_idx_slots[i] = malloc(total_vecs * sizeof(int));
    }
    /*int *valid_row_indices = table.row_idx_slots[0];
    int *valid_col_indices = table.col_idx_slots[0];
    for(int i = 0; i < total_vecs; i++){
        valid_row_indices[i] = i;
        valid_col_indices[i] = i;
    }*/
    table.rows.numvecs = 0;
    table.cols.numvecs = 0;
    table.rows.num_valid = total_vecs;
    table.cols.num_valid = total_vecs;
    table.rows.valid = valid_row_indices;
    table.cols.valid = valid_col_indices;
    table.num_unmatched = 0;
    int num_searched = 0;
    table.num_searched = &num_searched;
    table.unmatched = bitset_create(64 * r.bitarrays[0].size);
    
    search_aux(r, table, inters, NONE, 0);
    printf("num searched: %d\n", *table.num_searched);
    free(inters);
    free(r.vecs);
    free(r.label_to_elt);
    for(int i = 0; i < r.num_vecs; i++){
        bitset_free(r.bitarrays[i]);
    }
    free(r.bitarrays);
    free(valid_row_indices);
    free(valid_col_indices);
    bitset_free(table.unmatched);
    for(int i = 0; i < 3*VEC_SIZE; i++){
        free(table.row_idx_slots[i]);// = malloc(total_vecs * sizeof(int));
        free(table.col_idx_slots[i]);// = malloc(total_vecs * sizeof(int));
    }
}

void search(char *filename, int sum, int min_sum, int max_sum) {
    vecgroup group = read_vecs(filename);
    printf("read\n");
    double before = get_wall_time();
    if (sum != -1) {
        search_sum(group, sum);
    } else {
        if (min_sum == -1)
            min_sum = vec_sum(group.vecs[0]);
        if (max_sum == -1)
            max_sum = group.numsums + 1;
        for (int i = min_sum; i < max_sum; i++) {
            search_sum(group, i);
        }
    }
    double difference = get_wall_time() - before;
    printf("completed in %.5f secs\n", (double)difference);
    free(group.vecs);
    free(group.infos);
}

int main(int argc, char *argv[]) {
    char *filename;
    int sum = -1;
    int min_sum = -1;
    int max_sum = -1;
    int opt;

    while (1) {
        static struct option long_options[] = {
                {"file", required_argument, 0, 'f'},
                {"sum", required_argument, 0, 's'},
                {"min-sum", required_argument, 0, 'n'},
                {"max-sum", required_argument, 0, 'x'},
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
        case '?':
            break;
        default:
            abort();
        }
    }

    search(filename, sum, min_sum, max_sum);
}
