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
    if(table.numrows * table.numcols == VEC_SIZE * VEC_SIZE){
      //return;
        printf("solution found\n");
        for(int i = 0; i < VEC_SIZE; i++){
            printf("%d    ", table.rows[i]);
            for(int j = 0; j < VEC_SIZE; j++){
                printf("%d ", r.label_to_elt[r.vecs[table.rows[i]].elts[j]]);
            }
            printf("\n");
        }
        printf("\n");
        for(int i = 0; i < VEC_SIZE; i++){
            printf("%d    ", table.cols[i]);
            for(int j = 0; j < VEC_SIZE; j++){
                printf("%d ", r.label_to_elt[r.vecs[table.cols[i]].elts[j]]);
            }
            printf("\n");
        }
    }

}

int fill_valids(int *valid_rows, int last_row, int inter_val, unsigned char** inters, relabeling r, int numvalids, int *oldrows, int minvec){
    int row_ptr = 0;
    int last_max = r.vecs[last_row].elts[0];
    for(int i = 0; i < numvalids; i++){
        int test_row = oldrows[i];
        bool is_valid_row = (inters[last_row][test_row] == inter_val) && test_row >= minvec;
        int is_valid = (int) is_valid_row;
        valid_rows[row_ptr] = test_row * is_valid + valid_rows[row_ptr] * (1 - is_valid);
        row_ptr += is_valid;
        continue;
        if(is_valid_row){
            valid_rows[row_ptr] = test_row;
            row_ptr++;
        }
    }
    return row_ptr;
}


void search_aux(relabeling r, search_table table, unsigned char** inters, int last_appended, int minvec){
    search_table prev_table = table;
    if(!((table.numrows + table.num_valid_rows >= 6) && (table.numcols + table.num_valid_cols >= 6))){
        return;
    }
    record(table, r);
    int last_row = 0;
    if(table.numrows > 0){
        last_row = table.rows[table.numrows - 1];
    }
    int last_col = 0;
    if(table.numcols > 0){
        last_col = table.cols[table.numcols - 1];
    }
    int *valid_rows = table.row_idx_slots[table.numrows+table.numcols+1];//malloc(sizeof(int) * table.num_valid_rows);
    int *valid_cols = table.col_idx_slots[table.numrows+table.numcols+1];//malloc(sizeof(int) * table.num_valid_cols);

    if(table.numrows > 1 || table.numcols > 1){
      minvec = 0;
    }
    if(last_appended == ROW || last_appended == ROWCOL){

      
        table.num_valid_rows = fill_valids(valid_rows, last_row, 0, inters, r, table.num_valid_rows, table.valid_rows, minvec);//row_ptr;

        table.num_valid_cols = fill_valids(valid_cols, last_row, 1, inters, r, table.num_valid_cols, table.valid_cols, minvec);//col_ptr;
        table.valid_rows = valid_rows;
        table.valid_cols = valid_cols;
    }
    if(last_appended == COL || last_appended == ROWCOL){

        table.num_valid_rows = fill_valids(valid_rows, last_col, 1, inters, r, table.num_valid_rows, table.valid_rows, minvec);//row_ptr;

        table.num_valid_cols = fill_valids(valid_cols, last_col, 0, inters, r, table.num_valid_cols, table.valid_cols, minvec);//col_ptr;
        table.valid_rows = valid_rows;
        table.valid_cols = valid_cols;
    }
    if(!((table.numrows + table.num_valid_rows >= 6) && (table.numcols + table.num_valid_cols >= 6))){
        return;
    }


    int max_unmatched = bitset_maximum(table.unmatched);
   

    int row_ptr = 0;
    int col_ptr = 0;
    int prev_max_elt = 0;
    int col_start = 0;
    
    
    for(row_ptr = 0; row_ptr < table.num_valid_rows; row_ptr++){
      break;
      if(table.numrows + table.numcols < 2)
        break;
        
        int max_elt = r.vecs[table.valid_rows[row_ptr]].elts[0];
        if(max_elt <= max_unmatched){
            break;
        }
        if(max_elt != prev_max_elt){
            prev_max_elt = max_elt;
            while(col_ptr < table.num_valid_cols && r.vecs[table.valid_cols[col_ptr]].elts[0] > max_elt){
                col_ptr++;
            }
            if(col_ptr == table.num_valid_cols) break;
            col_start = col_ptr;
            if(r.vecs[table.valid_cols[col_ptr]].elts[0] != max_elt){
                continue;
            }
            while(col_ptr < table.num_valid_cols && r.vecs[table.valid_cols[col_ptr]].elts[0] == max_elt){
                col_ptr++;
            }

        }
        int roww = table.valid_rows[row_ptr];

        
        if(col_ptr == col_start){
            continue;
        }

        int z = 0;
        if(!any_match(r.vecs[roww], max_unmatched)){
          for(int col_idx = col_start; col_idx < col_ptr; col_idx++){
            int coll = table.valid_cols[col_idx];
            if(any_match(r.vecs[coll], max_unmatched)) z++;
          }
          if(z==0) continue;
        }
        search_table new_table = table;

        new_table.rows[table.numrows] = roww;
        new_table.numrows++;
        /*new_table.valid_rows += (row_ptr + 1);
        new_table.num_valid_rows -= (row_ptr + 1);
        new_table.valid_cols += (col_start);
        new_table.num_valid_cols -= (col_start);*/


        int *new_valid_rows = table.row_idx_slots[new_table.numrows+new_table.numcols+1];
        int *new_valid_cols = table.col_idx_slots[new_table.numrows+new_table.numcols+1];
        
        new_table.num_valid_rows = fill_valids(new_valid_rows, roww, 0, inters, r, new_table.num_valid_rows, new_table.valid_rows, minvec);
        new_table.num_valid_cols = fill_valids(new_valid_cols, roww, 1, inters, r, new_table.num_valid_cols, new_table.valid_cols, minvec);
        new_table.valid_rows = new_valid_rows;
        new_table.valid_cols = new_valid_cols;
        if(!((new_table.numrows + new_table.num_valid_rows >= 6) && (new_table.numcols + new_table.num_valid_cols >= 6))){
            continue;
        }
        

        bitset_inplace_xor(new_table.unmatched, r.bitarrays[roww]);
        
        
        for(int col_idx = col_start; col_idx < col_ptr; col_idx++){
            int coll = table.valid_cols[col_idx];

            if(table.numrows == 0 && table.numcols == 0 && coll < roww){
                continue;
            }
            if(any_match(r.vecs[coll], max_unmatched)) z++;
            if((!any_match(r.vecs[roww], max_unmatched) && !any_match(r.vecs[coll], max_unmatched)) || inters[roww][coll] != 1){
                continue;
            }
            
            search_table new_table_2 = new_table;
            new_table_2.cols[table.numcols] = coll;
            new_table_2.numcols++;

            bitset_inplace_xor(new_table_2.unmatched, r.bitarrays[coll]);
            search_aux(r, new_table_2, inters, COL, max(roww, coll));
            bitset_inplace_xor(new_table_2.unmatched, r.bitarrays[coll]);

        }
        //if(!any_match(r.vecs[roww], max_unmatched)) 
        //  printf("%d\n", z);
        bitset_inplace_xor(new_table.unmatched, r.bitarrays[roww]);

        
    }

    int smallest_col = 0;
    int d = 0;
    int e = 0;
    for(int i = row_ptr; i < table.num_valid_rows; i++){
        if(table.numrows == VEC_SIZE) break;
        int max_elt = r.vecs[table.valid_rows[i]].elts[0];
        if(max_elt < max_unmatched){
            break;
        }
        
        else{
         /* e++;
          int z = 0;
        for(int j = 0; j < 6; j++){
          if(r.vecs[table.valid_rows[i]].elts[j]==max_unmatched){
            d++;
            z++;
            
          }
        }*/
            if(!any_match(r.vecs[table.valid_rows[i]], max_unmatched) && table.numrows > 0){
              continue;
            }
            /*if(max_elt == max_unmatched){
              printf("y\n");
            }
            else{
              printf("n\n");
            }*/
            e++;
            search_table new_table = table;
            new_table.rows[table.numrows] = table.valid_rows[i];
            new_table.numrows++;
            
            int roww = table.valid_rows[i];
            /*new_table.valid_rows += (i + 1);
            new_table.num_valid_rows -= (i + 1);
            while(smallest_col < table.num_valid_cols && table.valid_cols[smallest_col] <= roww){
                smallest_col++;
            }
            new_table.valid_cols += smallest_col;
            new_table.num_valid_cols -= smallest_col;*/
            bitset_inplace_xor(new_table.unmatched, r.bitarrays[roww]);
            search_aux(r, new_table, inters, ROW, table.valid_rows[i]);
            bitset_inplace_xor(new_table.unmatched, r.bitarrays[roww]);
            
        }
    }
    int smallest_row = 0;
    for(int i = col_ptr; i < table.num_valid_cols; i++){
        if(table.numcols == VEC_SIZE) break;
        if(table.numrows == 0 && table.numcols == 0){
            break;
        }
        int max_elt = r.vecs[table.valid_cols[i]].elts[0];
        //if(max_elt > max_unmatched) continue;
        if(max_elt < max_unmatched){
            break;
        }
        else{
        

            if(!any_match(r.vecs[table.valid_cols[i]], max_unmatched)){
              continue;
            }
            /*if(max_elt == max_unmatched){
              printf("y\n");
            }
            else{
              printf("n\n");
            }*/
            e++;
            search_table new_table = table;
            new_table.cols[table.numcols] = table.valid_cols[i];
            new_table.numcols++;
            
            int coll = table.valid_cols[i];
            /*new_table.valid_cols += (i + 1);
            new_table.num_valid_cols -= (i + 1);
            while(smallest_row < table.num_valid_rows && table.valid_rows[smallest_row] <= coll){
                smallest_row++;
            }
            new_table.valid_rows += smallest_row;
            new_table.num_valid_rows -= smallest_row;*/
            
            bitset_inplace_xor(new_table.unmatched, r.bitarrays[coll]);
            search_aux(r, new_table, inters, COL, table.valid_cols[i]);
            bitset_inplace_xor(new_table.unmatched, r.bitarrays[coll]);
            //table.numcols--;

        }
    }
    //free(valid_rows);
    //free(valid_cols);
     /*if(table.numrows>=1 || table.numcols >= 1){
      printf("%d\n", e);
    }*/
    return;
}

void search_sum(vecgroup group, int sum){
    vec *start = group.infos[sum].start;
    int total_vecs = group.infos[sum].length;
    if(total_vecs == 0){
        return;
    }
    relabeling r = elt_relabeling(group, sum);
    unsigned char** inters = intersections(r, total_vecs);
    printf("inters calculated\n");
    
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
    table.numrows = 0;
    table.numcols = 0;
    table.num_valid_rows = total_vecs;
    table.num_valid_cols = total_vecs;
    table.valid_rows = valid_row_indices;
    table.valid_cols = valid_col_indices;
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
