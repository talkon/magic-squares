#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "arrangement.h"


int sum(vec r){
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
    vec* vec_arr = malloc(sizeof(vec) * numvecs);
    fp = fopen(filename,"r");
    int veccount = 0;
    int maxsum;
    while(fgets(buf, 1024, fp) != NULL){
        vec slot;
        sscanf(buf, "%d %d %d %d %d %d", &slot.elts[0], &slot.elts[1], &slot.elts[2], &slot.elts[3], &slot.elts[4], &slot.elts[5]);
        vec_arr[veccount] = slot;
        veccount++;
        maxsum = sum(slot);
    }
    int cursum = 0;
    // all infos are initialized to 0 length and null starts
    vecgroupinfo* infos = calloc(maxsum + 1, sizeof(vecgroupinfo));
    for(int i = 0; i < numvecs; i++){
        vec curvec = vec_arr[i];
        if(sum(curvec) != cursum){
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
    vec* vec_arr = group.infos[sum].start;
    int total_vecs = group.infos[sum].length;
    relabeling r;
    indexed_int* elt_freqs = calloc(sum, sizeof(indexed_int));
    for(int i = 0; i < sum; i++){
        elt_freqs[i].elt = i;
    }
    for(int i = 0; i < total_vecs; i++){
        for(int j = 0; j < VEC_SIZE; j++){
            elt_freqs[vec_arr[i].elts[j]].freq++;
        }
    }
    qsort(elt_freqs, sum, sizeof(indexed_int), index_cmp_reverse);
    int* elt_to_label = malloc(sum * sizeof(int));
    int* label_to_elt = malloc(sum * sizeof(int));
    for(int i = 0; i < sum; i++){
        elt_to_label[elt_freqs[i].elt] = i;
        label_to_elt[i] = elt_freqs[i].elt;
    }
    vec* new_vecs = malloc(total_vecs * sizeof(vec));
    for(int i = 0; i < total_vecs; i++){
        for(int j = 0; j < VEC_SIZE; j++){
            new_vecs[i].elts[j] = elt_to_label[vec_arr[i].elts[j]];
        }
        qsort(&(new_vecs[i].elts[0]), VEC_SIZE, sizeof(int), int_reverse);
    }
    qsort(new_vecs, total_vecs, sizeof(vec), vec_reverse);
    r.vecs = new_vecs;
    r.label_to_elt = label_to_elt;
    free(elt_to_label);
    free(elt_freqs);
    return r;
}

int num_inters(vec a, vec b){
    int p1 = 0;
    int p2 = 0;
    int tot = 0;
    while(p1 < VEC_SIZE && p2 < VEC_SIZE){
        if(a.elts[p1] > a.elts[p2]){
            p1++;
        }
        else{
            if(a.elts[p1] == a.elts[p2]){
                p1++;
                tot++;
            }
            p2++;
        }
    }
    return tot;
}

unsigned char** intersections(vec* vecs, int numvecs){
    unsigned char** inters = malloc(numvecs * sizeof(char*));
    for(int i = 0; i < numvecs; i++){
        inters[i] = malloc(numvecs * sizeof(char));
        for(int j = 0; j < numvecs; j++){
            inters[i][j] = num_inters(vecs[i], vecs[j]);
        }
    }
    return inters;
}

void search_sum(vecgroup group, int sum){
    vec* start = group.infos[sum].start;
    int total_vecs = group.infos[sum].length;
    if(total_vecs == 0){
        return;
    }
    relabeling r = elt_relabeling(group, sum);
    unsigned char** inters = intersections(r.vecs, total_vecs);
    int* valid_indices = malloc(total_vecs * sizeof(int));
    for(int i = 0; i < total_vecs; i++){
        valid_indices[i] = i;
    }
    int num_valid_vecs = total_vecs;
    free(inters);
    free(r.vecs);
    free(r.label_to_elt);
    free(valid_indices);
}

void test_print_vecs(char *filename){
    vecgroup group = read_vecs(filename);
    for(int i = 0; i < group.numsums; i++){
        printf("%d\n", group.infos[i].length);
    }
    search_sum(group, 327);
    free(group.vecs);
    free(group.infos);
}



int main(int argc, char *argv[]){
    int opt;
    char filename[1024];
    while ((opt = getopt(argc, argv, "f:")) != -1){
        test_print_vecs(optarg);
    }
}