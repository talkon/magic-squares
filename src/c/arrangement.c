#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "arrangement.h"


int sum(row r){
    int s = 0;
    for(int i = 0; i < ROW_SIZE; i++){
        s += r.els[i];
    }
    return s;
}

// IMPORTANT: the arrays pointed to by the returned rowgroup must be freed after use
rowgroup read_rows(char *filename){
    FILE *fp = fopen(filename,"r");
    int numrows = 0;
    char buf[1024];
    while(fgets(buf, 1024, fp) != NULL){
        numrows++;
    }
    row* row_arr = malloc(sizeof(row) * numrows);
    fp = fopen(filename,"r");
    int rowcount = 0;
    int maxsum;
    while(fgets(buf, 1024, fp) != NULL){
        row slot;
        sscanf(buf, "%d %d %d %d %d %d", &slot.els[0], &slot.els[1], &slot.els[2], &slot.els[3], &slot.els[4], &slot.els[5]);
        row_arr[rowcount] = slot;
        rowcount++;
        maxsum = sum(slot);
    }
    int cursum = 0;
    // all infos are initialized to 0 length and null starts
    rowgroupinfo* infos = calloc(maxsum + 1, sizeof(rowgroupinfo));
    for(int i = 0; i < numrows; i++){
        row currow = row_arr[i];
        if(sum(currow) != cursum){
            infos[sum(currow)].start = row_arr + i;
        }
        infos[sum(currow)].length++;
        cursum = sum(currow);
    }
    rowgroup r;
    r.rows = row_arr;
    r.infos = infos;
    r.numrows = numrows;
    r.numsums = maxsum + 1;
    return r;

}

void test_print_rows(char *filename){
    rowgroup group = read_rows(filename);
    for(int i = 0; i < group.numsums; i++){
        printf("%d\n", group.infos[i].length);
    }
    free(group.rows);
    free(group.infos);
}

int main(int argc, char *argv[]){
    int opt;
    char filename[1024];
    while ((opt = getopt (argc, argv, "f:")) != -1){
        test_print_rows(optarg);
    }
}