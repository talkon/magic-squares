#include <stdio.h>
#include <stdlib.h>

#include "new_inters.h"

inter_grid create_inter_grid(global_t g, row_table rows, row_table cols, search_table table){

    // need special case for when zero rows have been placed?
    

    int* counts = calloc(g.num_labels, sizeof(int));
    int* num_rows = calloc(g.num_labels, sizeof(int));
    for (int to_add = ROW; to_add <= COL; to_add++) {
        row_table active = to_add == ROW ? rows : cols;
        for (int i = 0; i < active.num_valid; i++) {
            vec_t v = g.vecs[active.valid[i]];
            for(int j = 0; j < VEC_SIZE; j++){
                counts[v.elts[j]]++;
                if(to_add == ROW){
                    num_rows[v.elts[j]]++;
                }
            }
        }
    }
    int* partial_counts = calloc(g.num_labels, sizeof(int));
    int** vecs_by_elt = malloc((g.num_labels)*sizeof(int*));
    for (int elt = 0; elt < g.num_labels; elt++){
        vecs_by_elt[elt] = malloc(counts[elt] * sizeof(int) + 1);
    }
    // second loop is not redundant; needed to know how much to malloc in vecs_by_elt
    for (int to_add = ROW; to_add <= COL; to_add++) {
        row_table active = to_add == ROW ? rows : cols;
        for (int i = 0; i < active.num_valid; i++) {
            vec_t v = g.vecs[active.valid[i]];
            for(int j = 0; j < VEC_SIZE; j++){
                vecs_by_elt[v.elts[j]][partial_counts[v.elts[j]]] = active.valid[i];
                partial_counts[v.elts[j]]++;
            }
        }
    }
    // Note: after one vector has been placed, it is impossible for a vector to be both a valid row and a valid col
    inter_grid grid;
    grid.rows = malloc((rows.num_valid+cols.num_valid)*sizeof(inter_row));
    int max_row = rows.valid[rows.num_valid - 1];
    int max_col = cols.valid[cols.num_valid - 1];
    grid.vec_to_idx = malloc((max_row > max_col ? max_row : max_col)*sizeof(int));
    for(int i = 0; i < rows.num_valid; i++){
        grid.vec_to_idx[rows.valid[i]] = i;
    }
    for(int i = 0; i < cols.num_valid; i++){
        grid.vec_to_idx[cols.valid[i]] = i + rows.num_valid;
    }

    for (int elt = 0; elt < g.num_labels; elt++){
        if(counts[elt] == 0){
            continue;
        }
        int size = counts[elt];
        int num_rows = num_rows[elt];
        uint64_t inters_0 = calloc((1 + (size - 1) / 64) * (rows.num_valid+cols.num_valid), sizeof(uint64_t));
        uint64_t inters_1 = calloc((1 + (size - 1) / 64) * (rows.num_valid+cols.num_valid), sizeof(uint64_t));
        uint64_t scratch = calloc((1 + (size - 1) / 64), sizeof(uint64_t));
        int* vecs = vecs_by_elt[elt];
        inter_row row;
        row.inters_0 = inters_0;
        row.inters_1 = inters_1;
        row.scratch = scratch;
        row.size = size;
        row.num_rows = num_rows;
        row.vecs = vecs;
        for(int v1 = 0; v1 < (rows.num_valid+cols.num_valid); v1++){
            int vec;
            if(v1 < rows.num_valid){
                vec = rows.valid[v1];
            }
            else{
                vec = cols.valid[v1 - rows.num_valid];
            }
            for(int v2 = 0; v2 < size; v2++){
                int vec2 = vecs_by_elt[elt][v2];
                if(g.inters[vec][vec2] == 0){
                    inter_set(row, vec, vec2, 0);
                }
                if(g.inters[vec][vec2] == 1){
                    inter_set(row, vec, vec2, 1);
                }
            }
        }
        grid.rows[elt] = row;
    }
    return grid;
    free(counts);
    free(partial_counts);
    free(num_rows);
    free(vecs_by_elt);
}

void free_inter_grid(inter_grid grid){
    
}
