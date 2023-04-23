#include <stdio.h>
#include <stdlib.h>

#include "search.h"

void record(search_table table, row_table rows, row_table cols, global_t g) {
  *table.num_searched += 1;
  if (rows.num_vecs * cols.num_vecs == VEC_SIZE * VEC_SIZE) {
    return;
    printf("solution found\n");
    for (int i = 0; i < rows.num_vecs; i++) {
      printf("%ld    ", rows.vecs[i]);
      for (int j = 0; j < VEC_SIZE; j++) {
        printf("%ld ", g.label_to_elt[g.vecs[rows.vecs[i]][j]]);
      }
      printf("\n");
    }
    printf("\n");
    for (int i = 0; i < cols.num_vecs; i++) {
      printf("%ld    ", cols.vecs[i]);
      for (int j = 0; j < VEC_SIZE; j++) {
        printf("%ld ", g.label_to_elt[g.vecs[cols.vecs[i]][j]]);
      }
      printf("\n");
    }
  }
}

bool any_match(vec_t vec, size_t target) {
  for (int i = 0; i < VEC_SIZE; i++) {
    if (vec[i] == target) {
      return true;
    }
  }
  return false;
}

bool any_row_match(row_table rows, size_t target){
  for(int i = 0; i < rows.num_vecs; i++){
    if(any_match(rows.vecs[i], target)){
      return true;
    }
  }
  return false;
}



void filter_valids(row_table zeros, row_table ones, inter_grid grid, int elt){
  inter_row row = grid.rows[elt];
  for(int i = 0; i < 1 + (row.size - 1 / 64); i++){
    row.scratch[i] = -1;
  }
  int filt_size = 1 + (row.size - 1 / 64);
  for(int r = zeros.filt_start; r < zeros.num_vecs; r++){
    for(int i = 0; i < filt_size; i++){
      row.scratch[i] &= row.inters_0[grid.vec_to_idx[zeros.vecs[r]] * filt_size + i];
    }
  }
  for(int r = ones.filt_start; r < ones.num_vecs; r++){
    for(int i = 0; i < filt_size; i++){
      row.scratch[i] &= row.inters_1[grid.vec_to_idx[ones.vecs[r]] * filt_size + i];
    }
  }
}

/* new_rows = {row for row in old_rows if inters[old_row] == inter_val} 
  This filters old_rows to include only rows/columns that are compatible to the newly-added vector*/
void fill_valids(row_table *old_rows, row_table *new_rows,
                 unsigned char *inters, unsigned char inter_val,
                 size_t minvec) {
  new_rows->num_valid = 0;
  /* ensure everything is >= minvec */
  size_t i = 0;
  while (old_rows->valid[i] < minvec)
    i++;
  /* keep only things that have inters == inter_val */
  for (; i < old_rows->num_valid; i++) {
    size_t old_row = old_rows->valid[i];
    int is_valid = inters[old_row] == inter_val;
    /* if is_valid, add old_row, else nothing (but branchfree) */
    new_rows->valid[new_rows->num_valid] =
        old_row * is_valid +
        new_rows->valid[new_rows->num_valid] * (1 - is_valid);
    new_rows->num_valid += is_valid;
  }
}


// need to filter before starting? but then first filter doesn't matter?\

void search_aux_grid(global_t g, row_table rows, row_table cols,
                search_table table, inter_grid grid) {
  record(table, rows, cols, g);
  if(rows.num_vecs + cols.num_vecs == 12){
    return;
  }
  int max_unmatched = bitset_maximum(table.unmatched);
  bool is_row = any_row_match(rows, max_unmatched);
  if(is_row){
    filter_valids(cols, rows, grid, max_unmatched);
  }
  else{
    filter_valids(rows, cols, grid, max_unmatched);
  }
  int min_vec = 0;
  // issue: suppose we have added multiple columns before starting search_aux_grid. We need to know what number cols.vecs[0] is
  // doesn't matter if we already filtered
  // need to add as filter before starting search
  if(cols.filt_start == 0 && cols.num_vecs == 1){
    min_vec = cols.vecs[0];
  }
  else if(rows.filt_start == 0 && rows.num_vecs == 1){
    min_vec = rows.vecs[0];
  }
  for(int x = grid.rows[max_unmatched].size - 1; x >= 0; x--){
    if(grid.rows[max_unmatched].scratch[x] == 0){
      continue;
    }
    int idx = 64 * x +  (63 - __builtin_clzll(grid.rows[max_unmatched].scratch[x]));
    grid.rows[max_unmatched].scratch[x] ^= ((uint64_t)1) << (63 - __builtin_clzll(grid.rows[max_unmatched].scratch[x]));
    int vec_num = grid.rows[max_unmatched].vecs[idx];
    if(vec_num < min_vec){
      continue;
    }
    row_table *new_active = is_row ? &cols : &rows;
    new_active->vecs[new_active->num_vecs] = vec_num;
    new_active->num_vecs++;
    bitset_inplace_xor(table.unmatched, g.bitarrays[vec_num]);
    search_aux_grid(g, rows, cols, table, grid);
    bitset_inplace_xor(table.unmatched, g.bitarrays[vec_num]);
    new_active->num_vecs--;

  }
}

/* tries to fill in a row or col */
void search_aux(global_t g, row_table rows, row_table cols,
                search_table table) {
  record(table, rows, cols, g);

  /* if we're done, return */
  if (!((rows.num_vecs + rows.num_valid >= 6) &&
        (cols.num_vecs + cols.num_valid >= 6))) {
    return;
  }

  int max_unmatched = bitset_maximum(table.unmatched);
  /*if(rows.num_vecs + cols.num_vecs >= 3){
    bitset_unset(table.unmatched, max_unmatched);

    int max_unmatched_2 = bitset_maximum(table.unmatched);

    bitset_set(table.unmatched, max_unmatched);
    int c1 = 0;
    int c2 = 0;
    for (int to_add = ROW; to_add <= COL; to_add++) {
      row_table active = to_add == ROW ? rows : cols;
      for (int i = 0; i < active.num_valid; i++) {
        int max_elt = g.vecs[active.valid[i]][0];
        if (max_elt < max_unmatched_2)
          break;
        if (any_match(g.vecs[active.valid[i]], max_unmatched)){
          c1++;
        }
        if (any_match(g.vecs[active.valid[i]], max_unmatched_2)){
          c2++;
        }
      }
    }
    //printf("\n%d %d %d %d\n", max_unmatched, max_unmatched_2, rows.num_vecs,cols.num_vecs);
    if(c1 == 0 || c2 == 0){
      return;
    }
    /*if(c1 > c2){
      max_unmatched = max_unmatched_2;
    }*//*
    int c3 = 0;
    for (int to_add = ROW; to_add <= COL; to_add++) {
      row_table active = to_add == ROW ? rows : cols;
      row_table other = to_add == COL ? rows : cols;
      for (int i = 0; i < active.num_valid; i++) {
        int max_elt = g.vecs[active.valid[i]][0];
        //printf("%d ", max_elt);
        if(bitset_get(table.unmatched, max_elt) == 1){
          continue;
        }
        c3++;
        for(int j = 0; j < other.num_valid; j++){
          int max_elt_2 = g.vecs[other.valid[j]][0];
          if (max_elt_2 < max_elt)
            break;
          if (any_match(g.vecs[other.valid[j]], max_elt)){
            c3--;
            break;
          }
          
        }

      }
    }
    //printf("%d %d\n", c3, rows.num_valid+cols.num_valid);
  }*/





  for (int to_add = ROW; to_add <= COL; to_add++) {
    /* don't add a col if we don't have a row */
    if (to_add == COL && rows.num_vecs == 0)
      break;
    /* we're trying to_add a vector to rows or cols */
    row_table active = to_add == ROW ? rows : cols;
    if (active.num_vecs == VEC_SIZE)
      continue;

    /* for each valid vector, try adding it: */
    for (int i = 0; i < active.num_valid; i++) {
      /* force this vector to match max_unmatched (unless we're adding our
       * first row, in which case there's nothing to match */
      int max_elt = g.vecs[active.valid[i]][0];
      if (max_elt < max_unmatched)
        break;
      if (!any_match(g.vecs[active.valid[i]], max_unmatched) &&
          (to_add == COL || rows.num_vecs > 0))
        continue;

      /* add the new vector to rows or cols */
      row_table new_rows = rows;
      row_table new_cols = cols;
      size_t new_vec = active.valid[i];
      row_table *new_active = to_add == ROW ? &new_rows : &new_cols;
      new_active->vecs[active.num_vecs] = new_vec;
      new_active->num_vecs++;

      /* update new_rows.valid and new_cols.valid */
      size_t new_minvec = new_vec;
      if (new_rows.num_vecs > 1 || new_cols.num_vecs > 1) {
        new_minvec = 0;
      }
      new_rows.valid = g.row_idx_slots[rows.num_vecs + cols.num_vecs + 1];
      new_cols.valid = g.col_idx_slots[rows.num_vecs + cols.num_vecs + 1];
      fill_valids(&rows, &new_rows, g.inters[new_vec], to_add, new_minvec);
      fill_valids(&cols, &new_cols, g.inters[new_vec], 1 - to_add, new_minvec);

      /* update table.unmatched and recurse */
      bitset_inplace_xor(table.unmatched, g.bitarrays[new_vec]);
      search_aux(g, new_rows, new_cols, table);
      bitset_inplace_xor(table.unmatched, g.bitarrays[new_vec]);
    }
  }

  return;
}
