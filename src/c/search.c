#include <stdio.h>
#include <stdlib.h>

#include "search.h"

void record(search_table table, row_table rows, row_table cols, global_t g) {
  *table.num_searched += 1;
  if (rows.num_vecs * cols.num_vecs == VEC_SIZE * VEC_SIZE) {
    // return;
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

// new_vecs = vec in old_vecs such that:
//   0. ensure g.inters[new_vec][vec] == inter_val
//   1. ensure vec >= min_vec
//   2. ensure max(vec) >= max_unmatched
//   2. ensure max_unmatched in vec
void fill_valids(const global_t g, const row_table *old_vecs,
                 row_table *new_vecs, const search_table table,
                 unsigned char inter_val, size_t new_vec, size_t min_vec) {
  new_vecs->num_valid = 0;
  int max_unmatched = bitset_maximum(table.unmatched);
  // optimization 1: ensure everything is >= minvec
  size_t i = 0;
  while (old_vecs->valid[i] < min_vec) {
    i++;
  }
  // optimization 2: ensure max(vec) >= max_unmatched
  size_t max_i = i;
  while (g.vecs[old_vecs->valid[max_i]][0] >= max_unmatched) {
    max_i++;
  }
  // for each vec in old_vecs...
  for (; i < max_i; i++) {
    size_t vec = old_vecs->valid[i];
    // requirement 0: ensure g.inters[new_vec][vec] == inter_val
    // optimization 3: max_unmatched in vec
    int is_valid = g.inters[new_vec][vec] == inter_val &&
                   any_match(g.vecs[vec], max_unmatched);
    // if is_valid, add old_row, else nothing (but branchfree)
    new_vecs->valid[new_vecs->num_valid] =
        vec * is_valid + new_vecs->valid[new_vecs->num_valid] * (1 - is_valid);
    new_vecs->num_valid += is_valid;
  }
}

// tries to fill in a row or col
void search_aux(global_t g, row_table rows, row_table cols,
                search_table table) {
  record(table, rows, cols, g);

  // if we're done, return
  if (!((rows.num_vecs + rows.num_valid >= 6) &&
        (cols.num_vecs + cols.num_valid >= 6))) {
    return;
  }

  for (int to_add = ROW; to_add <= COL; to_add++) {
    // don't add a col if we don't have a row
    if (to_add == COL && rows.num_vecs == 0) {
      break;
    }
    // we're trying to_add a vector to rows or cols
    row_table active = to_add == ROW ? rows : cols;
    if (active.num_vecs == VEC_SIZE) {
      continue;
    }

    // for each valid vector, try adding it:
    for (int i = 0; i < active.num_valid; i++) {
      row_table new_rows = rows;
      row_table new_cols = cols;
      size_t new_vec = active.valid[i];
      row_table *new_active = to_add == ROW ? &new_rows : &new_cols;
      new_active->vecs[active.num_vecs] = new_vec;
      new_active->num_vecs++;

      // update new_rows.valid and new_cols.valid
      bitset_inplace_xor(table.unmatched, g.bitarrays[new_vec]);
      size_t min_vec = new_vec;
      if (new_rows.num_vecs > 1 || new_cols.num_vecs > 1) {
        min_vec = 0;
      }

      new_rows.valid = g.row_idx_slots[rows.num_vecs + cols.num_vecs + 1];
      new_cols.valid = g.col_idx_slots[rows.num_vecs + cols.num_vecs + 1];
      fill_valids(g, &rows, &new_rows, table, to_add, new_vec, min_vec);
      fill_valids(g, &cols, &new_cols, table, 1 - to_add, new_vec, min_vec);

      // recurse
      search_aux(g, new_rows, new_cols, table);
      bitset_inplace_xor(table.unmatched, g.bitarrays[new_vec]);
    }
  }

  return;
}
