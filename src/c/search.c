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
    fflush(stdout);
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

/// new_vecs = vec in old_vecs such that:
//   - g.inters[new_vec][vec] == inter_val
//   - vec >= min_vec
void fill_valids(const global_t g, const row_table *old_vecs,
                 row_table *new_vecs, const search_table table,
                 unsigned char inter_val, size_t new_vec, size_t min_vec) {
  new_vecs->num_valid = 0;
  // ensure everything is >= minvec
  size_t i = 0;
  while (old_vecs->valid[i] < min_vec)
    i++;
  for (; i < old_vecs->num_valid; i++) {
    size_t vec = old_vecs->valid[i];
    // keep only things that have inters == inter_val
    int is_valid = g.inters[new_vec][vec] == inter_val;
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
  if (!((rows.num_vecs + rows.num_valid >= VEC_SIZE) &&
        (cols.num_vecs + cols.num_valid >= VEC_SIZE))) {
    return;
  }

  size_t max_unmatched = bitset_maximum(table.unmatched);

  // add to either ROW or COL
  for (size_t to_add = ROW; to_add <= COL; to_add++) {
    row_table active = to_add == ROW ? rows : cols;
    // don't add a col if we don't have a row
    if (to_add == COL && rows.num_vecs == 0) {
      break;
    }
    // don't add if we're done
    if (active.num_vecs == VEC_SIZE) {
      continue;
    }

    // for each valid vector, try adding it:
    for (size_t i = 0; i < active.num_valid; i++) {
      // want this vec to have max_unmatched
      // (unless we haven't added anything yet)
      size_t max_elt = g.vecs[active.valid[i]][0];
      if (max_elt < max_unmatched) {
        break;
      }
      if (!any_match(g.vecs[active.valid[i]], max_unmatched) &&
          (to_add == COL || rows.num_vecs > 0)) {
        continue;
      }

      // add the new vector to rows or cols, update table unmatched
      row_table new_rows = rows;
      row_table new_cols = cols;
      size_t new_vec = active.valid[i];
      row_table *new_active = to_add == ROW ? &new_rows : &new_cols;
      new_active->vecs[active.num_vecs] = new_vec;
      new_active->num_vecs++;
      bitset_inplace_xor(table.unmatched, g.bitarrays[new_vec]);

      // update new_rows.valid and new_cols.valid
      size_t min_vec = new_vec;
      if (new_rows.num_vecs > 1 || new_cols.num_vecs > 1) {
        min_vec = 0;
      }
      new_rows.valid = g.row_idx_slots[rows.num_vecs + cols.num_vecs + 1];
      new_cols.valid = g.col_idx_slots[rows.num_vecs + cols.num_vecs + 1];
      fill_valids(g, &rows, &new_rows, table, to_add, new_vec, min_vec);
      fill_valids(g, &cols, &new_cols, table, 1 - to_add, new_vec, min_vec);

      // recurse and reset
      search_aux(g, new_rows, new_cols, table);
      bitset_inplace_xor(table.unmatched, g.bitarrays[new_vec]);
    }
  }

  return;
}
