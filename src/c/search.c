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

/* new_rows = {row for row in old_rows if inters[old_row] == inter_val} */
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
