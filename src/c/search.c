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

bool in_any_row(global_t g, row_table *rows, size_t target) {
  for (size_t i = 0; i < rows->num_vecs; i++) {
    if (any_match(g.vecs[rows->vecs[i]], target)) {
      return true;
    }
  }
  return false;
}

/// new_vecs = vec in old_vecs such that:
//   - g.inters[new_vec][vec] == inter_val
//   - vec >= min_vec
void fill_valids(const global_t g, const row_table *old_vecs,
                 row_table *new_vecs, unsigned char inter_val, size_t new_vec,
                 size_t min_vec) {
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
void search_aux(global_t g, row_table rows_init, row_table cols_init,
                search_table table, size_t to_add_init) {
  typedef struct {
    row_table rows;
    row_table cols;
    size_t to_add;
    size_t *i;
    size_t *max_i;
  } frame;

  frame stack[NUM_VALID_SLOTS];
  frame *s = stack;

  stack[0].rows = rows_init;
  stack[0].cols = cols_init;
  stack[0].to_add = to_add_init;
  stack[0].i = NULL;
  stack[0].max_i = NULL;

  while (s >= stack) {
    row_table active = s->to_add == ROW ? s->rows : s->cols;

    // set i and max_i if not set
    if (s->i == NULL) {
      s->i = active.valid;
      s->max_i = s->i + active.num_valid - 1;
      while (s->max_i >= s->i && g.vecs[*s->max_i][0] < table.max_unmatched) {
        s->max_i--;
      }
      s->max_i = s->max_i;
    }

  loop_head:
    // cleanup from last time
    if (s->i > active.valid) {
      toggle_unmatched(&table, g.bitarrays[*(s->i - 1)]);
    }

    // want this vec to have max_unmatched
    // (unless we haven't added anything yet)
    if (s->rows.num_vecs > 0) {
      while (s->i <= s->max_i &&
             !any_match(g.vecs[*s->i], table.max_unmatched)) {
        s->i++;
      }
    }
    if (s->i > s->max_i) {
      s--;
      continue;
    }

    // add the new vector to rows or cols, update table unmatched
    row_table new_rows = s->rows;
    row_table new_cols = s->cols;
    row_table *new_active = s->to_add == ROW ? &new_rows : &new_cols;
    new_active->vecs[active.num_vecs] = *s->i;
    new_active->num_vecs++;
    toggle_unmatched(&table, g.bitarrays[*s->i]);

    // update new_rows.valid and new_cols.valid
    size_t min_vec = new_rows.num_vecs > 1 || new_cols.num_vecs > 1 ? 0 : *s->i;
    new_rows.valid = g.row_idx_slots[s->rows.num_vecs + s->cols.num_vecs + 1];
    new_cols.valid = g.col_idx_slots[s->rows.num_vecs + s->cols.num_vecs + 1];
    fill_valids(g, &s->rows, &new_rows, s->to_add, *s->i, min_vec);
    fill_valids(g, &s->cols, &new_cols, 1 - s->to_add, *s->i, min_vec);

    // recurse and reset
    record(table, new_rows, new_cols, g);

    if ((new_rows.num_vecs + new_rows.num_valid >= VEC_SIZE) &&
        (new_cols.num_vecs + new_cols.num_valid >= VEC_SIZE) &&
        !(new_rows.num_vecs == VEC_SIZE && new_cols.num_vecs == VEC_SIZE)) {
      size_t new_to_add =
          in_any_row(g, &new_rows, table.max_unmatched) ? COL : ROW;
      // freeze current stack
      s->i = s->i + 1;
      // add stack call
      s++;
      s->rows = new_rows;
      s->cols = new_cols;
      s->to_add = new_to_add;
      s->i = NULL;
      s->max_i = NULL;
      continue;
    }

    s->i++;
    goto loop_head;
  }

  return;
}
