#include <stdio.h>
#include <stdlib.h>

#include "search.h"
#include "simde/x86/avx512.h"

SIMDE_FUNCTION_ATTRIBUTES
simde__m512i simde_mm512_i32gather_epi32(simde__m512i vindex, void *base_addr,
                                         const int32_t scale)
    SIMDE_REQUIRE_CONSTANT(scale) HEDLEY_REQUIRE_MSG(
        (scale && scale <= 8 && !(scale & (scale - 1))),
        "`scale' must be a power of two less than or equal to 8") {
#if defined(SIMDE_X86_AVX512F_NATIVE)
  return _mm512_i32gather_epi32(vindex, base_addr, scale);
#else
  simde__m512i_private vindex_ = simde__m512i_to_private(vindex), r_;
  const uint8_t *addr = HEDLEY_REINTERPRET_CAST(const uint8_t *, base_addr);

  SIMDE_VECTORIZE
  for (size_t i = 0; i < (sizeof(vindex_.i32) / sizeof(vindex_.i32[0])); i++) {
    const uint8_t *src = addr + (HEDLEY_STATIC_CAST(size_t, vindex_.i32[i]) *
                                 HEDLEY_STATIC_CAST(size_t, scale));
    int32_t dst;
    simde_memcpy(&dst, src, sizeof(dst));
    r_.i32[i] = dst;
  }

  return simde__m512i_from_private(r_);
#endif
}

void record(search_table table, row_table rows, row_table cols, global_t g) {
  *table.num_searched += 1;
  if (rows.num_vecs * cols.num_vecs == VEC_SIZE * VEC_SIZE) {
    // return;
    printf("solution found\n");
    for (size_t i = 0; i < rows.num_vecs; i++) {
      printf("%ld    ", rows.vecs[i]);
      for (size_t j = 0; j < VEC_SIZE; j++) {
        printf("%ld ", g.label_to_elt[g.vecs[rows.vecs[i]][j]]);
      }
      printf("\n");
    }
    printf("\n");
    for (size_t i = 0; i < cols.num_vecs; i++) {
      printf("%ld    ", cols.vecs[i]);
      for (size_t j = 0; j < VEC_SIZE; j++) {
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

void fill_valids_vectorized(const global_t g, const row_table *old_vecs,
                            row_table *new_vecs, const search_table table,
                            unsigned char inter_val, size_t new_vec,
                            size_t min_vec) {
  int adjusted_num_valid = old_vecs->num_valid - old_vecs->num_valid % 16;
  int offset = 0;

  simde__m512i inter_vals = simde_mm512_set1_epi32(inter_val);
  simde__m512i min_vecs = simde_mm512_set1_epi32((uint32_t)min_vec);
  simde__m512i mask = simde_mm512_set1_epi32(255);

  for (int i = 0; i < adjusted_num_valid / 16; i++) {
    simde__m512i vecs = simde_mm512_loadu_si512(old_vecs->valid + i * 16);
    simde__m512i inters = simde_mm512_i32gather_epi32(
        vecs, (void *)g.inters[new_vec], sizeof(unsigned char));
    inters &= mask;
    simde__mmask16 size_mask = simde_mm512_cmpge_epi32_mask(vecs, min_vecs);
    if (size_mask == 0) {
      continue;
    }
    simde__mmask16 inter_mask =
        simde_mm512_cmpeq_epi32_mask(inters, inter_vals);
    simde_mm512_mask_compressstoreu_epi32(new_vecs->valid + offset,
                                          inter_mask & size_mask, vecs);
    int num_vals = __builtin_popcountll((uint64_t)(inter_mask & size_mask));
    offset += num_vals;
    new_vecs->num_valid += num_vals;
  }
}

/// new_vecs = vec in old_vecs such that:
//   - g.inters[new_vec][vec] == inter_val
//   - vec >= min_vec
void fill_valids(const global_t g, const row_table *old_vecs,
                 row_table *new_vecs, const search_table table,
                 unsigned char inter_val, size_t new_vec, size_t min_vec) {
  new_vecs->num_valid = 0;

  // vectorized
  fill_valids_vectorized(g, old_vecs, new_vecs, table, inter_val, new_vec,
                         min_vec);

  // remaining
  size_t i = old_vecs->num_valid - old_vecs->num_valid % 16;
  for (; i < old_vecs->num_valid; i++) {
    uint32_t vec = old_vecs->valid[i];
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
  size_t to_add = ROW;
  if (rows.num_vecs > 0) {
    to_add = in_any_row(g, &rows, max_unmatched) ? COL : ROW;
  }
  row_table active = to_add == ROW ? rows : cols;
  // don't add if we're done
  if (active.num_vecs == VEC_SIZE) {
    return;
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

  return;
}
