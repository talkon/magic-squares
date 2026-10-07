/*
 * Recursive search, specialized for a label-bitset width of W 64-bit words.
 * Included from arrange.c once per width.
 *
 * Invariant on entering SEARCH_REC(s, d): the candidate lists valid[d][a] and
 * their bit-sliced label counts cnt[d][a] are filled in, every unmatched cell
 * has at least one candidate, and each axis has enough candidates left. The
 * parent does this work while creating the child (see TRY_CHILD), so that
 * children which would be pruned are discarded as cheaply as possible.
 */
#ifndef W
#error "define W before including arrange_core.h"
#endif

#define CAT_(a, b) a##_##b
#define CAT(a, b) CAT_(a, b)
#define SEARCH_REC CAT(search_rec, W)
#define SEARCH_ROOT CAT(search_root, W)
#define FILTER_COUNT CAT(filter_count, W)
#define TRY_CHILD CAT(try_child, W)

/*
 * out = { u in in[0..cnt) : u >= min_v and bit u of row is set }, and the
 * bit-sliced counts (>= 1, >= 2, ..., >= NSLICE vectors) of each label over
 * out. Returns |out|.
 */
static inline uint32_t FILTER_COUNT(const sstate_t *s, const uint32_t *in,
                                    uint32_t cnt, const uint64_t *row,
                                    uint32_t min_v, uint32_t *out,
                                    uint64_t *g /* [NSLICE][W] */) {
  uint32_t k = filter_list(in, cnt, row, min_v, out);
#ifdef HAVE_SIMD_COUNT
  if (W == 2) {
    count_slices_simd_w2(s->bits, out, k, g);
    return k;
  }
  if (W == 4) {
    count_slices_simd_w4(s->bits, out, k, g);
    return k;
  }
  if (W == 8) {
    count_slices_simd_w8(s->bits, out, k, g);
    return k;
  }
#endif
  /* thermometer code: gg[t] = labels in more than t of the kept vectors */
  uint64_t gg[NSLICE][W];
  for (int t = 0; t < NSLICE; t++)
    for (int w = 0; w < W; w++)
      gg[t][w] = 0;
  for (uint32_t i = 0; i < k; i++) {
    const uint64_t *m = s->bits + (uint64_t)out[i] * W;
    for (int w = 0; w < W; w++) {
      for (int t = NSLICE - 1; t > 0; t--)
        gg[t][w] |= gg[t - 1][w] & m[w];
      gg[0][w] |= m[w];
    }
  }
  for (int t = 0; t < NSLICE; t++)
    for (int w = 0; w < W; w++)
      g[t * W + w] = gg[t][w];
  return k;
}

static void SEARCH_REC(sstate_t *s, int d);

/*
 * Create the child of depth d+1 obtained by placing v on axis b, and search
 * it unless it is pruned. min_v restricts candidates to indices >= min_v.
 */
static inline void TRY_CHILD(sstate_t *s, int d, int b, uint32_t v,
                             uint32_t min_v) {
  const int n = s->n;
  s->nodes++;
  if (s->opts.node_limit && s->nodes > s->opts.node_limit) {
    s->stop = 2;
    return;
  }
  const uint64_t *m = s->bits + (uint64_t)v * W;
  const uint64_t *rc = s->cells[d][ROW], *cc = s->cells[d][COL];
  uint64_t *nrc = s->cells[d + 1][ROW], *ncc = s->cells[d + 1][COL];
  for (int w = 0; w < W; w++) {
    nrc[w] = rc[w] | (b == ROW ? m[w] : 0);
    ncc[w] = cc[w] | (b == COL ? m[w] : 0);
  }
  s->placed[b][s->np[b]++] = v;
  if (s->np[ROW] == n && s->np[COL] == n) {
    report(s);
    s->np[b]--;
    return;
  }
  const uint64_t *ncell_b = b == ROW ? nrc : ncc;
  const uint64_t *ncell_o = b == ROW ? ncc : nrc;

  /* candidates on the other axis must meet v exactly once; they cover the
   * unmatched cells on axis b (including v's own new cells), which is where
   * most dead ends show up, so check those first */
  const int o = 1 - b;
  uint64_t *go = s->cnt[d + 1][o];
  uint32_t ko = FILTER_COUNT(s, s->valid[d][o], s->nvalid[d][o],
                             s->inters1 + (uint64_t)v * s->IW, min_v,
                             s->valid[d + 1][o], go);
  s->nvalid[d + 1][o] = ko;
  uint64_t dead = 0;
  if (s->np[o] + (int)ko < n)
    dead = 1;
  if (s->opts.forward_check)
    for (int w = 0; w < W; w++)
      dead |= ncell_b[w] & ~ncell_o[w] & ~go[w];
  if (!dead) {
    /* candidates on axis b must be disjoint from v; they cover the unmatched
     * cells on the other axis */
    uint64_t *gb = s->cnt[d + 1][b];
    uint32_t kb = FILTER_COUNT(s, s->valid[d][b], s->nvalid[d][b],
                               s->inters0 + (uint64_t)v * s->IW, min_v,
                               s->valid[d + 1][b], gb);
    s->nvalid[d + 1][b] = kb;
    if (s->np[b] + (int)kb < n)
      dead = 1;
    if (s->opts.forward_check)
      for (int w = 0; w < W; w++)
        dead |= ncell_o[w] & ~ncell_b[w] & ~gb[w];
    if (!dead)
      SEARCH_REC(s, d + 1);
  }
#ifdef PROFILE
  prof_created[d + 1]++;
  prof_kids[d]++;
  prof_pruned[d + 1] += dead != 0;
#endif
  s->np[b]--;
}

static void SEARCH_REC(sstate_t *s, int d) {
#ifdef PROFILE
  prof_vr[d] += s->nvalid[d][ROW];
  prof_vc[d] += s->nvalid[d][COL];
#endif
  if (s->stop)
    return;
  const uint64_t *rc = s->cells[d][ROW], *cc = s->cells[d][COL];
  /* candidates on axis a cover unmatched cells of axis 1-a, so cells of
   * placed rows are classified by the counts over candidate cols */
  const uint64_t *gr = s->cnt[d][COL], *gc = s->cnt[d][ROW];

  /* branch on the unmatched cell with the fewest candidates (exactly 1, 2,
   * ..., or >= NSLICE), ties broken by the largest label */
  uint64_t urow[W], ucol[W];
  for (int w = 0; w < W; w++) {
    urow[w] = rc[w] & ~cc[w];
    ucol[w] = cc[w] & ~rc[w];
  }
  int x = -1;
  for (int c = s->opts.mrv ? 0 : NSLICE - 1; c < NSLICE && x < 0; c++) {
    for (int w = W - 1; w >= 0; w--) {
      uint64_t cls;
      if (!s->opts.mrv) {
        cls = urow[w] | ucol[w];
      } else {
        uint64_t cr = urow[w] & gr[c * W + w], cl = ucol[w] & gc[c * W + w];
        if (c + 1 < NSLICE) {
          cr &= ~gr[(c + 1) * W + w];
          cl &= ~gc[(c + 1) * W + w];
        }
        cls = cr | cl;
      }
      if (cls) {
        x = w * 64 + 63 - __builtin_clzll(cls);
        break;
      }
    }
  }
  if (x < 0)
    return;
  const int xw = x / 64;
  const uint64_t xb = (uint64_t)1 << (x % 64);
  /* a row-unmatched cell is covered by a col, and vice versa */
  const int b = (urow[xw] & xb) ? COL : ROW;

  const uint32_t *val = s->valid[d][b];
  uint32_t cnt = s->nvalid[d][b];
  for (uint32_t i = 0; i < cnt && !s->stop; i++) {
    uint32_t v = val[i];
    if (s->bits[(uint64_t)v * W + xw] & xb)
      TRY_CHILD(s, d, b, v, 0);
  }
}

/* the first row r1 is the lowest-index vector of the square, so it contains
 * the square's largest label; every other vector has index > r1 */
static void SEARCH_ROOT(sstate_t *s) {
  s->nodes = 1;
  uint32_t count = s->N;
  for (int a = 0; a < 2; a++) {
    for (uint32_t i = 0; i < count; i++)
      s->valid[0][a][i] = i;
    s->nvalid[0][a] = count;
    for (int w = 0; w < W; w++)
      s->cells[0][a][w] = 0;
  }
  for (uint32_t r1 = 0; r1 < count && !s->stop; r1++)
    TRY_CHILD(s, 0, ROW, r1, r1 + 1);
}

#undef SEARCH_REC
#undef SEARCH_ROOT
#undef FILTER_COUNT
#undef TRY_CHILD
#undef CAT
#undef CAT_
