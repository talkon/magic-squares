/*
 * Recursive search, specialized for a label-bitset width of W 64-bit words.
 * Included from arrange.c once per width.
 *
 * Invariant on entering SEARCH_REC(s, d): the candidate lists valid[d][a] and
 * their label counts cnt[d][a] are filled in, every unmatched cell
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
#define COVERED CAT(covered, W)
#define IN_CLASS CAT(in_class, W)

/* the kind of label counters for this width (see arrange.c), and the number
 * of classes of counts used to pick the most constrained cell: exactly 1, 2,
 * ..., NCLASS - 1 candidates, or >= NCLASS */
#if W <= COUNT_BYTES_MAX_W
#define COUNT_BYTES 1
#define NCLASS 255
#else
#define NCLASS NSLICE
#endif

/*
 * Reading the label counts g = cnt[d][a] (CNT_WORDS * W words), for the
 * labels of word w:
 * COVERED: labels in at least one candidate;
 * IN_CLASS: labels in exactly c + 1 candidates, or, for the last class
 * c = NCLASS - 1, in at least NCLASS candidates.
 * Byte counters: byte i of words 8 * w ... 8 * w + 7 counts label 64 * w + i.
 * Bit-sliced counters: bit i of word t * W + w is set if label 64 * w + i is
 * in more than t candidates.
 */
static inline uint64_t COVERED(const uint64_t *g, int w) {
#ifdef COUNT_BYTES
  __m512i x = _mm512_loadu_si512(g + 8 * w);
  return _mm512_test_epi8_mask(x, x);
#else
  return g[w];
#endif
}

static inline uint64_t IN_CLASS(const uint64_t *g, int w, int c) {
#ifdef COUNT_BYTES
  __m512i x = _mm512_loadu_si512(g + 8 * w);
  if (c + 1 < NCLASS)
    return _mm512_cmpeq_epu8_mask(x, _mm512_set1_epi8((char)(c + 1)));
  return _mm512_cmpge_epu8_mask(x, _mm512_set1_epi8((char)NCLASS));
#else
  uint64_t r = g[c * W + w];
  if (c + 1 < NCLASS)
    r &= ~g[(c + 1) * W + w];
  return r;
#endif
}

/*
 * out = { u in in[0..cnt) : u >= min_v and bit u of row is set }, and the
 * counts g (see COVERED / IN_CLASS) of each label over out. Returns |out|.
 */
static inline uint32_t FILTER_COUNT(const sstate_t *s, const uint32_t *in,
                                    uint32_t cnt, const uint64_t *row,
                                    uint32_t min_v, uint32_t *out,
                                    uint64_t *g /* [CNT_WORDS * W] */) {
  uint32_t k = filter_list(in, cnt, row, min_v, out);
#ifdef COUNT_BYTES
  /* 4 vectors per iteration (fewer iterations, so fewer mispredicted loop
   * exits), padding the list with the empty vector N */
  const __m512i one = _mm512_set1_epi8(1);
  __m512i acc[W];
  for (int w = 0; w < W; w++)
    acc[w] = _mm512_setzero_si512();
  for (int j = 0; j < 4; j++)
    out[k + j] = s->N;
  for (uint32_t i = 0; i < k; i += 4)
    for (int j = 0; j < 4; j++) {
      const uint64_t *m = s->bits + (uint64_t)out[i + j] * W;
      for (int w = 0; w < W; w++)
        acc[w] = _mm512_mask_adds_epu8(acc[w], _cvtu64_mask64(m[w]), acc[w],
                                       one);
    }
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(g + 8 * w, acc[w]);
  return k;
#else
#ifdef HAVE_SIMD_COUNT_W8
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
#endif
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
      dead |= ncell_b[w] & ~ncell_o[w] & ~COVERED(go, w);
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
        dead |= ncell_o[w] & ~ncell_b[w] & ~COVERED(gb, w);
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
   * ..., or >= NCLASS), ties broken by the largest label */
  uint64_t urow[W], ucol[W];
  for (int w = 0; w < W; w++) {
    urow[w] = rc[w] & ~cc[w];
    ucol[w] = cc[w] & ~rc[w];
  }
  int x = -1;
  for (int c = s->opts.mrv ? 0 : NCLASS - 1; c < NCLASS && x < 0; c++) {
    for (int w = W - 1; w >= 0; w--) {
      uint64_t cls;
      if (!s->opts.mrv) {
        cls = urow[w] | ucol[w];
      } else {
        cls = (urow[w] & IN_CLASS(gr, w, c)) | (ucol[w] & IN_CLASS(gc, w, c));
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

  /* the children are the candidates on axis b through cell x; select them
   * all first (vectorized, no branch per candidate), then search them */
  uint32_t *kids = s->kids[d];
  uint32_t nk = filter_list(s->valid[d][b], s->nvalid[d][b],
                            s->has_label + (uint64_t)x * s->IW, 0, kids);
  for (uint32_t i = 0; i < nk && !s->stop; i++)
    TRY_CHILD(s, d, b, kids[i], 0);
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
#undef COVERED
#undef IN_CLASS
#undef COUNT_BYTES
#undef NCLASS
#undef CAT
#undef CAT_
