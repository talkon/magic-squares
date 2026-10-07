/*
 * Recursive search, specialized for a label-bitset width of W 64-bit words.
 * Included from arrange.c once per width.
 *
 * Invariant on entering SEARCH_REC(s, d): the candidate lists of depth d
 * (valid[d][a], or vw[d][a] when they carry their bitsets) and their label
 * counts cnt[d][a] are filled in, every unmatched cell has at least one
 * candidate, each axis has enough candidates left, and the lists are closed
 * under the support filter (see SUPPORT). From LAZY_DEPTH on (matrix path),
 * only the union of each list is filled in, and SEARCH_REC computes the
 * full counts itself; when the lists carry their bitsets, the unions are in
 * uni[d][a] and the counts are computed only for the children that survive
 * all the filters. The parent does this work while creating the child (see
 * TRY_CHILD), so that children which would be pruned are discarded as
 * cheaply as possible.
 */
#ifndef W
#error "define W before including arrange_core.h"
#endif

#define CAT_(a, b) a##_##b
#define CAT(a, b) CAT_(a, b)
#define SEARCH_REC CAT(search_rec, W)
#define SEARCH_ROOT CAT(search_root, W)
#define COUNT_LIST CAT(count_list, W)
#define UNION_LIST CAT(union_list, W)
#define STORE_UNION CAT(store_union, W)
#define KEEP_COUNT CAT(keep_count, W)
#define SUPPORT CAT(support, W)
#define FILTER_CARRY CAT(filter_carry, W)
#define KEEP_CARRY CAT(keep_carry, W)
#define COUNT_CARRY CAT(count_carry, W)
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
/* whether the candidate lists carry their bitsets (see vw in arrange.c) */
#if W <= CARRY_MAX_W
#define CARRY 1
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

#ifdef CARRY
/*
 * Candidate lists that carry their bitsets: word w of entry i of a list is
 * at list[w * cap + i].
 *
 * FILTER_CARRY: out = the entries u of in[0..cnt) with |u & v| == once
 * (once = 0 or 1), where vb is the bitset of v, and, for once = 0 and excl
 * not NULL, also disjoint from excl; un = the union of the kept entries.
 * Returns |out|.
 *
 * This replaces the N x N intersection matrices: |u & v| is computed from
 * the bitsets, 8 entries at a time, word w of the 8 entries being a single
 * contiguous load (no gather, no matrix row to fetch from L2 or beyond),
 * and the kept entries are compressed into the output with their bitsets.
 * The labels are not counted here: most children are pruned by the forward
 * check and the support filter right after this, which only need to know
 * which labels are covered, an OR of the kept entries. The counts are
 * computed afterwards for the children that survive (COUNT_CARRY).
 *
 * The support filter's first pass over the disjoint axis (excl = the labels
 * neither placed nor covered by the other axis) is free here: "disjoint
 * from v and from excl" is one test against v | excl.
 */
static inline uint32_t FILTER_CARRY(const sstate_t *s, const uint64_t *in,
                                    uint32_t cnt, const uint64_t *vb,
                                    const int once, const uint64_t *excl,
                                    uint64_t *out, uint64_t *un /* [W] */) {
  const size_t cap = s->cap;
  __m512i V[W], acc[W];
  /* with more than 2 words, skip the words where the test mask is empty
   * (v has at most 8 labels; with 2 words the test costs more than it
   * saves) */
  bool tw[W];
  for (int w = 0; w < W; w++) {
    uint64_t t = vb[w] | (!once && excl ? excl[w] : 0);
    V[w] = _mm512_set1_epi64((long long)t);
    acc[w] = _mm512_setzero_si512();
    tw[w] = W <= 2 || t;
  }
  uint32_t k = 0;
  for (uint32_t i = 0; i < cnt; i += 8) {
    __mmask8 live = cnt - i >= 8 ? 0xff : (__mmask8)((1u << (cnt - i)) - 1);
    __m512i x[W];
    for (int w = 0; w < W; w++)
      x[w] = _mm512_maskz_loadu_epi64(live, in + w * cap + i);
    __mmask8 keep;
    if (once) {
#ifdef __AVX512VPOPCNTDQ__
      __m512i c = _mm512_setzero_si512();
      for (int w = 0; w < W; w++)
        if (tw[w])
          c = _mm512_add_epi64(
              c, _mm512_popcnt_epi64(_mm512_and_si512(x[w], V[w])));
      keep = _mm512_mask_cmpeq_epi64_mask(live, c, _mm512_set1_epi64(1));
#else
      /* exactly one bit in all the words a = u & v: some bit, no word with
       * two bits (a & (a - 1) != 0), and no two nonzero words (their
       * minimum is nonzero) */
      __m512i any = _mm512_setzero_si512(), bad = _mm512_setzero_si512();
      for (int w = 0; w < W; w++)
        if (tw[w]) {
          __m512i a = _mm512_and_si512(x[w], V[w]);
          __m512i a1 = _mm512_sub_epi64(a, _mm512_set1_epi64(1));
          bad = _mm512_ternarylogic_epi64(bad, a, a1, 0xF8); /* A | B & C */
          bad = _mm512_or_si512(bad, _mm512_min_epu64(any, a));
          any = _mm512_or_si512(any, a);
        }
      keep = _mm512_mask_testn_epi64_mask(
          _mm512_mask_test_epi64_mask(live, any, any), bad, bad);
#endif
    } else {
      keep = live;
      for (int w = 0; w < W; w++)
        if (tw[w])
          keep = _mm512_mask_testn_epi64_mask(keep, x[w], V[w]);
    }
    for (int w = 0; w < W; w++) {
      /* compress in a register and store all 8 lanes (the lists have
       * room) */
      __m512i y = _mm512_maskz_compress_epi64(keep, x[w]);
      _mm512_storeu_si512(out + w * cap + k, y);
      acc[w] = _mm512_or_si512(acc[w], y);
    }
    k += (uint32_t)__builtin_popcount(keep);
  }
  for (int w = 0; w < W; w++)
    un[w] = (uint64_t)_mm512_reduce_or_epi64(acc[w]);
  return k;
}

/*
 * Keep the entries of list[0..cnt) disjoint from bad (compressed to the
 * front of the list, in place: the 8 lanes stored at position k <= i
 * overwrite only entries already loaded), set un to their union, and return
 * their number. This is the support filter's pass on carried bitsets: 8
 * entries per register per word, from contiguous loads, where a list of
 * indices loads the bitset of every entry from the vector table (KEEP_COUNT,
 * one scalar test per entry).
 */
static inline uint32_t KEEP_CARRY(const sstate_t *s, uint64_t *list,
                                  uint32_t cnt, const uint64_t *bad,
                                  uint64_t *un /* [W] */) {
  const size_t cap = s->cap;
  __m512i B[W], acc[W];
  bool tw[W];
  for (int w = 0; w < W; w++) {
    B[w] = _mm512_set1_epi64((long long)bad[w]);
    acc[w] = _mm512_setzero_si512();
    tw[w] = bad[w] != 0;
  }
  uint32_t k = 0;
  for (uint32_t i = 0; i < cnt; i += 8) {
    __mmask8 live = cnt - i >= 8 ? 0xff : (__mmask8)((1u << (cnt - i)) - 1);
    __m512i x[W];
    for (int w = 0; w < W; w++)
      x[w] = _mm512_maskz_loadu_epi64(live, list + w * cap + i);
    __mmask8 keep = live;
    for (int w = 0; w < W; w++)
      if (tw[w])
        keep = _mm512_mask_testn_epi64_mask(keep, x[w], B[w]);
    for (int w = 0; w < W; w++) {
      __m512i y = _mm512_maskz_compress_epi64(keep, x[w]);
      _mm512_storeu_si512(list + w * cap + k, y);
      acc[w] = _mm512_or_si512(acc[w], y);
    }
    k += (uint32_t)__builtin_popcount(keep);
  }
  for (int w = 0; w < W; w++)
    un[w] = (uint64_t)_mm512_reduce_or_epi64(acc[w]);
  return k;
}

/*
 * The label counts g (byte counters, see COVERED) over the entries of
 * list[0..cnt), for the words w with need[w] != 0 (the others are left at
 * 0: only the counts of unmatched cells are ever read).
 *
 * With VBMI, GFNI and BITALG, 8 entries at a time: the 8 x 64 bit matrix of
 * word w (one row per entry) is transposed by bytes with a vpermb, then
 * each 8 x 8 bit block by a gf2p8affineqb against the identity, after which
 * byte j holds the 8 entries' bits of label 64 w + j and a vpopcntb gives
 * its count: 4 instructions per word per 8 entries. Otherwise one masked
 * byte add per word per entry, with the word as the mask.
 */
static inline void COUNT_CARRY(const sstate_t *s, uint64_t *list,
                               uint32_t cnt, const uint64_t *need,
                               uint64_t *g /* [8 * W] */) {
  const size_t cap = s->cap;
#if defined(__AVX512VBMI__) && defined(__GFNI__) && defined(__AVX512BITALG__)
  __m512i acc[W];
  for (int w = 0; w < W; w++)
    acc[w] = _mm512_setzero_si512();
  /* byte 8 j + c of the transpose is byte 8 c + j of the input */
  const __m512i tr = _mm512_set_epi64(
      0x3f372f271f170f07, 0x3e362e261e160e06, 0x3d352d251d150d05,
      0x3c342c241c140c04, 0x3b332b231b130b03, 0x3a322a221a120a02,
      0x3931292119110901, 0x3830282018100800);
  /* byte j = 1 << j: gf2p8affineqb then gathers bit j of the 8 bytes of
   * each qword of the matrix operand into byte j */
  const __m512i id = _mm512_set1_epi64(0x8040201008040201);
  for (uint32_t i = 0; i < cnt; i += 8) {
    __mmask8 live = cnt - i >= 8 ? 0xff : (__mmask8)((1u << (cnt - i)) - 1);
    for (int w = 0; w < W; w++)
      if (need[w]) {
        __m512i y = _mm512_maskz_loadu_epi64(live, list + w * cap + i);
        y = _mm512_gf2p8affine_epi64_epi8(id, _mm512_permutexvar_epi8(tr, y),
                                          0);
        acc[w] = _mm512_adds_epu8(acc[w], _mm512_popcnt_epi8(y));
      }
  }
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(g + 8 * w, acc[w]);
#else
  /* 4 entries per iteration, padding the list with empty sets; even and
   * odd entries go to separate counters, to halve the chain of adds */
  const __m512i one = _mm512_set1_epi8(1);
  __m512i acc[2][W];
  for (int w = 0; w < W; w++) {
    acc[0][w] = acc[1][w] = _mm512_setzero_si512();
    for (int j = 0; j < 4; j++)
      list[w * cap + cnt + j] = 0;
  }
  for (uint32_t i = 0; i < cnt; i += 4)
    for (int j = 0; j < 4; j++)
      for (int w = 0; w < W; w++)
        if (need[w])
          acc[j % 2][w] =
              _mm512_mask_adds_epu8(acc[j % 2][w],
                                    _cvtu64_mask64(list[w * cap + i + j]),
                                    acc[j % 2][w], one);
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(g + 8 * w, _mm512_adds_epu8(acc[0][w], acc[1][w]));
#endif
}

#else /* !CARRY */

/*
 * The counts g (see COVERED / IN_CLASS) of each label over list[0..k). The
 * list must have room for 4 entries past its end.
 */
static inline void COUNT_LIST(const sstate_t *s, uint32_t *list, uint32_t k,
                              uint64_t *g /* [CNT_WORDS * W] */) {
#ifdef COUNT_BYTES
  /* 4 vectors per iteration (fewer iterations, so fewer mispredicted loop
   * exits), padding the list with the empty vector N. With W = 2, even and
   * odd entries go to separate counters, which halves the chain of dependent
   * adds (5% faster; no gain with more words) */
  const int nacc = W <= 2 ? 2 : 1;
  const __m512i one = _mm512_set1_epi8(1);
  __m512i acc[2][W];
  for (int w = 0; w < W; w++)
    acc[0][w] = acc[1][w] = _mm512_setzero_si512();
  for (int j = 0; j < 4; j++)
    list[k + j] = s->N;
  for (uint32_t i = 0; i < k; i += 4)
    for (int j = 0; j < 4; j++) {
      const uint64_t *m = s->bits + (uint64_t)list[i + j] * W;
      __m512i *a = acc[j % nacc];
      for (int w = 0; w < W; w++)
        a[w] = _mm512_mask_adds_epu8(a[w], _cvtu64_mask64(m[w]), a[w], one);
    }
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(g + 8 * w, nacc == 2
                                       ? _mm512_adds_epu8(acc[0][w], acc[1][w])
                                       : acc[0][w]);
#else
#ifdef HAVE_SIMD_COUNT_W8
  if (W == 8) {
    count_slices_simd_w8(s->bits, list, k, g);
    return;
  }
#endif
  /* thermometer code: gg[t] = labels in more than t of the vectors */
  uint64_t gg[NSLICE][W];
  for (int t = 0; t < NSLICE; t++)
    for (int w = 0; w < W; w++)
      gg[t][w] = 0;
  for (uint32_t i = 0; i < k; i++) {
    const uint64_t *m = s->bits + (uint64_t)list[i] * W;
    for (int w = 0; w < W; w++) {
      for (int t = NSLICE - 1; t > 0; t--)
        gg[t][w] |= gg[t - 1][w] & m[w];
      gg[0][w] |= m[w];
    }
  }
  for (int t = 0; t < NSLICE; t++)
    for (int w = 0; w < W; w++)
      g[t * W + w] = gg[t][w];
#endif
}

/*
 * Store the union u (W words) of a list as its counts g, for COVERED only:
 * slice 0 of bit-sliced counters, or byte counters of 0 or 255.
 */
static inline void STORE_UNION(uint64_t *g, const uint64_t *u) {
#ifdef COUNT_BYTES
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(g + 8 * w, _mm512_movm_epi8(_cvtu64_mask64(u[w])));
#else
  for (int w = 0; w < W; w++)
    g[w] = u[w];
#endif
}

/* only the union of list[0..k), as counts readable by COVERED */
static inline void UNION_LIST(const sstate_t *s, const uint32_t *list,
                              uint32_t k, uint64_t *g) {
#ifdef HAVE_SIMD_COUNT_W8
  if (W == 8) {
    union_simd_w8(s->bits, list, k, g);
    return;
  }
#endif
  uint64_t acc[W];
  for (int w = 0; w < W; w++)
    acc[w] = 0;
  for (uint32_t i = 0; i < k; i++) {
    const uint64_t *m = s->bits + (uint64_t)list[i] * W;
    for (int w = 0; w < W; w++)
      acc[w] |= m[w];
  }
  STORE_UNION(g, acc);
}

/*
 * Keep the vectors of list[0..k) disjoint from bad (moved to the front of
 * list, in place), compute their counts g (with bit-sliced counters and
 * union_only, only slice 0), and return their number. A vector meeting bad
 * is counted as the empty set, so there is no branch per vector.
 */
static inline uint32_t KEEP_COUNT(const sstate_t *s, uint32_t *list,
                                  uint32_t k, const uint64_t *bad,
                                  uint64_t *g, int union_only) {
  uint32_t kk = 0;
#ifdef COUNT_BYTES
  if (union_only) {
    uint64_t acc[W];
    for (int w = 0; w < W; w++)
      acc[w] = 0;
    for (uint32_t i = 0; i < k; i++) {
      uint32_t u = list[i];
      const uint64_t *m = s->bits + (uint64_t)u * W;
      uint64_t hit = 0;
      for (int w = 0; w < W; w++)
        hit |= m[w] & bad[w];
      const uint64_t keep = (uint64_t)0 - (hit == 0);
      for (int w = 0; w < W; w++)
        acc[w] |= m[w] & keep;
      list[kk] = u;
      kk += (uint32_t)(keep & 1);
    }
    STORE_UNION(g, acc);
    return kk;
  }
  const int nacc = W <= 2 ? 2 : 1;
  const __m512i one = _mm512_set1_epi8(1);
  __m512i acc[2][W];
  for (int w = 0; w < W; w++)
    acc[0][w] = acc[1][w] = _mm512_setzero_si512();
  for (uint32_t i = 0; i < k; i++) {
    uint32_t u = list[i];
    const uint64_t *m = s->bits + (uint64_t)u * W;
    uint64_t hit = 0;
    for (int w = 0; w < W; w++)
      hit |= m[w] & bad[w];
    const uint64_t keep = (uint64_t)0 - (hit == 0);
    __m512i *a = acc[i % nacc];
    for (int w = 0; w < W; w++)
      a[w] = _mm512_mask_adds_epu8(a[w], _cvtu64_mask64(m[w] & keep), a[w],
                                   one);
    list[kk] = u;
    kk += (uint32_t)(keep & 1);
  }
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(g + 8 * w, nacc == 2
                                       ? _mm512_adds_epu8(acc[0][w], acc[1][w])
                                       : acc[0][w]);
#else
#ifdef HAVE_SIMD_COUNT_W8
  if (W == 8)
    return union_only ? keep_union_simd_w8(s->bits, list, k, bad, g)
                      : keep_count_simd_w8(s->bits, list, k, bad, g);
#endif
  const int nsl = union_only ? 1 : NSLICE;
  uint64_t gg[NSLICE][W];
  for (int t = 0; t < nsl; t++)
    for (int w = 0; w < W; w++)
      gg[t][w] = 0;
  for (uint32_t i = 0; i < k; i++) {
    uint32_t u = list[i];
    const uint64_t *m = s->bits + (uint64_t)u * W;
    uint64_t hit = 0;
    for (int w = 0; w < W; w++)
      hit |= m[w] & bad[w];
    const uint64_t keep = (uint64_t)0 - (hit == 0);
    for (int w = 0; w < W; w++) {
      const uint64_t mw = m[w] & keep;
      for (int t = nsl - 1; t > 0; t--)
        gg[t][w] |= gg[t - 1][w] & mw;
      gg[0][w] |= mw;
    }
    list[kk] = u;
    kk += (uint32_t)(keep & 1);
  }
  for (int t = 0; t < nsl; t++)
    for (int w = 0; w < W; w++)
      g[t * W + w] = gg[t][w];
#endif
  return kk;
}
#endif /* CARRY */

/*
 * Support filter. In any completion of a node, every cell of a remaining col
 * lies in some row of the square: in a placed row, or in a remaining row,
 * which is one of the row candidates. So a candidate col with a cell that is
 * neither in a placed row nor in any candidate row can be dropped, and vice
 * versa. Forward checking asks this only of the unmatched cells; asking it
 * of every cell of every candidate (in particular of the cells in no placed
 * vector) is much stronger: about 85% of the nodes with two rows and two
 * cols placed that pass forward checking are dead by this test, and nearly
 * all of those one level deeper.
 *
 * The union of a list (COVERED) tells whether it needs filtering at all:
 * only when the union has a cell outside the other axis' support, and then
 * at least one candidate goes. Dropping candidates shrinks the support of
 * the other axis, so the two axes are filtered in turn, axis a first, until
 * neither changes. Filters valid[d][*] in place, recounts cnt[d][*] (only
 * the unions if union_only), and returns nonzero if the node is dead (an
 * unmatched cell without a candidate, or too few candidates). When the
 * lists carry their bitsets, filters vw[d][*] and updates only the unions
 * uni[d][*] (the counts are computed afterwards, if the node survives).
 */
static int SUPPORT(sstate_t *s, int d, int a, int union_only) {
  const int n = s->n;
  int quiet = 0; /* axes in a row found closed */
#ifdef CARRY
  /* the unions are uni[d][*] (the counts are computed afterwards) */
  (void)union_only;
#define UNION_OF(g, w) (g)[w]
#else
#define UNION_OF(g, w) COVERED(g, w)
#endif
  while (quiet < 2) {
    const uint64_t *C = s->cells[d][1 - a];
#ifdef CARRY
    const uint64_t *G = s->uni[d][1 - a];
    uint64_t *g = s->uni[d][a];
#else
    const uint64_t *G = s->cnt[d][1 - a];
    uint64_t *g = s->cnt[d][a];
#endif
    uint64_t bad[W], any = 0;
    for (int w = 0; w < W; w++) {
      bad[w] = UNION_OF(g, w) & ~(C[w] | UNION_OF(G, w));
      any |= bad[w];
    }
    if (!any) {
      quiet++;
    } else {
      quiet = 1;
#ifdef CARRY
      uint32_t k = KEEP_CARRY(s, s->vw[d][a], s->nvalid[d][a], bad, g);
#else
      uint32_t k =
          KEEP_COUNT(s, s->valid[d][a], s->nvalid[d][a], bad, g, union_only);
#endif
      s->nvalid[d][a] = k;
      if (s->np[a] + (int)k < n)
        return 1;
      /* forward checking: the unmatched cells of the other axis */
      const uint64_t *Ca = s->cells[d][a];
      uint64_t dead = 0;
      for (int w = 0; w < W; w++)
        dead |= C[w] & ~Ca[w] & ~UNION_OF(g, w);
      if (dead)
        return 1;
    }
    a = 1 - a;
  }
#undef UNION_OF
  return 0;
}

static void SEARCH_REC(sstate_t *s, int d);

/*
 * Create the child of depth d+1 obtained by placing v (index v, bitset m)
 * on axis b, and search it unless it is pruned. min_v restricts candidates
 * to indices >= min_v. When the lists carry their bitsets, the vectors are
 * identified by their bitsets alone, and v is unused.
 */
static inline void TRY_CHILD(sstate_t *s, int d, int b, const uint64_t *m,
                             uint32_t v, uint32_t min_v) {
  const int n = s->n;
  s->nodes++;
  if (s->opts.node_limit && s->nodes > s->opts.node_limit) {
    s->stop = 2;
    return;
  }
  const uint64_t *rc = s->cells[d][ROW], *cc = s->cells[d][COL];
  uint64_t *nrc = s->cells[d + 1][ROW], *ncc = s->cells[d + 1][COL];
  for (int w = 0; w < W; w++) {
    nrc[w] = rc[w] | (b == ROW ? m[w] : 0);
    ncc[w] = cc[w] | (b == COL ? m[w] : 0);
  }
#ifdef CARRY
  (void)v;
  for (int w = 0; w < W; w++)
    s->placedw[b][s->np[b]][w] = m[w];
  s->np[b]++;
  if (s->np[ROW] == n && s->np[COL] == n) {
    report_bits(s, W);
    s->np[b]--;
    return;
  }
#else
  s->placed[b][s->np[b]++] = v;
  if (s->np[ROW] == n && s->np[COL] == n) {
    report(s);
    s->np[b]--;
    return;
  }
#endif
  const uint64_t *ncell_b = b == ROW ? nrc : ncc;
  const uint64_t *ncell_o = b == ROW ? ncc : nrc;
  const int sup = s->opts.support && s->opts.forward_check;
  const int o = 1 - b;
  uint64_t dead = 0;

  /* candidates on the other axis must meet v exactly once; they cover the
   * unmatched cells on axis b (including v's own new cells), which is where
   * most dead ends show up, so check those first. Then the candidates on
   * axis b must be disjoint from v; they cover the unmatched cells on the
   * other axis. With the support filter, they must also lie within the
   * placed vectors and the candidates of the other axis, which costs
   * nothing extra there (one test against v | excl when the lists carry
   * their bitsets, folded into the counting pass otherwise). */
#ifdef CARRY
  /* the lists at depth 0 are all the vectors in order, so the entries
   * >= min_v are the ones from position min_v on (min_v = 0 below the
   * root) */
  uint64_t *uo = s->uni[d + 1][o], *ub = s->uni[d + 1][b];
  uint32_t ko = FILTER_CARRY(s, s->vw[d][o] + min_v, s->nvalid[d][o] - min_v,
                             m, 1, NULL, s->vw[d + 1][o], uo);
  s->nvalid[d + 1][o] = ko;
  if (s->np[o] + (int)ko < n)
    dead = 1;
  if (s->opts.forward_check)
    for (int w = 0; w < W; w++)
      dead |= ncell_b[w] & ~ncell_o[w] & ~uo[w];
  if (!dead) {
    uint64_t excl[W];
    for (int w = 0; w < W; w++)
      excl[w] = ~(ncell_o[w] | uo[w]);
    uint32_t kb = FILTER_CARRY(s, s->vw[d][b] + min_v,
                               s->nvalid[d][b] - min_v, m, 0,
                               sup ? excl : NULL, s->vw[d + 1][b], ub);
    s->nvalid[d + 1][b] = kb;
    if (s->np[b] + (int)kb < n)
      dead = 1;
    if (s->opts.forward_check)
      for (int w = 0; w < W; w++)
        dead |= ncell_o[w] & ~ncell_b[w] & ~ub[w];
    /* the b candidates are closed already, the o candidates may not be */
    if (!dead && sup)
      dead = SUPPORT(s, d + 1, o, 1);
    if (!dead) {
      /* the child survives: count the labels of its lists, for the cells
       * each list covers (the unmatched cells of the other axis); only
       * MRV reads the counts */
      if (s->opts.mrv) {
        uint64_t need_o[W], need_b[W];
        for (int w = 0; w < W; w++) {
          need_o[w] = ncell_b[w] & ~ncell_o[w];
          need_b[w] = ncell_o[w] & ~ncell_b[w];
        }
        COUNT_CARRY(s, s->vw[d + 1][o], s->nvalid[d + 1][o], need_o,
                    s->cnt[d + 1][o]);
        COUNT_CARRY(s, s->vw[d + 1][b], s->nvalid[d + 1][b], need_b,
                    s->cnt[d + 1][b]);
      }
      SEARCH_REC(s, d + 1);
    }
  }
#else
  /* deep children are mostly pruned by the support filter, which needs only
   * the unions: leave the full counts to SEARCH_REC */
  const int lazy = sup && d + 1 >= LAZY_DEPTH;
  uint64_t *go = s->cnt[d + 1][o];
  uint32_t *lo = s->valid[d + 1][o];
  uint32_t ko = filter_list(s->valid[d][o], s->nvalid[d][o],
                            s->inters1 + (uint64_t)v * s->IW, min_v, s->nseg,
                            lo);
  if (lazy)
    UNION_LIST(s, lo, ko, go);
  else
    COUNT_LIST(s, lo, ko, go);
  s->nvalid[d + 1][o] = ko;
  if (s->np[o] + (int)ko < n)
    dead = 1;
  if (s->opts.forward_check)
    for (int w = 0; w < W; w++)
      dead |= ncell_b[w] & ~ncell_o[w] & ~COVERED(go, w);
  if (!dead) {
    uint64_t *gb = s->cnt[d + 1][b];
    uint32_t *lb = s->valid[d + 1][b];
    uint32_t kb = filter_list(s->valid[d][b], s->nvalid[d][b],
                              s->inters0 + (uint64_t)v * s->IW, min_v, s->nseg,
                              lb);
    if (sup) {
      uint64_t bad[W];
      for (int w = 0; w < W; w++)
        bad[w] = ~(ncell_o[w] | COVERED(go, w));
      kb = KEEP_COUNT(s, lb, kb, bad, gb, lazy);
    } else {
      COUNT_LIST(s, lb, kb, gb);
    }
    s->nvalid[d + 1][b] = kb;
    if (s->np[b] + (int)kb < n)
      dead = 1;
    if (s->opts.forward_check)
      for (int w = 0; w < W; w++)
        dead |= ncell_o[w] & ~ncell_b[w] & ~COVERED(gb, w);
    /* the b candidates are closed already, the o candidates may not be */
    if (!dead && sup)
      dead = SUPPORT(s, d + 1, o, lazy);
    if (!dead)
      SEARCH_REC(s, d + 1);
  }
#endif
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
#ifndef CARRY
  if (d >= LAZY_DEPTH && s->opts.support && s->opts.forward_check &&
      s->opts.mrv) {
    COUNT_LIST(s, s->valid[d][ROW], s->nvalid[d][ROW], s->cnt[d][ROW]);
    COUNT_LIST(s, s->valid[d][COL], s->nvalid[d][COL], s->cnt[d][COL]);
  }
#endif
  const uint64_t *rc = s->cells[d][ROW], *cc = s->cells[d][COL];
  /* candidates on axis a cover unmatched cells of axis 1-a, so cells of
   * placed rows are classified by the counts over candidate cols */
  const uint64_t *gr = s->cnt[d][COL], *gc = s->cnt[d][ROW];

  /* branch on the unmatched cell with the fewest candidates (exactly 1, 2,
   * ..., or >= NCLASS), ties broken by the smallest label (a common number;
   * ~0.5% fewer nodes than the largest with labels by plain frequency, the
   * same with those of assign_labels). Among cells with >= NCLASS
   * candidates, and without MRV, by the largest label instead, which is a
   * much better stand-in for the fewest candidates there (with 8 bit-sliced
   * classes, the smallest label gives 4x more nodes). */
  uint64_t urow[W], ucol[W];
  for (int w = 0; w < W; w++) {
    urow[w] = rc[w] & ~cc[w];
    ucol[w] = cc[w] & ~rc[w];
  }
  int x = -1;
  if (s->opts.mrv)
    for (int c = 0; c < NCLASS - 1 && x < 0; c++)
      for (int w = 0; w < W; w++) {
        uint64_t cls =
            (urow[w] & IN_CLASS(gr, w, c)) | (ucol[w] & IN_CLASS(gc, w, c));
        if (cls) {
          x = w * 64 + __builtin_ctzll(cls);
          break;
        }
      }
  for (int w = W - 1; w >= 0 && x < 0; w--) {
    uint64_t cls = urow[w] | ucol[w];
    if (s->opts.mrv)
      cls = (urow[w] & IN_CLASS(gr, w, NCLASS - 1)) |
            (ucol[w] & IN_CLASS(gc, w, NCLASS - 1));
    if (cls)
      x = w * 64 + 63 - __builtin_clzll(cls);
  }
  if (x < 0)
    return;
  const int xw = x / 64;
  const uint64_t xb = (uint64_t)1 << (x % 64);
  /* a row-unmatched cell is covered by a col, and vice versa */
  const int b = (urow[xw] & xb) ? COL : ROW;

  /* the children are the candidates on axis b through cell x; select them
   * all first (vectorized, no branch per candidate), then search them */
#ifdef CARRY
  /* the entries whose word xw has bit x, with their bitsets */
  uint64_t *kids = s->kidw[d];
  const size_t cap = s->cap;
  const uint32_t cnt = s->nvalid[d][b];
  const uint64_t *in = s->vw[d][b];
  const __m512i xm = _mm512_set1_epi64((long long)xb);
  uint32_t nk = 0;
  for (uint32_t i = 0; i < cnt; i += 8) {
    __mmask8 live = cnt - i >= 8 ? 0xff : (__mmask8)((1u << (cnt - i)) - 1);
    __mmask8 keep = _mm512_mask_test_epi64_mask(
        live, _mm512_maskz_loadu_epi64(live, in + xw * cap + i), xm);
    /* (unmasked loads, which need not wait for keep: they may read up to 7
     * entries past the list, inside its array, which are then dropped) */
    for (int w = 0; w < W; w++)
      _mm512_storeu_si512(kids + w * cap + nk,
                          _mm512_maskz_compress_epi64(
                              keep, _mm512_loadu_si512(in + w * cap + i)));
    nk += (uint32_t)__builtin_popcount(keep);
  }
  for (uint32_t i = 0; i < nk && !s->stop; i++) {
    uint64_t m[W];
    for (int w = 0; w < W; w++)
      m[w] = kids[w * cap + i];
    TRY_CHILD(s, d, b, m, 0, 0);
  }
#else
  uint32_t *kids = s->kids[d];
  uint32_t nk = filter_list(s->valid[d][b], s->nvalid[d][b],
                            s->has_label + (uint64_t)x * s->IW, 0, s->nseg,
                            kids);
  for (uint32_t i = 0; i < nk && !s->stop; i++)
    TRY_CHILD(s, d, b, s->bits + (uint64_t)kids[i] * W, kids[i], 0);
#endif
}

/* the first row r1 is the lowest-index vector of the square, so it contains
 * the square's largest label x (as does c1, the col through x); every other
 * vector has index > r1, so labels < x (see assign_labels in arrange.c).
 * MRV branches on x at depth 1 for ~95% of the r1 that survive; always
 * branching on x there (choosing c1 right after r1) is no better (+0.1%
 * nodes; +2% with labels by plain frequency, where MRV picks x ~75% of the
 * time). */
static void SEARCH_ROOT(sstate_t *s) {
  s->nodes = 1;
  uint32_t count = s->N;
  for (int a = 0; a < 2; a++) {
#ifdef CARRY
    for (uint32_t i = 0; i < count; i++)
      for (int w = 0; w < W; w++)
        s->vw[0][a][w * s->cap + i] = s->bits[(uint64_t)i * W + w];
#else
    for (uint32_t i = 0; i < count; i++)
      s->valid[0][a][i] = i;
#endif
    s->nvalid[0][a] = count;
    for (int w = 0; w < W; w++)
      s->cells[0][a][w] = 0;
  }
  for (uint32_t r1 = 0; r1 < count && !s->stop; r1++)
    TRY_CHILD(s, 0, ROW, s->bits + (uint64_t)r1 * W, r1, r1 + 1);
}

#undef SEARCH_REC
#undef SEARCH_ROOT
#undef COUNT_LIST
#undef UNION_LIST
#undef STORE_UNION
#undef KEEP_COUNT
#undef SUPPORT
#undef FILTER_CARRY
#undef KEEP_CARRY
#undef COUNT_CARRY
#undef TRY_CHILD
#undef COVERED
#undef IN_CLASS
#undef COUNT_BYTES
#undef CARRY
#undef NCLASS
#undef CAT
#undef CAT_
