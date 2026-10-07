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
#define FILTER_STEP CAT(filter_step, W)
#define KEEP_CARRY CAT(keep_carry, W)
#define COUNT_CARRY CAT(count_carry, W)
#define TRY_CHILD CAT(try_child, W)
#define COVERED CAT(covered, W)
#define IN_CLASS CAT(in_class, W)
#define CROSS_AXIS CAT(cross_axis, W)
#define CROSS_TEST CAT(cross_test, W)
#define CROSS CAT(cross, W)
#define SUPPORT_AGAIN CAT(support_again, W)

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

/* on the matrix path, the cross support filter (CROSS) pays off with the
 * AVX-512 code below and byte counters (up to 256 labels), for n <= 6; in
 * plain C, or with wider bitsets and bit-sliced counters, it cost ~10% more
 * time than it saved, and for n = 7 (where the search dies further down)
 * 10-25% more, so there it only runs when forced (opts.cross = 2, to test
 * that code). With carried bitsets it runs for n <= 6 (n = 7: 7% fewer
 * nodes, 33% more time). */
#if defined(COUNT_BYTES) && defined(__BMI2__) && !defined(CARRY)
#define CROSS_SIMD 1
#endif
/* the most cells y in a cross pass (2 (n - 2) for the (2,2) nodes) */
#define CROSS_MAXY 16
/* TRY_CHILD is inlined into SEARCH_REC on the carried path (3-6% faster on
 * bench/quick.txt and full.txt, with or without cross support; GCC stopped
 * inlining it for W = 3 when cross support was added, 5% slower there), not
 * with the matrices (7% slower with W = 2) */
#ifdef CARRY
#define TRY_CHILD_INLINE __attribute__((always_inline))
#else
#define TRY_CHILD_INLINE
#endif

/*
 * Reading the label counts g = cnt[d][a] (CNT_WORDS * W words), for the
 * labels of word w:
 * COVERED: labels in at least one candidate;
 * IN_CLASS (bit-sliced counters; SEARCH_REC reads byte counters directly):
 * labels in exactly c + 1 candidates, or, for the last class
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

#ifndef COUNT_BYTES
static inline uint64_t IN_CLASS(const uint64_t *g, int w, int c) {
  uint64_t r = g[c * W + w];
  if (c + 1 < NCLASS)
    r &= ~g[(c + 1) * W + w];
  return r;
}
#endif

#ifdef CARRY
/*
 * Candidate lists that carry their bitsets: word w of entry i of a list is
 * at list[w * cap + i].
 *
 * FILTER_STEP: the test and compression of FILTER_CARRY and KEEP_CARRY for
 * the 8 entries of in from position i on, of which those in live exist;
 * appends the kept ones (word w at out + w * cap + k) and returns the new
 * k. full (a compile-time constant) says that all 8 exist: the callers run
 * it on the full groups of 8 entries, then once with a lane mask on the
 * last one, rather than computing a lane mask in every iteration (~8
 * scalar instructions and masked loads: the loop split of gcc's
 * profile-guided build, 4-5% less time on bench/prod.txt, see
 * research/ideas.md). For once = 1, V[w] and R[w] are word w of v and its
 * rotation (see FILTER_CARRY); otherwise V[w] is word w of the bitset that
 * the kept entries miss (R is unused).
 */
static inline __attribute__((always_inline)) uint32_t
FILTER_STEP(const uint64_t *in, size_t cap, uint32_t i, const int full,
            __mmask8 live, const __m512i *V, const __m512i *R,
            const bool *tw, const int once, uint64_t *out, uint32_t k,
            __m512i *acc) {
  __m512i x[W];
  for (int w = 0; w < W; w++)
    x[w] = full ? _mm512_loadu_si512(in + w * cap + i)
                : _mm512_maskz_loadu_epi64(live, in + w * cap + i);
  __mmask8 keep;
  if (once) {
    /* c = the words of u & v, each rotated by R[w], or-ed into one word:
     * the rotated words of v are disjoint, so |c| = |u & v| */
    __m512i c = _mm512_and_si512(x[0], V[0]);
    for (int w = 1; w < W; w++)
      if (tw[w]) /* c |= rot(u) & rot(v) */
        c = _mm512_ternarylogic_epi64(c, _mm512_rolv_epi64(x[w], R[w]), V[w],
                                      0xF8);
#ifdef __AVX512VPOPCNTDQ__
    keep = _mm512_mask_cmpeq_epi64_mask(live, _mm512_popcnt_epi64(c),
                                        _mm512_set1_epi64(1));
#else
    /* exactly one bit: c != 0 and c & (c - 1) == 0 */
    keep = _mm512_mask_testn_epi64_mask(
        _mm512_mask_test_epi64_mask(live, c, c), c,
        _mm512_sub_epi64(c, _mm512_set1_epi64(1)));
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
  return k + (uint32_t)__builtin_popcount(keep);
}

/*
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
 *
 * For |u & v| == 1, word w of v is rotated by some r_w such that the
 * rotated words are disjoint: then |u & v| is the number of bits of one
 * word, the or of the words of u & v rotated the same way, so the test is
 * one popcount, or one exactly-one-bit test without VPOPCNTDQ, instead of
 * one per word and a sum or cross-word test (W = 2: 5 vector instructions
 * per 8 entries instead of 6 with VPOPCNTDQ, 6 instead of 11 without it,
 * which made the -march=cascadelake build ~5% faster on bench/prod.txt).
 * Such an r_w exists: v has n <= 8 labels in all, so at most 16 of the 64
 * rotations of a word hit the bits taken by the words before it.
 */
static inline uint32_t FILTER_CARRY(const sstate_t *s, const uint64_t *in,
                                    uint32_t cnt, const uint64_t *vb,
                                    const int once, const uint64_t *excl,
                                    uint64_t *out, uint64_t *un /* [W] */) {
  const size_t cap = s->cap;
  __m512i V[W], R[W], acc[W];
  /* with more than 2 words, skip the words where the test mask is empty
   * (v has at most 8 labels; with 2 words the test costs more than it
   * saves) */
  bool tw[W];
  uint64_t used = 0; /* the bits of the rotated words of v so far */
  for (int w = 0; w < W; w++) {
    uint64_t t = vb[w] | (!once && excl ? excl[w] : 0);
    int r = 0;
    if (once && w > 0)
      while (rotl64(t, r) & used)
        r++;
    t = rotl64(t, r);
    used |= t;
    V[w] = _mm512_set1_epi64((long long)t);
    R[w] = _mm512_set1_epi64(r);
    acc[w] = _mm512_setzero_si512();
    tw[w] = W <= 2 || t;
  }
  /* the full groups of 8 entries, then the last one (see FILTER_STEP) */
  uint32_t k = 0, i = 0;
  for (; i + 8 <= cnt; i += 8)
    k = FILTER_STEP(in, cap, i, 1, 0xff, V, R, tw, once, out, k, acc);
  if (i < cnt)
    k = FILTER_STEP(in, cap, i, 0, (__mmask8)((1u << (cnt - i)) - 1), V, R,
                    tw, once, out, k, acc);
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
  /* the test of FILTER_CARRY with once = 0, full groups of 8 entries
   * first (see FILTER_STEP) */
  uint32_t k = 0, i = 0;
  for (; i + 8 <= cnt; i += 8)
    k = FILTER_STEP(list, cap, i, 1, 0xff, B, NULL, tw, 0, list, k, acc);
  if (i < cnt)
    k = FILTER_STEP(list, cap, i, 0, (__mmask8)((1u << (cnt - i)) - 1), B,
                    NULL, tw, 0, list, k, acc);
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
  /* full groups of 8 entries, then the last one with a lane mask (see
   * FILTER_STEP) */
  uint32_t i = 0;
  for (; i + 8 <= cnt; i += 8)
    for (int w = 0; w < W; w++)
      if (need[w]) {
        __m512i y = _mm512_loadu_si512(list + w * cap + i);
        y = _mm512_gf2p8affine_epi64_epi8(id, _mm512_permutexvar_epi8(tr, y),
                                          0);
        acc[w] = _mm512_adds_epu8(acc[w], _mm512_popcnt_epi8(y));
      }
  if (i < cnt)
    for (int w = 0; w < W; w++)
      if (need[w]) {
        __m512i y = _mm512_maskz_loadu_epi64(
            (__mmask8)((1u << (cnt - i)) - 1), list + w * cap + i);
        y = _mm512_gf2p8affine_epi64_epi8(id, _mm512_permutexvar_epi8(tr, y),
                                          0);
        acc[w] = _mm512_adds_epu8(acc[w], _mm512_popcnt_epi8(y));
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

/*
 * Cross support. Let y be an unmatched cell of a placed vector of axis f, a
 * cell in no placed vector of the other axis h. In any completion of the
 * node, y lies in one remaining vector of axis h, which is one of the
 * candidates of axis h through y, and every remaining vector of axis f
 * meets it. So a candidate of axis f that misses U_y, the union of the
 * candidates of axis h through y, is in no completion.
 *
 * The support filter asks every cell of a candidate to be covered by a
 * candidate of the other axis; this asks every candidate to cross each of
 * the remaining vectors of the other axis, which the unmatched cells pin
 * down to a few candidates each. It is a relaxation of "some candidate
 * through y meets the candidate exactly once", which kills hardly more but
 * needs a pairwise test. See TRY_CHILD for where it runs.
 */
#ifdef CARRY
/*
 * The test of a cross pass on carried lists: keeps the entries of
 * fl[0..kf) that meet U_y for every cell y and are disjoint from bad
 * (compressed in place, as KEEP_CARRY), sets un to their union and returns
 * their number. Ub holds the U_y of ng groups of 8 cells, word w of U_y of
 * cell 8 g + j at Ub[8 W g + 8 w + j]. 8 entries at a time: t_y = u & U_y
 * (an and per word), and the entry is kept if the minimum of the t_y is
 * nonzero.
 */
static inline __attribute__((always_inline)) uint32_t
CROSS_TEST(const sstate_t *s, uint64_t *fl, uint32_t kf, const uint64_t *Ub,
           const int ng, const uint64_t *bad, uint64_t *un /* [W] */) {
  const size_t cap = s->cap;
  __m512i acc[W], Bad[W];
  for (int w = 0; w < W; w++) {
    acc[w] = _mm512_setzero_si512();
    Bad[w] = _mm512_set1_epi64((long long)bad[w]);
  }
  /* the list padded with empty sets, which fail every U_y (the lists have
   * room), so that no lane mask is needed */
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(fl + w * cap + kf, _mm512_setzero_si512());
  uint32_t k = 0;
  for (uint32_t i = 0; i < kf; i += 8) {
    __m512i x[W];
    for (int w = 0; w < W; w++)
      x[w] = _mm512_loadu_si512(fl + w * cap + i);
    __mmask8 keep = _mm512_testn_epi64_mask(x[0], Bad[0]);
    for (int w = 1; w < W; w++)
      keep = _mm512_mask_testn_epi64_mask(keep, x[w], Bad[w]);
    for (int g = 0; g < ng; g++) {
      __m512i t[8];
      for (int j = 0; j < 8; j++) {
        const uint64_t *u = Ub + 8 * W * g + j;
        t[j] = _mm512_and_si512(x[0], _mm512_set1_epi64((long long)u[0]));
        for (int w = 1; w < W; w++) /* t |= x & U */
          t[j] = _mm512_ternarylogic_epi64(
              t[j], x[w], _mm512_set1_epi64((long long)u[8 * w]), 0xF8);
      }
      /* (a tree of minimums, rather than a chain of masked tests: 2-5%
       * faster) */
      for (int j = 0; j < 4; j++)
        t[j] = _mm512_min_epu64(t[j], t[j + 4]);
      t[0] = _mm512_min_epu64(_mm512_min_epu64(t[0], t[2]),
                              _mm512_min_epu64(t[1], t[3]));
      keep = _mm512_mask_test_epi64_mask(keep, t[0], t[0]);
    }
    for (int w = 0; w < W; w++) {
      __m512i y = _mm512_maskz_compress_epi64(keep, x[w]);
      _mm512_storeu_si512(fl + w * cap + k, y);
      acc[w] = _mm512_or_si512(acc[w], y);
    }
    k += (uint32_t)__builtin_popcount(keep);
  }
  for (int w = 0; w < W; w++)
    un[w] = (uint64_t)_mm512_reduce_or_epi64(acc[w]);
  return k;
}

/*
 * One cross pass on carried lists: drops the candidates of axis f that
 * miss U_y for some unmatched cell y of the placed vectors of axis f, and
 * those meeting bad (the support filter's test, for free), in place; sets
 * uni[d][f] and nvalid[d][f], and returns the number of candidates left.
 *
 * The U_y are built 8 cells at a time: for 8 entries of the list of axis
 * h = 1 - f, one test per cell y gives the entries through y, which are
 * or-ed into that cell's accumulators (one per word); the 8 accumulators
 * of a word are then reduced at once, lane j getting U_y of the j-th cell
 * (or_lanes8). Then CROSS_TEST.
 */
static inline uint32_t CROSS_AXIS(sstate_t *s, int d, int f,
                                  const uint64_t *bad) {
  const int h = 1 - f;
  const size_t cap = s->cap;
  const uint64_t *Cf = s->cells[d][f], *Ch = s->cells[d][h];
  uint64_t Y[W]; /* the cells y */
  int ny = 0, nw[W]; /* nw[w]: cells in the words before w */
  for (int w = 0; w < W; w++) {
    Y[w] = Cf[w] & ~Ch[w];
    nw[w] = ny;
    ny += __builtin_popcountll(Y[w]);
  }
  const uint32_t kf = s->nvalid[d][f];
  if (ny == 0 || ny > CROSS_MAXY)
    return kf;
  const int ng = (ny + 7) / 8;
  /* all ones past the cells, which then pass */
  uint64_t Ub[CROSS_MAXY * W] __attribute__((aligned(64)));
  uint64_t *hl = s->vw[d][h];
  const uint32_t kh = s->nvalid[d][h];
  /* padded with empty sets, which are in no cell */
  for (int w = 0; w < W; w++)
    _mm512_storeu_si512(hl + w * cap + kh, _mm512_setzero_si512());
  for (int g = 0; g < ng; g++) {
    const uint64_t *base[8];
    __m512i B[8], acc[8][W];
#pragma GCC unroll 8
    for (int j = 0; j < 8; j++) {
      /* cell 8 g + j: its word w and bit (with pdep, no branch) */
      const int c = 8 * g + j;
      int w = 0;
      for (int v = 1; v < W; v++)
        w = c >= nw[v] ? v : w;
      base[j] = hl + w * cap;
      const uint64_t bit = c < ny ? nth_bit(Y[w], c - nw[w]) : 0;
      B[j] = _mm512_set1_epi64((long long)bit);
#pragma GCC unroll 8
      for (int v = 0; v < W; v++)
        acc[j][v] = _mm512_setzero_si512();
    }
    for (uint32_t i = 0; i < kh; i += 8) {
      __m512i x[W];
      for (int w = 0; w < W; w++)
        x[w] = _mm512_loadu_si512(hl + w * cap + i);
      for (int j = 0; j < 8; j++) {
        const __mmask8 m =
            _mm512_test_epi64_mask(_mm512_loadu_si512(base[j] + i), B[j]);
        for (int w = 0; w < W; w++)
          acc[j][w] = _mm512_mask_or_epi64(acc[j][w], m, acc[j][w], x[w]);
      }
    }
    const __mmask8 pad =
        ny - 8 * g >= 8 ? 0 : (__mmask8)(0xff << (ny - 8 * g));
    for (int w = 0; w < W; w++) {
      __m512i a[8];
      for (int j = 0; j < 8; j++)
        a[j] = acc[j][w];
      _mm512_store_si512(Ub + 8 * W * g + 8 * w,
                         _mm512_mask_set1_epi64(or_lanes8(a), pad, -1));
    }
  }
  uint64_t *fl = s->vw[d][f];
  uint64_t *un = s->uni[d][f];
  /* (specialized for the common single group, whose U_y stay in
   * registers) */
  const uint32_t k = ng == 1 ? CROSS_TEST(s, fl, kf, Ub, 1, bad, un)
                             : CROSS_TEST(s, fl, kf, Ub, ng, bad, un);
  s->nvalid[d][f] = k;
  return k;
}

/*
 * Cross support on both axes, axis `first` first, each pass with the
 * support filter's test against the current union of the other axis, and
 * forward checking after each. Returns nonzero if the node is dead.
 */
static int CROSS(sstate_t *s, int d, int first) {
  const int n = s->n;
  for (int p = 0; p < 2; p++) {
    const int f = p ? 1 - first : first;
    const int h = 1 - f;
    const uint64_t *Cf = s->cells[d][f], *Ch = s->cells[d][h];
    uint64_t bad[W];
    for (int w = 0; w < W; w++)
      bad[w] = ~(Ch[w] | s->uni[d][h][w]);
    const uint32_t k0 = s->nvalid[d][f];
    const uint32_t k = CROSS_AXIS(s, d, f, bad);
    if (k == k0)
      continue;
    if (s->np[f] + (int)k < n)
      return 1;
    /* forward checking: the cells that the f candidates cover */
    uint64_t dead = 0;
    for (int w = 0; w < W; w++)
      dead |= Ch[w] & ~Cf[w] & ~s->uni[d][f][w];
    if (dead)
      return 1;
  }
  return 0;
}

#else /* !CARRY */

/*
 * One pass on lists of indices: drops the candidates of axis o = 1 - b that
 * miss some U_y (y: the unmatched cells of the placed vectors of axis o),
 * in place (counts not updated); returns nonzero if any was dropped.
 */
static int CROSS_AXIS(sstate_t *s, int d, int b) {
  const int o = 1 - b;
  const uint64_t *Co = s->cells[d][o], *Cb = s->cells[d][b];
  uint64_t Y[W]; /* the cells y */
  int ny = 0, off[W];
  for (int w = 0; w < W; w++) {
    Y[w] = Co[w] & ~Cb[w];
    off[w] = ny; /* index of the first cell of word w among the cells */
    ny += __builtin_popcountll(Y[w]);
  }
  if (ny == 0)
    return 0;
  const uint32_t *Lb = s->valid[d][b];
  const uint32_t kb = s->nvalid[d][b];
  uint32_t *Lo = s->valid[d][o];
  const uint32_t ko = s->nvalid[d][o];
  uint32_t kk = 0;
#ifdef CROSS_SIMD
  if (ny <= 16) {
    /* lane j of U[w][z] is word w of U_y for the (8z + j)-th cell y. Each
     * candidate of axis b is or-ed into the lanes of its cells of Y, which
     * pext extracts from its bitset as a mask of cell indices; then each
     * candidate of axis o is tested against all the U_y at once. */
    __m512i U[W][2];
    for (int w = 0; w < W; w++)
      U[w][0] = U[w][1] = _mm512_setzero_si512();
    for (uint32_t i = 0; i < kb; i++) {
      const uint64_t *vb = s->bits + (uint64_t)Lb[i] * W;
      uint32_t cm = 0;
      for (int w = 0; w < W; w++)
        cm |= (uint32_t)_pext_u64(vb[w], Y[w]) << off[w];
      for (int w = 0; w < W; w++) {
        /* broadcast from memory: a load, not a shuffle */
        const __m512i x = _mm512_broadcastq_epi64(
            _mm_loadl_epi64((const __m128i *)(vb + w)));
        U[w][0] = _mm512_mask_or_epi64(U[w][0], (__mmask8)cm, U[w][0], x);
        if (ny > 8)
          U[w][1] =
              _mm512_mask_or_epi64(U[w][1], (__mmask8)(cm >> 8), U[w][1], x);
      }
    }
    const __mmask8 all0 = ny >= 8 ? 0xff : (__mmask8)((1u << ny) - 1);
    const __mmask8 all1 = ny > 8 ? (__mmask8)((1u << (ny - 8)) - 1) : 0;
    for (uint32_t i = 0; i < ko; i++) {
      const uint32_t u = Lo[i];
      const uint64_t *ub = s->bits + (uint64_t)u * W;
      __m512i a0 =
          _mm512_and_si512(U[0][0], _mm512_set1_epi64((long long)ub[0]));
      for (int w = 1; w < W; w++) /* a0 |= U & u */
        a0 = _mm512_ternarylogic_epi64(
            U[w][0], _mm512_set1_epi64((long long)ub[w]), a0, 0xEA);
      __mmask8 miss = _mm512_mask_testn_epi64_mask(all0, a0, a0);
      if (ny > 8) {
        __m512i a1 =
            _mm512_and_si512(U[0][1], _mm512_set1_epi64((long long)ub[0]));
        for (int w = 1; w < W; w++)
          a1 = _mm512_ternarylogic_epi64(
              U[w][1], _mm512_set1_epi64((long long)ub[w]), a1, 0xEA);
        miss |= _mm512_mask_testn_epi64_mask(all1, a1, a1);
      }
      Lo[kk] = u;
      kk += miss == 0;
    }
  } else
#endif
  {
    uint64_t U[ny][W];
    for (int j = 0; j < ny; j++)
      for (int w = 0; w < W; w++)
        U[j][w] = 0;
    for (uint32_t i = 0; i < kb; i++) {
      const uint64_t *vb = s->bits + (uint64_t)Lb[i] * W;
      for (int w = 0; w < W; w++)
        for (uint64_t m = vb[w] & Y[w]; m; m &= m - 1) {
          /* index of the cell: number of cells before it */
          const uint64_t below = (m & -m) - 1;
          uint64_t *Uy = U[off[w] + __builtin_popcountll(Y[w] & below)];
          for (int w2 = 0; w2 < W; w2++)
            Uy[w2] |= vb[w2];
        }
    }
    for (uint32_t i = 0; i < ko; i++) {
      const uint32_t u = Lo[i];
      const uint64_t *ub = s->bits + (uint64_t)u * W;
      int ok = 1;
      for (int j = 0; j < ny; j++) {
        uint64_t hit = 0;
        for (int w = 0; w < W; w++)
          hit |= ub[w] & U[j][w];
        ok &= hit != 0;
      }
      Lo[kk] = u;
      kk += (uint32_t)ok;
    }
  }
  s->nvalid[d][o] = kk;
  return kk < ko;
}

/*
 * Cross support on both axes (the col candidates first, which was slightly
 * faster), with forward checking after each; recounts the lists that
 * changed (only their unions if union_only). Returns nonzero if the node is
 * dead.
 */
static int CROSS(sstate_t *s, int d, int union_only) {
  const int n = s->n;
  for (int b = ROW; b <= COL; b++) {
    const int o = 1 - b;
    if (!CROSS_AXIS(s, d, b))
      continue;
    if (s->np[o] + (int)s->nvalid[d][o] < n)
      return 1;
    uint64_t *go = s->cnt[d][o];
    if (union_only)
      UNION_LIST(s, s->valid[d][o], s->nvalid[d][o], go);
    else
      COUNT_LIST(s, s->valid[d][o], s->nvalid[d][o], go);
    /* forward checking: the cells that the o candidates cover */
    const uint64_t *Co = s->cells[d][o], *Cb = s->cells[d][b];
    uint64_t dead = 0;
    for (int w = 0; w < W; w++)
      dead |= Cb[w] & ~Co[w] & ~COVERED(go, w);
    if (dead)
      return 1;
  }
  return 0;
}
#endif /* CARRY */

#ifdef CARRY
/* (a separate function, so that the main call of SUPPORT in TRY_CHILD is
 * inlined as before) */
static __attribute__((noinline)) int SUPPORT_AGAIN(sstate_t *s, int d,
                                                   int a) {
  return SUPPORT(s, d, a, 1);
}
#endif

static void SEARCH_REC(sstate_t *s, int d);

/*
 * Create the child of depth d+1 obtained by placing v (index v, bitset m)
 * on axis b, and search it unless it is pruned. min_v restricts candidates
 * to indices >= min_v. When the lists carry their bitsets, the vectors are
 * identified by their bitsets alone, and v is unused.
 */
static inline TRY_CHILD_INLINE void TRY_CHILD(sstate_t *s, int d, int b,
                                            const uint64_t *m, uint32_t v,
                                            uint32_t min_v) {
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
    /* with two rows and two cols placed, most children are dead, and cross
     * support (with the support filter's test folded in) finds it faster
     * than the support filter's cascade: on the o candidates first, which
     * the support filter would also check first. That is, where the support
     * filter alone lets many of them through: where it kills nearly all of
     * them, cross support is cheaper after it, on the few that it lets
     * through (4-8% less time than before it on P = 16 5 4 2, S = 1200,
     * where the support filter kills 98% of these children, but ~8% more
     * on production-sized instances, where it kills about half). So the
     * search estimates that kill rate as it goes: one child in 16 (all of
     * them in the "after" mode) runs the support filter first, and the
     * mode is "after" when it killed more than 3/4 of the last ~256 of
     * those. Elsewhere than at the (2,2) children, cross support cost more
     * than it saved (see research/ideas.md). */
    const int cross = sup && s->np[ROW] == 2 && s->np[COL] == 2 &&
                      (s->opts.cross > 1 || (s->opts.cross && n <= 6));
    int after = 0; /* cross support after the support filter */
    if (!dead && cross) {
      after = s->cross_after || (++s->xs_tick & 15) == 0;
      if (!after)
        dead = CROSS(s, d + 1, o);
    }
    /* the b candidates are closed already, the o candidates may not be */
    if (!dead && sup)
      dead = SUPPORT(s, d + 1, o, 1);
    if (after) {
      s->xs_k += dead != 0;
      if (++s->xs_n == 256) {
        s->cross_after = s->xs_k > 192;
        s->xs_n = 128; /* (decay) */
        s->xs_k /= 2;
      }
      if (!dead) {
        const uint32_t k0 = s->nvalid[d + 1][ROW] + s->nvalid[d + 1][COL];
        dead = CROSS(s, d + 1, o);
        if (!dead && s->nvalid[d + 1][ROW] + s->nvalid[d + 1][COL] < k0)
          dead = SUPPORT_AGAIN(s, d + 1, o);
      }
    }
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
    /* with two rows and two cols placed, most children are dead, and cross
     * support finds it faster than the support filter, which needs several
     * passes for it (and nearly always where the support filter does not);
     * at the other nodes it cost more than it saved */
#ifdef CROSS_SIMD
    const int cross = s->opts.cross > 1 || (s->opts.cross && n <= 6);
#else
    const int cross = s->opts.cross > 1;
#endif
    if (!dead && sup && cross && s->np[ROW] == 2 && s->np[COL] == 2)
      dead = CROSS(s, d + 1, lazy);
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
#ifdef COUNT_BYTES
  /* with byte counters, the fewest candidates directly: t = count - 1 for
   * the unmatched cells (wrapping, so that a count of 0, which only
   * --no-fc leaves, becomes 255 like the other cells), its minimum over
   * all the cells, then the cells with that minimum. One pass instead of
   * one per class up to the minimum (~7 classes on average at the searched
   * nodes of bench/full.txt, with 4 loads and compares each, and a loop
   * exit that is hard to predict): 2-3% less time on bench/prod.txt, 5% on
   * bench/full.txt. */
  if (s->opts.mrv) {
    const __m512i one = _mm512_set1_epi8(1);
    __m512i t[W], lo = _mm512_set1_epi8(-1);
    for (int w = 0; w < W; w++) {
      t[w] = _mm512_mask_sub_epi8(lo, _cvtu64_mask64(urow[w]),
                                  _mm512_loadu_si512(gr + 8 * w), one);
      t[w] = _mm512_mask_sub_epi8(t[w], _cvtu64_mask64(ucol[w]),
                                  _mm512_loadu_si512(gc + 8 * w), one);
      lo = _mm512_min_epu8(lo, t[w]);
    }
    const int mn = hmin_epu8(lo);
    if (mn == 255) /* no unmatched cell with a candidate */
      return;
    const __m512i M = _mm512_set1_epi8((char)mn);
    if (mn < NCLASS - 1) {
      for (int w = 0; w < W; w++) {
        const uint64_t k = _cvtmask64_u64(_mm512_cmpeq_epi8_mask(t[w], M));
        if (k) {
          x = w * 64 + __builtin_ctzll(k);
          break;
        }
      }
    } else {
      /* >= NCLASS candidates (saturated): the largest label, see below */
      for (int w = W - 1; w >= 0; w--) {
        const uint64_t k = _cvtmask64_u64(_mm512_cmpeq_epi8_mask(t[w], M));
        if (k) {
          x = w * 64 + 63 - __builtin_clzll(k);
          break;
        }
      }
    }
  }
#else
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
#endif
  for (int w = W - 1; w >= 0 && x < 0; w--) {
    uint64_t cls = urow[w] | ucol[w];
#ifndef COUNT_BYTES
    if (s->opts.mrv)
      cls = (urow[w] & IN_CLASS(gr, w, NCLASS - 1)) |
            (ucol[w] & IN_CLASS(gc, w, NCLASS - 1));
#endif
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
    /* (a lane mask only for the last group, see FILTER_STEP) */
    __mmask8 keep =
        cnt - i >= 8
            ? _mm512_test_epi64_mask(_mm512_loadu_si512(in + xw * cap + i), xm)
            : _mm512_mask_test_epi64_mask(
                  (__mmask8)((1u << (cnt - i)) - 1),
                  _mm512_maskz_loadu_epi64((__mmask8)((1u << (cnt - i)) - 1),
                                           in + xw * cap + i),
                  xm);
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
#undef FILTER_STEP
#undef KEEP_CARRY
#undef COUNT_CARRY
#undef TRY_CHILD
#undef COVERED
#undef IN_CLASS
#undef CROSS_AXIS
#undef CROSS_TEST
#undef CROSS
#undef SUPPORT_AGAIN
#undef TRY_CHILD_INLINE
#undef CROSS_SIMD
#undef COUNT_BYTES
#undef CARRY
#undef NCLASS
#undef CAT
#undef CAT_
