/*
 * Arrangement search, label-bitset version.
 *
 * Every vector is a bitset over labels (labels number the distinct numbers
 * in the vectors, roughly the most frequent first: see assign_labels, which
 * also explains how the search is split between the labels). For one
 * (P, S) there are typically only 50-200 labels, so a vector fits in a few
 * machine words, and all of the bookkeeping of the search is done with a
 * handful of bitwise operations:
 *
 *   rowcells / colcells   = union of the placed rows / cols
 *   row-unmatched cells   = rowcells & ~colcells   (need a col through them)
 *   col-unmatched cells   = colcells & ~rowcells   (need a row through them)
 *
 * At every node, we count for every label how many of the remaining
 * candidate rows/cols contain it (saturating byte counters with AVX-512BW,
 * bit-sliced counters otherwise), which tells us for every unmatched cell
 * whether it can be covered by 0, 1, 2 or more candidates.
 * A cell with 0 candidates prunes the node (forward checking), and so does
 * a candidate list that loses all its candidates through some cell under
 * the support filter or cross support (see SUPPORT and CROSS in
 * arrange_core.h); otherwise we branch on the cell with the fewest
 * candidates (ties: smallest label).
 *
 * The axis whose cells are most likely to lose all candidates is filtered
 * first, so that most dead children are discarded after a single pass, and
 * the forward check and the support filter only need the union of each
 * list: the full counts are computed only for the children that survive
 * them (from LAZY_DEPTH on with the N x N intersection matrices, always
 * when the candidate lists carry their bitsets, see CARRY_MAX_W below).
 */
/* clock_gettime and the CPU-time clocks also under a strict -std=c17 */
#define _POSIX_C_SOURCE 200809L
#include "arrange.h"

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#ifdef __AVX512F__
#include <immintrin.h>
#endif

static double wall_time(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + 1e-6 * (double)t.tv_usec;
}

void search_opts_default(search_opts_t *o) {
  o->node_limit = 0;
  o->forward_check = 1;
  o->mrv = 1;
  o->min_words = 0;
  o->gather = 0;
  o->support = 1;
  o->cross = 1;
  o->r1_width = 1;
  o->pretest_min = 5;
  o->r1_stride = 0;
  o->r1_offset = 0;
  o->r1_nstrata = 0;
  for (int h = 0; h < 8; h++)
    o->r1_sstride[h] = o->r1_soffset[h] = 0;
  o->r1_log = NULL;
  o->top_numbers = NULL;
  o->n_top = 0;
  o->top_root_only = 0;
  o->r1_list = NULL;
  o->r1_nlist = 0;
  o->class_support = 0;
}

/*
 * How candidates are counted per label (read with COVERED / IN_CLASS in
 * arrange_core.h):
 *
 * - with AVX-512BW and W <= COUNT_BYTES_MAX_W, and always when the lists
 *   carry their bitsets (W <= CARRY_MAX_W, up to 8 words), one saturating
 *   byte counter per label: adding a vector is one masked add per 64
 *   labels, with the
 *   bitset word loaded straight into a mask register, and the W
 *   accumulators are independent chains of 1-cycle adds (when the lists
 *   carry their bitsets, 8 vectors are added at a time, with a bit
 *   transpose: see COUNT_CARRY). The counts are exact up to 255, so the
 *   most constrained cell is chosen exactly (8% fewer nodes than with 8
 *   slices on bench/quick.txt).
 * - otherwise (the matrix path), NSLICE bit-sliced thermometer counters
 *   (labels in >= 1, 2, ..., NSLICE candidates), vectorized for W = 8 with
 *   AVX-512 (one register per slice; byte counters were 35% slower there),
 *   plain C otherwise.
 *   More slices give a more exact choice (fewer nodes) at a higher cost per
 *   node, which pays off when the counting is vectorized.
 *
 * (Vectorized bit-sliced counters for W = 2 and 4 keep the slices in the
 * lanes of a register, and adding a vector is a lane shift and a ternary
 * op in one dependency chain, ~4 cycles: counting was then a quarter of the
 * search time.)
 */
#ifdef __AVX512BW__
#define COUNT_BYTES_MAX_W 4
#else
#define COUNT_BYTES_MAX_W 0
#endif
#ifndef NSLICE
#ifdef __AVX512F__
#define NSLICE 8
#else
#define NSLICE 4
#endif
#endif
/* with the N x N intersection matrices, nodes from this depth on (two rows
 * and two cols placed, for n = 6) are mostly pruned by the support filter,
 * which needs only the union of each candidate list, so the full counts for
 * choosing the branching cell are computed only for the nodes that are
 * actually searched; above it, counting while filtering is cheaper (2-4
 * are within noise of each other, 5 or never is 15-25% slower, see
 * research/ideas.md). When the lists carry their bitsets, the counts are
 * always computed only for the nodes that are searched (COUNT_CARRY). */
#ifndef LAZY_DEPTH
#define LAZY_DEPTH 4
#endif
/* words of cnt[d][a] per word of labels, for either kind of counters */
#define CNT_WORDS (NSLICE > 8 ? NSLICE : 8)

/*
 * How candidate lists are filtered by the vector v just placed (keeping the
 * vectors u with |u & v| = 0 on v's axis, 1 on the other):
 *
 * - with AVX-512BW and up to CARRY_MAX_W words of labels (<= 512 labels,
 *   nearly all searches), the lists carry the label bitsets of their
 *   entries (vw below) and |u & v| is computed from them, 8 entries at a
 *   time (FILTER_CARRY in arrange_core.h); the support filter tests the
 *   carried words in the same way. No N x N matrices: setup is 4-8x faster
 *   and the search touches only the lists, not a matrix row of N bits per
 *   filter and a random bitset load per candidate (N = 1000-5000);
 * - otherwise, with precomputed N x N intersection bit matrices, looking up
 *   bit u of row v for every entry u (filter_list below).
 *
 * The carried bitsets are fastest with VPOPCNTDQ (|u & v| == 1 in
 * FILTER_CARRY) and VBMI, GFNI and BITALG (the bit transpose of
 * COUNT_CARRY): Ice Lake and later, Zen 4. The other AVX-512BW CPUs
 * (Skylake-X, Cascade Lake) use fallbacks for those, and are still much
 * faster than with the matrices since the support filter (built with
 * -march=cascadelake: 40% less time on bench/full.txt; without the support
 * filter they were 5% slower, see research/ideas.md).
 * Compile with -DCARRY_MAX_W=0 to always use the matrices.
 *
 * The carried path uses exactly the words it needs: W = ceil(L / 64) (at
 * least 2), and each first row r1 is searched with only the words of the
 * labels <= its largest label x, the only ones in its subproblem (see
 * search_root and opts.r1_width). Above 256 labels this replaces the
 * matrices, which cost x2.7 time and x2.8 nodes on the same r1 at N = 23k
 * (8-slice MRV, no cross support); below, each extra word cost x1.38 time.
 */
#ifdef __AVX512BW__
#ifndef CARRY_MAX_W
#define CARRY_MAX_W 8
#endif
#else
#undef CARRY_MAX_W
#define CARRY_MAX_W 0
#endif
#if CARRY_MAX_W > 8
#error "CARRY_MAX_W is at most 8 (see CARRY_WORDS)"
#endif
/* words of the per-depth unions and placed bitsets of the carried path */
#define CARRY_WORDS 8

#define ROW 0
#define COL 1
/* per-depth state for d = 0 .. 2n + 1 (see the allocation in search_vectors:
 * cells up to depth 2n, lists one past it) */
#define MAXD (2 * SQ_MAX_N + 2)

/* state shared by all bitset widths */
typedef struct {
  int n;
  uint32_t N, L;
  uint16_t (*lab)[8];   /* lab[v][p], descending */
  uint64_t *label_val;  /* label -> number */
  uint64_t IW;          /* words per row of the bit matrices below */
  uint32_t nseg;        /* 1024-bit segments covering a row (N bits), or 0
                           to always gather */
  uint64_t *inters0;    /* bit u of row v: u, v disjoint */
  uint64_t *inters1;    /* bit u of row v: |u & v| == 1 */
  uint64_t *has_label;  /* bit v of row x: vector v contains label x */
  uint64_t *bits;       /* bits[v * bw ...]: label bitset of vector v, for
                           v < N, and an empty set for v = N */
  int bw;               /* words per vector in bits: the width W of the
                           search, or on the carried path the widest W of
                           the r1 (the searches at a narrower width read
                           the first W words, which hold all the labels of
                           their vectors) */

  /* candidate lists, as vector indices ... */
  uint32_t *valid[MAXD][2];
  uint32_t nvalid[MAXD][2];
  uint32_t *kids[MAXD];     /* candidates through the branching cell */
  /* ... or, if they carry their bitsets (W <= CARRY_MAX_W), as bitsets
   * only: word w of the bitset of entry i of the list of axis a at depth d
   * is vw[d][a][w * cap + i] (valid, kids and placed are then unused, the
   * vectors being identified by their bitsets), and uni[d][a] is the union
   * of the list, which the forward check and the support filter read
   * instead of the counts */
  uint32_t cap;
  uint64_t *vw[MAXD][2];
  uint64_t *kidw[MAXD];
  uint64_t uni[MAXD][2][CARRY_WORDS];
  uint64_t placedw[2][SQ_MAX_N][CARRY_WORDS];
  uint64_t *cells[MAXD][2]; /* W words each: union of placed rows / cols */
  uint64_t *cnt[MAXD][2];   /* CNT_WORDS * W words each: number of
                               candidates valid[d][a] through each label
                               (see COVERED / IN_CLASS in arrange_core.h) */

  uint32_t placed[2][SQ_MAX_N];
  int np[2];

  search_opts_t opts;
  square_cb cb;
  void *ctx;
  uint64_t nodes, squares;
  /* where cross support runs (see TRY_CHILD in arrange_core.h): after the
   * support filter (1) or before it (0); xs_n of the children where the
   * support filter ran first (a decaying count), xs_k of them killed by
   * it, and a counter to sample one child in 16 */
  int cross_after;
  uint32_t xs_tick, xs_n, xs_k;
  /* the children with at least this many vectors placed run the pretest
   * (opts.pretest_min; INT_MAX when off), which keeps the masks of the
   * entries its two filters keep in ptmask[0] (the other axis) and
   * ptmask[1] (the axis of the new vector), one byte per 8 entries */
  int pretest_min;
  uint8_t *ptmask[2];
  int stop;
  /* the root loop takes r1 < root_limit (top_root_only, see
   * search_opts_t) */
  uint32_t root_limit;
  /* r1 sampling (see search_opts_t and r1_next): ns strata (0: every r1,
   * the default loop), stratum h = root list [lo, hi) taken from lo + off
   * in steps of stride, or the r1 of r1_list (ascending, one stratum of
   * stride 1), and per stratum the sampled count and the sums (and sums of
   * squares) of the per-r1 CPU seconds, squares and nodes. The r1 are
   * global indices into the root list, whatever width each is searched at
   * (see search_root). */
  int r1_ns;
  struct {
    uint32_t lo, hi, stride, off;
    uint64_t n;
    double t, t2, q, q2, nodes;
  } r1s[8];
  uint32_t *r1_list, r1_nlist;
  /* class support (opts.class_support; CLASS_SUP in arrange_core.h): the
   * layers it runs at (bit r + c: the children with r rows and c cols
   * placed; 0 = off), the number of classes and their labels (the top
   * ncls) */
  int cls_np, ncls;
  uint64_t clsmask[CARRY_WORDS];
  /* the words of the class labels (cls_whi: of label L - 1, cls_wlo: of the
   * lowest), the lzcnt of label L - 1 in its word, and the number of class
   * labels in word cls_whi when they span two words (see CLS_AT), and
   * per axis the keep masks of the class support's passes */
  int cls_whi, cls_wlo;
  long long cls_lz0, cls_nhi;
  uint8_t *clskm[2];
} sstate_t;

static inline bool bit_get(const uint64_t *row, uint32_t u) {
  return (row[u >> 6] >> (u & 63)) & 1;
}

#ifdef __AVX512BW__
/* the minimum of the 64 unsigned bytes of v: halving to 16 bytes, then
 * phminposuw on the minima of their pairs (the high byte of each 16-bit
 * lane is then 0) */
static inline int hmin_epu8(__m512i v) {
  __m256i a = _mm256_min_epu8(_mm512_castsi512_si256(v),
                              _mm512_extracti64x4_epi64(v, 1));
  __m128i b = _mm_min_epu8(_mm256_castsi256_si128(a),
                           _mm256_extracti128_si256(a, 1));
  b = _mm_min_epu8(b, _mm_srli_epi16(b, 8));
  return _mm_cvtsi128_si32(_mm_minpos_epu16(b)) & 0xffff;
}
#endif

/*
 * Filtering a candidate list by a row of one of the bit matrices (inters0,
 * inters1, has_label): with AVX-512, 16 entries at a time, looking up the
 * dword holding bit u of the row for every entry u. For N <= 3072 (rows of
 * at most PERM_SEGS 1024-bit segments, which covers nearly all production
 * searches) the row is loaded into registers once and the dwords are picked
 * with one permute per segment: 10-15% faster overall than a gather, whose
 * long latency is exposed after every mispredicted branch. Beyond that the
 * permutes cost more than the gather, which is used instead.
 */
#ifndef PERM_SEGS
#define PERM_SEGS 3 /* at most 3; 0 always gathers (for testing) */
#endif

#ifdef __AVX512F__
static inline __attribute__((always_inline)) uint32_t
filter_list_seg(const uint32_t *in, uint32_t cnt, const uint64_t *row,
                uint32_t min_v, uint32_t *out,
                const int nseg /* 1..3, or 0 for a gather */) {
  __m512i T[6];
  for (int t = 0; t < 2 * nseg; t++)
    T[t] = _mm512_loadu_si512(row + 8 * t);
  const __m512i minv = _mm512_set1_epi32((int)min_v);
  const __m512i low5 = _mm512_set1_epi32(31);
  const __m512i one = _mm512_set1_epi32(1);
  uint32_t k = 0;
  for (uint32_t i = 0; i < cnt; i += 16) {
    __mmask16 live = cnt - i >= 16 ? 0xffff : (__mmask16)((1u << (cnt - i)) - 1);
    __m512i u = _mm512_maskz_loadu_epi32(live, in + i);
    if (min_v)
      live &= _mm512_cmpge_epu32_mask(u, minv);
    /* dword containing bit u, then the bit itself */
    __m512i widx = _mm512_srli_epi32(u, 5), words;
    if (nseg == 0) {
      words = _mm512_mask_i32gather_epi32(_mm512_setzero_si512(), live, widx,
                                          row, 4);
    } else {
      /* a permute picks among 32 dwords (index bits 0-4) */
      words = _mm512_permutex2var_epi32(T[0], widx, T[1]);
      for (int sg = 1; sg < nseg; sg++)
        words = _mm512_mask_mov_epi32(
            words, _mm512_cmpge_epu32_mask(widx, _mm512_set1_epi32(32 * sg)),
            _mm512_permutex2var_epi32(T[2 * sg], widx, T[2 * sg + 1]));
    }
    __m512i bits = _mm512_srlv_epi32(words, _mm512_and_si512(u, low5));
    __mmask16 keep = _mm512_mask_test_epi32_mask(live, bits, one);
    /* compress in a register and store all 16 lanes (out has room): unlike
     * a compressing store, a plain store can be forwarded to the loads of
     * the counting loop that follows */
    _mm512_storeu_si512(out + k, _mm512_maskz_compress_epi32(keep, u));
    k += (uint32_t)__builtin_popcount(keep);
  }
  return k;
}
#endif

/* out = { u in in[0..cnt) : u >= min_v and bit u of row is set };
 * returns |out|. nseg = number of 1024-bit segments of the rows. out must
 * have room for 16 entries past the result (the lists have 64). */
static inline uint32_t filter_list(const uint32_t *in, uint32_t cnt,
                                   const uint64_t *row, uint32_t min_v,
                                   uint32_t nseg, uint32_t *out) {
#ifdef __AVX512F__
  if (nseg == 1 && PERM_SEGS >= 1)
    return filter_list_seg(in, cnt, row, min_v, out, 1);
  if (nseg == 2 && PERM_SEGS >= 2)
    return filter_list_seg(in, cnt, row, min_v, out, 2);
  if (nseg == 3 && PERM_SEGS >= 3)
    return filter_list_seg(in, cnt, row, min_v, out, 3);
  return filter_list_seg(in, cnt, row, min_v, out, 0);
#else
  (void)nseg;
  uint32_t k = 0;
  for (uint32_t i = 0; i < cnt; i++) {
    uint32_t u = in[i];
    out[k] = u;
    k += (u >= min_v) & bit_get(row, u);
  }
  return k;
#endif
}

#if defined(__AVX512F__) && NSLICE == 8
static inline void count_slices_simd_w8(const uint64_t *bits,
                                        const uint32_t *list, uint32_t k,
                                        uint64_t *g) {
  /* 512-bit label sets: one register per slice */
  __m512i A[NSLICE];
  for (int t = 0; t < NSLICE; t++)
    A[t] = _mm512_setzero_si512();
  for (uint32_t i = 0; i < k; i++) {
    __m512i m = _mm512_loadu_si512(bits + 8 * (uint64_t)list[i]);
    for (int t = NSLICE - 1; t > 0; t--)
      A[t] = _mm512_ternarylogic_epi64(A[t], A[t - 1], m, 0xF8);
    A[0] = _mm512_or_si512(A[0], m);
  }
  for (int t = 0; t < NSLICE; t++)
    _mm512_storeu_si512(g + 8 * t, A[t]);
}

/* the same, over the vectors of list[0..k) disjoint from bad, which are moved
 * to the front of list (in place); returns their number. A vector meeting
 * bad is counted as the empty set. */
static inline uint32_t keep_count_simd_w8(const uint64_t *bits, uint32_t *list,
                                          uint32_t k, const uint64_t *bad,
                                          uint64_t *g) {
  const __m512i b = _mm512_loadu_si512(bad);
  __m512i A[NSLICE];
  for (int t = 0; t < NSLICE; t++)
    A[t] = _mm512_setzero_si512();
  uint32_t kk = 0;
  for (uint32_t i = 0; i < k; i++) {
    uint32_t u = list[i];
    __m512i x = _mm512_loadu_si512(bits + 8 * (uint64_t)u);
    int keep = _mm512_test_epi64_mask(x, b) == 0;
    __m512i m = _mm512_maskz_mov_epi64((__mmask8)-keep, x);
    for (int t = NSLICE - 1; t > 0; t--)
      A[t] = _mm512_ternarylogic_epi64(A[t], A[t - 1], m, 0xF8);
    A[0] = _mm512_or_si512(A[0], m);
    list[kk] = u;
    kk += (uint32_t)keep;
  }
  for (int t = 0; t < NSLICE; t++)
    _mm512_storeu_si512(g + 8 * t, A[t]);
  return kk;
}

/* only slice 0 (the union), of all of list[0..k) or of the vectors disjoint
 * from bad (in place, as above) */
static inline void union_simd_w8(const uint64_t *bits, const uint32_t *list,
                                 uint32_t k, uint64_t *g) {
  __m512i acc = _mm512_setzero_si512();
  for (uint32_t i = 0; i < k; i++)
    acc = _mm512_or_si512(acc,
                          _mm512_loadu_si512(bits + 8 * (uint64_t)list[i]));
  _mm512_storeu_si512(g, acc);
}

static inline uint32_t keep_union_simd_w8(const uint64_t *bits, uint32_t *list,
                                          uint32_t k, const uint64_t *bad,
                                          uint64_t *g) {
  const __m512i b = _mm512_loadu_si512(bad);
  __m512i acc = _mm512_setzero_si512();
  uint32_t kk = 0;
  for (uint32_t i = 0; i < k; i++) {
    uint32_t u = list[i];
    __m512i x = _mm512_loadu_si512(bits + 8 * (uint64_t)u);
    int keep = _mm512_test_epi64_mask(x, b) == 0;
    acc = _mm512_mask_or_epi64(acc, (__mmask8)-keep, acc, x);
    list[kk] = u;
    kk += (uint32_t)keep;
  }
  _mm512_storeu_si512(g, acc);
  return kk;
}
#define HAVE_SIMD_COUNT_W8 1
#endif

/* report a square, given the labels (descending) of its rows and cols */
static void report_labels(sstate_t *s, const uint16_t *lab[2][SQ_MAX_N]) {
  s->squares++;
  if (!s->cb)
    return;
  square_t sq;
  sq.n = s->n;
  for (int i = 0; i < s->n; i++)
    for (int j = 0; j < s->n; j++) {
      sq.rows[i][j] = s->label_val[lab[ROW][i][j]];
      sq.cols[i][j] = s->label_val[lab[COL][i][j]];
    }
  if (s->cb(&sq, s->ctx))
    s->stop = 1;
}

/* the square of the placed vectors (by index) */
static void report(sstate_t *s) {
  const uint16_t *lab[2][SQ_MAX_N];
  for (int a = 0; a < 2; a++)
    for (int i = 0; i < s->n; i++)
      lab[a][i] = s->lab[s->placed[a][i]];
  report_labels(s, lab);
}

#if CARRY_MAX_W > 0
/* the k-th lowest set bit of y (k < popcount(y)) */
static inline uint64_t nth_bit(uint64_t y, int k) {
#ifdef __BMI2__
  return _pdep_u64((uint64_t)1 << k, y);
#else
  for (; k > 0; k--)
    y &= y - 1;
  return y & -y;
#endif
}

/* y rotated left by r (0 <= r < 64) */
static inline uint64_t rotl64(uint64_t y, int r) {
  return r ? y << r | y >> (64 - r) : y;
}

/* lane j of the result: the or of the 8 lanes of a[j] */
static inline __m512i or_lanes8(const __m512i a[8]) {
  /* t[k], 128-bit lane l: the or of lanes 2l, 2l + 1 of a[2k], a[2k + 1] */
  __m512i t[4], u[2];
  for (int k = 0; k < 4; k++)
    t[k] = _mm512_or_si512(_mm512_unpacklo_epi64(a[2 * k], a[2 * k + 1]),
                           _mm512_unpackhi_epi64(a[2 * k], a[2 * k + 1]));
  /* u[k], 128-bit lane l: the or of 128-bit lanes 2l, 2l + 1 of t[2k]
   * (l < 2) or of t[2k + 1] (l >= 2) */
  for (int k = 0; k < 2; k++)
    u[k] = _mm512_or_si512(_mm512_shuffle_i64x2(t[2 * k], t[2 * k + 1], 0x88),
                           _mm512_shuffle_i64x2(t[2 * k], t[2 * k + 1], 0xDD));
  return _mm512_or_si512(_mm512_shuffle_i64x2(u[0], u[1], 0x88),
                         _mm512_shuffle_i64x2(u[0], u[1], 0xDD));
}

/* the square of the placed vectors (by bitset, of W words) */
static void report_bits(sstate_t *s, int W) {
  uint16_t buf[2][SQ_MAX_N][8];
  const uint16_t *lab[2][SQ_MAX_N];
  for (int a = 0; a < 2; a++)
    for (int i = 0; i < s->n; i++) {
      int p = 0;
      for (int w = W - 1; w >= 0; w--)
        for (uint64_t x = s->placedw[a][i][w]; x && p < 8;
             x &= ~((uint64_t)1 << (63 - __builtin_clzll(x))))
          buf[a][i][p++] = (uint16_t)(64 * w + 63 - __builtin_clzll(x));
      lab[a][i] = buf[a][i];
    }
  report_labels(s, lab);
}
#endif

#ifdef PROFILE
/* per depth: children created, children pruned before being searched, and
 * for the searched ones, children generated and candidates on each axis */
uint64_t prof_created[MAXD], prof_pruned[MAXD], prof_kids[MAXD],
    prof_vr[MAXD], prof_vc[MAXD];
void prof_print(void) {
  for (int d = 0; d < MAXD; d++) {
    uint64_t searched = prof_created[d] - prof_pruned[d];
    if (prof_created[d])
      fprintf(stderr,
              "depth %2d created %10lu pruned %5.1f%% | searched: kids %5.2f "
              "candidate rows %7.1f cols %7.1f\n",
              d, prof_created[d], 100.0 * prof_pruned[d] / prof_created[d],
              searched ? (double)prof_kids[d] / searched : 0,
              searched ? (double)prof_vr[d] / searched : 0,
              searched ? (double)prof_vc[d] / searched : 0);
  }
  memset(prof_created, 0, sizeof(prof_created));
  memset(prof_pruned, 0, sizeof(prof_pruned));
  memset(prof_kids, 0, sizeof(prof_kids));
  memset(prof_vr, 0, sizeof(prof_vr));
  memset(prof_vc, 0, sizeof(prof_vc));
}
#endif

#ifdef CHILD_PROF
/*
 * Profile of the creation of children on the carried path (build with
 * -DCHILD_PROF; printed on exit): by the layer of the child (rows, cols
 * placed), TSC cycles per phase of TRY_CHILD (and for choosing the kids in
 * SEARCH_REC), how the children die, and the list sizes. The rdtsc pairs
 * cost ~50% more time, which they add to the short phases most: for the
 * share of a phase, timing builds that stop after it was more reliable
 * (see research/ideas.md, "The cost of creating children").
 */
#include <x86intrin.h>
enum { CP_FO, CP_FB, CP_CROSS, CP_SUP, CP_AFTER, CP_COUNT, CP_SEL, CP_PT,
       CP_NPH };
enum { CD_CNT_O, CD_FC_O, CD_CNT_B, CD_FC_B, CD_CROSS, CD_SUP, CD_AFTER,
       CD_LIVE, CD_N };
#define CP_LAYERS (9 * 9)
uint64_t cp_cyc[CP_LAYERS][CP_NPH], cp_calls[CP_LAYERS][CP_NPH],
    cp_die[CP_LAYERS][CD_N], cp_in[CP_LAYERS][2], cp_out[CP_LAYERS][2],
    cp_pt_kill[CP_LAYERS], cp_total;
static void __attribute__((destructor)) cp_print(void) {
  static const char *ph[CP_NPH] = {"filter o", "filter b", "cross",
                                   "support",  "after",    "count",
                                   "kids",     "pretest"};
  static const char *dn[CD_N] = {"count o", "fc o",  "count b", "fc b",
                                 "cross",   "support", "after", "live"};
  fprintf(stderr, "search: %.3fG cycles\n", cp_total * 1e-9);
  for (int l = 0; l < CP_LAYERS; l++) {
    uint64_t made = 0;
    for (int k = 0; k < CD_N; k++)
      made += cp_die[l][k];
    if (!made)
      continue;
    fprintf(stderr,
            "(%d,%d): %llu children, lists in %.1f / %.1f, after the "
            "filters %.1f / %.1f (o / b)\n",
            l / 9, l % 9, (unsigned long long)made,
            (double)cp_in[l][0] / made, (double)cp_in[l][1] / made,
            cp_calls[l][CP_FO] ? (double)cp_out[l][0] / cp_calls[l][CP_FO] : 0,
            cp_calls[l][CP_FB] ? (double)cp_out[l][1] / cp_calls[l][CP_FB] : 0);
    for (int p = 0; p < CP_NPH; p++)
      if (cp_calls[l][p])
        fprintf(stderr, "  %-8s %11llu calls %7.1f cycles %5.1f%%\n", ph[p],
                (unsigned long long)cp_calls[l][p],
                (double)cp_cyc[l][p] / cp_calls[l][p],
                100.0 * cp_cyc[l][p] / cp_total);
    fprintf(stderr, "  died:");
    for (int k = 0; k < CD_N; k++)
      fprintf(stderr, " %s %.1f%%", dn[k], 100.0 * cp_die[l][k] / made);
    if (cp_pt_kill[l])
      fprintf(stderr, " (pretest %.1f%%)", 100.0 * cp_pt_kill[l] / made);
    fprintf(stderr, "\n");
  }
}
#define CP_START() uint64_t cp_t = __rdtsc()
#define CP_RESET() (cp_t = __rdtsc())
#define CP_PHASE(l, p)                                                         \
  do {                                                                         \
    const uint64_t cp_u = __rdtsc();                                           \
    cp_cyc[l][p] += cp_u - cp_t;                                               \
    cp_calls[l][p]++;                                                          \
    cp_t = cp_u;                                                               \
  } while (0)
#else
#define CP_START()
#define CP_RESET()
#define CP_PHASE(l, p)
#endif

#ifdef CLASS_PROF
/* the cost of the class support (build with -DMAGIC_DEFS=CLASS_PROF;
 * printed on exit): TSC cycles in CLASS_SUP and in the whole search, its
 * calls, the dead ones, its passes, and the entries in and dropped (from
 * the nodes that live) */
#include <x86intrin.h>
enum { CLP_ALL, CLP_SEARCH, CLP_N };
uint64_t clp_cyc[CLP_N], clp_calls, clp_dead, clp_passes, clp_removed, clp_in,
    clp_npass_hist[2][16];
static void __attribute__((destructor)) clp_print(void) {
  if (!clp_calls)
    return;
  const double t = (double)clp_cyc[CLP_SEARCH];
  fprintf(stderr,
          "class support: %lu calls, %.1f%% dead, %.2f passes, in %.1f "
          "removed %.1f entries per call; search %.3fG cycles, class %.2f%%, "
          "%.0f cycles per call\n",
          (unsigned long)clp_calls, 100.0 * clp_dead / clp_calls,
          (double)clp_passes / clp_calls, (double)clp_in / clp_calls,
          (double)clp_removed / (clp_calls - clp_dead + 1e-9), t * 1e-9,
          100.0 * clp_cyc[CLP_ALL] / t, (double)clp_cyc[CLP_ALL] / clp_calls);
  fprintf(stderr, "  passes of the live / dead calls:");
  for (int k = 0; k < 2; k++) {
    fprintf(stderr, k ? " | dead" : " live");
    for (int i = 0; i < 16; i++)
      fprintf(stderr, " %lu", (unsigned long)clp_npass_hist[k][i]);
  }
  fprintf(stderr, "\n");
}
#endif

/*
 * r1 sampling (search_opts_t r1_*), for timing sums too large to search
 * whole: the root loop runs only the sampled r1, by their global index in
 * the root list [0, root_limit), whatever width each is searched at (see
 * search_root): stratum h is [lo, hi) from lo + off in steps of stride
 * (one stratum [0, root_limit) with r1_stride / r1_offset), or the r1 of
 * r1_list. stride times the sampled totals of a stratum is unbiased for it.
 * The r1 log gets one line per sampled r1: r1, stratum, stride, squares,
 * nodes, thread CPU seconds, its largest label, the width it was searched
 * at; and a "# run" line per width run.
 */
static double thread_cpu(void) {
  struct timespec t;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

/* the first sampled r1 >= r (UINT32_MAX if none), and its stratum *hp */
static uint32_t r1_next(const sstate_t *s, uint32_t r, int *hp) {
  *hp = 0;
  if (s->r1_list) {
    uint32_t lo = 0, hi = s->r1_nlist;
    while (lo < hi) {
      const uint32_t mid = (lo + hi) / 2;
      if (s->r1_list[mid] < r)
        lo = mid + 1;
      else
        hi = mid;
    }
    return lo < s->r1_nlist ? s->r1_list[lo] : UINT32_MAX;
  }
  for (int h = 0; h < s->r1_ns; h++) {
    const uint64_t first = (uint64_t)s->r1s[h].lo + s->r1s[h].off,
                   k = s->r1s[h].stride, hi = s->r1s[h].hi;
    if (r >= hi)
      continue;
    const uint64_t x = r <= first ? first : first + (r - first + k - 1) / k * k;
    if (x < hi) {
      *hp = h;
      return (uint32_t)x;
    }
  }
  return UINT32_MAX;
}

static void r1_account(sstate_t *s, int h, uint32_t r1, uint64_t squares,
                       uint64_t nodes, double sec, int w) {
  s->r1s[h].n++;
  s->r1s[h].t += sec;
  s->r1s[h].t2 += sec * sec;
  s->r1s[h].q += (double)squares;
  s->r1s[h].q2 += (double)squares * (double)squares;
  s->r1s[h].nodes += (double)nodes;
  if (s->opts.r1_log)
    fprintf(s->opts.r1_log, "%u %d %u %lu %lu %.6f %u %d\n", r1, h,
            s->r1s[h].stride, (unsigned long)squares, (unsigned long)nodes,
            sec, (unsigned)s->lab[r1][0], w);
}

/* instantiate the recursive search for each bitset width (the carried path
 * also for 5-7 words) */
#define W 2
#include "arrange_core.h"
#undef W
#define W 3
#include "arrange_core.h"
#undef W
#define W 4
#include "arrange_core.h"
#undef W
#if CARRY_MAX_W >= 5
#define W 5
#include "arrange_core.h"
#undef W
#endif
#if CARRY_MAX_W >= 6
#define W 6
#include "arrange_core.h"
#undef W
#endif
#if CARRY_MAX_W >= 7
#define W 7
#include "arrange_core.h"
#undef W
#endif
#define W 8
#include "arrange_core.h"
#undef W
#define W 16
#include "arrange_core.h"
#undef W
#define W 64
#include "arrange_core.h"
#undef W

/* the root loop over r1 in [i0, i1) at width w (see SEARCH_ROOT) */
static void search_root_w(sstate_t *s, int w, uint32_t i0, uint32_t i1) {
  switch (w) {
  case 2:
    search_root_2(s, i0, i1);
    break;
  case 3:
    search_root_3(s, i0, i1);
    break;
  case 4:
    search_root_4(s, i0, i1);
    break;
#if CARRY_MAX_W >= 5
  case 5:
    search_root_5(s, i0, i1);
    break;
#endif
#if CARRY_MAX_W >= 6
  case 6:
    search_root_6(s, i0, i1);
    break;
#endif
#if CARRY_MAX_W >= 7
  case 7:
    search_root_7(s, i0, i1);
    break;
#endif
  case 8:
    search_root_8(s, i0, i1);
    break;
  case 16:
    search_root_16(s, i0, i1);
    break;
  default:
    search_root_64(s, i0, i1);
    break;
  }
}

/* the width of the subproblem of r1 = v: the words up to its largest label
 * (at least 2, and at least opts.min_words) */
static inline int r1_words(const sstate_t *s, uint32_t v) {
  const int w = s->lab[v][0] / 64 + 1, wmin = s->opts.min_words;
  return w >= wmin && w >= 2 ? w : wmin > 2 ? wmin : 2;
}

/*
 * The search, r1 by r1, at width W_ (the width of all the labels). On the
 * carried path with opts.r1_width, the r1 are searched in runs of equal
 * width instead: the subproblem of r1 has only the vectors of index > r1
 * (see SEARCH_ROOT), whose labels are all <= x, the largest label of r1
 * (the vectors are sorted by their labels, descending), so it needs only
 * the words up to x, and x falls with r1. Each run fills the lists of depth
 * 0 (the vectors from its first r1 on) with that many words of their
 * bitsets. The adaptive cross support's state (cross_after, xs_*) carries
 * over from one run to the next, so the search makes the same choices,
 * and visits the same nodes, as at a single width.
 *
 * The roots are the r1 < root_limit (top_root_only: those through the top
 * label, a prefix of the sorted list), partitioned into the width runs;
 * with r1 sampling, each run searches only the sampled r1 in it (r1_next,
 * by global index), so every sampled r1 is searched exactly once, at the
 * width of its run.
 */
static void search_root(sstate_t *s, int W_, int carry) {
  s->nodes = 1;
  /* the roots: every r1, or with top_root_only those through the top
   * label (a prefix of the list) */
  const uint32_t lim = s->root_limit < s->N ? s->root_limit : s->N;
  if (!carry || !s->opts.r1_width) {
    search_root_w(s, W_, 0, lim);
    return;
  }
  for (uint32_t i0 = 0, i1; i0 < lim && !s->stop; i0 = i1) {
    const int w = r1_words(s, i0);
    for (i1 = i0 + 1; i1 < lim && r1_words(s, i1) == w; i1++)
      ;
    if (s->r1_ns && s->opts.r1_log)
      fprintf(s->opts.r1_log, "# run W %d r1 %u .. %u\n", w, i0, i1);
#ifdef ARRANGE_DEBUG
    /* poison the words that a narrower run must not read (they hold the
     * lists of the wider runs before it) */
    if (w < W_)
      for (int d = 0; d <= 2 * s->n + 1; d++) {
        const size_t off = (size_t)w * s->cap,
                     len = (size_t)(W_ - w) * s->cap * sizeof(uint64_t);
        memset(s->vw[d][0] + off, 0xff, len);
        memset(s->vw[d][1] + off, 0xff, len);
        memset(s->kidw[d] + off, 0xff, len);
      }
#endif
    search_root_w(s, w, i0, i1);
  }
}

/* ---------------------------------------------------------------------- */
/* setup                                                                   */

/*
 * The labels (the order of the numbers). The search finds every square from
 * its largest label x: r1 and c1 are the row and col through x, and all the
 * other vectors have only labels < x (see SEARCH_ROOT). So it is one
 * subproblem per number x, over the vectors whose numbers all have labels
 * <= x, and the order decides how the squares (and the dead ends) are split
 * between the subproblems. A subproblem is cheap when x is in few of its
 * vectors: few choices of r1, few cols through x for each (MRV then nearly
 * always branches on x first), so the rarest numbers get the largest labels.
 *
 * Rarest where it matters, though: in the subproblem of x, the vectors
 * through larger labels are gone. So the labels are assigned from the
 * largest down, each to the number in the fewest of the vectors not yet
 * removed, and then its vectors are removed (a degeneracy order). The
 * rarest numbers overall still come first, but a number whose vectors
 * mostly go with rarer numbers moves down. Ties: rarer overall, then the
 * larger number, get the larger label (as in the plain frequency order this
 * replaces). 2-8% fewer nodes than the plain frequency order on every
 * instance of bench/quick.txt and bench/full.txt (6% and 3% in total, 2% on
 * the largest instance). Other orders, also relative to plain frequency on
 * bench/quick.txt: by the number of distinct numbers sharing a vector with
 * x, +6% nodes (as a degeneracy order +3%); by the pairs of vectors meeting
 * exactly at x, +0.3% (as a degeneracy order -5%); frequency ascending,
 * +43%; random, +29%. A local search over orders, running the search for
 * each, found only 3-4% fewer nodes than this one.
 *
 * vals[0..L) are the distinct numbers (ascending) and rk[v * n + p] is the
 * rank in vals of number p of vector v; sets lab_of[rank]. O(L^2 + n count).
 */
static void assign_labels(int n, size_t count, size_t L, const uint32_t *rk,
                          uint32_t *lab_of, const uint32_t *top, int ntop) {
  /* the vectors through each number (CSR) */
  uint32_t *first = calloc(L + 1, sizeof(uint32_t));
  for (size_t i = 0; i < count * n; i++)
    first[rk[i] + 1]++;
  for (size_t x = 0; x < L; x++)
    first[x + 1] += first[x];
  uint32_t *occ = malloc(count * n * sizeof(uint32_t));
  uint32_t *live = malloc(L * sizeof(uint32_t)); /* fill pointers for now */
  memcpy(live, first, L * sizeof(uint32_t));
  for (size_t v = 0; v < count; v++)
    for (int p = 0; p < n; p++)
      occ[live[rk[v * n + p]]++] = (uint32_t)v;
  /* freq: vectors through the number; live: those not yet removed */
  uint32_t *freq = malloc(L * sizeof(uint32_t));
  for (size_t x = 0; x < L; x++)
    freq[x] = live[x] = first[x + 1] - first[x];
  uint8_t *removed = calloc(count, 1), *done = calloc(L, 1);
  /* the given top numbers (ranks) first, without removing their vectors */
  size_t lab = L;
  for (int t = 0; t < ntop; t++)
    if (!done[top[t]]) {
      done[top[t]] = 1;
      lab_of[top[t]] = (uint32_t)--lab;
    }
  while (lab-- > 0) {
    size_t b = L; /* the rarest number left; ranks ascend with the numbers */
    for (size_t x = 0; x < L; x++)
      if (!done[x] && (b == L || live[x] < live[b] ||
                       (live[x] == live[b] && freq[x] <= freq[b])))
        b = x;
    done[b] = 1;
    lab_of[b] = (uint32_t)lab;
    for (uint32_t k = first[b]; k < first[b + 1]; k++) {
      uint32_t v = occ[k];
      if (removed[v])
        continue;
      removed[v] = 1;
      for (int p = 0; p < n; p++)
        live[rk[v * n + p]]--;
    }
  }
  free(first);
  free(occ);
  free(live);
  free(freq);
  free(removed);
  free(done);
}

typedef struct {
  uint64_t val;
  uint32_t id; /* dense id of val (see dense_ids) */
} vid_t;

static int vid_cmp(const void *a, const void *b) {
  const vid_t *x = a, *y = b;
  return (x->val > y->val) - (x->val < y->val);
}

static int lab_cmp_n;
static int lab_cmp(const void *a, const void *b) {
  const uint16_t *x = a, *y = b;
  for (int p = 0; p < lab_cmp_n; p++)
    if (x[p] != y[p])
      return (x[p] < y[p]) - (x[p] > y[p]); /* descending */
  return 0;
}

int search_carries(uint32_t labels, const search_opts_t *opts) {
  /* as the choice of W_ in search_vectors */
  const size_t Lw = labels > 64 * (size_t)opts->min_words
                        ? labels
                        : 64 * (size_t)opts->min_words;
  int W_ = Lw <= 128    ? 2
           : Lw <= 192  ? 3
           : Lw <= 256  ? 4
           : Lw <= 512  ? 8
           : Lw <= 1024 ? 16
                        : 64;
  if ((Lw + 63) / 64 <= CARRY_MAX_W && W_ > 4)
    W_ = (int)((Lw + 63) / 64);
  return W_ <= CARRY_MAX_W;
}

search_stats_t search_vectors(const vec_list_t *l, size_t start, size_t count,
                              const search_opts_t *opts, square_cb cb,
                              void *ctx) {
  search_stats_t st;
  memset(&st, 0, sizeof(st));
  int n = l->n;
  if (count < (size_t)(2 * n) || n > SQ_MAX_N)
    return st;
  double t0 = wall_time();

  sstate_t s;
  memset(&s, 0, sizeof(s));
  s.n = n;
  s.N = (uint32_t)count;
  s.opts = *opts;
  s.cb = cb;
  s.ctx = ctx;
  /* (the pretest repeats the forward check, so it needs it to be exact) */
  s.pretest_min = opts->pretest_min > 0 && opts->forward_check
                      ? opts->pretest_min
                      : INT_MAX;
  /* every r1, or with top_root_only those through the top label (set after
   * the sort below) */
  s.root_limit = (uint32_t)count;
  /* distinct values (dense ids, see dense_ids: hashing rather than sorting
   * all count * n elements, which was half of the setup), and the rank of
   * every element among them in ascending order (sorting only the L
   * distinct values) */
  size_t ne = count * n;
  uint64_t *vals = malloc(ne * sizeof(uint64_t));
  uint32_t *ids = malloc(ne * sizeof(uint32_t));
  size_t L = dense_ids(vec_list_get(l, start), ne, ids, vals);
  if (L > 64 * 64) {
    /* would need wider bitsets; not expected for realistic inputs */
    free(vals);
    free(ids);
    st.num_labels = (int)L;
    st.truncated = 1;
    return st;
  }
  vid_t *vi = malloc(L * sizeof(vid_t));
  for (size_t i = 0; i < L; i++) {
    vi[i].val = vals[i];
    vi[i].id = (uint32_t)i;
  }
  qsort(vi, L, sizeof(vid_t), vid_cmp);
  uint32_t *rank_of_id = malloc(L * sizeof(uint32_t));
  for (size_t i = 0; i < L; i++)
    rank_of_id[vi[i].id] = (uint32_t)i;
  uint32_t *rk = malloc(ne * sizeof(uint32_t));
  for (size_t i = 0; i < ne; i++)
    rk[i] = rank_of_id[ids[i]];
  uint32_t *lab_of = malloc(L * sizeof(uint32_t)); /* by rank */
  /* numbers forced onto the top labels (top_numbers: the ranks of those
   * that occur, in the given order) */
  uint32_t top_rank[SQ_MAX_N * 8];
  int ntop = 0;
  for (int t = 0; t < opts->n_top && ntop < SQ_MAX_N * 8; t++)
    for (size_t i = 0; i < L; i++)
      if (vi[i].val == opts->top_numbers[t])
        top_rank[ntop++] = (uint32_t)i;
  assign_labels(n, count, L, rk, lab_of, top_rank, ntop);
  s.L = (uint32_t)L;
  s.label_val = malloc((L + 1) * sizeof(uint64_t));
  for (size_t i = 0; i < L; i++)
    s.label_val[lab_of[i]] = vi[i].val;
  s.label_val[L] = 0;

  s.lab = malloc(count * sizeof(*s.lab));
  for (size_t v = 0; v < count; v++) {
    for (int p = 0; p < 8; p++)
      s.lab[v][p] = (uint16_t)L;
    for (int p = 0; p < n; p++) {
      uint16_t x = (uint16_t)lab_of[rk[v * n + p]];
      int b = p; /* insertion sort, descending */
      for (; b > 0 && s.lab[v][b - 1] < x; b--)
        s.lab[v][b] = s.lab[v][b - 1];
      s.lab[v][b] = x;
    }
  }
  free(rk);
  free(lab_of);
  free(rank_of_id);
  free(vi);
  free(ids);
  free(vals);
  /* sort by descending labels: the vector with the largest label first */
  int kb = 1; /* bits per label in a packed key */
  while (((size_t)1 << kb) < L)
    kb++;
  if (n * kb <= 64) {
    /* the labels of a vector packed into one integer, largest first, sort
     * the same way as the labels and determine the vector: radix sort the
     * keys, then unpack them (a qsort of the label arrays with a
     * comparator was most of the setup time left after dense_ids) */
    uint64_t *key = malloc(count * sizeof(uint64_t));
    uint64_t *tmp = malloc(count * sizeof(uint64_t));
    for (size_t v = 0; v < count; v++) {
      uint64_t k = 0;
      for (int p = 0; p < n; p++)
        k = k << kb | s.lab[v][p];
      key[v] = k;
    }
    for (int shift = 0; shift < n * kb; shift += 8) {
      size_t pos[256] = {0};
      for (size_t v = 0; v < count; v++)
        pos[(key[v] >> shift) & 255]++;
      for (size_t c = 0, t = 0; c < 256; c++) {
        size_t m = pos[c];
        pos[c] = t;
        t += m;
      }
      for (size_t v = 0; v < count; v++)
        tmp[pos[(key[v] >> shift) & 255]++] = key[v];
      uint64_t *sw = key;
      key = tmp;
      tmp = sw;
    }
    const uint64_t kmask = ((uint64_t)1 << kb) - 1;
    for (size_t v = 0; v < count; v++) {
      uint64_t k = key[count - 1 - v]; /* descending */
      for (int p = n - 1; p >= 0; p--, k >>= kb)
        s.lab[v][p] = (uint16_t)(k & kmask);
    }
    free(key);
    free(tmp);
  } else {
    lab_cmp_n = n;
    qsort(s.lab, count, sizeof(*s.lab), lab_cmp);
  }

  if (opts->top_root_only && ntop > 0 && count > 0 && s.lab[0][0] == L - 1) {
    /* the vectors through the top label come first */
    uint32_t k = 0;
    while (k < count && s.lab[k][0] == L - 1)
      k++;
    s.root_limit = k;
  }
  /* r1 sampling plan, over the roots [0, root_limit) */
  if (opts->r1_nlist > 0 && opts->r1_list) {
    /* the listed r1 (a sorted copy, without repeats), as one stratum of
     * stride 1 */
    s.r1_ns = 1;
    s.r1s[0].lo = 0;
    s.r1s[0].hi = s.root_limit;
    s.r1s[0].stride = 1;
    s.r1_list = malloc(opts->r1_nlist * sizeof(uint32_t));
    for (uint32_t i = 0; i < opts->r1_nlist; i++) {
      uint32_t j = i;
      for (; j > 0 && s.r1_list[j - 1] > opts->r1_list[i]; j--)
        s.r1_list[j] = s.r1_list[j - 1];
      s.r1_list[j] = opts->r1_list[i];
    }
    for (uint32_t i = 0; i < opts->r1_nlist; i++)
      if (!s.r1_nlist || s.r1_list[i] != s.r1_list[s.r1_nlist - 1])
        s.r1_list[s.r1_nlist++] = s.r1_list[i];
  } else if (opts->r1_nstrata > 0 || opts->r1_stride >= 1 || opts->r1_log) {
    /* (r1_stride 1, which msearch's --calib-r1-stride 1 sets, is a plan
     * with every r1: exact "estimates"; the default 0 is no plan) */
    const int K = opts->r1_nstrata > 0
                      ? (opts->r1_nstrata < 8 ? opts->r1_nstrata : 8)
                      : 1;
    const uint64_t nr = s.root_limit;
    s.r1_ns = K;
    for (int h = 0; h < K; h++) {
      s.r1s[h].lo = (uint32_t)(nr * h / K);
      s.r1s[h].hi = (uint32_t)(nr * (h + 1) / K);
      uint32_t k = opts->r1_nstrata > 0 ? opts->r1_sstride[h] : opts->r1_stride;
      uint32_t o = opts->r1_nstrata > 0 ? opts->r1_soffset[h] : opts->r1_offset;
      s.r1s[h].stride = k < 1 ? 1 : k;
      s.r1s[h].off = o % s.r1s[h].stride;
    }
  }

  /* bitset width: on the carried path the words needed, otherwise the next
   * of the instantiated widths of the matrix path */
  size_t Lw = L > 64 * (size_t)opts->min_words ? L : 64 * (size_t)opts->min_words;
  int W_ = Lw <= 128    ? 2
           : Lw <= 192  ? 3
           : Lw <= 256  ? 4
           : Lw <= 512  ? 8
           : Lw <= 1024 ? 16
                        : 64;
  if ((Lw + 63) / 64 <= CARRY_MAX_W && W_ > 4)
    W_ = (int)((Lw + 63) / 64);
  s.bw = W_;
  s.bits = calloc((count + 1) * W_ + 8, sizeof(uint64_t));
  for (size_t v = 0; v < count; v++)
    for (int p = 0; p < n; p++)
      s.bits[v * W_ + s.lab[v][p] / 64] |= (uint64_t)1 << (s.lab[v][p] % 64);

  /* intersection bit matrices, unless the lists carry their bitsets */
  const int carry = W_ <= CARRY_MAX_W;
  s.IW = (count + 63) / 64 + 1;
  s.nseg = opts->gather ? 0 : (uint32_t)((count + 1023) / 1024);
  if (!carry) {
    /* + 16: filter_list may read whole segments past the end of a row */
    s.inters0 = calloc(count * s.IW + 16, sizeof(uint64_t));
    s.inters1 = calloc(count * s.IW + 16, sizeof(uint64_t));
    uint32_t *deg = calloc(L + 2, sizeof(uint32_t));
    for (size_t v = 0; v < count; v++)
      for (int p = 0; p < n; p++)
        deg[s.lab[v][p] + 1]++;
    for (size_t x = 0; x < L; x++)
      deg[x + 1] += deg[x];
    uint32_t *occ = malloc(ne * sizeof(uint32_t));
    uint32_t *fill = malloc((L + 1) * sizeof(uint32_t));
    memcpy(fill, deg, (L + 1) * sizeof(uint32_t));
    for (size_t v = 0; v < count; v++)
      for (int p = 0; p < n; p++)
        occ[fill[s.lab[v][p]]++] = (uint32_t)v;
    uint8_t *inter = calloc(count, 1);
    uint32_t *touched = malloc(count * sizeof(uint32_t));
    for (size_t v = 0; v < count; v++) {
      size_t nt = 0;
      for (int p = 0; p < n; p++) {
        uint32_t x = s.lab[v][p];
        for (uint32_t k = deg[x]; k < deg[x + 1]; k++) {
          uint32_t w = occ[k];
          if (inter[w]++ == 0)
            touched[nt++] = w;
        }
      }
      uint64_t *r0 = s.inters0 + v * s.IW, *r1 = s.inters1 + v * s.IW;
      memset(r0, 0xff, (count / 64) * sizeof(uint64_t));
      if (count % 64)
        r0[count / 64] = ((uint64_t)1 << (count % 64)) - 1;
      for (size_t t = 0; t < nt; t++) {
        uint32_t w = touched[t];
        r0[w >> 6] &= ~((uint64_t)1 << (w & 63));
        if (inter[w] == 1)
          r1[w >> 6] |= (uint64_t)1 << (w & 63);
        inter[w] = 0;
      }
    }
    free(deg);
    free(occ);
    free(fill);
    free(inter);
    free(touched);
  }
  /* the same layout for label -> vectors, so that the children of a node
   * (the candidates through the branching cell) are selected with the same
   * vectorized filter as the candidate lists, instead of a scalar loop with
   * an unpredictable branch per candidate */
  if (!carry) {
    s.has_label = calloc((L + 1) * s.IW + 16, sizeof(uint64_t));
    for (size_t v = 0; v < count; v++)
      for (int p = 0; p < n; p++)
        s.has_label[s.lab[v][p] * s.IW + v / 64] |= (uint64_t)1 << (v % 64);
  }

  int maxd = 2 * n + 1;
  _Static_assert(2 * SQ_MAX_N + 1 < MAXD, "MAXD too small for n = SQ_MAX_N");
  /* room for 64 entries past the end of a list (filters store whole
   * vectors), and 64-byte aligned word arrays */
  s.cap = (uint32_t)((count + 64 + 7) & ~(size_t)7);
  const size_t wbytes = (size_t)W_ * s.cap * sizeof(uint64_t);
  for (int a = 0; a < 2; a++)
    for (int d = 0; d <= maxd; d++) {
      if (carry)
        s.vw[d][a] = aligned_alloc(64, wbytes);
      else
        s.valid[d][a] = malloc((count + 64) * sizeof(uint32_t));
      s.cells[d][a] = calloc(W_, sizeof(uint64_t));
      s.cnt[d][a] = calloc(CNT_WORDS * W_, sizeof(uint64_t));
    }
  for (int d = 0; d <= maxd; d++) {
    if (carry)
      s.kidw[d] = aligned_alloc(64, wbytes);
    else
      s.kids[d] = malloc((count + 64) * sizeof(uint32_t));
  }
  if (s.r1_ns && s.opts.r1_log)
    fprintf(s.opts.r1_log,
            "# N %u L %u W %d carry %d r1_width %d roots %u strata %d list %u\n",
            s.N, s.L, W_, carry, s.opts.r1_width, s.root_limit, s.r1_ns,
            s.r1_nlist);
  if (carry)
    for (int a = 0; a < 2; a++)
      s.ptmask[a] = malloc(s.cap / 8 + 1);
  /* class support: the top numbers are the classes when every vector has
   * exactly one of them, i.e. its largest label is a top one and its
   * second largest is not (checked here; V_d guarantees it). With fewer
   * than n of them occurring there is no square, and the classes that do
   * occur are as good as any. Carried path only, with the support filter
   * (which it extends). */
  if (carry && opts->class_support > 0 && opts->top_numbers &&
      opts->n_top == n && ntop > 0 && opts->support && opts->forward_check) {
    const size_t c0 = L - (size_t)ntop; /* the lowest class label */
    int ok = 1;
    for (size_t v = 0; v < count && ok; v++)
      ok = s.lab[v][0] >= c0 && s.lab[v][1] < c0;
    if (ok) {
      s.ncls = ntop;
      s.cls_np = 1 << 2; /* the (1,1) children */
      for (size_t x = c0; x < L; x++)
        s.clsmask[x >> 6] |= (uint64_t)1 << (x & 63);
      s.cls_whi = (int)((L - 1) >> 6);
      s.cls_wlo = (int)(c0 >> 6);
      s.cls_lz0 = 63 - (long long)((L - 1) & 63);
      s.cls_nhi = (long long)((L - 1) & 63) + 1;
      for (int a = 0; a < 2; a++)
        s.clskm[a] = malloc(s.cap / 8 + 2);
    }
  }
  double t1 = wall_time();

#ifdef CLASS_PROF
  const uint64_t clp_t0 = __rdtsc();
#endif
  search_root(&s, W_, carry);
#ifdef CLASS_PROF
  clp_cyc[CLP_SEARCH] += __rdtsc() - clp_t0;
#endif

  free(s.r1_list);
  st.nodes = s.nodes;
  st.squares = s.squares;
  st.truncated = s.stop == 2;
  st.class_used = s.ncls > 0;
  st.num_labels = (int)L;
  st.setup_seconds = t1 - t0;
  if (s.r1_ns) {
    double vt = 0, vq = 0;
    for (int h = 0; h < s.r1_ns; h++) {
      const double k = s.r1s[h].stride, m = (double)s.r1s[h].n;
      const double Nh = (double)(s.r1s[h].hi - s.r1s[h].lo);
      st.r1_sampled += s.r1s[h].n;
      st.est_seconds += k * s.r1s[h].t;
      st.est_squares += k * s.r1s[h].q;
      st.est_nodes += k * s.r1s[h].nodes;
      if (s.r1_list) {
        /* (no standard error for a list of r1: the estimates are then the
         * sampled totals) */
      } else if (m >= 2) {
        /* N_h^2 (1 - m / N_h) s_h^2 / m */
        const double f = Nh * Nh * (1 - m / Nh) / m / (m - 1);
        vt += f * (s.r1s[h].t2 - s.r1s[h].t * s.r1s[h].t / m);
        vq += f * (s.r1s[h].q2 - s.r1s[h].q * s.r1s[h].q / m);
      } else if (m < Nh) {
        /* fewer than 2 of the stratum's r1 sampled (and not all of them):
         * its variance cannot be estimated, so the standard errors leave it
         * out (they are too small), and with none sampled the estimates
         * have nothing of it either; counted, so callers can flag it */
        st.r1_strata_nose++;
      }
    }
    st.est_nodes += 1; /* the root */
    st.se_seconds = vt > 0 ? sqrt(vt) : 0;
    st.se_squares = vq > 0 ? sqrt(vq) : 0;
  }
  st.seconds = wall_time() - t1;

  for (int a = 0; a < 2; a++)
    for (int d = 0; d <= maxd; d++) {
      free(s.valid[d][a]);
      free(s.cells[d][a]);
      free(s.cnt[d][a]);
      free(s.vw[d][a]);
    }
  for (int d = 0; d <= maxd; d++) {
    free(s.kids[d]);
    free(s.kidw[d]);
  }
  for (int a = 0; a < 2; a++) {
    free(s.ptmask[a]);
    free(s.clskm[a]);
  }
  free(s.inters0);
  free(s.inters1);
  free(s.has_label);
  free(s.bits);
  free(s.lab);
  free(s.label_val);
  return st;
}
