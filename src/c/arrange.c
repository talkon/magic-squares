/*
 * Arrangement search, label-bitset version.
 *
 * Every vector is a bitset over labels (labels are the distinct numbers in
 * the vectors, most frequent first). For one (P, S) there are typically only
 * 50-200 labels, so a vector fits in a few machine words, and all of the
 * bookkeeping of the search is done with a handful of bitwise operations:
 *
 *   rowcells / colcells   = union of the placed rows / cols
 *   row-unmatched cells   = rowcells & ~colcells   (need a col through them)
 *   col-unmatched cells   = colcells & ~rowcells   (need a row through them)
 *
 * At every node, we count for every label how many of the remaining
 * candidate rows/cols contain it (saturating byte counters with AVX-512BW,
 * bit-sliced counters otherwise), which tells us for every unmatched cell
 * whether it can be covered by 0, 1, 2 or more candidates.
 * A cell with 0 candidates prunes the node (forward checking); otherwise we
 * branch on the cell with the fewest candidates (ties: largest label).
 *
 * The counting is fused into the filtering of the candidate lists, and the
 * axis whose cells are most likely to lose all candidates is filtered first,
 * so that most dead children are discarded after a single pass.
 */
#include "arrange.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
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
}

/*
 * How candidates are counted per label (read with COVERED / IN_CLASS in
 * arrange_core.h):
 *
 * - with AVX-512BW and W <= COUNT_BYTES_MAX_W, one saturating byte counter
 *   per label: adding a vector is one masked add per 64 labels, with the
 *   bitset word loaded straight into a mask register, and the W
 *   accumulators are independent chains of 1-cycle adds. The counts are
 *   exact up to 255, so the most constrained cell is chosen exactly (8%
 *   fewer nodes than with 8 slices on bench/quick.txt).
 * - otherwise, NSLICE bit-sliced thermometer counters (labels in >= 1, 2,
 *   ..., NSLICE candidates), vectorized for W = 8 with AVX-512 (one register
 *   per slice; byte counters were 35% slower there), plain C otherwise.
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
/* words of cnt[d][a] per word of labels, for either kind of counters */
#define CNT_WORDS (NSLICE > 8 ? NSLICE : 8)

#define ROW 0
#define COL 1
#define MAXD (2 * SQ_MAX_N + 1)

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
  uint64_t *bits;       /* bits[v * W ...]: label bitset of vector v, for
                           v < N, and an empty set for v = N */

  uint32_t *valid[MAXD][2];
  uint32_t nvalid[MAXD][2];
  uint32_t *kids[MAXD];     /* candidates through the branching cell */
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
  int stop;
} sstate_t;

static inline bool bit_get(const uint64_t *row, uint32_t u) {
  return (row[u >> 6] >> (u & 63)) & 1;
}

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
#define HAVE_SIMD_COUNT_W8 1
#endif

static void report(sstate_t *s) {
  s->squares++;
  if (!s->cb)
    return;
  square_t sq;
  sq.n = s->n;
  for (int i = 0; i < s->n; i++)
    for (int j = 0; j < s->n; j++) {
      sq.rows[i][j] = s->label_val[s->lab[s->placed[ROW][i]][j]];
      sq.cols[i][j] = s->label_val[s->lab[s->placed[COL][i]][j]];
    }
  if (s->cb(&sq, s->ctx))
    s->stop = 1;
}

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

/* instantiate the recursive search for each bitset width */
#define W 2
#include "arrange_core.h"
#undef W
#define W 3
#include "arrange_core.h"
#undef W
#define W 4
#include "arrange_core.h"
#undef W
#define W 8
#include "arrange_core.h"
#undef W
#define W 16
#include "arrange_core.h"
#undef W
#define W 64
#include "arrange_core.h"
#undef W

/* ---------------------------------------------------------------------- */
/* setup                                                                   */

typedef struct {
  uint64_t val;
  uint32_t freq;
} vf_t;

static int vf_cmp(const void *a, const void *b) {
  const vf_t *x = a, *y = b;
  if (x->freq != y->freq)
    return (x->freq < y->freq) - (x->freq > y->freq); /* most frequent first */
  return (x->val > y->val) - (x->val < y->val);
}

static int u64_cmp(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}

static int lab_cmp_n;
static int lab_cmp(const void *a, const void *b) {
  const uint16_t *x = a, *y = b;
  for (int p = 0; p < lab_cmp_n; p++)
    if (x[p] != y[p])
      return (x[p] < y[p]) - (x[p] > y[p]); /* descending */
  return 0;
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

  /* distinct values, by decreasing frequency */
  size_t ne = count * n;
  uint64_t *vals = malloc(ne * sizeof(uint64_t));
  memcpy(vals, vec_list_get(l, start), ne * sizeof(uint64_t));
  qsort(vals, ne, sizeof(uint64_t), u64_cmp);
  vf_t *vf = malloc(ne * sizeof(vf_t));
  size_t L = 0;
  for (size_t i = 0; i < ne; i++) {
    if (i == 0 || vals[i] != vals[i - 1]) {
      vf[L].val = vals[i];
      vf[L].freq = 0;
      L++;
    }
    vf[L - 1].freq++;
  }
  /* vals[0..L) = distinct values ascending, vf[0..L) same order for now */
  for (size_t i = 0; i < L; i++)
    vals[i] = vf[i].val;
  qsort(vf, L, sizeof(vf_t), vf_cmp);
  s.L = (uint32_t)L;
  s.label_val = malloc((L + 1) * sizeof(uint64_t));
  uint32_t *val_to_lab = malloc(L * sizeof(uint32_t)); /* by rank in vals */
  for (size_t i = 0; i < L; i++) {
    s.label_val[i] = vf[i].val;
    size_t lo = 0, hi = L;
    while (lo < hi) {
      size_t mid = (lo + hi) / 2;
      if (vals[mid] < vf[i].val)
        lo = mid + 1;
      else
        hi = mid;
    }
    val_to_lab[lo] = (uint32_t)i;
  }
  s.label_val[L] = 0;

  s.lab = malloc(count * sizeof(*s.lab));
  for (size_t v = 0; v < count; v++) {
    const uint64_t *vec = vec_list_get(l, start + v);
    for (int p = 0; p < 8; p++)
      s.lab[v][p] = (uint16_t)L;
    for (int p = 0; p < n; p++) {
      size_t lo = 0, hi = L;
      while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (vals[mid] < vec[p])
          lo = mid + 1;
        else
          hi = mid;
      }
      uint16_t x = (uint16_t)val_to_lab[lo];
      int b = p; /* insertion sort, descending */
      for (; b > 0 && s.lab[v][b - 1] < x; b--)
        s.lab[v][b] = s.lab[v][b - 1];
      s.lab[v][b] = x;
    }
  }
  lab_cmp_n = n;
  /* sort by descending labels: the vector with the largest label first */
  qsort(s.lab, count, sizeof(*s.lab), lab_cmp);

  /* bitset width */
  size_t Lw = L > 64 * (size_t)opts->min_words ? L : 64 * (size_t)opts->min_words;
  int W_ = Lw <= 128    ? 2
           : Lw <= 192  ? 3
           : Lw <= 256  ? 4
           : Lw <= 512  ? 8
           : Lw <= 1024 ? 16
                        : 64;
  if (L > 64 * 64) {
    /* would need wider bitsets; not expected for realistic inputs */
    free(vals);
    free(vf);
    free(val_to_lab);
    free(s.lab);
    free(s.label_val);
    st.num_labels = (int)L;
    st.truncated = 1;
    return st;
  }
  s.bits = calloc((count + 1) * W_ + 8, sizeof(uint64_t));
  for (size_t v = 0; v < count; v++)
    for (int p = 0; p < n; p++)
      s.bits[v * W_ + s.lab[v][p] / 64] |= (uint64_t)1 << (s.lab[v][p] % 64);

  /* intersection bit matrices */
  s.IW = (count + 63) / 64 + 1;
  s.nseg = opts->gather ? 0 : (uint32_t)((count + 1023) / 1024);
  /* + 16: filter_list may read whole segments past the end of a row */
  s.inters0 = calloc(count * s.IW + 16, sizeof(uint64_t));
  s.inters1 = calloc(count * s.IW + 16, sizeof(uint64_t));
  {
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
  s.has_label = calloc((L + 1) * s.IW + 16, sizeof(uint64_t));
  for (size_t v = 0; v < count; v++)
    for (int p = 0; p < n; p++)
      s.has_label[s.lab[v][p] * s.IW + v / 64] |= (uint64_t)1 << (v % 64);

  int maxd = 2 * n + 1;
  for (int a = 0; a < 2; a++)
    for (int d = 0; d <= maxd; d++) {
      s.valid[d][a] = malloc((count + 64) * sizeof(uint32_t));
      s.cells[d][a] = calloc(W_, sizeof(uint64_t));
      s.cnt[d][a] = calloc(CNT_WORDS * W_, sizeof(uint64_t));
    }
  for (int d = 0; d <= maxd; d++)
    s.kids[d] = malloc((count + 64) * sizeof(uint32_t));
  double t1 = wall_time();

  switch (W_) {
  case 2:
    search_root_2(&s);
    break;
  case 3:
    search_root_3(&s);
    break;
  case 4:
    search_root_4(&s);
    break;
  case 8:
    search_root_8(&s);
    break;
  case 16:
    search_root_16(&s);
    break;
  default:
    search_root_64(&s);
    break;
  }

  st.nodes = s.nodes;
  st.squares = s.squares;
  st.truncated = s.stop == 2;
  st.num_labels = (int)L;
  st.setup_seconds = t1 - t0;
  st.seconds = wall_time() - t1;

  for (int a = 0; a < 2; a++)
    for (int d = 0; d <= maxd; d++) {
      free(s.valid[d][a]);
      free(s.cells[d][a]);
      free(s.cnt[d][a]);
    }
  for (int d = 0; d <= maxd; d++)
    free(s.kids[d]);
  free(s.inters0);
  free(s.inters1);
  free(s.has_label);
  free(s.bits);
  free(s.lab);
  free(s.label_val);
  free(val_to_lab);
  free(vals);
  free(vf);
  return st;
}
