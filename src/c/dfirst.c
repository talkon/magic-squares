/*
 * Diagonal-first search: see dfirst.h.
 */
/* clock_gettime and the CPU-time clocks also under a strict -std=c17 */
#define _POSIX_C_SOURCE 200809L
#include "dfirst.h"

#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "square.h"

static double wall_time(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + 1e-6 * (double)t.tv_usec;
}

static double thread_cpu(void) {
  struct timespec t;
  clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t);
  return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

/* the involutions of gen_inv: for even n the perfect matchings, (n - 1)!!,
 * for odd n those with one fixed point, n (n - 2)!! (105 for n = 7 and 8) */
#define DF_DFACT(k)                                                            \
  ((k) <= 1 ? 1 : (k) == 3 ? 3 : (k) == 5 ? 15 : (k) == 7 ? 105 : (k) == 9 ? 945 \
                                                                 : 10395)
#define DF_NINV(n) ((n) % 2 ? (n) * DF_DFACT((n) - 2) : DF_DFACT((n) - 1))
#define DF_MAX_INV 128
_Static_assert(SQ_MAX_N <= 12, "DF_DFACT covers n <= 12");
_Static_assert(DF_NINV(SQ_MAX_N) <= DF_MAX_INV &&
                   DF_NINV(SQ_MAX_N - 1) <= DF_MAX_INV,
               "dfirst_t inv[] too small for the involutions of n = SQ_MAX_N");

struct dfirst_s {
  int n;
  const vec_list_t *vecs;
  size_t vstart, vcount;
  const vec_list_t *diags;
  size_t dstart, dcount;
  /* number -> dense id over the numbers of vecs (open addressing; hval is
   * id + 1, 0 = empty) */
  int hbits;
  uint64_t *hkey;
  uint32_t *hval;
  uint32_t L;
  /* postings: the vectors through number id, ascending, are
   * post[first[id] .. first[id + 1]) */
  uint32_t *first, *post;
  /* per-d scratch: hit counters, the touched vectors, and a bitmap of the
   * vectors hit exactly once */
  uint8_t *hit;
  uint32_t *touched;
  uint64_t *mark;
  /* the diagonal set: sorted (descending) copies of its n-sets, and an
   * open-addressing table of their indices + 1 */
  uint64_t *dsorted;
  int sbits;
  uint32_t *sidx;
  /* the involutions that relate the two diagonals (no fixed point for even
   * n, one for odd n), as in square.c */
  int ninv;
  int inv[DF_MAX_INV][SQ_MAX_N];
  int top_root; /* see dfirst_set_top_root */
  /* per-d scratch of vd_labels: a stamp per number id */
  uint32_t *lstamp, lcur;
};

void dfirst_set_top_root(dfirst_t *df, int on) { df->top_root = on; }

static inline uint32_t hash_u64(uint64_t x, int bits) {
  return (uint32_t)((x * 0x9E3779B97F4A7C15ull) >> (64 - bits));
}

static uint64_t hash_tuple(const uint64_t *t, int n) {
  uint64_t h = 0x243F6A8885A308D3ull;
  for (int i = 0; i < n; i++) {
    h ^= t[i] + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
    h *= 0xBF58476D1CE4E5B9ull;
  }
  return h ^ (h >> 31);
}

static void sort_desc(uint64_t *t, int n) {
  for (int i = 1; i < n; i++) {
    uint64_t x = t[i];
    int j = i;
    for (; j > 0 && t[j - 1] < x; j--)
      t[j] = t[j - 1];
    t[j] = x;
  }
}

static int64_t num_id(const dfirst_t *df, uint64_t x) {
  const uint32_t mask = (1u << df->hbits) - 1;
  for (uint32_t h = hash_u64(x, df->hbits);; h = (h + 1) & mask) {
    if (!df->hval[h])
      return -1;
    if (df->hkey[h] == x)
      return (int64_t)df->hval[h] - 1;
  }
}

static void gen_inv(dfirst_t *df, int *cur, int fixed_left) {
  const int n = df->n;
  int a = 0;
  while (a < n && cur[a] >= 0)
    a++;
  if (a == n) {
    if (df->ninv < DF_MAX_INV) /* (always, see the static assert) */
      memcpy(df->inv[df->ninv++], cur, SQ_MAX_N * sizeof(int));
    return;
  }
  if (fixed_left) {
    cur[a] = a;
    gen_inv(df, cur, 0);
    cur[a] = -1;
  }
  for (int b = a + 1; b < n; b++) {
    if (cur[b] >= 0)
      continue;
    cur[a] = b;
    cur[b] = a;
    gen_inv(df, cur, fixed_left);
    cur[a] = cur[b] = -1;
  }
}

dfirst_t *dfirst_new(const vec_list_t *vecs, size_t vstart, size_t vcount,
                     const vec_list_t *diags, size_t dstart, size_t dcount) {
  dfirst_t *df = calloc(1, sizeof(*df));
  const int n = vecs->n;
  df->n = n;
  df->vecs = vecs;
  df->vstart = vstart;
  df->vcount = vcount;
  df->diags = diags;
  df->dstart = dstart;
  df->dcount = dcount;

  /* number ids */
  const size_t ne = vcount * (size_t)n;
  df->hbits = 4;
  while (((size_t)1 << df->hbits) < 2 * ne + 16)
    df->hbits++;
  df->hkey = malloc(((size_t)1 << df->hbits) * sizeof(uint64_t));
  df->hval = calloc((size_t)1 << df->hbits, sizeof(uint32_t));
  uint32_t *ids = malloc((ne + 1) * sizeof(uint32_t));
  const uint32_t mask = (1u << df->hbits) - 1;
  for (size_t e = 0; e < ne; e++) {
    const uint64_t x = vecs->elts[vstart * n + e];
    uint32_t h = hash_u64(x, df->hbits);
    while (df->hval[h] && df->hkey[h] != x)
      h = (h + 1) & mask;
    if (!df->hval[h]) {
      df->hkey[h] = x;
      df->hval[h] = ++df->L;
    }
    ids[e] = df->hval[h] - 1;
  }
  /* postings (CSR, vectors ascending) */
  df->first = calloc(df->L + 1, sizeof(uint32_t));
  for (size_t e = 0; e < ne; e++)
    df->first[ids[e] + 1]++;
  for (uint32_t x = 0; x < df->L; x++)
    df->first[x + 1] += df->first[x];
  df->post = malloc((ne + 1) * sizeof(uint32_t));
  uint32_t *fill = malloc((df->L + 1) * sizeof(uint32_t));
  memcpy(fill, df->first, (df->L + 1) * sizeof(uint32_t));
  for (size_t v = 0; v < vcount; v++)
    for (int p = 0; p < n; p++)
      df->post[fill[ids[v * n + p]]++] = (uint32_t)v;
  free(fill);
  free(ids);
  df->hit = calloc(vcount + 1, 1);
  df->touched = malloc((ne + 1) * sizeof(uint32_t));
  df->mark = calloc(vcount / 64 + 2, sizeof(uint64_t));

  /* the diagonal set */
  df->dsorted = malloc((dcount * (size_t)n + 1) * sizeof(uint64_t));
  df->sbits = 4;
  while (((size_t)1 << df->sbits) < 2 * dcount + 16)
    df->sbits++;
  df->sidx = calloc((size_t)1 << df->sbits, sizeof(uint32_t));
  const uint32_t smask = (1u << df->sbits) - 1;
  for (size_t i = 0; i < dcount; i++) {
    uint64_t *t = df->dsorted + i * n;
    memcpy(t, vec_list_get(diags, dstart + i), n * sizeof(uint64_t));
    sort_desc(t, n);
    uint32_t h = (uint32_t)(hash_tuple(t, n) >> (64 - df->sbits));
    int dup = 0;
    while (df->sidx[h]) {
      if (!memcmp(df->dsorted + (size_t)(df->sidx[h] - 1) * n, t,
                  n * sizeof(uint64_t))) {
        dup = 1; /* a repeated diagonal: the set keeps the first */
        break;
      }
      h = (h + 1) & smask;
    }
    if (!dup)
      df->sidx[h] = (uint32_t)i + 1;
  }

  int cur[SQ_MAX_N];
  for (int i = 0; i < SQ_MAX_N; i++)
    cur[i] = -1;
  gen_inv(df, cur, n % 2);
  df->top_root = 1;
  return df;
}

void dfirst_free(dfirst_t *df) {
  if (!df)
    return;
  free(df->hkey);
  free(df->hval);
  free(df->first);
  free(df->post);
  free(df->hit);
  free(df->touched);
  free(df->mark);
  free(df->dsorted);
  free(df->sidx);
  free(df->lstamp);
  free(df);
}

uint32_t dfirst_num_labels(const dfirst_t *df) { return df->L; }

int dfirst_in_set(const dfirst_t *df, const uint64_t *t0) {
  const int n = df->n;
  uint64_t t[SQ_MAX_N];
  memcpy(t, t0, n * sizeof(uint64_t));
  sort_desc(t, n);
  const uint32_t smask = (1u << df->sbits) - 1;
  for (uint32_t h = (uint32_t)(hash_tuple(t, n) >> (64 - df->sbits));;
       h = (h + 1) & smask) {
    if (!df->sidx[h])
      return 0;
    if (!memcmp(df->dsorted + (size_t)(df->sidx[h] - 1) * n, t,
                n * sizeof(uint64_t)))
      return 1;
  }
}

size_t dfirst_vd(dfirst_t *df, size_t i, vec_list_t *out) {
  const int n = df->n;
  const uint64_t *d = vec_list_get(df->diags, df->dstart + i);
  size_t nt = 0;
  for (int p = 0; p < n; p++) {
    const int64_t x = num_id(df, d[p]);
    if (x < 0)
      continue;
    for (uint32_t k = df->first[x]; k < df->first[x + 1]; k++) {
      const uint32_t v = df->post[k];
      if (df->hit[v]++ == 0)
        df->touched[nt++] = v;
    }
  }
  /* the vectors hit once, in list order: a bitmap over the touched words */
  size_t wlo = (size_t)-1, whi = 0;
  for (size_t t = 0; t < nt; t++) {
    const uint32_t v = df->touched[t];
    if (df->hit[v] == 1) {
      df->mark[v >> 6] |= (uint64_t)1 << (v & 63);
      if ((v >> 6) < wlo)
        wlo = v >> 6;
      if ((v >> 6) > whi)
        whi = v >> 6;
    }
    df->hit[v] = 0;
  }
  size_t k = 0;
  for (size_t w = wlo; w <= whi && wlo != (size_t)-1; w++) {
    uint64_t m = df->mark[w];
    df->mark[w] = 0;
    while (m) {
      const size_t v = w * 64 + (size_t)__builtin_ctzll(m);
      m &= m - 1;
      vec_list_push(out, vec_list_get(df->vecs, df->vstart + v));
      k++;
    }
  }
  return k;
}

/* the top-label root pays with the carried bitsets and with the matrices
 * of up to 4 words (0.64x and 0.68x the CPU of the plain root at 159
 * labels), not with the matrices of 8 words (257-512 labels without
 * AVX-512BW, or with CARRY_MAX_W < 8: 1.27x, as the top labels then cost
 * more nodes there than on the carried path; research/ideas.md) */
static int top_root_pays(uint32_t labels, const search_opts_t *opts) {
  const uint32_t lw = labels > 64u * (uint32_t)opts->min_words
                          ? labels
                          : 64u * (uint32_t)opts->min_words;
  return lw <= 256 || search_carries(labels, opts);
}

/* the distinct numbers of V_d (only computed when the sum's labels do not
 * decide top_root_pays) */
static uint32_t vd_labels(dfirst_t *df, const vec_list_t *sub) {
  if (!df->lstamp)
    df->lstamp = calloc(df->L + 1, sizeof(uint32_t));
  if (++df->lcur == 0) {
    memset(df->lstamp, 0, (df->L + 1) * sizeof(uint32_t));
    df->lcur = 1;
  }
  uint32_t k = 0;
  for (size_t v = 0; v < sub->count; v++) {
    const uint64_t *e = vec_list_get(sub, v);
    for (int p = 0; p < df->n; p++) {
      const int64_t id = num_id(df, e[p]);
      if (id >= 0 && df->lstamp[id] != df->lcur) {
        df->lstamp[id] = df->lcur;
        k++;
      }
    }
  }
  return k;
}

/* the search's callback on V_d: one (square, d) pair */
typedef struct {
  dfirst_t *df;
  size_t d_index;
  const uint64_t *d;
  dsquare_cb cb;
  void *ctx;
  uint64_t pairs, partners;
  int stopped;
} inner_t;

static int next_perm(int *p, int n) {
  int i = n - 2;
  while (i >= 0 && p[i] >= p[i + 1])
    i--;
  if (i < 0)
    return 0;
  int j = n - 1;
  while (p[j] <= p[i])
    j--;
  int t = p[i];
  p[i] = p[j];
  p[j] = t;
  for (int a = i + 1, b = n - 1; a < b; a++, b--) {
    t = p[a];
    p[a] = p[b];
    p[b] = t;
  }
  return 1;
}

static int inner_cb(const square_t *sq, void *vctx) {
  inner_t *c = vctx;
  const dfirst_t *df = c->df;
  const int n = df->n;
  grid_t g;
  square_to_grid(sq, &g);
  /* sigma: the col of d's number in each row */
  int sigma[SQ_MAX_N];
  for (int i = 0; i < n; i++) {
    sigma[i] = -1;
    for (int j = 0; j < n; j++)
      for (int p = 0; p < n; p++)
        if (g.g[i][j] == c->d[p])
          sigma[i] = j;
  }
  dsquare_t ds;
  ds.sq = sq;
  ds.d = c->d;
  ds.d_index = c->d_index;
  /* traversals in the set */
  int perm[SQ_MAX_N];
  for (int i = 0; i < n; i++)
    perm[i] = i;
  ds.set_count = 0;
  do {
    uint64_t t[SQ_MAX_N];
    for (int i = 0; i < n; i++)
      t[i] = g.g[i][perm[i]];
    ds.set_count += dfirst_in_set(df, t);
  } while (next_perm(perm, n));
  /* a partner on the other diagonal (sigma is complete: every row of a
   * square of V_d meets d once) */
  ds.partner = 0;
  int complete = 1;
  for (int i = 0; i < n; i++)
    complete &= sigma[i] >= 0;
  for (int k = 0; k < df->ninv && !ds.partner && complete; k++) {
    uint64_t t[SQ_MAX_N];
    for (int i = 0; i < n; i++)
      t[i] = g.g[i][sigma[df->inv[k][i]]];
    ds.partner = dfirst_in_set(df, t);
  }
  c->pairs++;
  c->partners += ds.partner;
  if (c->cb && c->cb(&ds, c->ctx)) {
    c->stopped = 1;
    return 1;
  }
  return 0;
}

dfirst_stats_t dfirst_search(dfirst_t *df, size_t lo, size_t hi, size_t stride,
                             size_t off, const search_opts_t *opts,
                             dsquare_cb cb, void *ctx, FILE *dlog) {
  dfirst_stats_t st;
  memset(&st, 0, sizeof(st));
  const int n = df->n;
  if (stride < 1)
    stride = 1;
  if (hi > df->dcount)
    hi = df->dcount;
  vec_list_t sub;
  vec_list_init(&sub, n);
  const double c0 = thread_cpu();
  for (size_t i = lo + off % stride; i < hi; i += stride) {
    const double ci = thread_cpu(), w0 = wall_time();
    sub.count = 0;
    dfirst_vd(df, i, &sub);
    const double w1 = wall_time();
    st.nd++;
    st.vd_total += sub.count;
    st.vd_seconds += w1 - w0;
    inner_t ic = {df, i, vec_list_get(df->diags, df->dstart + i), cb, ctx,
                  0, 0, 0};
    search_stats_t ss;
    memset(&ss, 0, sizeof(ss));
    search_opts_t o = *opts;
    uint64_t top[SQ_MAX_N];
    /* the top-label root, where it pays (see dfirst_set_top_root): every
     * V_d of the sum, or else this V_d by its own labels */
    if (df->top_root && (top_root_pays(df->L, opts) ||
                         top_root_pays(vd_labels(df, &sub), opts))) {
      /* d's numbers by the number of vectors of V_d through them */
      const uint64_t *d = ic.d;
      size_t cnt[SQ_MAX_N] = {0};
      for (size_t v = 0; v < sub.count; v++) {
        const uint64_t *e = vec_list_get(&sub, v);
        for (int p = 0; p < n; p++)
          for (int q = 0; q < n; q++)
            cnt[q] += e[p] == d[q];
      }
      int ord[SQ_MAX_N];
      for (int q = 0; q < n; q++) {
        int b = q;
        for (; b > 0 && cnt[ord[b - 1]] > cnt[q]; b--)
          ord[b] = ord[b - 1];
        ord[b] = q;
      }
      for (int q = 0; q < n; q++)
        top[q] = d[ord[q]];
      o.top_numbers = top;
      o.n_top = n;
      o.top_root_only = 1;
    }
    if (sub.count >= (size_t)(2 * n))
      ss = search_vectors(&sub, 0, sub.count, &o, inner_cb, &ic);
    st.nodes += ss.nodes;
    st.pairs += ic.pairs;
    st.partners += ic.partners;
    st.setup_seconds += ss.setup_seconds;
    st.search_seconds += ss.seconds;
    st.truncated |= ss.truncated;
    const double cd = thread_cpu() - ci;
    st.d_cpu += cd;
    st.d_cpu2 += cd * cd;
    st.d_pairs2 += (double)ic.pairs * (double)ic.pairs;
    if (dlog)
      fprintf(dlog, "%zu %zu %d %lu %lu %.6f %.6f %.6f %.6f\n", i, sub.count,
              ss.num_labels, (unsigned long)ss.nodes, (unsigned long)ic.pairs,
              w1 - w0, ss.setup_seconds, ss.seconds, cd);
    if (ic.stopped) {
      st.stopped = 1;
      break;
    }
  }
  st.cpu_seconds = thread_cpu() - c0;
  vec_list_free(&sub);
  return st;
}
