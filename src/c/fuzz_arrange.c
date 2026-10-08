/*
 * Differential test of the arrangement search against an independent
 * brute-force oracle.
 *
 * usage: fuzz_arrange [-v] [--variants a,b,..] [--mode M] [--n N]
 *                     [--budget NODES] [--diff | --noracle] [--cross]
 *                     [--dfirst]
 *                     seed0 count
 *
 * For every seed in [seed0, seed0 + count): generate a synthetic family of
 * n-sets (planted squares with Latin-square alternates, random transversals,
 * mutants and noise; random dense families; overlapping grids; wide label
 * universes up to W = 16; large sparse families; hubs with more than 255
 * vectors through one number; saturated counts at every unmatched cell,
 * --mode 6), find all semi-magic squares with a simple
 * backtracking oracle (r1 -> one col per cell of r1 -> exact cover by rows),
 * canonicalize them (up to transposition and row/col order), and compare
 * with what search_vectors reports for each option set in variants[] (and
 * check that every reported square is valid). --cross also checks the oracle
 * against a second, clique-based one; --diff compares the variants with
 * each other where the oracle exceeds its budget (--noracle: always; the
 * first selected variant is the reference). --dfirst tests the
 * diagonal-first search instead (see dfirst_check). Exit status is nonzero on
 * any mismatch, and the failing family is written to fail_<seed>.txt.
 *
 * Written for the code review of October 2026 (3000 seeds x 16 variants x 5
 * builds, plus builds with other LAZY_DEPTH, PERM_SEGS, NSLICE, CARRY_MAX_W:
 * no mismatch); ctest runs a few hundred seeds.
 */
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static double now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static int verbose;

#include "arrange.h"
#include "enumerate.h"
#include "dfirst.h"

/* ------------------------------------------------------------------ rng */
static uint64_t rs;
static uint64_t rnd(void) {
  uint64_t z = (rs += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}
static int rint_(int lo, int hi) { /* inclusive */
  return lo + (int)(rnd() % (uint64_t)(hi - lo + 1));
}
static void shuffle_int(int *a, int k) {
  for (int i = k - 1; i > 0; i--) {
    int j = (int)(rnd() % (uint64_t)(i + 1));
    int t = a[i];
    a[i] = a[j];
    a[j] = t;
  }
}

/* -------------------------------------------------------------- family */
#define MAXN 8
static int n;
static int plant_nt_max = 24;
static int U;              /* universe ids 0..U-1 */
static uint64_t *uval;     /* id -> value */
static int nv, capv;
static int (*vec)[MAXN];   /* ids, any order */

static void add_vec(const int *ids) {
  /* require distinct */
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (ids[i] == ids[j])
        return;
  if (nv == capv) {
    capv = capv ? 2 * capv : 256;
    vec = realloc(vec, capv * sizeof(*vec));
  }
  for (int i = 0; i < n; i++)
    vec[nv][i] = ids[i];
  nv++;
}

static int cmp_int(const void *a, const void *b) {
  int x = *(const int *)a, y = *(const int *)b;
  return (x > y) - (x < y);
}
static int cmp_vec(const void *a, const void *b) {
  const int *x = a, *y = b;
  for (int i = 0; i < n; i++)
    if (x[i] != y[i])
      return (x[i] > y[i]) - (x[i] < y[i]);
  return 0;
}
static void dedupe(void) {
  for (int v = 0; v < nv; v++)
    qsort(vec[v], n, sizeof(int), cmp_int);
  qsort(vec, nv, sizeof(*vec), cmp_vec);
  int k = 0;
  for (int v = 0; v < nv; v++)
    if (k == 0 || cmp_vec(vec[k - 1], vec[v]) != 0)
      memcpy(vec[k++], vec[v], sizeof(*vec));
  nv = k;
  /* shuffle order */
  for (int i = nv - 1; i > 0; i--) {
    int j = (int)(rnd() % (uint64_t)(i + 1));
    int t[MAXN];
    memcpy(t, vec[i], sizeof(t));
    memcpy(vec[i], vec[j], sizeof(t));
    memcpy(vec[j], t, sizeof(t));
  }
}

static void random_subset(int *out, int k, int lo, int hi) {
  /* k distinct ids in [lo, hi) */
  for (int i = 0; i < k; i++) {
    int x;
    bool again;
    do {
      x = lo + (int)(rnd() % (uint64_t)(hi - lo));
      again = false;
      for (int j = 0; j < i; j++)
        again |= out[j] == x;
    } while (again);
    out[i] = x;
  }
}

/* random latin square L[i][j] in 0..n-1 */
static void random_latin(int L[MAXN][MAXN]) {
  int pr[MAXN], pc[MAXN], ps[MAXN];
  for (int i = 0; i < n; i++)
    pr[i] = pc[i] = ps[i] = i;
  shuffle_int(pr, n);
  shuffle_int(pc, n);
  shuffle_int(ps, n);
  int s = rint_(1, n - 1);
  /* cyclic with step s works when gcd(s, n) == 1; otherwise use step 1 */
  int g = n, t = s;
  while (t) {
    int r = g % t;
    g = t;
    t = r;
  }
  if (g != 1)
    s = 1;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      L[pr[i]][pc[j]] = ps[(i * s + j) % n];
}

/* plant a grid G of ids; add its rows/cols and variations */
static void plant(int G[MAXN][MAXN], int extra_lo, int extra_hi) {
  int v[MAXN];
  int addr = rnd() % 8 != 0, addc = rnd() % 8 != 0;
  if (addr)
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++)
        v[j] = G[i][j];
      add_vec(v);
    }
  if (addc)
    for (int j = 0; j < n; j++) {
      for (int i = 0; i < n; i++)
        v[i] = G[i][j];
      add_vec(v);
    }
  /* alternate column sets / row sets via latin squares */
  int nl = rint_(0, 2);
  for (int t = 0; t < nl; t++) {
    int L[MAXN][MAXN];
    random_latin(L);
    int asrows = rnd() & 1;
    for (int sym = 0; sym < n; sym++) {
      if (rnd() % 10 == 0)
        continue; /* sometimes incomplete */
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          if (L[i][j] == sym)
            v[i] = asrows ? G[j][i] : G[i][j];
      add_vec(v);
    }
  }
  /* random transversals */
  int nt = rint_(0, plant_nt_max);
  for (int t = 0; t < nt; t++) {
    int asrows = rnd() & 1;
    for (int i = 0; i < n; i++) {
      int j = rint_(0, n - 1);
      v[i] = asrows ? G[j][i] : G[i][j];
    }
    add_vec(v);
  }
  /* mutants of rows/cols */
  int nm = rint_(0, n);
  for (int t = 0; t < nm; t++) {
    int i = rint_(0, n - 1), axis = rnd() & 1;
    for (int j = 0; j < n; j++)
      v[j] = axis ? G[j][i] : G[i][j];
    int k = rint_(0, n - 1);
    if (extra_hi > extra_lo && rnd() % 2)
      v[k] = rint_(extra_lo, extra_hi - 1);
    else
      v[k] = G[rint_(0, n - 1)][rint_(0, n - 1)];
    add_vec(v);
  }
}

static int style;

static void make_values(void) {
  uval = realloc(uval, U * sizeof(uint64_t));
  for (int i = 0; i < U; i++) {
    uint64_t x;
    bool again;
    do {
      switch (style) {
      case 0:
        x = (uint64_t)i + 1;
        break;
      case 1:
        x = (rnd() >> 1) | 1;
        break;
      case 2:
        x = ((uint64_t)rint_(1, 4 * U) << 32); /* low bits all zero */
        break;
      default:
        x = (uint64_t)rint_(1, 1000000);
        break;
      }
      again = false;
      for (int j = 0; j < i && !again; j++)
        again = uval[j] == x;
    } while (again);
    uval[i] = x;
  }
}

/*
 * Mode 6: saturated counts at the node below r1. The n numbers of a vector
 * r1 are the rarest of the family, so that they take the top labels and r1
 * is the first row of the squares through it, and each of them is in more
 * than 255 vectors meeting r1 only there: with r1 placed, every unmatched
 * cell has a saturated byte counter (>= 255 candidates), so MRV picks the
 * cell among the saturated ones (by the largest label). Every square goes
 * through r1, and there is a planted one. Ids: r1 = 0 .. n-1, a core K of
 * 11 numbers, and E: n groups of n - 1 numbers. Vectors: a random ~60-90%
 * of the n-sets with one number of r1 and n - 1 of K (the cols through the
 * cells of r1: at least 262 per cell), all n-sets of K (rows for the
 * support filter, and to make K common), the planted square (r1, the cols
 * {r1_i} + group i, the rows taking the k-th number of every group), and
 * random n-sets of E, enough to make the numbers of E more common than
 * those of r1 (these also make more squares through r1, with the E rows
 * that take one number of each group). n = 5 (or 6 with --n 6), < 64
 * labels: the bug this catches (October 2026: the masked minimum of the
 * cell choice carried other words' counts into the lanes of cells that
 * are not unmatched, and the scan for the largest saturated label could
 * then pick one of those, dropping the whole subtree) shows with the words
 * past the labels, with every --min-words.
 */
static int sat_rows[MAXN][MAXN], sat_cols[MAXN][MAXN]; /* planted, ids */
static void gen_saturated(void) {
  const int K = 11, E = n * (n - 1);
  U = n + K + E;
  int v[MAXN];
  /* the cols through the cells of r1: for each cell i, target of the
   * C(11, n - 1) sets of n - 1 numbers of K, chosen at random */
  static int masks[512];
  int nm = 0;
  for (int m = 0; m < 1 << K; m++)
    if (__builtin_popcount(m) == n - 1)
      masks[nm++] = m;
  const int target = rint_(262, 300);
  for (int i = 0; i < n; i++) {
    shuffle_int(masks, nm);
    for (int t = 0; t < target; t++) {
      int k = 0;
      v[k++] = i;
      for (int b = 0; b < K; b++)
        if (masks[t] >> b & 1)
          v[k++] = n + b;
      add_vec(v);
    }
  }
  /* the n-sets of K */
  for (int m = 0; m < 1 << K; m++) {
    if (__builtin_popcount(m) != n)
      continue;
    int k = 0;
    for (int b = 0; b < K; b++)
      if (m >> b & 1)
        v[k++] = n + b;
    add_vec(v);
  }
  const int e0 = n + K;
  /* the planted square: r1, the cols {i} + group i, the rows k */
  for (int i = 0; i < n; i++)
    v[i] = sat_rows[0][i] = i;
  add_vec(v);
  for (int i = 0; i < n; i++) {
    v[0] = i;
    for (int k = 0; k < n - 1; k++)
      v[1 + k] = e0 + i * (n - 1) + k;
    memcpy(sat_cols[i], v, sizeof(v));
    add_vec(v);
  }
  for (int k = 0; k < n - 1; k++) {
    for (int i = 0; i < n; i++)
      v[i] = e0 + i * (n - 1) + k;
    memcpy(sat_rows[1 + k], v, sizeof(v));
    add_vec(v);
  }
  /* each number of E in ~ N n / E = N / (n - 1) of them (a few percent
   * fewer after dedupe): 100-160 more than a number of r1 has */
  const int N = (target + rint_(100, 160)) * (n - 1);
  for (int t = 0; t < N; t++) {
    random_subset(v, n, e0, U);
    add_vec(v);
  }
}

/* modes: 0 = one dense grid + transversals, 1 = many planted grids (wide),
 * 2 = random dense, 3 = a few grids overlapping in a small universe; only
 * with --mode: 4 = large sparse, 5 = hubs, 6 = saturated counts (see
 * gen_saturated), 7 = many planted grids over 257-500 labels (the carried
 * path with 5-8 words, its r1 searched at 2 to 8 words) */
static int gen(uint64_t seed, int force_mode, int force_n) {
  rs = seed * 0x2545F4914F6CDD1Dull + 12345;
  nv = 0;
  n = force_n ? force_n : rint_(4, 7);
  if (!force_n && rnd() % 10 == 0)
    n = rnd() % 2 ? 3 : 8;
  int mode = force_mode >= 0 ? force_mode : rint_(0, 3);
  plant_nt_max = mode == 1 || mode == 4 || mode == 7 ? n : 3 * n;
  style = rint_(0, 3);
  int G[MAXN][MAXN];
  static int perm[8192];
  if (mode == 0) {
    int extra = rint_(0, 2 * n);
    U = n * n + extra;
    for (int i = 0; i < U; i++)
      perm[i] = i;
    shuffle_int(perm, U);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        G[i][j] = perm[i * n + j];
    plant(G, n * n, U);
    int noise = rint_(0, 3 * n);
    for (int t = 0; t < noise; t++) {
      int v[MAXN];
      random_subset(v, n, 0, U);
      add_vec(v);
    }
  } else if (mode == 1 || mode == 7) {
    /* target label ranges: W=2, 3, 4, 8, 16 (mode 7: 257-500) */
    static const int lo_[] = {40, 129, 193, 257, 513};
    static const int hi_[] = {128, 192, 256, 512, 1100};
    int r = rint_(0, 4);
    U = mode == 7 ? rint_(257, 500) : rint_(lo_[r], hi_[r]);
    if (U < n * n)
      U = n * n;
    /* grids over random subsets; enough grids to cover most ids */
    int k = U / (n * n) + rint_(0, 3);
    for (int g = 0; g < k; g++) {
      for (int i = 0; i < U; i++)
        perm[i] = i;
      /* partial shuffle: first n^2 */
      for (int i = 0; i < n * n; i++) {
        int j = i + (int)(rnd() % (uint64_t)(U - i));
        int t = perm[i];
        perm[i] = perm[j];
        perm[j] = t;
      }
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          G[i][j] = perm[i * n + j];
      plant(G, 0, U);
    }
    /* cover every id at least once (so L == U, roughly) */
    static char used[4096];
    memset(used, 0, sizeof(used));
    for (int v = 0; v < nv; v++)
      for (int i = 0; i < n; i++)
        used[vec[v][i]] = 1;
    for (int x = 0; x < U; x++)
      if (!used[x] || rnd() % 8 == 0) {
        int v[MAXN];
        random_subset(v, n, 0, U);
        v[0] = x;
        bool dup = false;
        for (int j = 1; j < n; j++)
          dup |= v[j] == x;
        if (!dup)
          add_vec(v);
      }
  } else if (mode == 4) {
    /* large sparse: N in the thousands, L up to 4096 */
    U = rint_(150, 1000);
    int k = rint_(1100, 3600) / (2 * n + n / 2 + 4);
    for (int g = 0; g < k; g++) {
      for (int i = 0; i < n * n; i++) {
        bool again;
        do {
          perm[i] = rint_(0, U - 1);
          again = false;
          for (int j = 0; j < i; j++)
            again |= perm[j] == perm[i];
        } while (again);
      }
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          G[i][j] = perm[i * n + j];
      plant(G, 0, U);
    }
    int noise = rint_(0, 200);
    for (int t = 0; t < noise; t++) {
      int v[MAXN];
      random_subset(v, n, 0, U);
      add_vec(v);
    }
  } else if (mode == 5) {
    /* hubs: > 255 vectors through one number (saturating byte counters) */
    U = rint_(n * n + 20, 250);
    int k = rint_(1, 3);
    for (int g = 0; g < k; g++) {
      for (int i = 0; i < U; i++)
        perm[i] = i;
      for (int i = 0; i < n * n; i++) {
        int j = i + (int)(rnd() % (uint64_t)(U - i));
        int t = perm[i];
        perm[i] = perm[j];
        perm[j] = t;
      }
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          G[i][j] = perm[i * n + j];
      plant(G, 0, U);
    }
    int nh = rint_(1, 3);
    for (int h = 0; h < nh; h++) {
      int hub = G[rint_(0, n - 1)][rint_(0, n - 1)];
      int cntv = rint_(256, 600);
      for (int t = 0; t < cntv; t++) {
        int v[MAXN];
        random_subset(v, n, 0, U);
        bool dup = false;
        for (int j = 1; j < n; j++)
          dup |= v[j] == hub;
        if (dup)
          continue;
        v[0] = hub;
        add_vec(v);
      }
    }
  } else if (mode == 6) {
    n = force_n == 6 ? 6 : 5;
    gen_saturated();
  } else if (mode == 2) {
    U = n * n + rint_(0, 2 * n);
    int N = rint_(2 * n, n == 3 ? 60 : n == 4 ? 160 : 260);
    for (int t = 0; t < N; t++) {
      int v[MAXN];
      random_subset(v, n, 0, U);
      add_vec(v);
    }
  } else {
    /* a few grids sharing many numbers */
    int extra = rint_(0, 3 * n);
    U = n * n + extra;
    int k = rint_(2, 4);
    for (int g = 0; g < k; g++) {
      for (int i = 0; i < U; i++)
        perm[i] = i;
      for (int i = 0; i < n * n; i++) {
        int j = i + (int)(rnd() % (uint64_t)(U - i));
        int t = perm[i];
        perm[i] = perm[j];
        perm[j] = t;
      }
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          G[i][j] = perm[i * n + j];
      plant(G, 0, U);
    }
    int noise = rint_(0, 2 * n);
    for (int t = 0; t < noise; t++) {
      int v[MAXN];
      random_subset(v, n, 0, U);
      add_vec(v);
    }
  }
  make_values();
  dedupe();
  return mode;
}

/* -------------------------------------------------------- canonical forms */
typedef struct {
  uint64_t v[2 * MAXN * MAXN];
} canon_t;
static int canon_len;

static int cmp_u64(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}
static int cmp_line(const void *a, const void *b) {
  const uint64_t *x = a, *y = b;
  for (int i = 0; i < n; i++)
    if (x[i] != y[i])
      return (x[i] > y[i]) - (x[i] < y[i]);
  return 0;
}
static int cmp_canon(const void *a, const void *b) {
  const canon_t *x = a, *y = b;
  for (int i = 0; i < canon_len; i++)
    if (x->v[i] != y->v[i])
      return (x->v[i] > y->v[i]) - (x->v[i] < y->v[i]);
  return 0;
}
static void canonicalize(uint64_t rows[MAXN][MAXN], uint64_t cols[MAXN][MAXN],
                         canon_t *out) {
  uint64_t r[MAXN * MAXN], c[MAXN * MAXN];
  for (int i = 0; i < n; i++) {
    memcpy(r + i * n, rows[i], n * 8);
    memcpy(c + i * n, cols[i], n * 8);
    qsort(r + i * n, n, 8, cmp_u64);
    qsort(c + i * n, n, 8, cmp_u64);
  }
  qsort(r, n, n * 8, cmp_line);
  qsort(c, n, n * 8, cmp_line);
  int rc = memcmp(r, c, 0); /* placeholder */
  rc = 0;
  for (int i = 0; i < n * n && !rc; i++)
    rc = (r[i] > c[i]) - (r[i] < c[i]);
  memset(out, 0, sizeof(*out));
  memcpy(out->v, rc <= 0 ? r : c, n * n * 8);
  memcpy(out->v + n * n, rc <= 0 ? c : r, n * n * 8);
}

typedef struct {
  canon_t *a;
  int k, cap;
} clist_t;
static void cl_push(clist_t *l, const canon_t *c) {
  if (l->k == l->cap) {
    l->cap = l->cap ? 2 * l->cap : 64;
    l->a = realloc(l->a, l->cap * sizeof(canon_t));
  }
  l->a[l->k++] = *c;
}

/* ----------------------------------------------------------------- oracle */
static int UW;
static uint64_t *vb;   /* nv * UW bitsets over ids */
static int *occ_first, *occ; /* CSR id -> vectors */
static clist_t oracle_out;
static uint64_t oracle_work, oracle_budget;
static int oracle_blown;

static uint8_t *imat; /* nv x nv intersection sizes, or NULL */
static inline int inter_cnt(int a, int b) {
  if (imat)
    return imat[(size_t)a * nv + b];
  int c = 0;
  for (int w = 0; w < UW; w++)
    c += __builtin_popcountll(vb[a * UW + w] & vb[b * UW + w]);
  return c;
}

static int o_r1, o_cols[MAXN], o_rows[MAXN];
static int o_nrows;
static uint64_t o_union[64], o_cov[64];
static int *o_cand[MAXN], o_ncand[MAXN];
static int *o_rowcand, o_nrowcand;

static void o_emit(void) {
  uint64_t rows[MAXN][MAXN], cols[MAXN][MAXN];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      rows[i][j] = uval[vec[o_rows[i]][j]];
      cols[i][j] = uval[vec[o_cols[i]][j]];
    }
  canon_t c;
  canonicalize(rows, cols, &c);
  cl_push(&oracle_out, &c);
}

/* exact cover of union \ covered by rows from o_rowcand */
static void o_rows_rec(void) {
  if (++oracle_work > oracle_budget) {
    oracle_blown = 1;
    return;
  }
  if (o_nrows == n) {
    o_emit();
    return;
  }
  /* smallest uncovered id */
  int y = -1;
  for (int w = 0; w < UW && y < 0; w++) {
    uint64_t m = o_union[w] & ~o_cov[w];
    if (m)
      y = w * 64 + __builtin_ctzll(m);
  }
  if (y < 0)
    return;
  for (int t = 0; t < o_nrowcand; t++) {
    int r = o_rowcand[t];
    if (!((vb[r * UW + y / 64] >> (y % 64)) & 1))
      continue;
    bool ok = true;
    for (int w = 0; w < UW && ok; w++)
      ok = (vb[r * UW + w] & o_cov[w]) == 0;
    if (!ok)
      continue;
    for (int w = 0; w < UW; w++)
      o_cov[w] |= vb[r * UW + w];
    o_rows[o_nrows++] = r;
    o_rows_rec();
    o_nrows--;
    for (int w = 0; w < UW; w++)
      o_cov[w] &= ~vb[r * UW + w];
  }
}

static void o_cols_rec(int k) {
  if (oracle_blown)
    return;
  if (++oracle_work > oracle_budget) {
    oracle_blown = 1;
    return;
  }
  if (k == n) {
    /* union of cols */
    for (int w = 0; w < UW; w++) {
      o_union[w] = 0;
      for (int j = 0; j < n; j++)
        o_union[w] |= vb[o_cols[j] * UW + w];
    }
    /* row candidates: subsets of union, != r1, disjoint from r1, meeting
     * every col exactly once */
    o_nrowcand = 0;
    /* every row meets col 0 exactly once: scan the vectors through its
     * cells (each such vector meets col 0 in one cell, so it is seen once) */
    for (int i0 = 0; i0 < n; i0++)
    for (int t0 = occ_first[vec[o_cols[0]][i0]]; t0 < occ_first[vec[o_cols[0]][i0] + 1]; t0++) {
      int v = occ[t0];
      if (v == o_r1)
        continue;
      bool ok = true;
      for (int w = 0; w < UW && ok; w++)
        ok = (vb[v * UW + w] & ~o_union[w]) == 0 &&
             (vb[v * UW + w] & vb[o_r1 * UW + w]) == 0;
      for (int j = 0; j < n && ok; j++)
        ok = inter_cnt(v, o_cols[j]) == 1;
      if (ok)
        o_rowcand[o_nrowcand++] = v;
    }
    for (int w = 0; w < UW; w++)
      o_cov[w] = vb[o_r1 * UW + w];
    o_rows[0] = o_r1;
    o_nrows = 1;
    o_rows_rec();
    return;
  }
  for (int t = 0; t < o_ncand[k]; t++) {
    int c = o_cand[k][t];
    bool ok = true;
    for (int j = 0; j < k && ok; j++)
      ok = inter_cnt(c, o_cols[j]) == 0;
    if (!ok)
      continue;
    o_cols[k] = c;
    /* every cell y of c off r1 needs a possible row through it: disjoint
     * from r1 and meeting each chosen col exactly once */
    for (int i = 0; i < n && ok; i++) {
      int y = vec[c][i];
      if ((vb[o_r1 * UW + y / 64] >> (y % 64)) & 1)
        continue;
      bool any = false;
      for (int t2 = occ_first[y]; t2 < occ_first[y + 1] && !any; t2++) {
        int r = occ[t2];
        if (r == o_r1 || inter_cnt(r, o_r1) != 0)
          continue;
        bool g = true;
        for (int j = 0; j <= k && g; j++)
          g = inter_cnt(r, o_cols[j]) == 1;
        any = g;
      }
      ok = any;
    }
    if (!ok)
      continue;
    o_cols_rec(k + 1);
  }
}

/* returns 0 if ok, 1 if the budget was exceeded */
static int oracle(uint64_t budget) {
  UW = (U + 63) / 64;
  free(vb);
  vb = calloc((size_t)nv * UW, 8);
  for (int v = 0; v < nv; v++)
    for (int i = 0; i < n; i++)
      vb[v * UW + vec[v][i] / 64] |= 1ull << (vec[v][i] % 64);
  free(occ_first);
  free(occ);
  occ_first = calloc(U + 1, sizeof(int));
  occ = malloc((size_t)nv * n * sizeof(int));
  for (int v = 0; v < nv; v++)
    for (int i = 0; i < n; i++)
      occ_first[vec[v][i] + 1]++;
  for (int x = 0; x < U; x++)
    occ_first[x + 1] += occ_first[x];
  int *fill = malloc((U + 1) * sizeof(int));
  memcpy(fill, occ_first, (U + 1) * sizeof(int));
  for (int v = 0; v < nv; v++)
    for (int i = 0; i < n; i++)
      occ[fill[vec[v][i]]++] = v;
  free(fill);
  free(imat);
  imat = NULL;
  if (nv <= 6000) {
    uint8_t *m = calloc((size_t)nv * nv, 1);
    for (int x = 0; x < U; x++)
      for (int a = occ_first[x]; a < occ_first[x + 1]; a++)
        for (int b = occ_first[x]; b < occ_first[x + 1]; b++)
          m[(size_t)occ[a] * nv + occ[b]]++;
    imat = m;
  }
  for (int k = 0; k < MAXN; k++) {
    free(o_cand[k]);
    o_cand[k] = malloc(nv * sizeof(int));
  }
  free(o_rowcand);
  o_rowcand = malloc(nv * sizeof(int));
  oracle_out.k = 0;
  oracle_work = 0;
  oracle_budget = budget;
  oracle_blown = 0;
  for (o_r1 = 0; o_r1 < nv && !oracle_blown; o_r1++) {
    for (int k = 0; k < n; k++) {
      int x = vec[o_r1][k];
      o_ncand[k] = 0;
      for (int t = occ_first[x]; t < occ_first[x + 1]; t++) {
        int c = occ[t];
        if (c != o_r1 && inter_cnt(c, o_r1) == 1)
          o_cand[k][o_ncand[k]++] = c;
      }
    }
    o_cols_rec(0);
  }
  if (oracle_blown)
    return 1;
  canon_len = 2 * n * n;
  if (oracle_out.k) /* (a is NULL before the first square) */
    qsort(oracle_out.a, oracle_out.k, sizeof(canon_t), cmp_canon);
  int k = 0;
  for (int i = 0; i < oracle_out.k; i++)
    if (k == 0 || cmp_canon(&oracle_out.a[k - 1], &oracle_out.a[i]) != 0)
      oracle_out.a[k++] = oracle_out.a[i];
  oracle_out.k = k;
  return 0;
}

/* --------------------------------------------- second oracle: cliques */
static clist_t oracle2_out;
static int *cl_list; /* cliques: n ids each */
static int ncl, capcl;
static int cl_cur[MAXN];
static int clq_blown;
static void clq_rec(int k, int from, uint64_t *used) {
  if (clq_blown)
    return;
  if (k == n) {
    if (ncl >= 6000) {
      clq_blown = 1;
      return;
    }
    if (ncl == capcl) {
      capcl = capcl ? 2 * capcl : 1024;
      cl_list = realloc(cl_list, (size_t)capcl * MAXN * sizeof(int));
    }
    memcpy(cl_list + (size_t)ncl * n, cl_cur, n * sizeof(int));
    ncl++;
    return;
  }
  for (int v = from; v < nv; v++) {
    bool ok = true;
    for (int w = 0; w < UW && ok; w++)
      ok = (vb[v * UW + w] & used[w]) == 0;
    if (!ok)
      continue;
    for (int w = 0; w < UW; w++)
      used[w] |= vb[v * UW + w];
    cl_cur[k] = v;
    clq_rec(k + 1, v + 1, used);
    for (int w = 0; w < UW; w++)
      used[w] &= ~vb[v * UW + w];
  }
}
static int oracle2(void) {
  ncl = 0;
  clq_blown = 0;
  uint64_t used[64] = {0};
  clq_rec(0, 0, used);
  if (clq_blown)
    return 1;
  oracle2_out.k = 0;
  /* all pairs of cliques with every cross intersection == 1 */
  for (int a = 0; a < ncl; a++)
    for (int b = a + 1; b < ncl; b++) {
      bool ok = true;
      for (int i = 0; i < n && ok; i++)
        for (int j = 0; j < n && ok; j++)
          ok = inter_cnt(cl_list[(size_t)a * n + i], cl_list[(size_t)b * n + j]) == 1;
      if (!ok)
        continue;
      uint64_t rows[MAXN][MAXN], cols[MAXN][MAXN];
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          rows[i][j] = uval[vec[cl_list[(size_t)a * n + i]][j]];
          cols[i][j] = uval[vec[cl_list[(size_t)b * n + i]][j]];
        }
      canon_t c;
      canonicalize(rows, cols, &c);
      cl_push(&oracle2_out, &c);
    }
  canon_len = 2 * n * n;
  if (oracle2_out.k)
    qsort(oracle2_out.a, oracle2_out.k, sizeof(canon_t), cmp_canon);
  return 0;
}

/* ------------------------------------------------------------------ search */
static clist_t found;
static int crossed;
static uint64_t *inset; /* sorted value tuples of the input vectors */
static int invalid;
static int cmp_tuple(const void *a, const void *b) { return cmp_line(a, b); }
static void build_inset(void) {
  free(inset);
  inset = malloc((size_t)nv * n * 8);
  for (int v = 0; v < nv; v++) {
    for (int i = 0; i < n; i++)
      inset[(size_t)v * n + i] = uval[vec[v][i]];
    qsort(inset + (size_t)v * n, n, 8, cmp_u64);
  }
  qsort(inset, nv, n * 8, cmp_tuple);
}
static bool is_input(const uint64_t *t) {
  uint64_t x[MAXN];
  memcpy(x, t, n * 8);
  qsort(x, n, 8, cmp_u64);
  return bsearch(x, inset, nv, n * 8, cmp_tuple) != NULL;
}
static int meet(const uint64_t *a, const uint64_t *b) {
  int c = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      c += a[i] == b[j];
  return c;
}
static int cb(const square_t *sq, void *ctx) {
  (void)ctx;
  uint64_t rows[MAXN][MAXN], cols[MAXN][MAXN];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      rows[i][j] = sq->rows[i][j];
      cols[i][j] = sq->cols[i][j];
    }
  bool ok = sq->n == n;
  for (int i = 0; i < n && ok; i++) {
    ok = is_input(rows[i]) && is_input(cols[i]);
    for (int j = 0; j < n && ok; j++) {
      ok = meet(rows[i], cols[j]) == 1;
      if (ok && j != i)
        ok = meet(rows[i], rows[j]) == 0 && meet(cols[i], cols[j]) == 0;
    }
  }
  if (!ok)
    invalid++;
  canon_t c;
  canonicalize(rows, cols, &c);
  cl_push(&found, &c);
  return 0;
}

typedef struct {
  const char *name;
  int fc, mrv, support, min_words, gather, cross;
  int fixw;    /* 1: opts.r1_width = 0 (every r1 at the width of all labels) */
  int pretest; /* opts.pretest_min: 0 = the default, -1 = never */
} variant_t;
static const variant_t variants[] = {
    {"default", 1, 1, 1, 0, 0, 1},  {"nofc", 0, 1, 1, 0, 0, 1},
    {"nomrv", 1, 0, 1, 0, 0, 1},    {"nosup", 1, 1, 0, 0, 0, 1},
    {"nofc_nomrv", 0, 0, 1, 0, 0, 1}, {"w3", 1, 1, 1, 3, 0, 1},
    {"w4", 1, 1, 1, 4, 0, 1},       {"w8", 1, 1, 1, 8, 0, 1},
    {"w8g", 1, 1, 1, 8, 1, 1},      {"w16", 1, 1, 1, 16, 0, 1},
    {"w64", 1, 1, 1, 64, 0, 1},     {"w8_nosup", 1, 1, 0, 8, 0, 1},
    {"w16_nomrv", 1, 0, 1, 16, 0, 1}, {"w4_nosup", 1, 1, 0, 4, 0, 1},
    {"w3_nomrv", 1, 0, 1, 3, 0, 1}, {"w8_nofc", 0, 1, 1, 8, 0, 1},
    {"nocross", 1, 1, 1, 0, 0, 0},  {"xcross", 1, 1, 1, 0, 0, 2},
    {"xcross_nomrv", 1, 0, 1, 0, 0, 2}, {"xcross_w3", 1, 1, 1, 3, 0, 2},
    {"xcross_w4", 1, 1, 1, 4, 0, 2}, {"xcross_w8", 1, 1, 1, 8, 0, 2},
    {"xcross_w16", 1, 1, 1, 16, 0, 2},
    {"w5", 1, 1, 1, 5, 0, 1},       {"w6", 1, 1, 1, 6, 0, 1},
    {"w7", 1, 1, 1, 7, 0, 1},       {"fixw", 1, 1, 1, 0, 0, 1, 1},
    {"w5_fixw", 1, 1, 1, 5, 0, 1, 1}, {"w6_nomrv", 1, 0, 1, 6, 0, 1},
    {"xcross_w5", 1, 1, 1, 5, 0, 2},
    /* the pretest of the deep children from other depths (it must not
     * change the squares) */
    {"pt_off", 1, 1, 1, 0, 0, 1, 0, -1}, {"pt_1", 1, 1, 1, 0, 0, 1, 0, 1},
    {"pt_3", 1, 1, 1, 0, 0, 1, 0, 3},    {"pt_4", 1, 1, 1, 0, 0, 1, 0, 4},
    {"pt_1_w4", 1, 1, 1, 4, 0, 1, 0, 1}, {"pt_2_nosup", 1, 1, 0, 0, 0, 1, 0, 2},
    {"pt_2_nomrv", 1, 0, 1, 0, 0, 1, 0, 2}, {"pt_2_xcross", 1, 1, 1, 0, 0, 2, 0, 2},
};
#define NVAR (int)(sizeof(variants) / sizeof(variants[0]))

static void dump(const char *path, uint64_t seed) {
  FILE *f = fopen(path, "w");
  if (!f)
    return;
  fprintf(f, "# seed %lu n %d nv %d\n", (unsigned long)seed, n, nv);
  for (int v = 0; v < nv; v++) {
    for (int i = 0; i < n; i++)
      fprintf(f, "%lu%c", (unsigned long)uval[vec[v][i]], i + 1 < n ? ' ' : '\n');
  }
  fclose(f);
}

/* ------------------------------------------------------------ d-first */
/*
 * --dfirst: the diagonal-first search (dfirst.h) against brute force. A
 * random set of candidate diagonals (random n-sets, numbers outside the
 * family, traversals of the oracle's squares, sometimes the family's own
 * vectors, and a planted pair: a square's traversals t1, t2 that fit on
 * the diagonal and the anti-diagonal together); dfirst_search over the
 * whole set must report exactly the multiset of (square, traversal in the
 * set) pairs, with the number of the square's traversals in the set and
 * the partner flag (another traversal in the set on the other diagonal),
 * which brute force computes from all pairs of the square's traversals in
 * the set. The planted square must be reported twice and flagged both
 * times; a d-stride 3 split (offsets 0, 1, 2) must add up to the same.
 */
typedef struct {
  canon_t sq;
  uint64_t d[MAXN]; /* ascending */
  int partner, set_count;
} dpair_t;
typedef struct {
  dpair_t *a;
  int k, cap;
} dplist_t;
static void dp_push(dplist_t *l, const dpair_t *p) {
  if (l->k == l->cap) {
    l->cap = l->cap ? 2 * l->cap : 64;
    l->a = realloc(l->a, l->cap * sizeof(dpair_t));
  }
  l->a[l->k++] = *p;
}
static int cmp_dpair(const void *a, const void *b) {
  const dpair_t *x = a, *y = b;
  int c = cmp_canon(&x->sq, &y->sq);
  if (c)
    return c;
  for (int i = 0; i < n; i++)
    if (x->d[i] != y->d[i])
      return (x->d[i] > y->d[i]) - (x->d[i] < y->d[i]);
  if (x->partner != y->partner)
    return x->partner - y->partner;
  return x->set_count - y->set_count;
}
static uint64_t *dset; /* the diagonal set: ascending tuples, sorted */
static int ndset, capdset;
static bool in_dset(const uint64_t *t0) {
  uint64_t t[MAXN];
  memcpy(t, t0, n * 8);
  qsort(t, n, 8, cmp_u64);
  return ndset && bsearch(t, dset, ndset, n * 8, cmp_tuple) != NULL;
}
static void dset_add(const uint64_t *t0) {
  uint64_t t[MAXN];
  memcpy(t, t0, n * 8);
  qsort(t, n, 8, cmp_u64);
  for (int i = 0; i + 1 < n; i++)
    if (t[i] == t[i + 1])
      return; /* not an n-set */
  if ((size_t)(ndset + 1) * n > (size_t)capdset) { /* capdset: words */
    capdset = capdset ? 2 * capdset : 256 * MAXN;
    while ((size_t)(ndset + 1) * n > (size_t)capdset)
      capdset *= 2;
    dset = realloc(dset, (size_t)capdset * 8);
  }
  memcpy(dset + (size_t)ndset * n, t, n * 8);
  ndset++;
}
/* the grid of a canonical form: first block as rows, second as cols */
static void canon_grid(const canon_t *c, uint64_t g[MAXN][MAXN]) {
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      g[i][j] = 0;
      for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
          if (c->v[i * n + a] == c->v[n * n + j * n + b])
            g[i][j] = c->v[i * n + a];
    }
}
static int next_perm_i(int *p, int k) {
  int i = k - 2;
  while (i >= 0 && p[i] >= p[i + 1])
    i--;
  if (i < 0)
    return 0;
  int j = k - 1;
  while (p[j] <= p[i])
    j--;
  int t = p[i];
  p[i] = p[j];
  p[j] = t;
  for (int a = i + 1, b = k - 1; a < b; a++, b--) {
    t = p[a];
    p[a] = p[b];
    p[b] = t;
  }
  return 1;
}
/* can traversals p1, p2 (row -> col) be the diagonal and anti-diagonal
 * after permuting rows and cols: p1^-1 p2 an involution with n mod 2 fixed
 * points */
static bool diag_pair(const int *p1, const int *p2) {
  int inv1[MAXN], rel[MAXN], fixed = 0;
  for (int i = 0; i < n; i++)
    inv1[p1[i]] = i;
  for (int r = 0; r < n; r++)
    rel[r] = inv1[p2[r]];
  for (int r = 0; r < n; r++) {
    if (rel[rel[r]] != r)
      return false;
    fixed += rel[r] == r;
  }
  return fixed == n % 2;
}
static dplist_t d_brute, d_found;
static int d_invalid;
static void brute_pairs(void) {
  d_brute.k = 0;
  static int (*tp)[MAXN];
  static int captp;
  for (int s = 0; s < oracle_out.k; s++) {
    uint64_t g[MAXN][MAXN];
    canon_grid(&oracle_out.a[s], g);
    int perm[MAXN], k = 0;
    for (int i = 0; i < n; i++)
      perm[i] = i;
    do {
      uint64_t t[MAXN];
      for (int i = 0; i < n; i++)
        t[i] = g[i][perm[i]];
      if (in_dset(t)) {
        if (k == captp) {
          captp = captp ? 2 * captp : 64;
          tp = realloc(tp, captp * sizeof(*tp));
        }
        memcpy(tp[k++], perm, sizeof(perm));
      }
    } while (next_perm_i(perm, n));
    for (int a = 0; a < k; a++) {
      dpair_t dp;
      memset(&dp, 0, sizeof(dp));
      dp.sq = oracle_out.a[s];
      for (int i = 0; i < n; i++)
        dp.d[i] = g[i][tp[a][i]];
      qsort(dp.d, n, 8, cmp_u64);
      dp.set_count = k;
      for (int b = 0; b < k && !dp.partner; b++)
        dp.partner = b != a && diag_pair(tp[a], tp[b]);
      dp_push(&d_brute, &dp);
    }
  }
  canon_len = 2 * n * n;
  if (d_brute.k)
    qsort(d_brute.a, d_brute.k, sizeof(dpair_t), cmp_dpair);
}
static int d_cb(const dsquare_t *ds, void *ctx) {
  (void)ctx;
  const square_t *sq = ds->sq;
  uint64_t rows[MAXN][MAXN], cols[MAXN][MAXN];
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      rows[i][j] = sq->rows[i][j];
      cols[i][j] = sq->cols[i][j];
    }
  bool ok = sq->n == n && in_dset(ds->d);
  for (int i = 0; i < n && ok; i++) {
    ok = is_input(rows[i]) && is_input(cols[i]) &&
         meet(rows[i], ds->d) == 1 && meet(cols[i], ds->d) == 1;
    for (int j = 0; j < n && ok; j++) {
      ok = meet(rows[i], cols[j]) == 1;
      if (ok && j != i)
        ok = meet(rows[i], rows[j]) == 0 && meet(cols[i], cols[j]) == 0;
    }
  }
  if (!ok)
    d_invalid++;
  dpair_t dp;
  memset(&dp, 0, sizeof(dp));
  canonicalize(rows, cols, &dp.sq);
  memcpy(dp.d, ds->d, n * 8);
  qsort(dp.d, n, 8, cmp_u64);
  dp.partner = ds->partner;
  dp.set_count = ds->set_count;
  dp_push(&d_found, &dp);
  return 0;
}
static bool same_dpairs(void) {
  canon_len = 2 * n * n;
  if (d_found.k)
    qsort(d_found.a, d_found.k, sizeof(dpair_t), cmp_dpair);
  if (d_found.k != d_brute.k)
    return false;
  for (int i = 0; i < d_found.k; i++)
    if (cmp_dpair(&d_found.a[i], &d_brute.a[i]))
      return false;
  return true;
}
static uint64_t d_pairs_total, d_partner_total, d_planted;
/* returns the number of failures */
static int dfirst_check(uint64_t seed, const vec_list_t *l,
                        const search_opts_t *o, const char *vname,
                        int top_root) {
  /* the diagonal set */
  ndset = 0;
  const int nr = rint_(0, 2 * n);
  for (int t = 0; t < nr; t++) {
    int ids[MAXN];
    uint64_t v[MAXN];
    random_subset(ids, n, 0, U);
    for (int i = 0; i < n; i++)
      v[i] = uval[ids[i]];
    if (rnd() % 4 == 0)
      v[rint_(0, n - 1)] = rnd() | (1ull << 63); /* not in the family */
    dset_add(v);
  }
  if (rnd() % 3 == 0)
    for (int v = 0; v < nv; v++) {
      uint64_t t[MAXN];
      for (int i = 0; i < n; i++)
        t[i] = uval[vec[v][i]];
      dset_add(t);
    }
  int planted = -1;
  uint64_t pt[2][MAXN];
  for (int s = 0; s < oracle_out.k; s++) {
    uint64_t g[MAXN][MAXN];
    canon_grid(&oracle_out.a[s], g);
    const int k = rint_(0, 2);
    for (int t = 0; t < k; t++) {
      int perm[MAXN];
      for (int i = 0; i < n; i++)
        perm[i] = i;
      shuffle_int(perm, n);
      uint64_t v[MAXN];
      for (int i = 0; i < n; i++)
        v[i] = g[i][perm[i]];
      dset_add(v);
    }
  }
  if (oracle_out.k) {
    /* a planted diagonal pair: sigma and sigma o tau, tau an involution
     * with n mod 2 fixed points */
    planted = rint_(0, oracle_out.k - 1);
    uint64_t g[MAXN][MAXN];
    canon_grid(&oracle_out.a[planted], g);
    int sigma[MAXN], ord[MAXN], tau[MAXN];
    for (int i = 0; i < n; i++)
      sigma[i] = ord[i] = i;
    shuffle_int(sigma, n);
    shuffle_int(ord, n);
    for (int i = 0; i + 1 < n; i += 2) {
      tau[ord[i]] = ord[i + 1];
      tau[ord[i + 1]] = ord[i];
    }
    if (n % 2)
      tau[ord[n - 1]] = ord[n - 1];
    for (int i = 0; i < n; i++) {
      pt[0][i] = g[i][sigma[i]];
      pt[1][i] = g[i][sigma[tau[i]]];
    }
    dset_add(pt[0]);
    dset_add(pt[1]);
  }
  qsort(dset, ndset, n * 8, cmp_tuple);
  int k = 0;
  for (int i = 0; i < ndset; i++)
    if (!k || cmp_tuple(dset + (size_t)(k - 1) * n, dset + (size_t)i * n))
      memmove(dset + (size_t)k++ * n, dset + (size_t)i * n, n * 8);
  ndset = k;
  /* the list given to dfirst: shuffled, each descending */
  vec_list_t dl;
  vec_list_init(&dl, n);
  int *ord = malloc((ndset + 1) * sizeof(int));
  for (int i = 0; i < ndset; i++)
    ord[i] = i;
  shuffle_int(ord, ndset);
  for (int i = 0; i < ndset; i++) {
    uint64_t t[MAXN];
    for (int j = 0; j < n; j++)
      t[j] = dset[(size_t)ord[i] * n + n - 1 - j];
    vec_list_push(&dl, t);
  }
  free(ord);
  brute_pairs();
  int fails = 0;
  dfirst_t *df = dfirst_new(l, 0, l->count, &dl, 0, dl.count);
  dfirst_set_top_root(df, top_root);
  d_found.k = 0;
  d_invalid = 0;
  dfirst_stats_t st = dfirst_search(df, 0, dl.count, 1, 0, o, d_cb, NULL, NULL);
  bool same = same_dpairs() && !d_invalid && !st.truncated &&
              st.pairs == (uint64_t)d_brute.k;
  /* the planted square, twice, flagged */
  int seen = 0;
  if (planted >= 0)
    for (int i = 0; i < d_found.k; i++)
      if (!cmp_canon(&d_found.a[i].sq, &oracle_out.a[planted])) {
        for (int t = 0; t < 2; t++) {
          uint64_t x[MAXN];
          memcpy(x, pt[t], n * 8);
          qsort(x, n, 8, cmp_u64);
          if (!memcmp(x, d_found.a[i].d, n * 8) && d_found.a[i].partner)
            seen |= 1 << t;
        }
      }
  const bool planted_ok = planted < 0 || seen == 3;
  /* d-stride 3, offsets 0, 1, 2 */
  bool split_ok = true;
  if (seed % 4 == 0) {
    d_found.k = 0;
    for (int off = 0; off < 3; off++)
      dfirst_search(df, 0, dl.count, 3, off, o, d_cb, NULL, NULL);
    split_ok = same_dpairs() && !d_invalid;
  }
  dfirst_free(df);
  vec_list_free(&dl);
  if (!same || !planted_ok || !split_ok) {
    fails++;
    printf("DFIRST MISMATCH seed %lu n %d nv %d diags %d variant %s top_root "
           "%d: found %d brute %d invalid %d planted %s split %s\n",
           (unsigned long)seed, n, nv, ndset, vname, top_root, d_found.k,
           d_brute.k,
           d_invalid, planted_ok ? "ok" : "MISSING", split_ok ? "ok" : "BAD");
    char path[256];
    snprintf(path, sizeof(path), "fail_%lu.txt", (unsigned long)seed);
    dump(path, seed);
  }
  d_pairs_total += d_brute.k;
  for (int i = 0; i < d_brute.k; i++)
    d_partner_total += d_brute.a[i].partner;
  d_planted += planted >= 0;
  if (verbose)
    printf("seed %lu n %d nv %d diags %d variant %s: %d pairs (%lu nodes)\n",
           (unsigned long)seed, n, nv, ndset, vname, d_brute.k,
           (unsigned long)st.nodes);
  return fails;
}

int main(int argc, char **argv) {
  int ai = 1;
  setvbuf(stdout, NULL, _IOLBF, 0);
  int cross = 0, force_mode = -1, force_n = 0, diff = 0, noracle = 0, diffonly = 0;
  int dfirst = 0;
  uint64_t budget = 20000000;
  const char *varsel = NULL;
  while (ai < argc && argv[ai][0] == '-') {
    if (!strcmp(argv[ai], "-v"))
      verbose = 1;
    else if (!strcmp(argv[ai], "--diff"))
      diff = 1;
    else if (!strcmp(argv[ai], "--noracle"))
      noracle = diff = 1;
    else if (!strcmp(argv[ai], "--cross"))
      cross = 1;
    else if (!strcmp(argv[ai], "--dfirst"))
      dfirst = 1;
    else if (!strcmp(argv[ai], "--mode"))
      force_mode = atoi(argv[++ai]);
    else if (!strcmp(argv[ai], "--n"))
      force_n = atoi(argv[++ai]);
    else if (!strcmp(argv[ai], "--budget"))
      budget = strtoull(argv[++ai], NULL, 10);
    else if (!strcmp(argv[ai], "--variants"))
      varsel = argv[++ai];
    ai++;
  }
  uint64_t seed0 = strtoull(argv[ai], NULL, 10);
  uint64_t cnt = strtoull(argv[ai + 1], NULL, 10);
  int fails = 0, skipped = 0, ran = 0;
  uint64_t total_sq = 0, inst_with_sq = 0;
  int Lhist[5] = {0};
  for (uint64_t seed = seed0; seed < seed0 + cnt; seed++) {
    int mode = gen(seed, force_mode, force_n);
    if (nv < 2 * n) {
      skipped++;
      continue;
    }
    double t0 = now();
    build_inset();
    int have_oracle = 1;
    if (noracle || oracle(budget)) {
      if (verbose)
        printf("seed %lu oracle budget exceeded or disabled (%.2fs)\n", (unsigned long)seed, now() - t0);
      if (!diff) {
        skipped++;
        continue;
      }
      have_oracle = 0;
      diffonly++;
    }
    if (have_oracle && cross && !oracle2()) {
      crossed++;
      bool same = oracle2_out.k == oracle_out.k;
      for (int i = 0; i < oracle_out.k && same; i++)
        same = cmp_canon(&oracle_out.a[i], &oracle2_out.a[i]) == 0;
      if (!same) {
        printf("ORACLE MISMATCH seed %lu: o1 %d o2 %d\n", (unsigned long)seed,
               oracle_out.k, oracle2_out.k);
        fails++;
      }
    }
    ran++;
    vec_list_t l;
    vec_list_init(&l, n);
    for (int v = 0; v < nv; v++) {
      uint64_t e[MAXN];
      for (int i = 0; i < n; i++)
        e[i] = uval[vec[v][i]];
      qsort(e, n, 8, cmp_u64);
      for (int i = 0; i < n / 2; i++) { /* descending */
        uint64_t t = e[i];
        e[i] = e[n - 1 - i];
        e[n - 1 - i] = t;
      }
      vec_list_push(&l, e);
    }
    if (dfirst) {
      if (have_oracle != 1) {
        skipped++;
        ran--;
        vec_list_free(&l);
        continue;
      }
      /* the default options, and one other option set in turn */
      search_opts_t o;
      search_opts_default(&o);
      fails += dfirst_check(seed, &l, &o, "default", 1);
      const variant_t *V = &variants[seed % NVAR];
      o.forward_check = V->fc;
      o.mrv = V->mrv;
      o.support = V->support;
      o.min_words = V->min_words;
      o.gather = V->gather;
      o.cross = V->cross;
      fails += dfirst_check(seed, &l, &o, V->name, (int)(seed / NVAR % 2));
      total_sq += oracle_out.k;
      inst_with_sq += oracle_out.k > 0;
      vec_list_free(&l);
      continue;
    }
    /* mode 6 (where the oracle is too slow): every variant must find the
     * planted square, besides agreeing with the others */
    canon_t planted;
    if (mode == 6) {
      uint64_t rows[MAXN][MAXN], cols[MAXN][MAXN];
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          rows[i][j] = uval[sat_rows[i][j]];
          cols[i][j] = uval[sat_cols[i][j]];
        }
      canonicalize(rows, cols, &planted);
    }
    int labels = -1;
    for (int vi = 0; vi < NVAR; vi++) {
      const variant_t *V = &variants[vi];
      if (varsel && !strstr(varsel, V->name))
        continue;
      search_opts_t o;
      search_opts_default(&o);
      o.forward_check = V->fc;
      o.mrv = V->mrv;
      o.support = V->support;
      o.min_words = V->min_words;
      o.gather = V->gather;
      o.cross = V->cross;
      o.r1_width = !V->fixw;
      if (V->pretest)
        o.pretest_min = V->pretest > 0 ? V->pretest : 0;
      found.k = 0;
      invalid = 0;
      double t1 = now();
      search_stats_t st = search_vectors(&l, 0, l.count, &o, cb, NULL);
      if (verbose)
        printf("seed %lu mode %d n %d nv %d L %d variant %s: %d squares, %lu nodes, %.3fs (oracle %.3fs)\n",
               (unsigned long)seed, mode, n, nv, st.num_labels, V->name, found.k,
               (unsigned long)st.nodes, now() - t1, t1 - t0);
      labels = st.num_labels;
      canon_len = 2 * n * n;
      if (found.k)
        qsort(found.a, found.k, sizeof(canon_t), cmp_canon);
      if (!have_oracle) {
        /* the first variant's squares (deduplicated) become the reference */
        oracle_out.k = 0;
        for (int i = 0; i < found.k; i++)
          if (i == 0 || cmp_canon(&found.a[i - 1], &found.a[i]) != 0)
            cl_push(&oracle_out, &found.a[i]);
        have_oracle = 2;
      }
      bool same = found.k == oracle_out.k && !st.truncated && !invalid;
      for (int i = 0; i < found.k && same; i++)
        same = cmp_canon(&found.a[i], &oracle_out.a[i]) == 0;
      const bool lost =
          mode == 6 && (!found.k || !bsearch(&planted, found.a, found.k,
                                             sizeof(canon_t), cmp_canon));
      /* the pretest must not change the nodes either: the same search with
       * it off (no callback: only the counts are compared) */
      if (V->pretest > 0 && !st.truncated) {
        search_opts_t o0 = o;
        o0.pretest_min = 0;
        const search_stats_t st0 =
            search_vectors(&l, 0, l.count, &o0, NULL, NULL);
        if (st0.nodes != st.nodes || st0.squares != st.squares) {
          fails++;
          printf("NODES MISMATCH seed %lu variant %s: %lu nodes, %lu "
                 "without the pretest\n",
                 (unsigned long)seed, V->name, (unsigned long)st.nodes,
                 (unsigned long)st0.nodes);
        }
      }
      if (!same || lost) {
        fails++;
        printf("MISMATCH seed %lu mode %d n %d nv %d labels %d variant %s: "
               "search %d oracle %d truncated %d invalid %d ref %s%s\n",
               (unsigned long)seed, mode, n, nv, st.num_labels, V->name,
               found.k, oracle_out.k, st.truncated, invalid,
               have_oracle == 2 ? "diff" : "oracle",
               lost ? " (planted square missing)" : "");
        char path[256];
        snprintf(path, sizeof(path), "fail_%lu.txt", (unsigned long)seed);
        dump(path, seed);
        fflush(stdout);
      }
    }
    /* (the oracle's squares, or the reference variant's with --diff) */
    total_sq += oracle_out.k;
    inst_with_sq += oracle_out.k > 0;
    int b = labels <= 128 ? 0 : labels <= 192 ? 1 : labels <= 256 ? 2 : labels <= 512 ? 3 : 4;
    Lhist[b]++;
    vec_list_free(&l);
  }
  if (cross)
    printf("cross-checked %d\n", crossed);
  if (dfirst)
    printf("dfirst: (square, d) pairs %lu, flagged %lu, planted pairs %lu\n",
           (unsigned long)d_pairs_total, (unsigned long)d_partner_total,
           (unsigned long)d_planted);
  printf("diff-only %d\n", diffonly);
  printf("ran %d skipped %d fails %d squares %lu inst_with_squares %lu "
         "W2 %d W3 %d W4 %d W8 %d W16+ %d\n",
         ran, skipped, fails, (unsigned long)total_sq,
         (unsigned long)inst_with_sq, Lhist[0], Lhist[1], Lhist[2], Lhist[3],
         Lhist[4]);
  return fails != 0;
}
