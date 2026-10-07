/*
 * Differential test of the arrangement search against an independent
 * brute-force oracle.
 *
 * usage: fuzz_arrange [-v] [--variants a,b,..] [--mode M] [--n N]
 *                     [--budget NODES] [--diff | --noracle] [--cross]
 *                     seed0 count
 *
 * For every seed in [seed0, seed0 + count): generate a synthetic family of
 * n-sets (planted squares with Latin-square alternates, random transversals,
 * mutants and noise; random dense families; overlapping grids; wide label
 * universes up to W = 16; large sparse families; hubs with more than 255
 * vectors through one number), find all semi-magic squares with a simple
 * backtracking oracle (r1 -> one col per cell of r1 -> exact cover by rows),
 * canonicalize them (up to transposition and row/col order), and compare
 * with what search_vectors reports for each option set in variants[] (and
 * check that every reported square is valid). --cross also checks the oracle
 * against a second, clique-based one; --diff compares the variants with
 * each other where the oracle exceeds its budget. Exit status is nonzero on
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

/* modes: 0 = one dense grid + transversals, 1 = many planted grids (wide),
 * 2 = random dense, 3 = a few grids overlapping in a small universe */
static int gen(uint64_t seed, int force_mode, int force_n) {
  rs = seed * 0x2545F4914F6CDD1Dull + 12345;
  nv = 0;
  n = force_n ? force_n : rint_(4, 7);
  if (!force_n && rnd() % 10 == 0)
    n = rnd() % 2 ? 3 : 8;
  int mode = force_mode >= 0 ? force_mode : rint_(0, 3);
  plant_nt_max = mode == 1 || mode == 4 ? n : 3 * n;
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
  } else if (mode == 1) {
    /* target label ranges: W=2, 3, 4, 8, 16 */
    static const int lo_[] = {40, 129, 193, 257, 513};
    static const int hi_[] = {128, 192, 256, 512, 1100};
    int r = rint_(0, 4);
    U = rint_(lo_[r], hi_[r]);
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

int main(int argc, char **argv) {
  int ai = 1;
  setvbuf(stdout, NULL, _IOLBF, 0);
  int cross = 0, force_mode = -1, force_n = 0, diff = 0, noracle = 0, diffonly = 0;
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
    total_sq += oracle_out.k;
    inst_with_sq += oracle_out.k > 0;
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
      if (!same) {
        fails++;
        printf("MISMATCH seed %lu mode %d n %d nv %d labels %d variant %s: "
               "search %d oracle %d truncated %d invalid %d ref %s\n",
               (unsigned long)seed, mode, n, nv, st.num_labels, V->name,
               found.k, oracle_out.k, st.truncated, invalid, have_oracle == 2 ? "diff" : "oracle");
        char path[256];
        snprintf(path, sizeof(path), "fail_%lu.txt", (unsigned long)seed);
        dump(path, seed);
        fflush(stdout);
      }
    }
    int b = labels <= 128 ? 0 : labels <= 192 ? 1 : labels <= 256 ? 2 : labels <= 512 ? 3 : 4;
    Lhist[b]++;
    vec_list_free(&l);
  }
  if (cross)
    printf("cross-checked %d\n", crossed);
  printf("diff-only %d\n", diffonly);
  printf("ran %d skipped %d fails %d squares %lu inst_with_squares %lu "
         "W2 %d W3 %d W4 %d W8 %d W16+ %d\n",
         ran, skipped, fails, (unsigned long)total_sq,
         (unsigned long)inst_with_sq, Lhist[0], Lhist[1], Lhist[2], Lhist[3],
         Lhist[4]);
  return fails != 0;
}
