#include "square.h"

#include <stdlib.h>
#include <string.h>

void square_to_grid(const square_t *sq, grid_t *grid) {
  int n = sq->n;
  grid->n = n;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      grid->g[i][j] = 0;
      for (int a = 0; a < n; a++)
        for (int b = 0; b < n; b++)
          if (sq->rows[i][a] == sq->cols[j][b])
            grid->g[i][j] = sq->rows[i][a];
    }
}

static int u64_cmp(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}

static int line_n;
static int line_cmp(const void *a, const void *b) {
  const uint64_t *x = a, *y = b;
  for (int i = 0; i < line_n; i++)
    if (x[i] != y[i])
      return (x[i] > y[i]) - (x[i] < y[i]);
  return 0;
}

void square_canonical(const square_t *sq, uint64_t *out) {
  int n = sq->n;
  uint64_t r[SQ_MAX_N * SQ_MAX_N], c[SQ_MAX_N * SQ_MAX_N];
  for (int i = 0; i < n; i++) {
    memcpy(r + i * n, sq->rows[i], n * sizeof(uint64_t));
    memcpy(c + i * n, sq->cols[i], n * sizeof(uint64_t));
    qsort(r + i * n, n, sizeof(uint64_t), u64_cmp);
    qsort(c + i * n, n, sizeof(uint64_t), u64_cmp);
  }
  line_n = n;
  qsort(r, n, n * sizeof(uint64_t), line_cmp);
  qsort(c, n, n * sizeof(uint64_t), line_cmp);
  int rc = 0;
  for (int i = 0; i < n * n && !rc; i++)
    rc = (r[i] > c[i]) - (r[i] < c[i]);
  const uint64_t *first = rc <= 0 ? r : c, *second = rc <= 0 ? c : r;
  memcpy(out, first, n * n * sizeof(uint64_t));
  memcpy(out + n * n, second, n * n * sizeof(uint64_t));
}

uint64_t square_hash(const square_t *sq) {
  uint64_t canon[2 * SQ_MAX_N * SQ_MAX_N];
  square_canonical(sq, canon);
  uint64_t h = 1469598103934665603ULL; /* FNV-1a */
  for (int i = 0; i < 2 * sq->n * sq->n; i++) {
    uint64_t x = canon[i];
    for (int b = 0; b < 8; b++) {
      h ^= (x >> (8 * b)) & 0xff;
      h *= 1099511628211ULL;
    }
  }
  return h;
}

/* ---------------------------------------------------------------------- */
/* diagonals                                                               */

/* all involutions of {0..n-1} with no fixed points (n even) or exactly one
 * fixed point (n odd): these are the possible relations between the two
 * diagonals after permuting rows and cols */
static void gen_partners_rec(int n, int *cur, int fixed_left,
                             int (*out)[SQ_MAX_N], int *count) {
  int a = 0;
  while (a < n && cur[a] >= 0)
    a++;
  if (a == n) {
    memcpy(out[(*count)++], cur, SQ_MAX_N * sizeof(int));
    return;
  }
  if (fixed_left) {
    cur[a] = a;
    gen_partners_rec(n, cur, 0, out, count);
    cur[a] = -1;
  }
  for (int b = a + 1; b < n; b++) {
    if (cur[b] >= 0)
      continue;
    cur[a] = b;
    cur[b] = a;
    gen_partners_rec(n, cur, fixed_left, out, count);
    cur[a] = cur[b] = -1;
  }
}

static int gen_partners(int n, int (*out)[SQ_MAX_N]) {
  int cur[SQ_MAX_N];
  for (int i = 0; i < SQ_MAX_N; i++)
    cur[i] = -1;
  int count = 0;
  gen_partners_rec(n, cur, n % 2, out, &count);
  return count;
}

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

static int perm_index(const int *p, int n) {
  /* lexicographic rank */
  int idx = 0;
  for (int i = 0; i < n; i++) {
    int smaller = 0;
    for (int j = i + 1; j < n; j++)
      smaller += p[j] < p[i];
    idx = idx * (n - i) + smaller;
  }
  return idx;
}

void square_diag_stats(const square_t *sq, diag_stats_t *ds) {
  int n = sq->n;
  grid_t grid;
  square_to_grid(sq, &grid);
  memset(ds, 0, sizeof(*ds));
  uint64_t S = 0;
  __uint128_t P = 1;
  for (int j = 0; j < n; j++) {
    S += sq->rows[0][j];
    P *= sq->rows[0][j];
  }
  ds->sum = S;
  ds->product = (uint64_t)P;

  int nperm = 1;
  for (int i = 2; i <= n; i++)
    nperm *= i;
  uint8_t *score = malloc(nperm);
  int perm[SQ_MAX_N];
  for (int i = 0; i < n; i++)
    perm[i] = i;
  int idx = 0;
  do {
    uint64_t s = 0;
    __uint128_t p = 1;
    int p_ok = 1;
    for (int i = 0; i < n; i++) {
      uint64_t x = grid.g[i][perm[i]];
      s += x;
      if (p_ok) {
        p *= x;
        if (p > P)
          p_ok = 0;
      }
    }
    int is_s = s == S, is_p = p_ok && p == P;
    score[idx] = is_s && is_p ? DIAG_SP : is_p ? DIAG_P : is_s ? DIAG_S : 0;
    ds->s_count += is_s;
    ds->p_count += is_p;
    ds->sp_count += is_s && is_p;
    idx++;
  } while (next_perm(perm, n));

  /* best pair of diagonals */
  int partners[128][SQ_MAX_N];
  int npart = gen_partners(n, partners);
  int best = -1;
  int best_p1[SQ_MAX_N], best_tau = 0;
  for (int i = 0; i < n; i++)
    perm[i] = i;
  idx = 0;
  do {
    if (score[idx] || best < 0) {
      for (int t = 0; t < npart; t++) {
        int p2[SQ_MAX_N];
        for (int i = 0; i < n; i++)
          p2[i] = perm[partners[t][i]];
        int sc = score[idx] + score[perm_index(p2, n)];
        if (sc > best) {
          best = sc;
          best_tau = t;
          memcpy(best_p1, perm, sizeof(perm));
        }
      }
    }
    idx++;
  } while (next_perm(perm, n));
  ds->best_score = best;
  free(score);

  /* reorder: rows by pi, cols by rho, so that cell (r, best_p1[r]) is on the
   * main diagonal and (r, best_p1[tau[r]]) on the anti-diagonal */
  const int *tau = partners[best_tau];
  int pi[SQ_MAX_N], rho[SQ_MAX_N];
  int k = 0;
  for (int a = 0; a < n; a++) {
    int b = tau[a];
    if (b == a)
      pi[a] = n / 2;
    else if (b > a) {
      pi[a] = k;
      pi[b] = n - 1 - k;
      k++;
    }
  }
  for (int r = 0; r < n; r++)
    rho[best_p1[r]] = pi[r];
  ds->best.n = n;
  for (int r = 0; r < n; r++)
    for (int c = 0; c < n; c++)
      ds->best.g[pi[r]][rho[c]] = grid.g[r][c];
}

void grid_print_json(FILE *fp, const grid_t *g) {
  fputc('[', fp);
  for (int i = 0; i < g->n; i++) {
    fputs(i ? ",[" : "[", fp);
    for (int j = 0; j < g->n; j++)
      fprintf(fp, j ? ",%lu" : "%lu", (unsigned long)g->g[i][j]);
    fputc(']', fp);
  }
  fputc(']', fp);
}
