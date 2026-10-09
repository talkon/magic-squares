/*
 * The class support (search_opts_t class_support) through the API with the
 * plain root (top_root_only = 0), which the V_d searches of dfirst.c never
 * use: there the r1 without the top label are roots too, their largest
 * label x is still a class label, and the per-r1 widths run them at
 * ceil((x + 1) / 64) words, which is one word fewer than the class labels
 * when these straddle a word boundary (e.g. 130 labels: the classes are
 * labels 124..129, and the r1 with x = 124..127 run at 2 words).
 *
 * usage: class_support_test [seed0 count]
 *
 * Each family plants 1-3 squares of V_d (d = 6 numbers, each row and col
 * holding exactly one of them) among random vectors with exactly one of
 * d's numbers, over U = 128-260 numbers. Every combination of
 * class_support, top_root_only and r1_width must find the same squares,
 * and r1_width must not change the nodes. Built with ARRANGE_DEBUG, which
 * checks every list's class blocks. Exit status is nonzero on a mismatch.
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "arrange.h"
#include "enumerate.h"

static uint64_t rs;
static uint64_t rnd(void) {
  rs ^= rs << 13;
  rs ^= rs >> 7;
  rs ^= rs << 17;
  return rs;
}
static int rint_(int lo, int hi) { /* [lo, hi] */
  return lo + (int)(rnd() % (uint64_t)(hi - lo + 1));
}
static void shuffle(int *a, int k) {
  for (int i = k - 1; i > 0; i--) {
    const int j = rint_(0, i), t = a[i];
    a[i] = a[j];
    a[j] = t;
  }
}
static int cmp_desc(const void *a, const void *b) {
  const uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return x < y ? 1 : x > y ? -1 : 0;
}
static uint64_t mix(uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdull;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ull;
  x ^= x >> 33;
  return x;
}

/* the squares found: their number and an order-free hash (of the rows and
 * the cols, each as a set of sets: invariant under the row and col order
 * and under transposition, as the root decides which axis is rows) */
typedef struct {
  uint64_t k, h;
} found_t;
static uint64_t hash_lines(const uint64_t (*m)[SQ_MAX_N], int n) {
  uint64_t h = 0;
  for (int i = 0; i < n; i++) {
    uint64_t r = 0;
    for (int j = 0; j < n; j++)
      r += mix(m[i][j]);
    h += mix(r);
  }
  return mix(h ^ 0x9e3779b97f4a7c15ull);
}
static int cb(const square_t *sq, void *ctx) {
  found_t *f = ctx;
  f->k++;
  f->h += hash_lines(sq->rows, sq->n) + hash_lines(sq->cols, sq->n);
  return 0;
}

#define N_ 6
#define UMAX 260

int main(int argc, char **argv) {
  const uint64_t seed0 = argc > 2 ? strtoull(argv[1], NULL, 10) : 1;
  const int count = argc > 2 ? atoi(argv[2]) : 200;
  /* class labels straddling a word boundary, or not */
  static const int Us[] = {130, 131, 133, 129, 194, 196, 258, 260, 128, 200};
  const int nU = (int)(sizeof(Us) / sizeof(Us[0]));
  int fails = 0;
  for (int t = 0; t < count; t++) {
    const uint64_t seed = seed0 + (uint64_t)t;
    rs = 0x2545f4914f6cdd1dull ^ mix(seed);
    for (int i = 0; i < 4; i++)
      rnd();
    const int n = N_, U = Us[t % nU];
    vec_list_t l;
    vec_list_init(&l, n);
    int used[UMAX] = {0};
    /* planted squares: d_i at (i, perm[i]), the other cells distinct
     * numbers outside d */
    const int nsq = rint_(1, 3);
    for (int q = 0; q < nsq; q++) {
      int perm[N_], pool[UMAX], g[N_][N_], k = 0;
      for (int i = 0; i < n; i++)
        perm[i] = i;
      shuffle(perm, n);
      for (int x = n; x < U; x++)
        pool[x - n] = x;
      shuffle(pool, U - n);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          g[i][j] = j == perm[i] ? i : pool[k++];
      for (int i = 0; i < n; i++) {
        uint64_t r[N_], c[N_];
        for (int j = 0; j < n; j++) {
          r[j] = 1000 + (uint64_t)g[i][j];
          c[j] = 1000 + (uint64_t)g[j][i];
          used[g[i][j]] = 1;
        }
        qsort(r, n, 8, cmp_desc);
        qsort(c, n, 8, cmp_desc);
        vec_list_push(&l, r);
        vec_list_push(&l, c);
      }
    }
    /* noise, until every number is used (so that there are U labels) */
    const int extra = rint_(150, 450);
    for (int e = 0;; e++) {
      int all = 1;
      for (int x = 0; x < U; x++)
        all &= used[x];
      if (all && e >= extra)
        break;
      uint64_t v[N_];
      int m = 0;
      const int c = rint_(0, n - 1);
      v[m++] = 1000 + (uint64_t)c;
      used[c] = 1;
      while (m < n) {
        const int x = rint_(n, U - 1);
        int dup = 0;
        for (int j = 1; j < m; j++)
          dup |= v[j] == 1000 + (uint64_t)x;
        if (!dup) {
          v[m++] = 1000 + (uint64_t)x;
          used[x] = 1;
        }
      }
      qsort(v, n, 8, cmp_desc);
      vec_list_push(&l, v);
    }
    uint64_t top[N_];
    for (int i = 0; i < n; i++)
      top[i] = 1000 + (uint64_t)i;
    found_t ref = {0, 0};
    uint64_t nodes[4] = {0};
    /* cfg bit 0: class support, bit 1: top root only, bit 2: no r1 widths */
    for (int cfg = 0; cfg < 8; cfg++) {
      search_opts_t o;
      search_opts_default(&o);
      o.top_numbers = top;
      o.n_top = n;
      o.class_support = cfg & 1;
      o.top_root_only = cfg >> 1 & 1;
      o.r1_width = !(cfg >> 2 & 1);
      found_t f = {0, 0};
      const search_stats_t st = search_vectors(&l, 0, l.count, &o, cb, &f);
      if (cfg == 0)
        ref = f;
      if (f.k != ref.k || f.h != ref.h || f.k < (uint64_t)nsq) {
        printf("MISMATCH seed %lu U %d class_support %d top_root_only %d "
               "r1_width %d: %lu squares, %lu without class support "
               "(%d planted)\n",
               (unsigned long)seed, U, o.class_support, o.top_root_only,
               o.r1_width, (unsigned long)f.k, (unsigned long)ref.k, nsq);
        fails++;
      }
      if (cfg < 4) {
        nodes[cfg] = st.nodes;
      } else if (st.nodes != nodes[cfg & 3]) {
        printf("MISMATCH seed %lu U %d class_support %d top_root_only %d: "
               "%lu nodes with the per-r1 widths, %lu without\n",
               (unsigned long)seed, U, o.class_support, o.top_root_only,
               (unsigned long)nodes[cfg & 3], (unsigned long)st.nodes);
        fails++;
      }
    }
    vec_list_free(&l);
  }
  printf("%d families, %d mismatches\n", count, fails);
  return fails != 0;
}
