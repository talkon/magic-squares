#include "enumerate.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const uint64_t ENUM_PRIMES[ENUM_MAX_PRIMES] = {2,  3,  5,  7,  11,
                                               13, 17, 19, 23, 29};

/*
 * Divisors are represented by their value plus a packed exponent vector, with
 * 6 bits per prime: 5 value bits and a guard bit. This allows checking
 * divisibility and dividing with one subtraction (SWAR).
 */
#define FIELD_BITS 6
#define MAX_EXP 31

static uint64_t guard_mask(void) {
  uint64_t g = 0;
  for (int i = 0; i < ENUM_MAX_PRIMES; i++)
    g |= (uint64_t)1 << (FIELD_BITS * i + FIELD_BITS - 1);
  return g;
}

uint64_t pexp_value(const prime_exps_t *p) {
  uint64_t v = 1;
  for (int i = 0; i < p->num_primes; i++)
    for (int j = 0; j < p->exps[i]; j++) {
      if (v > UINT64_MAX / ENUM_PRIMES[i])
        return 0;
      v *= ENUM_PRIMES[i];
    }
  return v;
}

uint64_t pexp_num_divisors(const prime_exps_t *p) {
  uint64_t t = 1;
  for (int i = 0; i < p->num_primes; i++)
    t *= (uint64_t)(p->exps[i] + 1);
  return t;
}

int pexp_parse(prime_exps_t *p, int argc, char **argv) {
  memset(p, 0, sizeof(*p));
  if (argc > ENUM_MAX_PRIMES)
    return -1;
  for (int i = 0; i < argc; i++) {
    char *end;
    long e = strtol(argv[i], &end, 10);
    if (*end != '\0' || e < 0 || e > MAX_EXP)
      return -1;
    p->exps[i] = (int)e;
  }
  p->num_primes = argc;
  /* strip trailing zeros so "10 4 3 2 0" == "10 4 3 2" */
  while (p->num_primes > 0 && p->exps[p->num_primes - 1] == 0)
    p->num_primes--;
  return 0;
}

void pexp_to_str(const prime_exps_t *p, char *buf, const char *sep) {
  buf[0] = '\0';
  for (int i = 0; i < p->num_primes; i++) {
    char tmp[16];
    snprintf(tmp, sizeof(tmp), "%s%d", i ? sep : "", p->exps[i]);
    strcat(buf, tmp);
  }
}

void vec_list_init(vec_list_t *l, int n) {
  l->elts = NULL;
  l->count = 0;
  l->cap = 0;
  l->n = n;
}

void vec_list_free(vec_list_t *l) {
  free(l->elts);
  l->elts = NULL;
  l->count = l->cap = 0;
}

void vec_list_push(vec_list_t *l, const uint64_t *vec) {
  if (l->count == l->cap) {
    l->cap = l->cap ? 2 * l->cap : 1024;
    l->elts = realloc(l->elts, l->cap * (size_t)l->n * sizeof(uint64_t));
    if (!l->elts) {
      fprintf(stderr, "out of memory in vec_list_push\n");
      exit(1);
    }
  }
  memcpy(l->elts + l->count * (size_t)l->n, vec, (size_t)l->n * sizeof(uint64_t));
  l->count++;
}

/* ---------------------------------------------------------------------- */
/* divisor table                                                           */

typedef struct {
  uint64_t val;
  uint64_t packed;
} divisor_t;

static int divisor_cmp(const void *a, const void *b) {
  uint64_t x = ((const divisor_t *)a)->val, y = ((const divisor_t *)b)->val;
  return (x > y) - (x < y);
}

static divisor_t *make_divisors(const prime_exps_t *p, size_t *count) {
  size_t total = pexp_num_divisors(p);
  divisor_t *divs = malloc(total * sizeof(divisor_t));
  divs[0].val = 1;
  divs[0].packed = 0;
  size_t len = 1;
  for (int i = 0; i < p->num_primes; i++) {
    size_t prev = len;
    for (size_t j = 0; j < prev; j++) {
      uint64_t v = divs[j].val, pk = divs[j].packed;
      for (int e = 1; e <= p->exps[i]; e++) {
        v *= ENUM_PRIMES[i];
        pk += (uint64_t)1 << (FIELD_BITS * i);
        divs[len].val = v;
        divs[len].packed = pk;
        len++;
      }
    }
  }
  qsort(divs, len, sizeof(divisor_t), divisor_cmp);
  *count = len;
  return divs;
}

/* first index i with divs[i].val >= x */
static size_t lower_bound(const divisor_t *divs, size_t count, uint64_t x) {
  size_t lo = 0, hi = count;
  while (lo < hi) {
    size_t mid = (lo + hi) / 2;
    if (divs[mid].val < x)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}

/* ---------------------------------------------------------------------- */
/* depth-first enumeration                                                 */

typedef struct {
  const divisor_t *divs;
  size_t num_divs;
  uint64_t guard;
  uint64_t sum_min, sum_max;
  int n;
  uint64_t cur[ENUM_MAX_N];
  vec_list_t *out;
  uint64_t nodes;
} enum_state_t;

/* largest integer r with r^k <= x (k >= 1) */
static uint64_t iroot_floor(uint64_t x, int k) {
  if (k == 1)
    return x;
  uint64_t r = (uint64_t)floor(pow((double)x, 1.0 / k));
  /* fix floating point error */
  for (;;) {
    /* check (r+1)^k <= x */
    __uint128_t t = 1;
    bool over = false;
    for (int i = 0; i < k; i++) {
      t *= (r + 1);
      if (t > x) {
        over = true;
        break;
      }
    }
    if (over)
      break;
    r++;
  }
  for (;;) {
    __uint128_t t = 1;
    bool over = false;
    for (int i = 0; i < k; i++) {
      t *= r;
      if (t > x) {
        over = true;
        break;
      }
    }
    if (!over || r == 0)
      break;
    r--;
  }
  return r;
}

/*
 * Upper bound on the sum of k distinct integers < upper with product rem.
 * Relaxing to reals x_i in [1, upper - i], the sum is a convex function of
 * the logs, so it is maximized at a vertex: the largest caps taken in full,
 * one partial element, and the rest equal to 1.
 */
static uint64_t max_completion(uint64_t rem, int k, uint64_t upper) {
  uint64_t total = 0;
  double r = (double)rem;
  for (int i = 1; i <= k; i++) {
    double cap = upper == UINT64_MAX ? 1e300 : (double)(upper - i);
    if (cap < 1)
      return total; /* infeasible anyway */
    if (r <= cap) {
      double t = (double)total + r + (double)(k - i);
      return t > 1.8e19 ? UINT64_MAX : (uint64_t)t + 1;
    }
    total += (uint64_t)cap;
    r /= cap;
  }
  return total + 1;
}

/*
 * choose the next (largest remaining) element. k elements remain, with
 * product rem (packed exponents rem_pk), each < upper, current sum `sum`.
 */
static void enum_rec(enum_state_t *st, int k, uint64_t rem, uint64_t rem_pk,
                     uint64_t upper, uint64_t sum) {
  st->nodes++;
  int pos = st->n - k;
  if (k == 1) {
    if (rem < upper && sum + rem >= st->sum_min && sum + rem <= st->sum_max) {
      st->cur[pos] = rem;
      vec_list_push(st->out, st->cur);
    }
    return;
  }
  if (sum >= st->sum_max)
    return;
  uint64_t budget = st->sum_max - sum;
  if (st->sum_min > sum && max_completion(rem, k, upper) < st->sum_min - sum)
    return;

  /* largest element d satisfies d^k > rem (elements are distinct) */
  uint64_t lo = iroot_floor(rem, k) + 1;
  size_t i = lower_bound(st->divs, st->num_divs, lo);

  if (k == 2) {
    /* d > e, d * e = rem, d + e in [sum_min - sum, sum_max - sum];
     * d + rem / d is increasing for d > sqrt(rem) */
    uint64_t need = st->sum_min > sum ? st->sum_min - sum : 0;
    for (; i < st->num_divs; i++) {
      uint64_t d = st->divs[i].val;
      if (d >= upper || d > rem)
        break;
      uint64_t q = st->guard | rem_pk;
      uint64_t diff = q - st->divs[i].packed;
      if ((diff & st->guard) != st->guard)
        continue;
      uint64_t e = rem / d;
      uint64_t s2 = d + e;
      if (s2 > budget)
        break;
      if (s2 < need)
        continue;
      st->cur[pos] = d;
      st->cur[pos + 1] = e;
      vec_list_push(st->out, st->cur);
    }
    return;
  }

  double inv = 1.0 / (k - 1);
  for (; i < st->num_divs; i++) {
    uint64_t d = st->divs[i].val;
    if (d >= upper || d > budget)
      break;
    uint64_t q = st->guard | rem_pk;
    uint64_t diff = q - st->divs[i].packed;
    if ((diff & st->guard) != st->guard)
      continue;
    uint64_t r2 = rem / d;
    /* remaining k-1 elements have product r2, so (AM-GM) their sum is at
     * least (k-1) * r2^(1/(k-1)); this lower bound is increasing in d for
     * d^k >= rem, so we can break */
    double lb = (double)d + (k - 1) * pow((double)r2, inv) * (1.0 - 1e-12);
    if (lb > (double)budget)
      break;
    st->cur[pos] = d;
    enum_rec(st, k - 1, r2, diff & ~st->guard, d, sum + d);
  }
}

static int n_global_for_cmp;

static int vec_cmp_by_sum(const void *a, const void *b) {
  const uint64_t *x = a, *y = b;
  uint64_t sx = 0, sy = 0;
  for (int i = 0; i < n_global_for_cmp; i++) {
    sx += x[i];
    sy += y[i];
  }
  if (sx != sy)
    return (sx > sy) - (sx < sy);
  /* same order as enumeration.cpp: lexicographic on the elements in
   * ascending order (vectors are stored descending) */
  for (int i = n_global_for_cmp - 1; i >= 0; i--)
    if (x[i] != y[i])
      return (x[i] > y[i]) - (x[i] < y[i]);
  return 0;
}

/* sort vectors [from, count) by sum, then lexicographically descending */
static void sort_by_sum(vec_list_t *l, size_t from, uint64_t sum_min,
                        uint64_t sum_max) {
  int n = l->n;
  size_t m = l->count - from;
  if (m < 2)
    return;
  n_global_for_cmp = n;
  uint64_t *base = l->elts + from * (size_t)n;
  if (sum_max - sum_min > (uint64_t)(1 << 22)) {
    qsort(base, m, (size_t)n * sizeof(uint64_t), vec_cmp_by_sum);
    return;
  }
  size_t nb = sum_max - sum_min + 2;
  size_t *start = calloc(nb, sizeof(size_t));
  for (size_t i = 0; i < m; i++)
    start[vec_list_sum(l, from + i) - sum_min + 1]++;
  for (size_t b = 1; b < nb; b++)
    start[b] += start[b - 1];
  uint64_t *tmp = malloc(m * (size_t)n * sizeof(uint64_t));
  for (size_t i = 0; i < m; i++) {
    size_t dst = start[vec_list_sum(l, from + i) - sum_min]++;
    memcpy(tmp + dst * n, base + i * n, (size_t)n * sizeof(uint64_t));
  }
  /* start[b] now holds the end of bucket b */
  size_t lo = 0;
  for (size_t b = 0; b + 1 < nb; b++) {
    size_t hi = start[b];
    if (hi - lo > 1)
      qsort(tmp + lo * n, hi - lo, (size_t)n * sizeof(uint64_t),
            vec_cmp_by_sum);
    lo = hi;
  }
  memcpy(base, tmp, m * (size_t)n * sizeof(uint64_t));
  free(tmp);
  free(start);
}

uint64_t enum_vectors(const prime_exps_t *p, int n, uint64_t sum_min,
                      uint64_t sum_max, vec_list_t *out) {
  if (n < 1 || n > ENUM_MAX_N || pexp_value(p) == 0)
    return 0;
  size_t num_divs;
  divisor_t *divs = make_divisors(p, &num_divs);
  uint64_t pk = 0;
  for (int i = 0; i < p->num_primes; i++)
    pk += (uint64_t)p->exps[i] << (FIELD_BITS * i);

  enum_state_t st = {.divs = divs,
                     .num_divs = num_divs,
                     .guard = guard_mask(),
                     .sum_min = sum_min,
                     .sum_max = sum_max,
                     .n = n,
                     .out = out,
                     .nodes = 0};
  size_t before = out->count;
  enum_rec(&st, n, pexp_value(p), pk, UINT64_MAX, 0);
  free(divs);

  sort_by_sum(out, before, sum_min, sum_max);
  return st.nodes;
}

uint64_t enum_min_sum(const prime_exps_t *p, int n) {
  uint64_t P = pexp_value(p);
  if (P == 0)
    return 0;
  double base = n * pow((double)P, 1.0 / n);
  uint64_t hi = (uint64_t)(base * 1.02) + n * n;
  for (int attempt = 0; attempt < 60; attempt++) {
    vec_list_t l;
    vec_list_init(&l, n);
    enum_vectors(p, n, 0, hi, &l);
    if (l.count > 0) {
      uint64_t best = vec_list_sum(&l, 0); /* sorted by sum */
      vec_list_free(&l);
      return best;
    }
    vec_list_free(&l);
    hi = hi * 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */
/* reduction                                                               */

static int u64_cmp(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}

size_t reduce_vectors(const vec_list_t *l, size_t start, size_t count,
                      vec_list_t *out, int strong) {
  int n = l->n;
  if (count < (size_t)(2 * n))
    return 0;

  /* relabel elements densely */
  size_t ne = count * n;
  uint64_t *vals = malloc(ne * sizeof(uint64_t));
  memcpy(vals, vec_list_get(l, start), ne * sizeof(uint64_t));
  qsort(vals, ne, sizeof(uint64_t), u64_cmp);
  size_t nlabels = 0;
  for (size_t i = 0; i < ne; i++)
    if (i == 0 || vals[i] != vals[i - 1])
      vals[nlabels++] = vals[i];

  uint32_t *lab = malloc(ne * sizeof(uint32_t)); /* lab[v*n + j] */
  for (size_t v = 0; v < count; v++) {
    const uint64_t *vec = vec_list_get(l, start + v);
    for (int j = 0; j < n; j++) {
      size_t lo = 0, hi = nlabels;
      while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (vals[mid] < vec[j])
          lo = mid + 1;
        else
          hi = mid;
      }
      lab[v * n + j] = (uint32_t)lo;
    }
  }

  /* element -> list of vectors containing it (CSR) */
  uint32_t *deg = calloc(nlabels + 1, sizeof(uint32_t));
  for (size_t i = 0; i < ne; i++)
    deg[lab[i] + 1]++;
  for (size_t x = 0; x < nlabels; x++)
    deg[x + 1] += deg[x];
  uint32_t *occ = malloc(ne * sizeof(uint32_t));
  uint32_t *fill = malloc(nlabels * sizeof(uint32_t));
  memcpy(fill, deg, nlabels * sizeof(uint32_t));
  for (size_t v = 0; v < count; v++)
    for (int j = 0; j < n; j++)
      occ[fill[lab[v * n + j]]++] = (uint32_t)v;

  bool *alive = malloc(count);
  memset(alive, 1, count);
  size_t num_alive = count;
  uint32_t *alive_deg = malloc(nlabels * sizeof(uint32_t));
  for (size_t x = 0; x < nlabels; x++)
    alive_deg[x] = deg[x + 1] - deg[x];

  uint8_t *inter = calloc(count, 1);
  uint32_t *touched = malloc(count * sizeof(uint32_t));

  bool changed = true;
  while (changed && num_alive >= (size_t)(2 * n)) {
    changed = false;
    for (size_t v = 0; v < count; v++) {
      if (!alive[v])
        continue;
      bool ok = true;
      /* every element must be in some other vector */
      for (int j = 0; j < n && ok; j++)
        if (alive_deg[lab[v * n + j]] < 2)
          ok = false;

      if (ok && strong) {
        /* intersection sizes with all alive vectors touching v */
        size_t nt = 0;
        for (int j = 0; j < n; j++) {
          uint32_t x = lab[v * n + j];
          for (uint32_t k = deg[x]; k < deg[x + 1]; k++) {
            uint32_t w = occ[k];
            if (w == v || !alive[w])
              continue;
            if (inter[w]++ == 0)
              touched[nt++] = w;
          }
        }
        /* need >= n-1 disjoint vectors */
        if (num_alive - 1 - nt < (size_t)(n - 1))
          ok = false;
        /* each element needs a vector meeting v in exactly that element */
        for (int j = 0; j < n && ok; j++) {
          uint32_t x = lab[v * n + j];
          bool found = false;
          for (uint32_t k = deg[x]; k < deg[x + 1] && !found; k++) {
            uint32_t w = occ[k];
            if (w != v && alive[w] && inter[w] == 1)
              found = true;
          }
          if (!found)
            ok = false;
        }
        for (size_t t = 0; t < nt; t++)
          inter[touched[t]] = 0;
      }

      if (!ok) {
        alive[v] = false;
        num_alive--;
        for (int j = 0; j < n; j++)
          alive_deg[lab[v * n + j]]--;
        changed = true;
      }
    }
  }

  size_t kept = 0;
  if (num_alive >= (size_t)(2 * n)) {
    for (size_t v = 0; v < count; v++)
      if (alive[v]) {
        vec_list_push(out, vec_list_get(l, start + v));
        kept++;
      }
  }

  free(vals);
  free(lab);
  free(deg);
  free(occ);
  free(fill);
  free(alive);
  free(alive_deg);
  free(inter);
  free(touched);
  return kept;
}
