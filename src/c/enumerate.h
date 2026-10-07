#ifndef ENUMERATE_H
#define ENUMERATE_H

#include <stddef.h>
#include <stdint.h>

/*
 * Low-memory enumeration of "vectors": sets of `n` distinct positive integers
 * with product P and sum in [sum_min, sum_max].
 *
 * Unlike the old enumerators (which materialize every multiset of divisors with
 * product P, and need O(prod C(e_i + 5, 5)) memory), this is a depth-first
 * search over divisors in decreasing order, pruned by the sum bounds. Memory
 * use is proportional to the output, and the cost is roughly proportional to
 * the output size when sum_max is small relative to the full range of sums.
 */

#define ENUM_MAX_PRIMES 10
#define ENUM_MAX_N 8

extern const uint64_t ENUM_PRIMES[ENUM_MAX_PRIMES];

/* exponents of P over ENUM_PRIMES (2, 3, 5, ...) */
typedef struct {
  int exps[ENUM_MAX_PRIMES];
  int num_primes;
} prime_exps_t;

/* value of P, or 0 on overflow */
uint64_t pexp_value(const prime_exps_t *p);

/* number of divisors of P */
uint64_t pexp_num_divisors(const prime_exps_t *p);

/* parse "10 4 3 2" style exponent list from argv; returns 0 on success */
int pexp_parse(prime_exps_t *p, int argc, char **argv);

/* string like "10_4_3_2" (buf must hold >= 64 chars) */
void pexp_to_str(const prime_exps_t *p, char *buf, const char *sep);

/* a list of vectors, `n` elements each, stored contiguously, each vector
 * sorted descending */
typedef struct {
  uint64_t *elts;
  size_t count;
  size_t cap;
  int n;
} vec_list_t;

void vec_list_init(vec_list_t *l, int n);
void vec_list_free(vec_list_t *l);
void vec_list_push(vec_list_t *l, const uint64_t *vec);
static inline uint64_t *vec_list_get(const vec_list_t *l, size_t i) {
  return l->elts + i * (size_t)l->n;
}
static inline uint64_t vec_list_sum(const vec_list_t *l, size_t i) {
  uint64_t s = 0;
  for (int j = 0; j < l->n; j++)
    s += l->elts[i * (size_t)l->n + j];
  return s;
}

/* smallest possible sum of n distinct divisors of P with product P (or 0 if
 * there is none) -- found by the same search, so it is exact */
uint64_t enum_min_sum(const prime_exps_t *p, int n);

/*
 * Enumerate all vectors with sum in [sum_min, sum_max] and append them to
 * `out`, sorted by (sum, then lexicographically descending elements).
 * Returns number of DFS nodes visited (for profiling).
 */
uint64_t enum_vectors(const prime_exps_t *p, int n, uint64_t sum_min,
                      uint64_t sum_max, vec_list_t *out);

/*
 * Given the vectors with one fixed sum (contiguous in `l` from `start`,
 * `count` long), remove vectors that cannot appear in any n x n semi-magic
 * square, iterating to a fixed point:
 *   - every element of a vector must lie in another vector meeting it in
 *     exactly that one element (the vector through that cell on the other
 *     axis), and
 *   - a vector must be disjoint from at least n - 1 other vectors (the other
 *     rows/cols on its axis).
 * If fewer than 2n vectors remain, everything is removed. The surviving
 * vectors are written to `out` (appended). Returns number of survivors.
 */
size_t reduce_vectors(const vec_list_t *l, size_t start, size_t count,
                      vec_list_t *out, int strong);

#endif // !ENUMERATE_H
