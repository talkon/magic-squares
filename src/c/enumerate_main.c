/*
 * Standalone enumeration tool, a drop-in replacement for enumeration.cpp /
 * enumeration.py that writes the same file format (one vector per line,
 * elements descending, ascending order of sum) but uses far less memory and
 * can be restricted to a range of sums.
 *
 * usage: enumerate [--vec-size N] [--min-sum S] [--max-sum S]
 *                  [--reduce weak|strong|none] [--counts] --file F  e1 e2 ...
 *        enumerate [--vec-size N] --print-min-sum [--window F [--until C]]
 *                  e1 e2 ...
 *
 * --counts writes "S count" lines instead of the vectors, and
 * --print-min-sum prints the smallest sum S_min of any vector and exits;
 * with --window F it also prints "S count" lines (before reduction) for the
 * sums in [S_min, F * S_min], stopping early (at the end of a chunk of sums)
 * once some sum has at least C vectors if --until C is given.
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "enumerate.h"

int main(int argc, char *argv[]) {
  int n = 6;
  uint64_t sum_min = 0, sum_max = UINT64_MAX;
  const char *file = NULL;
  int reduce = 1; /* 0 = none, 1 = weak (same as old enumerators), 2 = strong */
  int counts = 0, print_min = 0;
  double window = 0;
  uint64_t until = UINT64_MAX;

  static struct option long_options[] = {
      {"vec-size", required_argument, 0, 'n'},
      {"min-sum", required_argument, 0, 'a'},
      {"max-sum", required_argument, 0, 'b'},
      {"file", required_argument, 0, 'f'},
      {"reduce", required_argument, 0, 'r'},
      {"counts", no_argument, 0, 'c'},
      {"print-min-sum", no_argument, 0, 'm'},
      {"window", required_argument, 0, 'w'},
      {"until", required_argument, 0, 'u'},
      {0, 0, 0, 0}};
  int opt;
  while ((opt = getopt_long(argc, argv, "n:a:b:f:r:cmw:u:", long_options, NULL)) !=
         -1) {
    switch (opt) {
    case 'n':
      n = atoi(optarg);
      break;
    case 'a':
      sum_min = strtoull(optarg, NULL, 10);
      break;
    case 'b':
      sum_max = strtoull(optarg, NULL, 10);
      break;
    case 'f':
      file = optarg;
      break;
    case 'r':
      reduce = !strcmp(optarg, "none") ? 0 : !strcmp(optarg, "strong") ? 2 : 1;
      break;
    case 'c':
      counts = 1;
      break;
    case 'm':
      print_min = 1;
      break;
    case 'w':
      window = atof(optarg);
      break;
    case 'u':
      until = strtoull(optarg, NULL, 10);
      break;
    default:
      fprintf(stderr, "bad arguments\n");
      return 1;
    }
  }
  prime_exps_t p;
  if (print_min) {
    if (pexp_parse(&p, argc - optind, argv + optind) != 0)
      return 1;
    uint64_t smin = enum_min_sum(&p, n);
    printf("%lu\n", (unsigned long)smin);
    if (window > 1 && smin > 0) {
      uint64_t end = (uint64_t)(window * smin);
      uint64_t chunk = smin / 50 > 8 ? smin / 50 : 8;
      int done = 0;
      for (uint64_t lo = smin; lo <= end && !done; lo += chunk) {
        uint64_t hi = lo + chunk - 1 < end ? lo + chunk - 1 : end;
        vec_list_t l;
        vec_list_init(&l, n);
        enum_vectors_grouped(&p, n, lo, hi, &l); /* counts only */
        for (size_t i = 0; i < l.count;) {
          uint64_t S = vec_list_sum(&l, i);
          size_t j = i;
          while (j < l.count && vec_list_sum(&l, j) == S)
            j++;
          printf("%lu %zu\n", (unsigned long)S, j - i);
          done |= j - i >= until;
          i = j;
        }
        if (done)
          printf("# counted to %lu\n", (unsigned long)hi);
        vec_list_free(&l);
      }
    }
    return 0;
  }
  if (!file || pexp_parse(&p, argc - optind, argv + optind) != 0) {
    fprintf(stderr, "usage: enumerate [--vec-size N] [--min-sum S] "
                    "[--max-sum S] [--reduce none|weak|strong] [--counts] "
                    "--file F e1 e2 ...\n"
                    "       enumerate [--vec-size N] --print-min-sum "
                    "[--window F [--until C]] e1 e2 ...\n");
    return 1;
  }

  vec_list_t all;
  vec_list_init(&all, n);
  /* the order within a sum only matters when writing the vectors */
  uint64_t nodes = counts ? enum_vectors_grouped(&p, n, sum_min, sum_max, &all)
                          : enum_vectors(&p, n, sum_min, sum_max, &all);
  fprintf(stderr, ">>> (enum) %zu %d-vecs with sum in range (%lu dfs nodes)\n",
          all.count, n, (unsigned long)nodes);

  FILE *fp = fopen(file, "w");
  if (!fp) {
    perror(file);
    return 1;
  }
  size_t written = 0;
  vec_list_t red;
  vec_list_init(&red, n);
  for (size_t i = 0; i < all.count;) {
    uint64_t s = vec_list_sum(&all, i);
    size_t j = i;
    while (j < all.count && vec_list_sum(&all, j) == s)
      j++;
    red.count = 0;
    if (reduce)
      reduce_vectors(&all, i, j - i, &red, reduce == 2);
    else
      for (size_t k = i; k < j; k++)
        vec_list_push(&red, vec_list_get(&all, k));
    if (counts) {
      if (red.count)
        fprintf(fp, "%lu %zu\n", (unsigned long)s, red.count);
    } else {
      for (size_t k = 0; k < red.count; k++) {
        uint64_t *v = vec_list_get(&red, k);
        for (int t = 0; t < n; t++)
          fprintf(fp, t ? " %lu" : "%lu", (unsigned long)v[t]);
        fputc('\n', fp);
      }
    }
    written += red.count;
    i = j;
  }
  fclose(fp);
  fprintf(stderr, ">>> (enum) wrote %zu rows to %s in ascending order of sum\n",
          written, file);
  vec_list_free(&all);
  vec_list_free(&red);
  return 0;
}
