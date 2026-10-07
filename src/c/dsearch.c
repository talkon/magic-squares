/*
 * dsearch: experiment driver for the "diagonal-first" search.
 *
 * A semi-magic square has a vector d as a traversal (a possible diagonal)
 * iff every row and every col meets d in exactly one number. So the squares
 * with d as a traversal are exactly the semi-magic squares that can be built
 * from V_d = { v : |v & d| = 1 }, and a magic square is a square with two
 * disjoint vectors as traversals.
 *
 * For each (P, S) given, this runs
 *   (a) the plain semi-magic search (all squares; counts the SP traversals
 *       of each with square_diag_stats), and
 *   (b) for every vector d (from the unreduced list), the semi-magic search
 *       on V_d (taken from the reduced list), which finds every
 *       (square, SP traversal) pair exactly once; the totals must agree
 *       with (a). Each square found in (b) is checked for a second SP
 *       traversal disjoint from d (a magic square).
 * and prints the time of each, to see how the ratio behaves with N.
 *
 * usage: dsearch [--skip-plain] --sums S1,S2,... e1 e2 ...
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "enumerate.h"
#include "arrange.h"
#include "square.h"

static double wall_time(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + 1e-6 * (double)t.tv_usec;
}

typedef struct {
  uint64_t squares, sp_pairs, magic;
  uint64_t s_trav, p_trav;
} plain_ctx_t;

static int plain_cb(const square_t *sq, void *vctx) {
  plain_ctx_t *c = vctx;
  diag_stats_t ds;
  square_diag_stats(sq, &ds);
  c->squares++;
  c->sp_pairs += ds.sp_count;
  c->s_trav += ds.s_count;
  c->p_trav += ds.p_count;
  if (ds.best_score >= 14)
    c->magic++;
  return 0;
}

typedef struct {
  uint64_t squares, magic;
  int n;
} d_ctx_t;

static int d_cb(const square_t *sq, void *vctx) {
  d_ctx_t *c = vctx;
  c->squares++;
  diag_stats_t ds;
  square_diag_stats(sq, &ds);
  if (ds.best_score >= 14) {
    c->magic++;
    fprintf(stderr, "!!! MAGIC SQUARE !!!\n");
    grid_print_json(stderr, &ds.best);
    fputc('\n', stderr);
  }
  return 0;
}

static int inter_count(const uint64_t *a, const uint64_t *b, int n) {
  int k = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      k += a[i] == b[j];
  return k;
}

int main(int argc, char *argv[]) {
  int skip_plain = 0;
  search_opts_t opts;
  search_opts_default(&opts);
  static struct option long_options[] = {
      {"skip-plain", no_argument, 0, 'p'}, {"sums", required_argument, 0, 's'},
      {0, 0, 0, 0}};
  int opt;
  char *sums = NULL;
  while ((opt = getopt_long(argc, argv, "ps:", long_options, NULL)) != -1) {
    switch (opt) {
    case 'p':
      skip_plain = 1;
      break;
    case 's':
      sums = optarg;
      break;
    default:
      return 2;
    }
  }
  prime_exps_t p;
  if (!sums || optind >= argc || pexp_parse(&p, argc - optind, argv + optind)) {
    fprintf(stderr, "usage: dsearch [--skip-plain] --sums S1,S2,... e1 e2 ...\n");
    return 2;
  }
  uint64_t sum_list[256];
  int nsums = 0;
  for (char *tok = strtok(sums, ","); tok && nsums < 256; tok = strtok(NULL, ","))
    sum_list[nsums++] = strtoull(tok, NULL, 10);
  const int n = 6;
  char pstr[64];
  pexp_to_str(&p, pstr, ".");
  printf("%-6s %6s %6s | %9s %10s %6s %6s | %9s %9s %6s %5s %7s %11s | %7s %7s\n", "S", "N",
         "Nred", "plain(s)", "nodes", "sq", "SPprs", "diag(s)", "search(s)", "found", "magic",
         "avg|V_d|", "nodes", "t-ratio", "n-ratio");
  for (int i = 0; i < nsums; i++) {
    uint64_t S = sum_list[i];
    vec_list_t all, red;
    vec_list_init(&all, n);
    vec_list_init(&red, n);
    enum_vectors(&p, n, S, S, &all);
    reduce_vectors(&all, 0, all.count, &red, 0);
    if (red.count < (size_t)(2 * n)) {
      printf("%-10lu %6zu %6zu | no squares possible\n", (unsigned long)S,
             all.count, red.count);
      vec_list_free(&all);
      vec_list_free(&red);
      continue;
    }

    plain_ctx_t pc = {0};
    double t_plain = 0;
    uint64_t plain_nodes = 0;
    if (!skip_plain) {
      search_stats_t st =
          search_vectors(&red, 0, red.count, &opts, plain_cb, &pc);
      t_plain = st.seconds + st.setup_seconds;
      plain_nodes = st.nodes;
    }

    /* (b): for each d among all vectors, search V_d */
    double t0 = wall_time();
    d_ctx_t dc = {0, 0, n};
    vec_list_t sub;
    vec_list_init(&sub, n);
    double sum_sub = 0, t_search = 0;
    uint64_t nodes = 0;
    for (size_t d = 0; d < all.count; d++) {
      const uint64_t *dv = vec_list_get(&all, d);
      sub.count = 0;
      for (size_t v = 0; v < red.count; v++) {
        const uint64_t *vv = vec_list_get(&red, v);
        if (inter_count(dv, vv, n) == 1)
          vec_list_push(&sub, vv);
      }
      sum_sub += (double)sub.count;
      if (sub.count < (size_t)(2 * n))
        continue;
      search_stats_t st =
          search_vectors(&sub, 0, sub.count, &opts, d_cb, &dc);
      nodes += st.nodes;
      t_search += st.seconds;
    }
    double t_diag = wall_time() - t0;
    vec_list_free(&sub);

    printf("%-6lu %6zu %6zu | %9.3f %10lu %6lu %6lu | %9.3f %9.3f %6lu %5lu %7.1f %11lu | %7.2f %7.2f\n",
           (unsigned long)S, all.count, red.count, t_plain,
           (unsigned long)plain_nodes, (unsigned long)pc.squares,
           (unsigned long)pc.sp_pairs, t_diag, t_search,
           (unsigned long)dc.squares, (unsigned long)dc.magic,
           sum_sub / (double)all.count, (unsigned long)nodes,
           t_plain > 0 ? t_search / t_plain : 0.0,
           plain_nodes ? (double)nodes / (double)plain_nodes : 0.0);
    fflush(stdout);
    vec_list_free(&all);
    vec_list_free(&red);
  }
  return 0;
}
