/*
 * msearch: enumerate + arrange + check diagonals, for one P and a range of S,
 * in a single process with no intermediate files.
 *
 * usage: msearch [options] e1 e2 ...      (P = 2^e1 3^e2 5^e3 ...)
 *
 *   --vec-size N      square size (default 6)
 *   --min-sum S       smallest sum to search (default: smallest possible)
 *   --max-sum S       largest sum to search (default: min-sum + 99)
 *   --sums S1,S2,...  search exactly these sums instead of a range
 *   --node-limit X    stop the search of one sum after X nodes (record is
 *                     marked truncated)
 *   --total-nodes X   stop after the sum during which X total nodes are hit
 *   --time-limit T    stop after the sum during which T seconds pass
 *   --out FILE        append JSON lines to FILE (default stdout)
 *   --format legacy   print in the old arrangement.c format instead, for
 *                     use with postprocess.py
 *   --reduce strong   use the stronger (slower) vector reduction
 *   --no-fc/--no-mrv  disable forward checking / most-constrained branching
 *   --no-support      disable the support filter (see arrange_core.h)
 *   --no-cross        disable the cross support filter (see arrange_core.h)
 *   --pretest-min K   pretest the children with K or more vectors placed
 *                     (default 5; 0 = never; see search_opts_t)
 *
 * Output (JSON lines, one record per searched sum, flushed immediately so a
 * killed run keeps all completed sums):
 *   {"type":"sum","n":6,"P":[13,6,3,2],"Pval":...,"S":506,"nvecs":831,
 *    "nvecs_raw":831,"labels":77,"nodes":890845,"squares":0,
 *    "time":0.23,"setup_time":..,"enum_time":0.001,"cpu":0.24,"truncated":0,
 *    "engine":2}
 * (time, setup_time, enum_time: wall seconds of the search, of the reduction
 * and setup, and this sum's share of the enumeration; cpu: the process CPU
 * seconds of all three, which the scheduler fits its time law to) and one
 * record per semi-magic square found, preceding its sum record (a run killed
 * in the middle of a sum leaves squares without a sum record: the scheduler
 * counts a square only once its sum record follows):
 *   {"type":"square","n":6,"P":[...],"S":...,"s_count":..,"p_count":..,
 *    "sp_count":..,"best_score":..,"grid":[[...],...]}
 * where grid has the best pair of diagonals on its main diagonals. A final
 *   {"type":"done","min_sum":..,"last_sum":..,"complete":0|1,...}
 * says every sum in [min_sum, last_sum] was searched (sums with fewer than 2n
 * vectors after reduction cannot have a square and get no "sum" record); it
 * is omitted with --sums.
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "enumerate.h"
#include "arrange.h"
#include "square.h"

/* version of the search, written into every "sum" record: bump it whenever
 * the search time changes materially, since the scheduler fits its model of
 * the time per sum only to records of the newest version (1 = the first
 * version of this pipeline, without the field; 2 = support and cross
 * filters, carried bitsets, October 2026) */
#define ENGINE_VERSION 2

static double wall_time(void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (double)t.tv_sec + 1e-6 * (double)t.tv_usec;
}

/* CPU time of this (single-threaded) process: the "cpu" field of a sum
 * record, which the scheduler learns its time law from (wall time includes
 * waiting on a loaded machine) */
static double cpu_time(void) {
  struct timespec t;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
  return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

typedef struct {
  FILE *out;
  int legacy;
  const prime_exps_t *p;
  uint64_t S;
  int square_index;
} out_ctx_t;

static void print_p_json(FILE *fp, const prime_exps_t *p) {
  fputc('[', fp);
  for (int i = 0; i < p->num_primes; i++)
    fprintf(fp, i ? ",%d" : "%d", p->exps[i]);
  fputc(']', fp);
}

static int square_found(const square_t *sq, void *vctx) {
  out_ctx_t *ctx = vctx;
  FILE *fp = ctx->out;
  if (ctx->legacy) {
    fprintf(fp, "solution found\n");
    for (int i = 0; i < sq->n; i++) {
      fprintf(fp, "%d   ", i);
      for (int j = 0; j < sq->n; j++)
        fprintf(fp, " %lu", (unsigned long)sq->rows[i][j]);
      fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    for (int i = 0; i < sq->n; i++) {
      fprintf(fp, "%d   ", i);
      for (int j = 0; j < sq->n; j++)
        fprintf(fp, " %lu", (unsigned long)sq->cols[i][j]);
      fprintf(fp, "\n");
    }
  } else {
    diag_stats_t ds;
    square_diag_stats(sq, &ds);
    fprintf(fp, "{\"type\":\"square\",\"n\":%d,\"P\":", sq->n);
    print_p_json(fp, ctx->p);
    fprintf(fp,
            ",\"S\":%lu,\"s_count\":%d,\"p_count\":%d,\"sp_count\":%d,"
            "\"best_score\":%d,\"hash\":\"%016lx\",\"grid\":",
            (unsigned long)ctx->S, ds.s_count, ds.p_count, ds.sp_count,
            ds.best_score, (unsigned long)square_hash(sq));
    grid_print_json(fp, &ds.best);
    fprintf(fp, "}\n");
    if (ds.best_score >= 14)
      fprintf(stderr, "!!! MAGIC SQUARE FOUND (S=%lu) !!!\n",
              (unsigned long)ctx->S);
  }
  fflush(fp);
  ctx->square_index++;
  return 0;
}

static int parse_sums(const char *s, uint64_t **out) {
  int cap = 16, cnt = 0;
  *out = malloc(cap * sizeof(uint64_t));
  const char *c = s;
  while (*c) {
    char *end;
    uint64_t v = strtoull(c, &end, 10);
    if (end == c)
      return -1;
    if (cnt == cap)
      *out = realloc(*out, (cap *= 2) * sizeof(uint64_t));
    (*out)[cnt++] = v;
    c = end;
    if (*c == ',')
      c++;
  }
  return cnt;
}

static int u64_cmp(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x > y) - (x < y);
}

int main(int argc, char *argv[]) {
  int n = 6;
  uint64_t min_sum = 0, max_sum = 0;
  int have_min = 0, have_max = 0;
  uint64_t *sum_list = NULL;
  int num_sum_list = -1;
  uint64_t total_node_limit = 0;
  double time_limit = 0;
  const char *out_file = NULL;
  int legacy = 0, strong = 0;
  search_opts_t opts;
  search_opts_default(&opts);

  static struct option long_options[] = {
      {"vec-size", required_argument, 0, 'n'},
      {"min-sum", required_argument, 0, 'a'},
      {"max-sum", required_argument, 0, 'b'},
      {"sums", required_argument, 0, 's'},
      {"node-limit", required_argument, 0, 'l'},
      {"total-nodes", required_argument, 0, 't'},
      {"time-limit", required_argument, 0, 'T'},
      {"out", required_argument, 0, 'o'},
      {"format", required_argument, 0, 'F'},
      {"reduce", required_argument, 0, 'r'},
      {"no-fc", no_argument, 0, 'f'},
      {"no-mrv", no_argument, 0, 'm'},
      {"no-support", no_argument, 0, 'S'},
      {"no-cross", no_argument, 0, 'X'},
      {"pretest-min", required_argument, 0, 'P'},
      {0, 0, 0, 0}};
  int opt;
  while ((opt = getopt_long(argc, argv, "n:a:b:s:l:t:T:o:F:r:fmSXP:",
                            long_options, NULL)) != -1) {
    switch (opt) {
    case 'n':
      n = atoi(optarg);
      break;
    case 'a':
      min_sum = strtoull(optarg, NULL, 10);
      have_min = 1;
      break;
    case 'b':
      max_sum = strtoull(optarg, NULL, 10);
      have_max = 1;
      break;
    case 's':
      num_sum_list = parse_sums(optarg, &sum_list);
      if (num_sum_list < 0) {
        fprintf(stderr, "bad --sums\n");
        return 2;
      }
      break;
    case 'l':
      opts.node_limit = strtoull(optarg, NULL, 10);
      break;
    case 't':
      total_node_limit = strtoull(optarg, NULL, 10);
      break;
    case 'T':
      time_limit = atof(optarg);
      break;
    case 'o':
      out_file = optarg;
      break;
    case 'F':
      legacy = !strcmp(optarg, "legacy");
      break;
    case 'r':
      strong = !strcmp(optarg, "strong");
      break;
    case 'f':
      opts.forward_check = 0;
      break;
    case 'm':
      opts.mrv = 0;
      break;
    case 'S':
      opts.support = 0;
      break;
    case 'X':
      opts.cross = 0;
      break;
    case 'P':
      opts.pretest_min = atoi(optarg);
      break;
    default:
      return 2;
    }
  }
  prime_exps_t p;
  if (n < 3 || n > SQ_MAX_N || optind >= argc ||
      pexp_parse(&p, argc - optind, argv + optind) || pexp_value(&p) == 0) {
    fprintf(stderr, "usage: msearch [--vec-size N] [--min-sum S] [--max-sum S] "
                    "[--sums S1,S2,..] [--node-limit X] [--total-nodes X] "
                    "[--time-limit T] [--out FILE] [--format json|legacy] "
                    "e1 e2 ...\n");
    return 2;
  }

  FILE *out = stdout;
  if (out_file) {
    out = fopen(out_file, "a");
    if (!out) {
      perror(out_file);
      return 2;
    }
  }

  if (num_sum_list >= 0) {
    qsort(sum_list, num_sum_list, sizeof(uint64_t), u64_cmp);
    if (num_sum_list == 0)
      return 0;
    min_sum = sum_list[0];
    max_sum = sum_list[num_sum_list - 1];
  } else {
    if (!have_min)
      min_sum = enum_min_sum(&p, n);
    if (!have_max)
      max_sum = min_sum + 99;
  }

  out_ctx_t ctx = {.out = out, .legacy = legacy, .p = &p};
  uint64_t total_nodes = 0;
  double t_start = wall_time();
  int stop = 0, stop_nodes = 0; /* stop_nodes: by --total-nodes */
  if (legacy)
    fprintf(out, "read\n");

  /* enumerate in windows of sums to bound memory */
  /* the window adapts so each holds roughly 10^5 - 10^6 vectors */
  uint64_t window = 32;
  int list_pos = 0;
  uint64_t last_sum = min_sum ? min_sum - 1 : 0;
  for (uint64_t lo = min_sum, hi; lo <= max_sum && !stop; lo = hi + 1) {
    hi = lo + window - 1 < max_sum ? lo + window - 1 : max_sum;
    if (num_sum_list >= 0) {
      while (list_pos < num_sum_list && sum_list[list_pos] < lo)
        list_pos++;
      if (list_pos >= num_sum_list)
        break;
      if (sum_list[list_pos] > hi)
        continue;
    }
    double te = wall_time(), ce = cpu_time();
    vec_list_t all, red;
    vec_list_init(&all, n);
    vec_list_init(&red, n);
    enum_vectors_grouped(&p, n, lo, hi, &all);
    double enum_time = wall_time() - te;
    double enum_cpu = cpu_time() - ce;
    size_t i = 0;
    for (uint64_t S = lo; S <= hi && !stop; S++) {
      last_sum = S;
      size_t j = i;
      while (j < all.count && vec_list_sum(&all, j) == S)
        j++;
      size_t raw = j - i;
      size_t start_i = i;
      i = j;
      if (num_sum_list >= 0) {
        int found = 0;
        for (int k = list_pos; k < num_sum_list && sum_list[k] <= hi; k++)
          found |= sum_list[k] == S;
        if (!found)
          continue;
      }
      double tr = wall_time(), cr = cpu_time();
      red.count = 0;
      reduce_vectors(&all, start_i, raw, &red, strong);
      double reduce_time = wall_time() - tr;
      if (red.count == 0)
        continue; /* no square possible: not recorded, like arrangement.c */
      ctx.S = S;
      if (legacy)
        fprintf(out, "sum %lu nvecs %zu\n", (unsigned long)S, red.count);
      search_stats_t st =
          search_vectors(&red, 0, red.count, &opts, square_found, &ctx);
      total_nodes += st.nodes;
      double share = (double)raw / (all.count ? all.count : 1);
      double cpu = cpu_time() - cr + enum_cpu * share;
      if (legacy) {
        fprintf(out, "num searched: %lu\n", (unsigned long)st.nodes);
      } else {
        char pstr[64];
        pexp_to_str(&p, pstr, ",");
        fprintf(out,
                "{\"type\":\"sum\",\"n\":%d,\"P\":[%s],\"Pval\":%lu,\"S\":%lu,"
                "\"nvecs\":%zu,\"nvecs_raw\":%zu,\"labels\":%d,\"nodes\":%lu,"
                "\"squares\":%lu,\"time\":%.6f,\"setup_time\":%.6f,"
                "\"enum_time\":%.6f,\"cpu\":%.6f,\"truncated\":%d,\"engine\":%d}\n",
                n, pstr, (unsigned long)pexp_value(&p), (unsigned long)S,
                red.count, raw, st.num_labels, (unsigned long)st.nodes,
                (unsigned long)st.squares, st.seconds,
                st.setup_seconds + reduce_time,
                enum_time * share, cpu, st.truncated, ENGINE_VERSION);
      }
      fflush(out);
      if (total_node_limit && total_nodes >= total_node_limit)
        stop = stop_nodes = 1;
      if (time_limit > 0 && wall_time() - t_start >= time_limit)
        stop = 1;
    }
    if (all.count < 100000 && window < 65536)
      window *= 2;
    else if (all.count > 1000000 && window > 1)
      window /= 2;
    vec_list_free(&all);
    vec_list_free(&red);
  }
  if (!legacy && num_sum_list < 0) {
    /* every sum in [min_sum, last_sum] has been searched (sums with too few
     * vectors produce no "sum" record) */
    char pstr[64];
    pexp_to_str(&p, pstr, ",");
    fprintf(out,
            "{\"type\":\"done\",\"n\":%d,\"P\":[%s],\"min_sum\":%lu,"
            "\"last_sum\":%lu,\"complete\":%d,\"time\":%.3f}\n",
            n, pstr, (unsigned long)min_sum, (unsigned long)last_sum,
            !stop || last_sum >= max_sum, wall_time() - t_start);
  }
  if (legacy) {
    if (stop_nodes)
      fprintf(out, "terminated: total count %lu exceeds cutoff %lu\n",
              (unsigned long)total_nodes, (unsigned long)total_node_limit);
    fprintf(out, "completed in %.5f secs\n", wall_time() - t_start);
  }
  if (out != stdout)
    fclose(out);
  free(sum_list);
  return 0;
}
