/*
 * Arrangement benchmark.
 *
 * usage: bench [--repeat K] [--only SUBSTR] [--no-fc] [--no-mrv] [--no-support]
 *              [--update] [--min-words K] [--gather] [--node-limit X]
 *              instances.txt
 *
 * Each instance line is
 *     n | e1 e2 ... | S | expected_squares | expected_hash
 * (blank lines and lines starting with '#' are ignored). Vectors are
 * enumerated in-process, then the search is timed (min over K repeats). The
 * squares found are checked against the expected count and an
 * order-independent hash (sum of per-square canonical hashes). Use "-" for
 * unknown expectations, and --update to print instance lines with the
 * observed values filled in.
 *
 * --no-fc, --no-mrv, --no-support, --min-words and --gather change the
 * search (to compare variants, or to test the code paths for wider label
 * bitsets or larger N on small instances). --node-limit X stops each search
 * after X nodes (to time instances too large to search completely; the
 * check is then meaningless, use "-" expectations).
 *
 * Exit status is nonzero if any instance does not match.
 */
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "enumerate.h"
#include "arrange.h"
#include "square.h"

typedef struct {
  uint64_t hash;
  uint64_t count;
} hash_ctx_t;

static int hash_cb(const square_t *sq, void *ctx) {
  hash_ctx_t *h = ctx;
  h->hash += square_hash(sq);
  h->count++;
  return 0;
}

static char *trim(char *s) {
  while (*s == ' ' || *s == '\t')
    s++;
  char *e = s + strlen(s);
  while (e > s && (e[-1] == ' ' || e[-1] == '\t' || e[-1] == '\n' ||
                   e[-1] == '\r'))
    *--e = '\0';
  return s;
}

int main(int argc, char *argv[]) {
  int repeat = 1, update = 0;
  const char *only = NULL;
  search_opts_t opts;
  search_opts_default(&opts);

  static struct option long_options[] = {{"repeat", required_argument, 0, 'r'},
                                         {"only", required_argument, 0, 'o'},
                                         {"no-fc", no_argument, 0, 'f'},
                                         {"no-mrv", no_argument, 0, 'm'},
                                         {"update", no_argument, 0, 'u'},
                                         {"min-words", required_argument, 0, 'w'},
                                         {"gather", no_argument, 0, 'g'},
                                         {"no-support", no_argument, 0, 's'},
                                         {"node-limit", required_argument, 0, 'L'},
                                         {0, 0, 0, 0}};
  int opt;
  while ((opt = getopt_long(argc, argv, "r:o:fmuw:gsL:", long_options, NULL)) !=
         -1) {
    switch (opt) {
    case 'r':
      repeat = atoi(optarg);
      break;
    case 'o':
      only = optarg;
      break;
    case 'f':
      opts.forward_check = 0;
      break;
    case 'm':
      opts.mrv = 0;
      break;
    case 'u':
      update = 1;
      break;
    case 'w':
      opts.min_words = atoi(optarg);
      break;
    case 'g':
      opts.gather = 1;
      break;
    case 's':
      opts.support = 0;
      break;
    case 'L':
      opts.node_limit = strtoull(optarg, NULL, 10);
      break;
    default:
      return 2;
    }
  }
  if (optind >= argc) {
    fprintf(stderr, "usage: bench [--repeat K] [--only SUBSTR] [--no-fc] "
                    "[--no-mrv] [--no-support] [--update] [--min-words K] "
                    "[--gather] [--node-limit X] instances.txt\n");
    return 2;
  }
  FILE *fp = fopen(argv[optind], "r");
  if (!fp) {
    perror(argv[optind]);
    return 2;
  }

  char line[1024];
  int failures = 0, ran = 0;
  double total_time = 0;
  uint64_t total_nodes = 0;
  if (!update)
    printf("%-22s %5s %5s %12s %4s %5s %9s %8s\n", "instance", "N", "labels",
           "nodes", "sq", "check", "time(s)", "Mnodes/s");
  while (fgets(line, sizeof(line), fp)) {
    char *l = trim(line);
    if (!*l || *l == '#')
      continue;
    char *fields[5];
    int nf = 0;
    char *save, *tok = strtok_r(l, "|", &save);
    while (tok && nf < 5) {
      fields[nf++] = trim(tok);
      tok = strtok_r(NULL, "|", &save);
    }
    if (nf < 3) {
      fprintf(stderr, "bad line\n");
      return 2;
    }
    int n = atoi(fields[0]);
    char *exps[ENUM_MAX_PRIMES];
    int ne = 0;
    char ebuf[256];
    snprintf(ebuf, sizeof(ebuf), "%s", fields[1]);
    for (char *e = strtok_r(ebuf, " ", &save); e && ne < ENUM_MAX_PRIMES;
         e = strtok_r(NULL, " ", &save))
      exps[ne++] = e;
    prime_exps_t p;
    if (pexp_parse(&p, ne, exps)) {
      fprintf(stderr, "bad exponents %s\n", fields[1]);
      return 2;
    }
    uint64_t S = strtoull(fields[2], NULL, 10);
    int have_count = nf > 3 && strcmp(fields[3], "-") != 0;
    int have_hash = nf > 4 && strcmp(fields[4], "-") != 0;
    uint64_t exp_count = have_count ? strtoull(fields[3], NULL, 10) : 0;
    uint64_t exp_hash = have_hash ? strtoull(fields[4], NULL, 16) : 0;

    char name[128], pstr[64];
    pexp_to_str(&p, pstr, ".");
    snprintf(name, sizeof(name), "%dx%d P=%s S=%lu", n, n, pstr,
             (unsigned long)S);
    if (only && !strstr(name, only))
      continue;

    vec_list_t all, red;
    vec_list_init(&all, n);
    vec_list_init(&red, n);
    enum_vectors(&p, n, S, S, &all);
    reduce_vectors(&all, 0, all.count, &red, 0);

    double best = 1e30;
    search_stats_t st = {0};
    hash_ctx_t h = {0, 0};
    for (int r = 0; r < repeat; r++) {
      h.hash = h.count = 0;
      st = search_vectors(&red, 0, red.count, &opts, hash_cb, &h);
      if (st.seconds < best)
        best = st.seconds;
    }
    int ok = (!have_count || exp_count == h.count) &&
             (!have_hash || exp_hash == h.hash);
    const char *check = !ok ? "FAIL" : (have_count || have_hash) ? "ok" : "-";
    failures += !ok;
    ran++;
    total_time += best;
    total_nodes += st.nodes;
    if (update)
      printf("%d | %s | %lu | %lu | %016lx\n", n, fields[1], (unsigned long)S,
             (unsigned long)h.count, (unsigned long)h.hash);
    else
      printf("%-22s %5zu %5d %12lu %4lu %5s %9.4f %8.2f\n", name, red.count,
             st.num_labels, (unsigned long)st.nodes, (unsigned long)h.count,
             check, best, st.nodes / best / 1e6);
    fflush(stdout);
    vec_list_free(&all);
    vec_list_free(&red);
  }
  fclose(fp);
  if (!update) {
    printf("%-22s %5s %5s %12lu %4s %5s %9.4f %8.2f\n", "TOTAL", "", "",
           (unsigned long)total_nodes, "", failures ? "FAIL" : "ok",
           total_time, total_nodes / total_time / 1e6);
    if (failures)
      printf("%d of %d instances FAILED\n", failures, ran);
  }
  return failures ? 1 : 0;
}
