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
 * Diagonal-first search (src/c/dfirst.h; research/ideas.md, "Diagonal-first
 * search in msearch"):
 *   --diag-first      search the sums with at least --diag-first-min-n
 *                     vectors (after reduction) diagonal-first: for every
 *                     vector d of the sum (unreduced), the semi-magic
 *                     search on V_d = {v : |v & d| = 1} finds every
 *                     (square, SP traversal) pair once, so every magic
 *                     square (twice, once per diagonal); the other sums
 *                     get the plain search
 *   --diag-first-min-n N0   default 5000, where d-first and the plain
 *                     search cost the same per sum with the integrated
 *                     search (1.32x at N = 4.1k, 0.78x at 5.9k, 0.63x at
 *                     6.7k; 0.3x at N = 15-23k on fdb77fc)
 *   --d-stride k, --d-offset o   search only d = lo + o, lo + o + k, ...
 *                     (an unbiased sample; o defaults to a random offset
 *                     from --sample-seed)
 *   --d-range lo:hi   only the d with index in [lo, hi) of the sum's
 *                     unreduced list (for splitting a sum into units); the
 *                     plain sums (below --diag-first-min-n) are searched
 *                     only by the unit with lo = 0, the others write a
 *                     "skip" record for them
 *   --d-chunk C       a "dchunk" checkpoint record every C indices of d
 *                     (default 256)
 *   --d-log FILE      one line per d (see dfirst_search)
 *   --d-plain-root    search all of V_d's first rows, not only those through
 *                     d's rarest number (as bin/dsearch; see
 *                     dfirst_set_top_root)
 *   --calib-r1-stride k   for each d-first sum, also run the plain search
 *                     on every k-th first row r1 (random offset), which
 *                     estimates the sum's semi-magic squares and plain time
 *
 * r1 sampling of the plain search (measurements; the squares and the sum
 * record become "csquare" / "csum" records with estimates):
 *   --r1-stride k, --r1-offset o, --r1-strata k1,k2,..   (see search_opts_t;
 *                     strata: equal index ranges of the root list)
 *   --r1-log FILE     one line per sampled r1
 *   --sample-seed X   seed of the default random offsets (default 1)
 *   (msearch ignores bench's SAMPLE_* environment variables: only these
 *   options sample)
 *
 * Output (JSON lines, one record per searched sum, flushed immediately so a
 * killed run keeps all completed sums):
 *   {"type":"sum","n":6,"P":[13,6,3,2],"Pval":...,"S":506,"nvecs":831,
 *    "nvecs_raw":831,"labels":77,"nodes":890845,"squares":0,
 *    "time":0.23,"setup_time":..,"enum_time":0.001,"cpu":0.24,"truncated":0,
 *    "engine":3}
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
 *
 * d-first sums write, instead of "square" and "sum":
 *   {"type":"dsquare",...,"d":i,"dvec":[...],"set_count":..,"sp_count":..,
 *    "best_score":..,"magic":0|1,"partner":0|1,"hash":..,"grid":..}
 * per (square, d) pair (a square appears once per SP traversal: dedupe by
 * hash; set_count = its traversals among the sum's vectors = sp_count;
 * partner: another of them fits on the other diagonal with d, i.e. magic),
 *   {"type":"dchunk","S":..,"d_lo":..,"d_hi":..,"nd":..,"pairs":..,...}
 * per completed chunk of d, and a final
 *   {"type":"dsum","mode":"dfirst",...,"nd":..,"d_stride":..,"pairs":..,
 *    "est_pairs":..,"time":..,"est_time":..,...,"cpu":..,"complete":0|1}
 * where pairs is the number of (square, SP traversal) pairs (what a magic
 * square needs twice), cpu the process CPU seconds as in a "sum" record
 * (the d loop, the reduction and the enumeration share; time: the d loop
 * only), and complete says that every d of the sum was
 * searched (not a --d-range part, a --d-stride sample, or a truncated
 * run). With --diag-first the "done" record has "mode":"dfirst" (and
 * diag_first_min_n): its range then also holds the d-first sums, which
 * have no "sum" record but are not empty, and with --d-range lo:hi, lo > 0,
 * the plain sums of the range get a "skip" record instead of a search
 * (the unit whose range starts at 0 searches them). Sampled plain searches
 * (--calib-r1-stride, --r1-*) write "csquare" records (the squares of the
 * sampled r1) and a "csum" record with est_squares, se_squares, est_time,
 * se_time (se_strata_missing: the strata with fewer than 2 sampled r1, left
 * out of the se_*; cpu: as in a "sum" record for --r1-*, the sampled search
 * only for --calib-r1-stride, whose "dsum" record has the rest). A run that
 * samples its plain sums (--r1-*) writes "r1_sample":1 (and
 * "mode":"sampled" without --diag-first) in its "done" record: its range
 * was not searched in full, and the scheduler skips the record.
 */
/* clock_gettime and the CPU-time clocks also under a strict -std=c17 */
#define _POSIX_C_SOURCE 200809L
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "enumerate.h"
#include "arrange.h"
#include "square.h"
#include "dfirst.h"

/* version of the search, written into every "sum" record: bump it whenever
 * the search time changes materially, since the scheduler fits its model of
 * the time per sum only to records of the newest version (1 = the first
 * version of this pipeline, without the field; 2 = support and cross
 * filters, carried bitsets, October 2026; 3 = per-r1 widths, carried
 * bitsets up to 512 labels, the pretest (cx/integrated, October 2026):
 * 0.66-0.89x engine 2's CPU at 129-256 labels, 0.26x above; keep
 * scripts/scheduler.py ENGINE in step) */
#define ENGINE_VERSION 3

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
  /* sampled plain search: "csquare" records with this weight, and the
   * traversal totals of its squares */
  int sampled;
  double weight;
  uint64_t sp_pairs, s_trav, p_trav;
  /* d-first: pairs with a magic best score */
  uint64_t magic_pairs;
} out_ctx_t;

static uint64_t splitmix(uint64_t x) {
  uint64_t z = x + 0x9E3779B97F4A7C15ull;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

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
    ctx->sp_pairs += ds.sp_count;
    ctx->s_trav += ds.s_count;
    ctx->p_trav += ds.p_count;
    fprintf(fp, "{\"type\":\"%s\",\"n\":%d,\"P\":",
            ctx->sampled ? "csquare" : "square", sq->n);
    print_p_json(fp, ctx->p);
    if (ctx->sampled)
      fprintf(fp, ",\"weight\":%g", ctx->weight);
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

static int dsquare_found(const dsquare_t *d, void *vctx) {
  out_ctx_t *ctx = vctx;
  FILE *fp = ctx->out;
  const square_t *sq = d->sq;
  diag_stats_t ds;
  square_diag_stats(sq, &ds);
  ctx->magic_pairs += ds.best_score >= 14;
  fprintf(fp, "{\"type\":\"dsquare\",\"n\":%d,\"P\":", sq->n);
  print_p_json(fp, ctx->p);
  fprintf(fp, ",\"S\":%lu,\"d\":%zu,\"dvec\":[", (unsigned long)ctx->S,
          d->d_index);
  for (int i = 0; i < sq->n; i++)
    fprintf(fp, i ? ",%lu" : "%lu", (unsigned long)d->d[i]);
  fprintf(fp,
          "],\"set_count\":%d,\"s_count\":%d,\"p_count\":%d,"
          "\"sp_count\":%d,\"best_score\":%d,\"magic\":%d,\"partner\":%d,"
          "\"hash\":\"%016lx\",\"grid\":",
          d->set_count, ds.s_count, ds.p_count, ds.sp_count, ds.best_score,
          ds.best_score >= 14, d->partner, (unsigned long)square_hash(sq));
  grid_print_json(fp, &ds.best);
  fprintf(fp, "}\n");
  if (ds.best_score >= 14 || d->partner)
    fprintf(stderr, "!!! MAGIC SQUARE FOUND (S=%lu, d-first) !!!\n",
            (unsigned long)ctx->S);
  fflush(fp);
  return 0;
}

/* a whole decimal argument (no trailing junk), or exit 2 */
static uint64_t arg_u64(const char *s, const char *name) {
  char *e;
  const uint64_t v = strtoull(s, &e, 10);
  if (!*s || *s == '-' || *e != '\0') {
    fprintf(stderr, "msearch: bad value '%s' for --%s\n", s, name);
    exit(2);
  }
  return v;
}

static double arg_double(const char *s, const char *name) {
  char *e;
  double v = strtod(s, &e);
  if (e == s || *e != '\0') {
    fprintf(stderr, "msearch: bad value '%s' for --%s\n", s, name);
    exit(2);
  }
  return v;
}

/* comma-separated unsigned list; returns the count, -1 if malformed or
 * longer than max */
static int parse_u32_list(const char *s, uint32_t *out, int max) {
  int k = 0;
  while (*s) {
    char *end;
    if (k == max || *s == '-')
      return -1;
    out[k++] = (uint32_t)strtoul(s, &end, 10);
    if (end == s || (*end && *end != ','))
      return -1;
    s = end;
    if (*s == ',')
      s++;
  }
  return k;
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


/* a sampled plain search of one sum: "csquare" records and a "csum" */
static uint64_t run_sampled(FILE *out, out_ctx_t *ctx, const char *pstr, int n,
                            const prime_exps_t *p, uint64_t S, size_t raw,
                            const vec_list_t *red, const search_opts_t *so,
                            const char *mode, double enum_share,
                            double reduce_time, double pre_cpu) {
  ctx->S = S;
  ctx->sampled = 1;
  /* the weight of a square: the stride of its r1 (with strata, see the
   * r1 log: 0 here) */
  ctx->weight = so->r1_nstrata > 0 ? 0 : (so->r1_stride > 1 ? so->r1_stride : 1);
  ctx->sp_pairs = ctx->s_trav = ctx->p_trav = 0;
  const double c0 = cpu_time();
  search_stats_t st =
      search_vectors(red, 0, red->count, so, square_found, ctx);
  const double cpu = cpu_time() - c0;
  ctx->sampled = 0;
  fprintf(out,
          "{\"type\":\"csum\",\"mode\":\"%s\",\"n\":%d,\"P\":[%s],\"Pval\":%lu,"
          "\"S\":%lu,\"nvecs\":%zu,\"nvecs_raw\":%zu,\"labels\":%d,",
          mode, n, pstr, (unsigned long)pexp_value(p), (unsigned long)S,
          red->count, raw, st.num_labels);
  if (so->r1_nstrata > 0) {
    fprintf(out, "\"r1_strata\":[");
    for (int h = 0; h < so->r1_nstrata; h++)
      fprintf(out, h ? ",%u" : "%u", so->r1_sstride[h]);
    fprintf(out, "],\"r1_offsets\":[");
    for (int h = 0; h < so->r1_nstrata; h++)
      fprintf(out, h ? ",%u" : "%u", so->r1_soffset[h]);
    fprintf(out, "],");
  } else {
    fprintf(out, "\"r1_stride\":%u,\"r1_offset\":%u,",
            so->r1_stride > 1 ? so->r1_stride : 1, so->r1_offset);
  }
  fprintf(out,
          "\"r1_sampled\":%lu,\"squares\":%lu,\"sp_pairs\":%lu,\"s_trav\":%lu,"
          "\"p_trav\":%lu,\"nodes\":%lu,\"time\":%.6f,\"cpu\":%.6f,"
          "\"setup_time\":%.6f,\"est_squares\":%.6g,\"se_squares\":%.6g,"
          "\"est_sp_pairs\":%.6g,\"est_nodes\":%.6g,\"est_time\":%.6g,"
          "\"se_time\":%.6g,\"se_strata_missing\":%d,\"enum_time\":%.6f,"
          "\"reduce_time\":%.6f,\"truncated\":%d,\"engine\":%d}\n",
          (unsigned long)st.r1_sampled, (unsigned long)st.squares,
          (unsigned long)ctx->sp_pairs, (unsigned long)ctx->s_trav,
          (unsigned long)ctx->p_trav, (unsigned long)st.nodes, st.seconds,
          cpu + pre_cpu, st.setup_seconds, st.est_squares, st.se_squares,
          ctx->weight > 0 ? ctx->weight * (double)ctx->sp_pairs : -1.0,
          st.est_nodes, st.est_seconds, st.se_seconds, st.r1_strata_nose,
          enum_share, reduce_time, st.truncated, ENGINE_VERSION);
  fflush(out);
  return st.nodes;
}

/* the d-first search of one sum (see the header comment) */
static uint64_t run_dfirst(FILE *out, out_ctx_t *ctx, const char *pstr, int n,
                           const prime_exps_t *p, uint64_t S,
                           const vec_list_t *all, size_t start_i, size_t raw,
                           const vec_list_t *red, const search_opts_t *opts,
                           size_t d_stride, int64_t d_offset, uint64_t seed,
                           size_t d_lo, size_t d_hi, size_t d_chunk,
                           FILE *d_log, uint32_t calib_stride, int top_root,
                           double enum_share, double reduce_time,
                           double pre_cpu, double deadline) {
  const double c0 = cpu_time();
  dfirst_t *df = dfirst_new(red, 0, red->count, all, start_i, raw);
  dfirst_set_top_root(df, top_root);
  const double t_index = cpu_time() - c0;
  const uint64_t sd = splitmix(seed ^ (S * 0x9E3779B97F4A7C15ull));
  const size_t off = d_offset >= 0 ? (size_t)d_offset % d_stride
                                   : (size_t)(splitmix(sd + 101) % d_stride);
  const size_t hi = d_hi < raw ? d_hi : raw, lo = d_lo < hi ? d_lo : hi;
  ctx->S = S;
  ctx->magic_pairs = 0;
  if (d_log)
    fprintf(d_log, "# S %lu N %zu Nred %zu stride %zu offset %zu\n",
            (unsigned long)S, raw, red->count, d_stride, off);
  dfirst_stats_t tot;
  memset(&tot, 0, sizeof(tot));
  const size_t base = lo + off;
  int stopped = 0;
  for (size_t c = lo; c < hi; c += d_chunk) {
    const size_t c2 = c + d_chunk < hi ? c + d_chunk : hi;
    /* the first sampled index >= c */
    const size_t i0 =
        c <= base ? base : base + (c - base + d_stride - 1) / d_stride * d_stride;
    dfirst_stats_t st;
    memset(&st, 0, sizeof(st));
    const double cc = cpu_time();
    if (i0 < c2)
      st = dfirst_search(df, i0, c2, d_stride, 0, opts, dsquare_found, ctx,
                         d_log);
    tot.nd += st.nd;
    tot.vd_total += st.vd_total;
    tot.nodes += st.nodes;
    tot.pairs += st.pairs;
    tot.partners += st.partners;
    tot.vd_seconds += st.vd_seconds;
    tot.setup_seconds += st.setup_seconds;
    tot.search_seconds += st.search_seconds;
    tot.d_cpu += st.d_cpu;
    tot.d_cpu2 += st.d_cpu2;
    tot.cpu_seconds += st.cpu_seconds;
    tot.d_pairs2 += st.d_pairs2;
    tot.truncated |= st.truncated;
    fprintf(out,
            "{\"type\":\"dchunk\",\"n\":%d,\"P\":[%s],\"S\":%lu,\"d_lo\":%zu,"
            "\"d_hi\":%zu,\"d_stride\":%zu,\"d_offset\":%zu,\"nd\":%lu,"
            "\"nodes\":%lu,\"pairs\":%lu,\"partners\":%lu,\"time\":%.6f}\n",
            n, pstr, (unsigned long)S, c, c2, d_stride, off,
            (unsigned long)st.nd, (unsigned long)st.nodes,
            (unsigned long)st.pairs, (unsigned long)st.partners,
            cpu_time() - cc);
    fflush(out);
    /* --time-limit: also between the chunks of a d-first sum, which can
     * take hours (the sum is then not complete) */
    if (st.stopped || (deadline > 0 && c2 < hi && wall_time() >= deadline)) {
      stopped = 1;
      break;
    }
  }
  /* every d of the sum searched: the sum's (square, SP traversal) pairs,
   * hence its magic squares, are all found (a d-range unit, a d sample,
   * a truncated or stopped run is not complete) */
  const int complete = lo == 0 && hi == raw && d_stride == 1 &&
                       !tot.truncated && !stopped;
  const uint32_t labels = dfirst_num_labels(df);
  dfirst_free(df);
  const double cpu = cpu_time() - c0;
  /* estimates for the d in [lo, hi): stride x the sampled totals, and the
   * simple-random-sampling standard errors over the sampled d */
  const double k = (double)d_stride, m = (double)tot.nd;
  const double Nd = (double)(hi - lo);
  double se_time = 0, se_pairs = 0;
  if (m >= 2) {
    const double f = Nd * Nd * (1 - m / Nd) / m / (m - 1);
    const double ct = tot.d_cpu; /* the sum over d of the per-d CPU */
    const double vt = f * (tot.d_cpu2 - ct * ct / m);
    const double vp =
        f * (tot.d_pairs2 - (double)tot.pairs * (double)tot.pairs / m);
    se_time = vt > 0 ? sqrt(vt) : 0;
    se_pairs = vp > 0 ? sqrt(vp) : 0;
  }
  fprintf(out,
          "{\"type\":\"dsum\",\"mode\":\"dfirst\",\"n\":%d,\"P\":[%s],"
          "\"Pval\":%lu,\"S\":%lu,"
          "\"nvecs\":%zu,\"nvecs_raw\":%zu,\"labels\":%u,\"d_lo\":%zu,"
          "\"d_hi\":%zu,\"d_stride\":%zu,\"d_offset\":%zu,\"nd\":%lu,"
          "\"avg_vd\":%.1f,\"nodes\":%lu,\"pairs\":%lu,\"partners\":%lu,"
          "\"magic_pairs\":%lu,\"est_pairs\":%.6g,\"se_pairs\":%.6g,"
          "\"est_nodes\":%.6g,\"time\":%.6f,\"index_time\":%.6f,"
          "\"vd_time\":%.6f,\"setup_time\":%.6f,\"search_time\":%.6f,"
          "\"est_time\":%.6g,\"se_time\":%.6g,\"enum_time\":%.6f,"
          "\"reduce_time\":%.6f,\"cpu\":%.6f,\"truncated\":%d,"
          "\"complete\":%d,\"engine\":%d}\n",
          n, pstr, (unsigned long)pexp_value(p), (unsigned long)S, red->count,
          raw, labels, lo, hi, d_stride, off, (unsigned long)tot.nd,
          tot.nd ? (double)tot.vd_total / (double)tot.nd : 0.0,
          (unsigned long)tot.nodes, (unsigned long)tot.pairs,
          (unsigned long)tot.partners, (unsigned long)ctx->magic_pairs,
          k * (double)tot.pairs, se_pairs, k * (double)tot.nodes, cpu, t_index,
          tot.vd_seconds, tot.setup_seconds, tot.search_seconds,
          t_index + k * (cpu - t_index), se_time, enum_share, reduce_time,
          cpu + pre_cpu, tot.truncated, complete, ENGINE_VERSION);
  fflush(out);
  uint64_t nodes = tot.nodes;
  if (calib_stride > 0) {
    search_opts_t so = *opts;
    so.r1_stride = calib_stride;
    so.r1_offset = (uint32_t)(splitmix(sd + 202) % calib_stride);
    /* (its "cpu" is the sampled search's only: the reduction and the
     * enumeration are in the dsum record's) */
    nodes += run_sampled(out, ctx, pstr, n, p, S, raw, red, &so, "calib",
                         enum_share, reduce_time, 0.0);
  }
  return nodes;
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
  /* d-first */
  int diag_first = 0;
  size_t dfirst_min_n = 5000, d_stride = 1, d_chunk = 256;
  int64_t d_offset = -1; /* -1: random from the seed */
  size_t d_range_lo = 0, d_range_hi = (size_t)-1;
  const char *d_log_file = NULL;
  uint32_t calib_stride = 0;
  int d_top_root = 1; /* see dfirst_set_top_root */
  const char *dfirst_opt = NULL; /* a d-first option given (for the check) */
  /* r1 sampling of the plain search */
  uint32_t r1_stride = 0, r1_strata[8];
  int64_t r1_offset = -1;
  int r1_nstrata = 0;
  const char *r1_log_file = NULL;
  uint64_t seed = 1;
  /* (no SAMPLE_* environment as in bench: a variable left exported would
   * silently turn every scheduler unit into a sampled run) */

  static struct option long_options[] = {
      {"diag-first", no_argument, 0, 1000},
      {"diag-first-min-n", required_argument, 0, 1001},
      {"d-stride", required_argument, 0, 1002},
      {"d-offset", required_argument, 0, 1003},
      {"d-range", required_argument, 0, 1004},
      {"d-chunk", required_argument, 0, 1005},
      {"d-log", required_argument, 0, 1006},
      {"calib-r1-stride", required_argument, 0, 1007},
      {"r1-stride", required_argument, 0, 1008},
      {"r1-offset", required_argument, 0, 1009},
      {"r1-strata", required_argument, 0, 1010},
      {"r1-log", required_argument, 0, 1011},
      {"sample-seed", required_argument, 0, 1012},
      {"d-plain-root", no_argument, 0, 1013},
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
      n = (int)arg_u64(optarg, "vec-size");
      break;
    case 'a':
      min_sum = arg_u64(optarg, "min-sum");
      have_min = 1;
      break;
    case 'b':
      max_sum = arg_u64(optarg, "max-sum");
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
      opts.node_limit = arg_u64(optarg, "node-limit");
      break;
    case 't':
      total_node_limit = arg_u64(optarg, "total-nodes");
      break;
    case 'T':
      time_limit = arg_double(optarg, "time-limit");
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
      opts.pretest_min = (int)arg_u64(optarg, "pretest-min");
      break;
    case 1000:
      diag_first = 1;
      break;
    case 1001:
      dfirst_min_n = arg_u64(optarg, "diag-first-min-n");
      break;
    case 1002:
      d_stride = arg_u64(optarg, "d-stride");
      if (d_stride < 1)
        d_stride = 1;
      dfirst_opt = "d-stride";
      break;
    case 1003:
      d_offset = (int64_t)arg_u64(optarg, "d-offset");
      dfirst_opt = "d-offset";
      break;
    case 1004: {
      /* lo:hi with lo < hi, or lo: (to the end) */
      char *c = strchr(optarg, ':');
      if (!c) {
        fprintf(stderr, "msearch: --d-range takes lo:hi or lo:\n");
        return 2;
      }
      *c = '\0';
      d_range_lo = arg_u64(optarg, "d-range");
      d_range_hi = c[1] ? arg_u64(c + 1, "d-range") : (size_t)-1;
      if (d_range_lo >= d_range_hi) {
        fprintf(stderr, "msearch: --d-range needs lo < hi\n");
        return 2;
      }
      dfirst_opt = "d-range";
      break;
    }
    case 1005:
      d_chunk = arg_u64(optarg, "d-chunk");
      if (d_chunk < 1)
        d_chunk = 1;
      dfirst_opt = "d-chunk";
      break;
    case 1006:
      d_log_file = optarg;
      dfirst_opt = "d-log";
      break;
    case 1007:
      calib_stride = (uint32_t)arg_u64(optarg, "calib-r1-stride");
      dfirst_opt = "calib-r1-stride";
      break;
    case 1008:
      r1_stride = (uint32_t)arg_u64(optarg, "r1-stride");
      break;
    case 1009:
      r1_offset = (int64_t)arg_u64(optarg, "r1-offset");
      break;
    case 1010:
      r1_nstrata = parse_u32_list(optarg, r1_strata, 8);
      if (r1_nstrata <= 0) {
        fprintf(stderr, "msearch: bad --r1-strata (1 to 8 strides)\n");
        return 2;
      }
      break;
    case 1011:
      r1_log_file = optarg;
      break;
    case 1012:
      seed = arg_u64(optarg, "sample-seed");
      break;
    case 1013:
      d_top_root = 0;
      dfirst_opt = "d-plain-root";
      break;
    default:
      return 2;
    }
  }
  prime_exps_t p;
  if (n < 3 || n > SQ_MAX_N || optind >= argc ||
      pexp_parse(&p, argc - optind, argv + optind) || pexp_value(&p) == 0) {
    fprintf(stderr,
            "usage: msearch [--vec-size N] [--min-sum S] [--max-sum S] "
            "[--sums S1,S2,..] [--node-limit X] [--total-nodes X] "
            "[--time-limit T] [--out FILE] [--format json|legacy] "
            "[--reduce strong] [--no-fc] [--no-mrv] [--no-support] "
            "[--no-cross] [--pretest-min K]\n"
            "  [--diag-first [--diag-first-min-n N0] [--d-stride k] "
            "[--d-offset o] [--d-range lo:hi] [--d-chunk C] [--d-log FILE] "
            "[--d-plain-root] [--calib-r1-stride k]]\n"
            "  [--r1-stride k] [--r1-offset o] [--r1-strata k1,k2,..] "
            "[--r1-log FILE] [--sample-seed X]\n"
            "  e1 e2 ...   (P = 2^e1 3^e2 5^e3 ...; see the header of "
            "src/c/msearch.c)\n");
    return 2;
  }
  if (dfirst_opt && !diag_first) {
    fprintf(stderr, "msearch: --%s needs --diag-first\n", dfirst_opt);
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
  FILE *d_log = d_log_file ? fopen(d_log_file, "a") : NULL;
  FILE *r1_log = r1_log_file ? fopen(r1_log_file, "a") : NULL;
  if ((d_log_file && !d_log) || (r1_log_file && !r1_log)) {
    perror(d_log_file && !d_log ? d_log_file : r1_log_file);
    return 2;
  }
  const int plain_sampled = r1_stride > 1 || r1_nstrata > 0 || r1_log;
  if (legacy && (diag_first || plain_sampled)) {
    fprintf(stderr, "--format legacy: no d-first or sampling\n");
    return 2;
  }
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
      char pstr[64];
      pexp_to_str(&p, pstr, ",");
      /* this sum's share of the enumeration, and the process CPU seconds of
       * its reduction and enumeration share (in every record's "cpu") */
      const double share = (double)raw / (all.count ? all.count : 1);
      const double pre_cpu = cpu_time() - cr + enum_cpu * share;
      if (diag_first && red.count < dfirst_min_n &&
          (d_range_lo > 0 || (d_offset >= 0 && (size_t)d_offset % d_stride != 0))) {
        /* a plain sum in a --d-range or --d-offset unit: the units of a
         * range (or of the offsets of a stride) split the d of its d-first
         * sums between them, but a plain sum is searched whole, so only by
         * the unit whose d-range starts at 0 (and whose offset is 0), not
         * once per unit */
        fprintf(out,
                "{\"type\":\"skip\",\"mode\":\"dfirst\",\"n\":%d,\"P\":[%s],"
                "\"S\":%lu,\"nvecs\":%zu,\"nvecs_raw\":%zu,\"d_lo\":%zu,"
                "\"d_offset\":%ld}\n",
                n, pstr, (unsigned long)S, red.count, raw, d_range_lo,
                (long)(d_offset >= 0 ? (size_t)d_offset % d_stride : 0));
        fflush(out);
        continue;
      }
      if (diag_first && red.count >= dfirst_min_n) {
        total_nodes += run_dfirst(out, &ctx, pstr, n, &p, S, &all, start_i,
                                  raw, &red, &opts, d_stride, d_offset, seed,
                                  d_range_lo, d_range_hi, d_chunk, d_log,
                                  calib_stride, d_top_root, enum_time * share,
                                  reduce_time, pre_cpu,
                                  time_limit > 0 ? t_start + time_limit : 0);
        if (total_node_limit && total_nodes >= total_node_limit)
          stop = stop_nodes = 1;
        if (time_limit > 0 && wall_time() - t_start >= time_limit)
          stop = 1;
        continue;
      }
      if (plain_sampled) {
        search_opts_t so = opts;
        uint64_t sd = splitmix(seed ^ (S * 0x9E3779B97F4A7C15ull));
        if (r1_nstrata > 0) {
          so.r1_nstrata = r1_nstrata;
          for (int h = 0; h < r1_nstrata; h++) {
            so.r1_sstride[h] = r1_strata[h] ? r1_strata[h] : 1;
            so.r1_soffset[h] = (uint32_t)(splitmix(sd + h) % so.r1_sstride[h]);
          }
        } else {
          so.r1_stride = r1_stride > 1 ? r1_stride : 1;
          so.r1_offset = r1_offset >= 0 ? (uint32_t)r1_offset
                                        : (uint32_t)(sd % so.r1_stride);
        }
        so.r1_log = r1_log;
        if (r1_log)
          fprintf(r1_log, "# S %lu N %zu\n", (unsigned long)S, red.count);
        total_nodes += run_sampled(out, &ctx, pstr, n, &p, S, raw, &red, &so,
                                   "sampled", enum_time * share, reduce_time,
                                   pre_cpu);
        if (total_node_limit && total_nodes >= total_node_limit)
          stop = stop_nodes = 1;
        if (time_limit > 0 && wall_time() - t_start >= time_limit)
          stop = 1;
        continue;
      }
      search_stats_t st =
          search_vectors(&red, 0, red.count, &opts, square_found, &ctx);
      total_nodes += st.nodes;
      double cpu = cpu_time() - cr + enum_cpu * share;
      if (legacy) {
        fprintf(out, "num searched: %lu\n", (unsigned long)st.nodes);
      } else {
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
            "\"last_sum\":%lu,\"complete\":%d,\"time\":%.3f",
            n, pstr, (unsigned long)min_sum, (unsigned long)last_sum,
            !stop || last_sum >= max_sum, wall_time() - t_start);
    /* with --diag-first the range is not all plain sums: the sums with a
     * "dsum" record were searched d-first (all of their d only if its
     * "complete" is 1), so a reader must not take them for searched plain
     * sums without squares */
    if (diag_first)
      fprintf(out, ",\"mode\":\"dfirst\",\"diag_first_min_n\":%zu",
              dfirst_min_n);
    else if (plain_sampled)
      fprintf(out, ",\"mode\":\"sampled\"");
    /* r1-sampled plain sums ("csum" records) were not searched in full:
     * a sampled run's done record says so, and the scheduler skips it (as
     * every record with "r1_sample") */
    if (plain_sampled)
      fprintf(out, ",\"r1_sample\":1");
    fprintf(out, "}\n");
  }
  if (legacy) {
    if (stop_nodes)
      fprintf(out, "terminated: total count %lu exceeds cutoff %lu\n",
              (unsigned long)total_nodes, (unsigned long)total_node_limit);
    fprintf(out, "completed in %.5f secs\n", wall_time() - t_start);
  }
  if (out != stdout)
    fclose(out);
  if (d_log)
    fclose(d_log);
  if (r1_log)
    fclose(r1_log);
  free(sum_list);
  return 0;
}
