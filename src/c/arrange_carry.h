/*
 * The carried path without AVX-512BW (ARM, AVX2, SSE; or -DCARRY_PORTABLE):
 * plain C on 64-bit words, included from arrange_core.h in three parts
 * (CARRY_PART 1: the filters, the counts and the selection of the children;
 * 2: cross support; 3: class support) in place of the AVX-512 kernels, for
 * each width W <= CARRY_MAX_W. The rest of the carried path (TRY_CHILD,
 * SUPPORT, CROSS, SEARCH_REC, SEARCH_ROOT) is shared.
 *
 * The same lists (word w of entry i at list[w * cap + i]), the same tests
 * and the same results as the AVX-512 kernels: the same entries kept, in the
 * same order, the same unions, the same label counts at the cells that are
 * read, the same keep masks of the class support, so the search makes the
 * same choices, visits the same nodes and finds the same squares in the same
 * order (checked: the node counts of bench and of fuzz_arrange -v against
 * the AVX-512 build). The AVX-512 kernels work on groups of 8 entries (one
 * zmm per word) with compresses; these work on blocks of up to 64 entries:
 * the tests of a block give a 64-bit mask of the kept entries, in a loop
 * over the entries without stores, which compilers vectorize (2 entries per
 * NEON register, 4 per AVX2 register), and then only the kept entries are
 * copied, by the set bits of the mask (most entries are dropped). No
 * padding is read (PAD_LIST is a no-op here). The exactly-one-bit test of
 * FILTER_CARRY needs no popcount (a vector instruction on ARM).
 *
 * The byte counters of COUNT_CARRY are bytes of the uint64_t words of
 * cnt[d][a] (byte i of words 8 w .. 8 w + 7 counts label 64 w + i, as with
 * AVX-512BW), read and written through uint8_t pointers, and only those of
 * the unmatched cells are computed: only those are read (by the choice of
 * the branching cell in SEARCH_REC, see mrv_eq in arrange.c).
 *
 * (CLASS_PROF is not implemented here; CHILD_PROF works on x86.)
 */
#if CARRY_PART == 1

/* the words V[w] of the test and their rotations R[w] (see FILTER_CARRY of
 * the AVX-512 code) */
static inline void CARRY_SETUP(const uint64_t *vb, const int once,
                               const uint64_t *excl, uint64_t *V, int *R) {
  uint64_t used = 0;
  for (int w = 0; w < W; w++) {
    uint64_t t = vb[w] | (!once && excl ? excl[w] : 0);
    int r = 0;
    if (once && w > 0)
      while (rotl64(t, r) & used)
        r++;
    t = rotl64(t, r);
    used |= t;
    V[w] = t;
    R[w] = r;
  }
}

/* the mask of the entries i < m (m <= 64) of in that pass the test of
 * FILTER_CARRY: |u & v| == 1 (once), as one exactly-one-bit test on the or
 * of the words of u & v rotated so that the words of v are disjoint (no
 * popcount, which is a vector instruction on ARM), or u disjoint from V (v
 * and excl) */
static inline __attribute__((always_inline)) uint64_t
CARRY_TEST(const uint64_t *in, size_t cap, uint32_t m, const int once,
           const uint64_t *V, const int *R) {
  KB_DECL;
  if (once) {
    for (uint32_t i = 0; i < m; i++) {
      uint64_t c = in[i] & V[0];
      for (int w = 1; w < W; w++)
        c |= rotl64v(in[w * cap + i], R[w]) & V[w];
      KB_SET(i, (c != 0) & ((c & (c - 1)) == 0));
    }
  } else {
    for (uint32_t i = 0; i < m; i++) {
      uint64_t t = 0;
      for (int w = 0; w < W; w++)
        t |= in[w * cap + i] & V[w];
      KB_SET(i, t == 0);
    }
  }
  return KB_MASK(m);
}

/* appends the entries of in with a bit in mask to out from position k
 * (also in place: out + k <= in), or-ing them into acc; returns the new
 * length */
static inline __attribute__((always_inline)) uint32_t
CARRY_COMPRESS(const uint64_t *in, size_t cap, uint64_t mask, uint64_t *out,
               uint32_t k, uint64_t *acc) {
  for (; mask; mask &= mask - 1) {
    const uint32_t j = (uint32_t)__builtin_ctzll(mask);
    for (int w = 0; w < W; w++) {
      const uint64_t x = in[w * cap + j];
      out[w * cap + k] = x;
      acc[w] |= x;
    }
    k++;
  }
  return k;
}

/* the block [i0, i0 + 64) of a list of cnt entries: its length */
#ifndef CARRY_BLEN
#define CARRY_BLEN(cnt, i0) ((cnt) - (i0) < 64 ? (cnt) - (i0) : 64u)
#endif

/*
 * FILTER_CARRY: out = the entries u of in[0..cnt) with |u & v| == once
 * (once = 0 or 1), and, for once = 0 and excl not NULL, also disjoint from
 * excl; un = the union of the kept entries; returns |out|.
 */
static inline uint32_t FILTER_CARRY(const sstate_t *s, const uint64_t *in,
                                    uint32_t cnt, const uint64_t *vb,
                                    const int once, const uint64_t *excl,
                                    uint64_t *out, uint64_t *un /* [W] */) {
  const size_t cap = s->cap;
  uint64_t V[W], acc[W];
  int R[W];
  CARRY_SETUP(vb, once, excl, V, R);
  for (int w = 0; w < W; w++)
    acc[w] = 0;
  uint32_t k = 0;
  for (uint32_t i0 = 0; i0 < cnt; i0 += 64) {
    const uint64_t mask =
        CARRY_TEST(in + i0, cap, CARRY_BLEN(cnt, i0), once, V, R);
    k = CARRY_COMPRESS(in + i0, cap, mask, out, k, acc);
  }
  for (int w = 0; w < W; w++)
    un[w] = acc[w];
  return k;
}

/* FILTER_CARRY without the output list (the pretest, see TRY_CHILD): the
 * number of entries kept, their union, and the keep masks km (bit j of
 * km[g]: entry 8 g + j kept; 0 past cnt), for FILTER_APPLY */
static inline uint32_t FILTER_COUNT(const sstate_t *s, const uint64_t *in,
                                    uint32_t cnt, const uint64_t *vb,
                                    const int once, const uint64_t *excl,
                                    uint8_t *km, uint64_t *un /* [W] */) {
  const size_t cap = s->cap;
  uint64_t V[W], acc[W];
  int R[W];
  CARRY_SETUP(vb, once, excl, V, R);
  for (int w = 0; w < W; w++)
    acc[w] = 0;
  uint32_t k = 0;
  for (uint32_t i0 = 0; i0 < cnt; i0 += 64) {
    const uint32_t m = CARRY_BLEN(cnt, i0);
    const uint64_t mask = CARRY_TEST(in + i0, cap, m, once, V, R);
    for (uint32_t g = 0; g < (m + 7) / 8; g++)
      km[(i0 >> 3) + g] = (uint8_t)(mask >> (8 * g));
    for (uint64_t t = mask; t; t &= t - 1) {
      const uint32_t j = (uint32_t)__builtin_ctzll(t);
      for (int w = 0; w < W; w++)
        acc[w] |= in[w * cap + i0 + j];
      k++;
    }
  }
  for (int w = 0; w < W; w++)
    un[w] = acc[w];
  return k;
}

/* the output of FILTER_CARRY from the keep masks km of FILTER_COUNT */
static inline void FILTER_APPLY(const sstate_t *s, const uint64_t *in,
                                uint32_t cnt, const uint8_t *km,
                                uint64_t *out) {
  const size_t cap = s->cap;
  uint64_t acc[W] = {0}; /* (unused) */
  uint32_t k = 0;
  for (uint32_t i0 = 0; i0 < cnt; i0 += 64) {
    const uint32_t m = CARRY_BLEN(cnt, i0);
    uint64_t mask = 0;
    for (uint32_t g = 0; g < (m + 7) / 8; g++)
      mask |= (uint64_t)km[(i0 >> 3) + g] << (8 * g);
    k = CARRY_COMPRESS(in + i0, cap, mask, out, k, acc);
  }
}

/* the lists of the child from the pretest's keep masks (see TRY_CHILD) */
static __attribute__((noinline)) void PT_APPLY(sstate_t *s, int d, int b,
                                               uint32_t min_v) {
  const int o = 1 - b;
  FILTER_APPLY(s, s->vw[d][o] + min_v, s->nvalid[d][o] - min_v, s->ptmask[0],
               s->vw[d + 1][o]);
  FILTER_APPLY(s, s->vw[d][b] + min_v, s->nvalid[d][b] - min_v, s->ptmask[1],
               s->vw[d + 1][b]);
}

/* (no padding: the loops here stop at the end of the list) */
static inline void PAD_LIST(const sstate_t *s, uint64_t *list, uint32_t cnt) {
  (void)s;
  (void)list;
  (void)cnt;
}

/* the support filter's pass: keep the entries of list[0..cnt) disjoint from
 * bad (in place), un = their union; returns their number */
static inline uint32_t KEEP_CARRY(const sstate_t *s, uint64_t *list,
                                  uint32_t cnt, const uint64_t *bad,
                                  uint64_t *un /* [W] */) {
  const size_t cap = s->cap;
  uint64_t acc[W];
  for (int w = 0; w < W; w++)
    acc[w] = 0;
  uint32_t k = 0;
  for (uint32_t i0 = 0; i0 < cnt; i0 += 64) {
    const uint32_t m = CARRY_BLEN(cnt, i0);
    KB_DECL;
    for (uint32_t i = 0; i < m; i++) {
      uint64_t t = 0;
      for (int w = 0; w < W; w++)
        t |= list[w * cap + i0 + i] & bad[w];
      KB_SET(i, t == 0);
    }
    k = CARRY_COMPRESS(list + i0, cap, KB_MASK(m), list, k, acc);
  }
  for (int w = 0; w < W; w++)
    un[w] = acc[w];
  return k;
}

/*
 * The byte counters g of the labels of need (the unmatched cells that the
 * list covers; the other bytes are 0, and never read) over list[0..cnt),
 * saturating at 255: one pass over the list per cell, a sum of one bit per
 * entry (a loop that compilers vectorize, without a branch or a chain of
 * increments through memory; a histogram of the cells of each entry, with
 * a counter per cell in memory, was 3x slower)
 */
static inline void COUNT_CARRY(const sstate_t *s, uint64_t *list,
                               uint32_t cnt, const uint64_t *need,
                               uint64_t *g /* [8 * W] */) {
  const size_t cap = s->cap;
  for (int w = 0; w < W; w++) {
    for (int j = 0; j < 8; j++)
      g[8 * w + j] = 0;
    uint8_t *gb = (uint8_t *)(g + 8 * w);
    const uint64_t *lw = list + w * cap;
    for (uint64_t m = need[w]; m; m &= m - 1) {
      const int b = __builtin_ctzll(m);
      uint64_t c = 0;
      for (uint32_t i = 0; i < cnt; i++)
        c += (lw[i] >> b) & 1;
      gb[b] = (uint8_t)(c < 255 ? c : 255);
    }
  }
}

/* the entries of in[0..cnt) with bit xb in word xw into out (the children
 * of a node, see SEARCH_REC); returns their number */
static inline uint32_t CARRY_SELECT(const sstate_t *s, const uint64_t *in,
                                    uint32_t cnt, int xw, uint64_t xb,
                                    uint64_t *out) {
  const size_t cap = s->cap;
  uint64_t acc[W] = {0}; /* (unused) */
  uint32_t k = 0;
  for (uint32_t i0 = 0; i0 < cnt; i0 += 64) {
    const uint32_t m = CARRY_BLEN(cnt, i0);
    KB_DECL;
    for (uint32_t i = 0; i < m; i++)
      KB_SET(i, (in[xw * cap + i0 + i] & xb) != 0);
    k = CARRY_COMPRESS(in + i0, cap, KB_MASK(m), out, k, acc);
  }
  return k;
}

#elif CARRY_PART == 2

/* the cross test of the entries of fl[0..kf) (in place, as KEEP_CARRY):
 * kept if disjoint from bad and meeting U_y (UT[w][j], word w of U_y of
 * cell j) for the NY cells j; acc gets their union; returns their number */
static inline __attribute__((always_inline)) uint32_t
CROSS_TEST(uint64_t *fl, size_t cap, uint32_t kf,
           uint64_t (*UT)[CROSS_MAXY], const int NY,
           const uint64_t *bad, uint64_t *acc) {
  uint32_t k = 0;
  for (uint32_t i0 = 0; i0 < kf; i0 += 64) {
    const uint32_t m = CARRY_BLEN(kf, i0);
    KB_DECL;
    for (uint32_t i = 0; i < m; i++) {
      uint64_t x[W], t = 0;
      for (int w = 0; w < W; w++) {
        x[w] = fl[w * cap + i0 + i];
        t |= x[w] & bad[w];
      }
      uint64_t miss = (uint64_t)(t != 0);
      for (int j = 0; j < NY; j++) {
        uint64_t u = 0;
        for (int w = 0; w < W; w++)
          u |= x[w] & UT[w][j];
        miss |= (uint64_t)(u == 0);
      }
      KB_SET(i, miss ^ 1);
    }
    k = CARRY_COMPRESS(fl + i0, cap, KB_MASK(m), fl, k, acc);
  }
  return k;
}

/*
 * One cross pass (see CROSS_AXIS of the AVX-512 code): drops the candidates
 * of axis f that miss U_y for some unmatched cell y of the placed vectors of
 * axis f, or that meet bad, in place; sets uni[d][f] and nvalid[d][f] and
 * returns the number of candidates left. As there, nothing is tested (and
 * kf returned) without such cells or with more than CROSS_MAXY of them.
 * U_y: every entry of axis h through y; an entry of h meets each placed
 * vector of f once, at such a cell, so it is or-ed into np[f] of the U_y.
 */
static inline uint32_t CROSS_AXIS(sstate_t *s, int d, int f,
                                  const uint64_t *bad) {
  const int h = 1 - f;
  const size_t cap = s->cap;
  const uint64_t *Cf = s->cells[d][f], *Ch = s->cells[d][h];
  uint64_t Y[W]; /* the cells y */
  int ny = 0;
  for (int w = 0; w < W; w++) {
    Y[w] = Cf[w] & ~Ch[w];
    ny += __builtin_popcountll(Y[w]);
  }
  const uint32_t kf = s->nvalid[d][f];
  if (ny == 0 || ny > CROSS_MAXY)
    return kf;
  /* the index of each cell (by word and bit) among the cells, and U_y by
   * word: UT[w][j] (cells j >= ny: all ones, which every entry meets) */
  uint8_t idx[W][64];
  uint64_t UT[W][CROSS_MAXY];
  {
    int j = 0;
    for (int w = 0; w < W; w++)
      for (uint64_t m = Y[w]; m; m &= m - 1)
        idx[w][__builtin_ctzll(m)] = (uint8_t)j++;
    for (int w = 0; w < W; w++)
      for (j = 0; j < CROSS_MAXY; j++)
        UT[w][j] = j < ny ? 0 : ~(uint64_t)0;
  }
  const uint64_t *hl = s->vw[d][h];
  const uint32_t kh = s->nvalid[d][h];
  for (uint32_t i = 0; i < kh; i++) {
    uint64_t x[W];
    for (int w = 0; w < W; w++)
      x[w] = hl[w * cap + i];
    for (int w = 0; w < W; w++)
      for (uint64_t m = x[w] & Y[w]; m; m &= m - 1) {
        const int j = idx[w][__builtin_ctzll(m)];
        for (int w2 = 0; w2 < W; w2++)
          UT[w2][j] |= x[w2];
      }
  }
  uint64_t *fl = s->vw[d][f];
  uint64_t acc[W];
  for (int w = 0; w < W; w++)
    acc[w] = 0;
  uint32_t k = 0;
  /* (the test over 8 cells, or 16, in fixed loops, which compilers
   * vectorize over the cells) */
  if (ny <= 8)
    k = CROSS_TEST(fl, cap, kf, UT, 8, bad, acc);
  else
    k = CROSS_TEST(fl, cap, kf, UT, CROSS_MAXY, bad, acc);
  for (int w = 0; w < W; w++)
    s->uni[d][f][w] = acc[w];
  s->nvalid[d][f] = k;
  return k;
}

#elif CARRY_PART == 3

/*
 * Class support (see CLASS_SUP of the AVX-512 code, whose registers of 8
 * lanes, word w in lane w, are arrays of W words here): the same passes,
 * the same tests per entry, the same keep masks (a bit per entry in
 * s->clskm), compacted once if the node lives. (No CLASS_PROF here.)
 */

/* the union of each block of list (G[c]), every entry alive */
static void CLS_UNIONS(const sstate_t *s, const uint64_t *list,
                       const uint32_t *b, uint64_t (*G)[W]) {
  const size_t cap = s->cap;
  for (int c = 0; c < s->ncls; c++)
    for (int w = 0; w < W; w++) {
      uint64_t acc = 0;
      for (uint32_t i = b[c]; i < b[c + 1]; i++)
        acc |= list[w * cap + i];
      G[c][w] = acc;
    }
}

/* drop the entries of list[0..cnt) outside the keep masks km, in place */
static void CLS_COMPACT(const sstate_t *s, uint64_t *list, uint32_t cnt,
                        const uint8_t *km, uint32_t kept) {
  const size_t cap = s->cap;
  uint32_t k = 0, i = 0;
  /* (the groups before the first drop stay where they are) */
  for (; i < cnt && km[i >> 3] == 0xff; i += 8)
    k += 8;
  for (; i < cnt; i++) {
    for (int w = 0; w < W; w++)
      list[w * cap + k] = list[w * cap + i];
    k += (uint32_t)(km[i >> 3] >> (i & 7)) & 1;
  }
#ifdef ARRANGE_DEBUG
  if (k != kept) {
    fprintf(stderr, "CLS_COMPACT: %u entries kept, %u counted\n", k, kept);
    abort();
  }
#else
  (void)kept;
#endif
}

/* the tests of block [b0, b1) of list (keep masks km), specialized for NT
 * classes to meet (gx) and HB (bad cells): acc gets the union of the
 * entries kept; returns the number of alive entries dropped. By groups of
 * 64 entries aligned in the list (one 64-bit word of km): the tests of the
 * group's entries, the alive ones (outside the block: none) or-ed into acc,
 * in a loop without stores (which compilers vectorize), then the keep
 * mask updated once */
static inline __attribute__((always_inline)) uint32_t
CLS_LOOP(const uint64_t *list, size_t cap, uint8_t *km, uint32_t b0,
         uint32_t b1, const uint64_t *bad, uint64_t (*gx)[W], const int NT,
         const int HB, uint64_t *acc) {
  for (int w = 0; w < W; w++)
    acc[w] = 0;
  uint32_t rem = 0;
  for (uint32_t g0 = b0 & ~63u; g0 < b1; g0 += 64) {
    const uint32_t lo = b0 > g0 ? b0 - g0 : 0;
    const uint32_t hi = b1 - g0 < 64 ? b1 - g0 : 64;
    const uint64_t lm = (hi == 64 ? ~(uint64_t)0 : ((uint64_t)1 << hi) - 1) &
                        ~(((uint64_t)1 << lo) - 1);
    const uint64_t old = km_load64(km + (g0 >> 3)), al = old & lm;
    if (!al)
      continue;
    const uint64_t *xl = list + g0;
    KB_DECL;
    for (uint32_t i = 0; i < hi; i++) {
      uint64_t pass = (al >> i) & 1;
      if (HB) {
        uint64_t t = 0;
        for (int w = 0; w < W; w++)
          t |= xl[w * cap + i] & bad[w];
        pass &= (uint64_t)(t == 0);
      }
      for (int q = 0; q < NT; q++) {
        uint64_t t = 0;
        for (int w = 0; w < W; w++)
          t |= xl[w * cap + i] & gx[q][w];
        pass &= (uint64_t)(t != 0);
      }
      const uint64_t m = (uint64_t)0 - pass;
      for (int w = 0; w < W; w++)
        acc[w] |= xl[w * cap + i] & m;
      KB_SET(i, pass);
    }
    const uint64_t drop = al & ~KB_MASK(hi);
    rem += (uint32_t)__builtin_popcountll(drop);
    km_store64(km + (g0 >> 3), old & ~drop);
  }
  return rem;
}

/* one pass over the list of axis a (see CLASS_PASS of the AVX-512 code) */
static uint32_t CLASS_PASS(sstate_t *s, int d, int a, int full, int uvalid,
                           int pcls_oa, const uint32_t *b,
                           uint64_t (*G)[SQ_MAX_N][W],
                           uint64_t (*Lost)[SQ_MAX_N][W]) {
  const int nc = s->ncls, oa = 1 - a;
  const size_t cap = s->cap;
  const uint64_t *list = s->vw[d][a];
  uint8_t *km = s->clskm[a];
  const uint64_t *Cz = s->cells[d][oa], *cmz = s->clsmask;
  /* per class c of oa: its cells to meet (GX), whether it has any, and the
   * unions of the other classes' G and Lost */
  uint64_t GX[SQ_MAX_N][W], LX[SQ_MAX_N][W], GO[SQ_MAX_N][W], LO[SQ_MAX_N][W];
  int anyc = 0;
  for (int w = 0; w < W; w++) {
    uint64_t pg = 0, pl = 0;
    for (int c = 0; c < nc; c++) {
      GO[c][w] = pg;
      LO[c][w] = pl;
      pg |= G[oa][c][w];
      pl |= Lost[oa][c][w];
    }
    pg = pl = 0;
    for (int c = nc - 1; c >= 0; c--) {
      GO[c][w] |= pg;
      LO[c][w] |= pl;
      pg |= G[oa][c][w];
      pl |= Lost[oa][c][w];
    }
  }
  for (int c = 0; c < nc; c++) {
    uint64_t any = 0;
    for (int w = 0; w < W; w++) {
      GX[c][w] = G[oa][c][w] & ~Cz[w];
      LX[c][w] = Lost[oa][c][w] & ~Cz[w];
      any |= GX[c][w];
    }
    anyc |= (any != 0) << c;
  }
  const int avail = ((1 << nc) - 1) & ~pcls_oa;
  uint32_t removed = 0;
  for (int j = 0; j < nc; j++) {
    const uint32_t b0 = b[j], b1 = b[j + 1];
    if (b0 == b1)
      continue;
    const int aj = avail & ~(1 << j);
    /* bad = ~(C | (G_j & cm) | G_other) */
    uint64_t bad[W], hbw = 0;
    for (int w = 0; w < W; w++)
      bad[w] = ~(Cz[w] | GO[j][w] | (G[oa][j][w] & cmz[w]));
    int tests = aj;
    if (uvalid) /* (only the bad cells the block's entries have) */
      for (int w = 0; w < W; w++)
        bad[w] &= G[a][j][w];
    if (!full) {
      /* the cells that became bad and the classes that lost cells, among
       * the cells of the block's entries U */
      const uint64_t *Uz = G[a][j];
      for (int w = 0; w < W; w++)
        bad[w] &= LO[j][w] | (Lost[oa][j][w] & cmz[w]);
      int lost = 0;
      for (int c = 0; c < nc; c++) {
        uint64_t t = 0;
        for (int w = 0; w < W; w++)
          t |= LX[c][w] & Uz[w];
        lost |= (t != 0) << c;
      }
      tests &= lost;
    }
    for (int w = 0; w < W; w++)
      hbw |= bad[w];
    const int hb = hbw != 0;
    uint32_t rem = 0;
    uint64_t gn[W], acc[W];
    for (int w = 0; w < W; w++)
      gn[w] = 0;
    if (aj & ~anyc) {
      /* some class to meet has no candidate: the block goes */
      for (uint32_t i = b0; i < b1; i++) {
        rem += (uint32_t)(km[i >> 3] >> (i & 7)) & 1;
        km[i >> 3] = (uint8_t)(km[i >> 3] & ~(1u << (i & 7)));
      }
    } else if (!full && !hb && !tests) {
      continue; /* nothing that block j's entries meet has changed */
    } else {
      uint64_t gx[SQ_MAX_N][W];
      int nt = 0;
      for (int c = 0; c < nc; c++)
        if (tests >> c & 1) {
          for (int w = 0; w < W; w++)
            gx[nt][w] = GX[c][w];
          nt++;
        }
#define CLS_CASE(k)                                                            \
  case k:                                                                      \
    rem = hb ? CLS_LOOP(list, cap, km, b0, b1, bad, gx, k, 1, acc)             \
             : CLS_LOOP(list, cap, km, b0, b1, bad, gx, k, 0, acc);            \
    break;
      switch (nt) {
        CLS_CASE(0)
        CLS_CASE(1)
        CLS_CASE(2)
        CLS_CASE(3)
        CLS_CASE(4)
        CLS_CASE(5)
        CLS_CASE(6)
      default:
        CLS_CASE(7)
      }
#undef CLS_CASE
      if (!full && !rem)
        continue;
      for (int w = 0; w < W; w++)
        gn[w] = acc[w];
    }
    for (int w = 0; w < W; w++) {
      Lost[a][j][w] |= G[a][j][w] & ~gn[w];
      G[a][j][w] = gn[w];
    }
    removed += rem;
  }
  return removed;
}

static __attribute__((noinline)) int CLASS_SUP(sstate_t *s, int d, int a) {
  const int n = s->n, nc = s->ncls;
  /* per axis: the class blocks, the entries alive, and per class the union
   * of the alive entries (G) and the labels it lost since the last pass
   * over the other axis (Lost) */
  uint64_t G[2][SQ_MAX_N][W], Lost[2][SQ_MAX_N][W];
  uint32_t bnd[2][SQ_MAX_N + 1], alive[2];
  int pcls[2] = {0, 0}, changed[2] = {0, 0};
  for (int ax = 0; ax < 2; ax++) {
    for (int q = 0; q < s->np[ax]; q++) {
      const uint64_t *pw = s->placedw[ax][q];
      for (int w = W - 1; w >= 0; w--)
        if (pw[w] & s->clsmask[w]) {
          pcls[ax] |= 1 << (int)(s->L - 1 -
                                 (64 * w + 63 -
                                  __builtin_clzll(pw[w] & s->clsmask[w])));
          break;
        }
    }
    const uint32_t cnt = s->nvalid[d][ax];
    CLS_BOUNDS(s, s->vw[d][ax], cnt, bnd[ax]);
    alive[ax] = cnt;
    uint8_t *km = s->clskm[ax];
    memset(km, 0xff, (cnt + 7) / 8);
    if (cnt & 7)
      km[cnt / 8] = (uint8_t)((1u << (cnt & 7)) - 1);
    /* (CLS_LOOP reads whole 64-bit words of km) */
    memset(km + (cnt + 7) / 8, 0, (cnt + 63) / 64 * 8 - (cnt + 7) / 8);
    for (int c = 0; c < nc; c++)
      for (int w = 0; w < W; w++)
        G[ax][c][w] = Lost[ax][c][w] = 0;
  }
  CLS_UNIONS(s, s->vw[d][1 - a], bnd[1 - a], G[1 - a]);
  const int a0 = a;
  int full[2] = {1, 1}, quiet = 0;
  while (quiet < 2) {
    const int oa = 1 - a;
    const uint32_t removed = CLASS_PASS(s, d, a, full[a], !full[a] || a != a0,
                                        pcls[oa], bnd[a], G, Lost);
    for (int c = 0; c < nc; c++)
      for (int w = 0; w < W; w++)
        Lost[oa][c][w] = 0;
    full[a] = 0;
    if (!removed) {
      quiet++;
    } else {
      quiet = 1;
      changed[a] = 1;
      alive[a] -= removed;
      /* too few candidates, or an unmatched cell of oa without one */
      uint64_t unm = 0;
      for (int w = 0; w < W; w++) {
        uint64_t un = 0;
        for (int c = 0; c < nc; c++)
          un |= G[a][c][w];
        unm |= s->cells[d][oa][w] & ~s->cells[d][a][w] & ~un;
      }
      if (s->np[a] + (int)alive[a] < n || unm)
        return 1;
    }
    a = oa;
  }
  /* the node lives: the lists without the dropped entries */
  for (int ax = 0; ax < 2; ax++)
    if (changed[ax]) {
      CLS_COMPACT(s, s->vw[d][ax], s->nvalid[d][ax], s->clskm[ax], alive[ax]);
      s->nvalid[d][ax] = alive[ax];
      for (int w = 0; w < W; w++) {
        uint64_t un = 0;
        for (int c = 0; c < nc; c++)
          un |= G[ax][c][w];
        s->uni[d][ax][w] = un;
      }
    }
  return 0;
}

#endif
