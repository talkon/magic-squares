# Retrospective: what each speed-up bought, in expected magic squares

Written 2026-10-08.

* Part 1 ([retrospective/timing.md](retrospective/timing.md)) timed 9 versions of the arrangement search on 33 common sums,
  N = 451-31,743. It gives each version's time ratio to the current code, R_v(N).
* This part turns those ratios into E(C), the expected number of magic squares found with C
  CPU-years. It uses the existence study's cost model (`research/existence.md` §4.1, §6) and
  separates code speed from search strategy.
* **The primary model** is the existence study's revised analytic frontier:
  * an ideal greedy order of all sums with S <= 2 S_min;
  * all exponent assignments;
  * time per sum t_v(N) = R_v(N) x the existence time model.

  For the current code it reproduces existence.md §6 exactly. The model lens is a secondary check.
* All CPU is single-core. 1 CPU-year is ~3 months of the 4-core machine.
* "Legacy" is the code of the first search (`arrangement_6`). "Current" is fdb77fc.
* Tables T1-T8 below are generated.

## 1. Bottom line

* **Code speed: the CPU needed to reach any given E shrank ~17-18x** from the first search's code
  to the current code (T2):
  * 16.4x at E = 0.03, 17.4x at 0.1, 18.0x at 0.3 and 18.2x at 1.
  * The model lens gives 18-19x for E = 0.03-0.3.
  * It falls to ~11-13x only where the frontier depends on sums with more than 256 distinct
    numbers (analytic E = 2, lens E = 1). Every version uses the N x N matrices there (the carried
    bitsets of 2d2bc6d and later need at most 256 labels; legacy and df352df-a72cef3 always use
    the matrices). The measured ratio there is 10.4x, but legacy's part rests on 2 sampled first
    rows (±50%; 8.3-9.0x when chained through df352df and a72cef3 on the same rows).

  This rests on the measured time ratios: 17x at N = 2k and 19x at 10k (T1). The ~21x at 30k is an
  extrapolation of the power fit: no sum with at most 256 labels was measured above N = 23k. The
  model enters only through the N of the sums the frontier uses. Their E-weighted N is 2k-7.5k up
  to 1,000 CPU-years (note under T2), the range where the ratios are measured on ~20 sums.
* **At a fixed budget, the gain in E is much smaller**, because E(C) rises only +0.2-0.5 per
  decade of CPU. E grows x3.1 at 1 CPU-year, x2.7 at 10, x2.3 at 100 and x2.0 at 1,000 (T2). These
  ratios do not depend on the model's overall level F, but they do depend on its slope.
  * This x3.1 is not the gain a further 17-18x would bring today. Legacy at 1 CPU-year collects
    what the current code collects in 1/17.6 CPU-year (0.057), where E(C) is steeper in relative
    terms. Another 17.6x on the current code would give x2.6 at 1 CPU-year (existence.md §6: x2.2
    for 10x).
* **Absolute numbers** (analytic central; the model factor F is x/÷2 at 68%; T1):

  | code | E(1 CPU-yr) | E(10) | E(100) | CPU-yr for E = 1: central (F = 2 to 0.5) | marginal CPU-yr per expected magic square at 1 / 100 CPU-yr |
  |---|---|---|---|---|---|
  | legacy (first search) | 0.051 | 0.13 | 0.28 | 12,000 (710 to 760,000) | 49 / 1,100 |
  | current (fdb77fc) | 0.16 | 0.34 | 0.64 | 650 (39 to 59,000) | 17 / 620 |

  * The CPU for E = 1 is undetermined over 2-3 decades for every version.
  * The time model's level matters less. At part 1's measured level (x0.58), current E(1) is 0.19
    instead of 0.16, and the CPU for E = 1 is 380 instead of 650 CPU-years (T4).
* **What each step bought** (T3):

  | step | CPU factor at fixed E = 0.1 / 1 | E gain at 1 CPU-yr |
  |---|---|---|
  | df352df: new arrange.c (forward checking, MRV, label bitsets, AVX-512) | 3.0 / 2.5 | x1.57 |
  | 8d07935: byte counters, permute filters, vectorized child choice | 1.66 / 1.45 | x1.21 |
  | a72cef3: support filter | 1.33 / 1.78 | x1.11 |
  | 2d2bc6d: candidate lists carry their bitsets | 1.91 / 1.75 | x1.28 |
  | f61e719: cross support (and label order, msearch fixes) | 1.29 / 1.52 | x1.11 |
  | fdb77fc: round-2 micro-optimizations | 1.06 / 1.06 | x1.02 |
  | **all six** | **17.4 / 18.2** | **x3.06** |

  * The constant-factor steps lose part of their effect as the frontier moves to larger N: the
    new arrange.c, and the byte counters with vectorized child choice.
  * The pruning steps gain: the support filter goes from 1.2x to 1.8x, and cross support from 1.2x
    to 1.5x.
  * The carried bitsets help only sums with at most 256 labels.
* **Search strategy, at fixed code** (analytic truth; T5). The first search's strategy:
  * P were picked by hand, and each was swept from S_min upwards; a P's sweep stopped once its
    node count passed a cutoff, after a complete sum. 56% of its P (642 of 1,136, mostly small
    ones, including the 324 searched over all S) were swept to 8 S_min or beyond, far past where
    P(magic | square) is highest. Those P took only 15% of its CPU, so the long sweeps explain
    little of the gap below; most of it is in which P and sums were chosen.
  * It took 1.13 CPU-years with the legacy code, for E = 0.0057 in the analytic model. The first
    search's own data give E = 0.0037-0.0041.

  With the same legacy code:

  | strategy | CPU-yr to reach the first search's E | E after 1.13 CPU-yr |
  |---|---|---|
  | first search | 1.13 (actual) | 0.0057 |
  | scheduler as fitted (near S_min first; its fitted model and pool) | 0.074 (15x less) | 0.019 |
  | ideal frontier, analytic universe (sorted exponents) | 0.026 (44x less) | 0.030 |
  | ideal frontier, all exponent assignments | 0.010 (110x less) | 0.054 |

  * **These factors are conservative by ~1.5-2x.** The first search's CPU is its actual CPU, while
    the other rows are modelled with R_legacy x the existence time model. That model over-predicts
    the first search's own CPU: x1.45 in total on the 479 P whose sweeps stayed below 3.3 S_min,
    and x2.0 overall (T7). On a common time level the first search to scheduler factor is ~20-30x
    rather than 15x: ~22x with the first search's CPU scaled by 1.45, 30x at part 1's time level
    (T8).
  * **The scheduler's gain over the first search holds in all three models.** At the first
    search's E it needs 14x less CPU in the analytic model, 9x in the model lens and 17x in its
    own fitted model (current code, `check_strategy.json`). Across the sensitivity runs it needs
    13-30x less (T8). The same 1.5-2x caveat applies to all of these.
* **Beyond that, the scheduler as fitted falls behind the ideal frontier** (current code):

  | | scheduler as fitted | same pool, true E/t | ideal, analytic universe | ideal, all assignments |
  |---|---|---|---|---|
  | E(1 CPU-yr) | 0.053 | 0.061 | 0.082 | 0.16 |
  | E(100) | 0.15 | 0.17 | 0.30 | 0.64 |
  | E available | 0.24 | 0.24 | 1.4 | 3.9 |

  * The gap is 3.0x at 1 CPU-year and 4.2x at 100.
  * At 1 CPU-year it splits as follows: x1.15 from the scheduler's fitted model, x1.35 from its
    pool against the sorted universe, and x1.9 from the other exponent assignments, which are in
    the pool only partly.
  * The "~10x less" in README.md ("The scheduler as fitted would collect about 10x less") and in
    research/forecast.md's header compares the scheduler's own forecast (0.015 at 1 CPU-year, its
    own model) with the analytic frontier (0.16). Scored by the analytic model, the scheduler's
    choices are worth 0.053, 3.5x its own forecast, so the gap is ~3x, not ~10x. existence.md §6
    itself gives no factor. This rests on the simulation of the scheduler (§3, §4).
  * With the model lens as the truth, the same gap is 3.2x at 1 CPU-year (T6).
* **Code and strategy multiply.** Reaching the first search's E:

  | code, strategy | CPU |
  |---|---|
  | legacy code, first search's strategy | 1.13 CPU-yr |
  | current code, same sums | 0.065 CPU-yr (x17.6) |
  | + scheduler as fitted | 0.0046 CPU-yr = 40 CPU-h (x14) |
  | + ideal frontier | 0.00074 CPU-yr = 6.5 CPU-h (x6) |

  * In all that is ~1,500x: ~18x from code and ~90x from strategy.
  * The model lens gives ~500x: 17.6 x 8.7 x 3.4.
  * Only the x17.6 is measured, as a ratio of actual CPU. The strategy factors are in model units
    at a tiny E (0.0057). The x14 mixes the first search's actual CPU with modelled CPU, so it is
    low by ~1.5-2x (see above): ~2,000-3,000x in all on a common time level. Treat these totals as
    an order of magnitude, not a forecast.
  * At fixed budgets the gains are much smaller. At 1 CPU-year the current code takes E from
    0.051 to 0.16 on the ideal frontier. Better scheduling then takes E from 0.053 (scheduler as
    fitted) to 0.16. Both are factors of ~3.

## 2. Which numbers rest on what

* **Measured speed only.** These are the speed-ups at fixed N (T1, column 3).
* **Measured speed, plus where the frontier is.** These are the CPU factors at fixed E (T2 right,
  T3).
  * The analytic model and the lens differ 2-5x in E at a fixed budget, yet give the same
    legacy/current factor (16-19x) wherever the frontier stays below N ~ 20k: analytic E <= 1,
    lens E <= 0.3 (monotone P).
  * The intermediate versions agree between the two models within 15-30%. Their ratios change with
    N, and the two frontiers sit at somewhat different N for the same E.
  * So "the CPU to reach any fixed E shrank ~18x" holds for any density model that puts the
    yield at N ~ 2k-20k.
  * Across the sensitivity runs (T8) this factor stays at 17.3-17.5 at E = 0.1. At E = 1 it
    ranges 14-19. The low end is the analytic model's low corner, where E = 1 needs sums with
    N > 20k.
* **Also the slope of E(C)**, which is a property of the model:
  * the gain in E at a fixed budget (x3.1 at 1 CPU-year down to x2.0 at 1,000). Across T8 it is
    x2.75-3.3 at 1 CPU-year and x2.1-2.4 at 100;
  * the per-step E gains.
* **Also the level F of the magic-density model.** This is uncertain x/÷2 at 68%
  (existence.md §4.1). It affects:
  * every absolute E;
  * the CPU for E = 1, which moves 10-100x for a factor 2 on E;
  * the marginal costs.
* **Also the time model's level.** Part 1 found the current code ~1.7x faster than the
  existence time model. Rescaling acts like a budget x1.72 for every version alike, so it leaves
  the ratios unchanged (T4).
  * The existence model is used here. It keeps the current-code numbers identical to
    existence.md §6, and it is a production time: msearch wall time, including enumeration.
* **The strategy comparison** is in model units for E: the first search's E and the scheduler's
  yield both come from the same model. Its CPU is not: the first search's is actual, the others'
  modelled, which makes the factors from the first search conservative by ~1.5-2x (§4). The factor
  from the first search to the scheduler is 9-17x across three models and 13-30x across T8. The gap from the scheduler to the
  ideal frontier at 1 CPU-year is 3x in both the analytic model and the lens, and 2.5-3.1x across
  T8.

## 3. Method

* **Time per sum for version v:** t_v(N) = R_v(N) x t_current(N).
  * t_current is the existence study's model (`frontier_rev.py`): the scheduler's model up to
    N = 4,000, then 51.3 s (N/4000)^3.92, steepening to slope 5.2 above 2x10^4.
  * R_v comes from part 1's per-N-bin knots, linear in ln N. Beyond the last knot (~18k) it
    follows each version's power-fit slope.
  * Sums with more than 256 labels use the measured > 256-label ratio. The analytic segments
    carry no labels, so they take a blend by the model lens's share of such sums at each N: 4% at
    N = 13k, 55% at 20k, and all of them above 32k.
  * 8d07935 was not measured above 256 labels and is interpolated. 7cd9bc5 and d023d32 (bench
    sums only) are a72cef3 and 2d2bc6d times their bench ratio.
* **Ideal frontier.** `existence/synth/rev/frontier_rev.py` is executed as is; its arguments can
  be overridden.
  * It supplies its segments: 2,809 sampled P with S <= 2 S_min, plus the P·p sums.
  * It supplies its revised factors and its all-assignments multiplier m = 2.7.
  * For each version the segments are re-sorted by E / t_v.
  * "Analytic universe" means the base segments only, with m = 1.
* **Scheduler as fitted.**
  * Its decisions come from its own model (`model_6.json`, as in research/forecast.md) on its
    candidate pool (`gen_candidates`: 10,075 P, sums up to 3.3 S_min).
  * Each P's sums are pooled in S order into blocks of non-increasing predicted density, and
    the blocks are taken in decreasing predicted density.
  * The blocks are scored with the analytic E(P, S): `gmodel.eval_P` per pool P, times the
    revised factors. They are timed with t_v.
  * Check: with its own predictions as the truth, this reproduces forecast.md's base curve
    (0.015 / 0.020 / 0.021 at 1 / 10 / 100 CPU-years against 0.015 / 0.020 / 0.022).
* **First search.**
  * Its 1,136 P, swept from S_min to max_S (capped at 8 S_min for E).
  * CPU: the actual 35.8M CPU-s. For another version, each P's CPU is multiplied by that
    version's time-weighted ratio to legacy over the P's sums.
  * E: the analytic model on its sums (0.0057). The data give 0.0037 (the model-free lower bound
    from its SP counts) to 0.0041.

## 4. Caveats

* **Sums with more than 256 labels decide the large-E end.**
  * Legacy's ratio there (10.4) rests on 2 sampled first rows, and 8d07935's is interpolated.
  * Nothing was measured between N = 23k and 31.7k with 256 labels or fewer. The 21x for
    legacy at N = 30k is extrapolated.
  * Up to E ~ 1 (analytic) the frontier stays below these sums.
* **The time model.**
  * It is msearch wall time with a build from before round 2, and it runs ~1.7x above the
    current code's measured CPU.
  * R_legacy x t_current over-predicts the first search's own CPU (T7): x1.45 in total on the
    479 P whose sweeps stayed below 3.3 S_min (per-P median x1.18, sd 0.52 in ln), and x2.0 over
    all its P (2.25 CPU-yr predicted for its sums up to 3.3 S_min alone, against 1.13 actual for
    all of them). The rest is mostly small, exhaustively searched P, where the model's times are
    far too high. Part of the x1.45 is the time model's level (x1.7 above this machine's measured
    CPU). Also, the first search ran on SuperCloud, not this machine, and its times are the
    program's own wall-clock timer per sweep (enumeration and file reading excluded).
  * Because the first search's CPU is actual and the other strategies' CPU is modelled, the
    strategy factors are understated by about that much (x1.45-2). On a common level the first
    search to scheduler factor is ~22x (first search scaled by 1.45) to 30x (part 1's time level,
    which under-predicts the first search instead; T8), rather than 15x.
  * Enumeration and msearch overhead are in t_current, but not in the measured ratios. They
    matter only at small N.
* **The ratios were measured near S_min** (x <= 0.5) and on N <= 31.7k.
* **The first search's strategy cannot be extended past its own budget**, because its P were
  picked by hand. It is compared at its own E and CPU only.
* **The scheduler simulation is a smooth approximation of `scheduler.py`.** It has no
  unit-drop, active set or minimum unit size. It is the scheduler fitted in October, not the
  scheduler v2 now being built on the analytic model; the "ideal frontier" rows are that
  scheduler's target.
* **The all-assignments multiplier m = 2.7** assumes the extra assignments have the same E/t
  distribution (existence.md §6).

## 5. Files and re-running (`forecast/` in `research/retrospective/scripts.tar.xz`)

The scripts, the timing data (`measure/times.json`, `measure/runs.jsonl`, `measure/ratio_model.json`)
and `forecast/results.json` are archived in `research/retrospective/scripts.tar.xz`. The
scripts expect the existence study's scripts (`research/existence/scripts.tar.xz`) and data
directories next to them; the paths in them are those of the original run.

| file | content |
|---|---|
| `retro.py` | the computation. It writes `results.json`, `tables.md` and `retro.md` (this file is `retro_head.md` + the tables + `sens.md`). |
| `retro_lib.py` | settings (CFG), the ratio model, the time model, cell builders (lens population, scheduler pool, first search), and the frontier and block arithmetic |
| `an_cells.py` | the analytic E(P, S) on the scheduler's pool and the first search's sums |
| `tables.py` | the tables T1-T7 (`python3 tables.py` rebuilds them from `results.json`) |
| `sens.py`, `sens/`, `sens.json`, `sens.md` | the sensitivity runs (T8); `python3 sens.py --run` |
| `check_strategy.py`, `check_strategy.json` | the strategy under three truths |
| `check_sched.py`, `legacy_E.py`, `legacy_time_check.py` | the checks quoted above |
| `cache/` | cell caches (~1.2 GB), rebuilt automatically when the settings change |
| `old/` | the first draft (model lens as primary) |

**Re-running with revised model parameters:**

```
python3 retro.py [cfg.json] [key=value ...]
```

* **Analytic model.** `analytic_rev_args='["SQ12=1.1","MULT=3"]'` passes arguments to
  `frontier_rev.py`. `analytic_rev_script=...` points to a revised script, which must define
  `w, Em, ns, N, PART, VAR, fac`.
* **Time model.** Set `t_scale`, `t_a`, `t_a2` and `t_n2`. `t_scale_alt` sets the second level
  shown in T4.
* **Ratios.** Set `ratio_model` (part 1's `ratio_model.json`), `beyond`, `label_split` and
  `R256_override`.
* **Model lens.** Set `lens_pop`, `hA`, `hB` and `m_all`. `do_lens=false` skips it.
* A full run takes ~9 minutes; `do_strategy=false` takes ~2.
* This reviewed copy (`retro/final/retro.md`) has hand edits in the narrative and in T1's header
  and first note. `retro.py` regenerates `forecast/retro.md` from `retro_head.md`, so copy
  `final/retro_head.md` over it, and the T1 edits into `tables.py`, before re-running.

## 6. Review (2026-10-08)

A skeptical review checked the builds, the timings, the fits and the arithmetic. Its files are in
`final/` of `research/retrospective/scripts.tar.xz`: `build/`, `bin/`, `indep_runs.py`, `runs/indep.jsonl` and
`recompute.py`.

* **The builds are the stated commits.** Every binary of part 1 is byte-identical
  to one rebuilt here (gcc 13.3 -O3 -DNDEBUG -march=native -flto):
  * from `git archive <commit>` plus `patch_new.py`, for all 8 new-pipeline versions;
  * unpatched, for df352df and fdb77fc;
  * for legacy, from HEAD's `src/c` plus `patch_legacy.py`, and unpatched.

  The legacy sources are byte-identical to the first search's commit (dadd783, May 2023), whose
  CMake build used the same flags. The patches change only the root loop over r1 and the clocks.
  The legacy patch re-implements `search_aux`'s root iteration exactly (same break test, min_vec
  and slots). Part 1's `src/legacy/` holds the discarded slow patch; the kept binary is built from
  `src/legacy2/`.
* **The inputs and conditions are comparable.**
  * Every sum has the same N and label count in every version's runs.
  * The 1-minute load at the start of the large-sum runs averaged 3.2 for legacy and 3.4-3.55
    for the others.
  * Node rates are smooth in N, with no load outliers.
  * Offsets were drawn at random per sum, before the runs.
  * Per-r1 times correlate across versions (ln-time correlation 0.70 for legacy vs fdb77fc and
    0.89 for legacy vs df352df on the full N = 2-3k runs), so the paired ratios are meaningful.
* **Independent re-timing** (own builds, ~10 CPU-min, load 2-4):

  | sum | this review | part 1 |
  |---|---|---|
  | 12 6 3 2 1 0 1, S = 900, N = 2,994, full | legacy 127.5 s, df352df 40.5, a72cef3 17.2, fdb77fc 5.86 / 6.01: ratios 21.2-21.8 / 6.7-6.9 / 2.9 | 21.2 / 6.65 / 2.81 |
  | 14 7 4 4 1 0 0 1, S = 3125, N = 4,956 | legacy 787 s (stride 8, new offset), fdb77fc 49.1 / 49.1 s (full): 16.0 | 17.2 (direct), 17.4 (paired) |
  | 12 6 3 2 1 1, S = 1080, N = 8,891 | legacy 15,600 s (stride 128, 69 r1), fdb77fc 685 / 681 s (stride 16, same offset): 22.8 direct, 23.0 paired | 18.9 (direct, used), 21.7 (paired) |

  Legacy found the expected 18 squares at S = 900 (S = 3125: 5 in the sample). Per sum the
  ratios agree within -7% to +21% (+6% against the paired value at S = 1080), inside the fits' scatter (sd 0.10 in ln). Pooling both
  measurements at S = 1080 would raise the 5.5k-10k knot for legacy from 18.75 to ~19.4. That moves
  the legacy/current CPU factor by ~+2-3% at E = 1 and leaves E at 1-100 CPU-years unchanged, so
  the tables were not re-run.
* **Independent recomputation** (`recompute.py`: frontier_rev.py's own segments and time model;
  ratios refitted from part 1's per-sum records, with own bin knots or one power fit, and no
  > 256-label regime):
  * current code: E = 0.156 / 0.338 / 0.641 / 1.095 / 1.613 at 1 / 10 / 100 / 1,000 / 10^4
    CPU-years, and 39 / 651 / 58,800 CPU-years for E = 1 at F = 2 / 1 / 0.5, as in existence.md §6;
  * legacy: E = 0.051 / 0.126 / 0.282 at 1 / 10 / 100 CPU-years (T1: 0.051 / 0.13 / 0.28);
  * legacy CPU for E = 1: 12,400-12,600 at F = 1 and 714-716 at F = 2 (T1: 12,000 / 710; the
    > 256-label blend accounts for the difference). Marginal at E = 1: 55,700 (T1: 53,000);
  * legacy marginal cost: 49 and 1,076 CPU-years per square at 1 and 100 CPU-years (§1: 49 / 1,100);
  * legacy/current CPU factor at E = 0.03 / 0.1 / 0.3 / 1: 16.9 / 17.5 / 18.0 / 19.1 (T2: 16.4 /
    17.4 / 18.0 / 18.2);
  * df352df, 8d07935, a72cef3, 2d2bc6d and f61e719 CPU for E = 1: 4,900 / 3,350 / 1,880 /
    1,020 / 690 (T1: 4,700 / 3,200 / 1,800 / 1,000 / 690).
* **What the review changed in the narrative.**
  * The strategy factors are now flagged as conservative by ~1.5-2x (actual vs modelled CPU),
    with the total at ~2,000-3,000x on a common time level.
  * The share of the first search's CPU in its long sweeps is now given (15%).
  * The "~10x" claim is attributed correctly: it is in README.md and forecast.md, not
    existence.md §6.
  * The 30k speed-ups are marked as extrapolated, and the > 256-label wording corrected.
  * A note now explains why x3.1 at 1 CPU-year is larger than what a further 17x would give today.

  No table value changed.

## 7. The integrated build: re-forecast (2026-10-09)

The integrated build (cx/integrated, e49934c) combines three prototypes:

* cx/wide: carried bitsets up to 512 labels, and per-r1 widths.
* cx/pretest: the count-only pretest.
* cx/dfirst: `msearch --diag-first`.

The later fixes (86ea201-ccb9328) leave the native search times unchanged. The top-label rule
changes only matrix builds, and the bench files run at 0.96-1.03x fdb77fc's time.

The build is added here as one more version: t_new(N) = R_new(N) x the existence time model, as
in §3. Everything else is unchanged: `frontier_rev.py`'s segments, its factors, m = 2.7 and the
ideal order. fdb77fc still reproduces T1 exactly (0.156 / 0.338 / 0.641 / 1.09, and 39 / 650 /
58,700 CPU-years).

### 7.1 Bottom line

* **At a fixed budget, E grows x1.11 / 1.14 / 1.19 / 1.24 at 1 / 10 / 100 / 1,000 CPU-years**
  (T10). These ratios do not depend on F.
  * The gain grows with the budget, because the build lowers the exponent of N rather than the
    level. Measured, the best mode goes like N^3.5 and fdb77fc like N^4.5.
  * From 100 to 1,000 CPU-years it is worth as much as existence.md §6's "t ~ N^3 beyond
    N = 4,000" (x1.21 / x1.25). At 100 CPU-years it equals a uniform 2x (x1.19).
* **The CPU needed to reach a given E shrinks 1.0x at E = 0.03, 1.2x at 0.1, 1.5x at 0.3, 2.3x
  at 1 and 5.6x at 2** (T10).
  * At E <= 0.03 the frontier sits at N ~ 2-3k. Nothing changed there: those sums have <= 128
    labels and run the plain search.
  * At E = 1 the factor is 1.9-4.7x across the analytic model's corners and m = 2.0 (T11). It
    depends on the density model only through the N where the frontier sits.
* **Absolute numbers** (central; F is x/÷2 at 68%):
  * E = 0.17 / 0.39 / 0.76 / 1.35 at 1 / 10 / 100 / 1,000 CPU-years.
  * E = 1 needs 23 / 280 / 10,000 CPU-years at F = 2 / 1 / 0.5, against 39 / 650 / 59,000 for
    fdb77fc. At 90% (F = 3.2 to 1/3.2) the range is 5 to 9x10^6, against 8 to 8x10^7.
  * E(C) now rises 0.21 / 0.38 / 0.59 / 0.64 per decade of CPU from 1 to 10^4 CPU-years
    (fdb77fc: 0.18 / 0.30 / 0.45 / 0.52).
* **From E ~ 0.3 upwards, most of the gain is the d-first mode.**
  * The plain search alone (cx/wide + cx/pretest) gives x1.05-1.10 in E and 1.2-1.3x in CPU
    at E = 0.1-1. Per sum it runs at 0.73-0.89x fdb77fc's time at N = 6.7-23k with <= 256
    labels, and at 0.26x above 256 labels.
  * d-first on top of that saves 1.0 / 1.2 / 1.8 / 2.3x CPU at E = 0.1 / 0.3 / 1 / 2.
* **Since the first search's code**, the CPU to reach E = 0.03 / 0.1 / 0.3 / 1 shrank
  17 / 21 / 27 / 43x (fdb77fc: 16-18x, T2). E at 1 CPU-year grew x3.4 (fdb77fc: x3.1).
* **The cost per expected magic square now grows ~N^2.8-3.2**, i.e. x7-9 per octave of N from
  2k to 32k. For fdb77fc it was N^3.1-4.6, i.e. x9-24 per octave (T12).
* **Where the frontier is.** The E-weighted N of the sums taken is 3.2k / 4.4k / 6.3k / 9.0k at
  1 / 10 / 100 / 1,000 CPU-years (fdb77fc: 3.0 / 4.1 / 5.6 / 7.5k).
  * At E = 1, sums with more than 256 labels hold 5% of the frontier's E and 13% of its CPU.
  * At E = 2 they hold 28% of the E and 62% of the CPU.

### 7.2 R_new(N)

* **Input.** The perf verifier's 13 sums (research/ideas.md, "Measurements on the integrated
  binary"):
  * process CPU per sum;
  * the plain search paired with fdb77fc on identical r1, the per-r1 minimum of 2 alternating
    runs;
  * d-first d-sampled;
  * relSE <= 8%.
* **Mode.** msearch's default policy: plain below N = 5,000 (`--diag-first-min-n`), d-first at
  N >= 5,000.
  * The crossover is at 3.8-4.9k.
  * At 4.1k the plain search is at 0.663x fdb77fc and d-first at 0.83x.
* **Calibration stream.** d-first writes no semi-magic squares, but the scheduler's models need
  them.
  * Assumed: every d-first sum also runs `--calib-r1-stride k`, which costs **10% of the sum's
    CPU**.
  * That is a stride k ~ 9 / (d-first / plain), i.e. k ~ 15-30. It yields ~1/20 of the sum's
    semi-magic squares, ~50 at the frontier's ~10^3 per sum.
  * The stream finds no extra magic squares, because d-first finds them all. It is pure
    overhead: t = t_dfirst / 0.9 at N >= 5,000.
  * Streams costing 0 / 5 / 20% are in T11.
* **Construction**, as in part 1:
  * At <= 256 labels: per-bin geometric-mean knots (1-2k, 2-3.1k, 3.1-5.5k, 5.5-10k, 10-25k),
    linear in ln N. The knot at N = 758 is the bench ratio, 0.99.
  * Beyond the last knot (16.9k): the ratio's power-fit slope over N >= 3k (-0.51 ± 0.08,
    residual sd 0.13), flat beyond 10^5.
  * Above 256 labels: the geometric mean of the three measured sums, 0.093 (0.084 / 0.114 /
    0.085 at N = 20.9k / 26.6k / 31.7k), held flat.
  * The two are blended by the share of sums with more than 256 labels at each N (§3).

| N | 2k | 3k | 4k | 5k | 7k | 10k | 15k | 20k | 25k | >= 32k |
|---|---|---|---|---|---|---|---|---|---|---|
| R_new, <= 256 labels (best mode + stream) | 1.01 | 0.98 | 0.69 | 0.68 | 0.59 | 0.49 | 0.40 | 0.34 | 0.31 | 0.27 |
| R_new, > 256 labels | 0.104 | | | | | | | | | 0.104 |
| share of sums with > 256 labels | 0 | 0 | 0 | 0 | 0 | 0.01 | 0.11 | 0.51 | 0.91 | 1 |
| **R_new used** (blend) | 1.01 | 0.98 | 0.69 | 0.68 | 0.59 | 0.49 | 0.37 | 0.22 | 0.12 | 0.104 |
| plain only, used (> 256 labels: 0.26) | 1.01 | 0.98 | 0.69 | 0.73 | 0.87 | 0.84 | 0.75 | 0.53 | 0.31 | 0.26 |

### 7.3 Tables

**T9. Key table** (analytic model, ideal frontier, all assignments, as T1)

| version | speed-up vs legacy at N = 2k / 10k / 30k extrapolated (> 256 labels) | E(1 CPU-yr) | E(10) | E(100) | E(1,000) | CPU-yr for E = 1 at F = 2 / **1** / 0.5 | marginal CPU-yr per magic square at 1 / 100 CPU-yr; at E = 1 |
|---|---|---|---|---|---|---|---|
| legacy | 1 / 1 / 1 (1) | 0.051 | 0.126 | 0.282 | 0.546 | 710 / **12,000** / 760,000 | 49 / 1,100; 53,000 |
| fdb77fc | 17.0 / 19.0 / 20.6 (10.4) | 0.156 | 0.338 | 0.641 | 1.095 | 39 / **650** / 59,000 | 17 / 620; 2,900 |
| integrated, plain only | 16.9 / 22.5 / 25.1 (40) | 0.169 | 0.360 | 0.676 | 1.169 | 32 / **500** / 24,000 | 16 / 580; 2,100 |
| **integrated, best mode** | 16.9 / 38.6 / 74 (100) | **0.173** | **0.387** | **0.763** | **1.353** | 23 / **280** / 10,000 | 15 / 470; 1,100 |

* The best mode's speed-ups include the 10% stream. The value in parentheses is the > 256-label
  ratio (for the integrated build, the geometric mean of its three such sums). The 30k column is
  extrapolated, as in T1.
* The CPU for E = 1 at F is the central curve's C(E = 1/F). At F = 2 it is therefore the
  E = 0.5 column of T10, and at F = 0.5 the E = 2 column.

**T10. Ratios that do not depend on F**

| | E_v / E_fdb77fc at C = 1 | 10 | 100 | 1,000 | 10^4 CPU-yr | C_fdb77fc / C_v at E = 0.03 | 0.1 | 0.3 | 0.5 | 1 | 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **integrated, best mode** | **1.11** | **1.14** | **1.19** | **1.24** | 1.23 | 1.02 | **1.21** | **1.49** | 1.69 | **2.34** | 5.6 |
| integrated, plain only | 1.08 | 1.07 | 1.05 | 1.07 | 1.10 | 1.02 | 1.19 | 1.24 | 1.22 | 1.31 | 2.5 |
| d-first on top (best / plain only) | 1.03 | 1.07 | 1.13 | 1.16 | 1.12 | 1.00 | 1.01 | 1.20 | 1.38 | 1.79 | 2.3 |
| legacy -> integrated (x legacy) | 3.40 | 3.07 | 2.71 | 2.48 | 2.07 | 16.6 | 21.0 | 26.8 | 31.0 | 42.7 | 73 |
| legacy -> fdb77fc (T2) | 3.06 | 2.68 | 2.28 | 2.01 | 1.67 | 16.4 | 17.4 | 18.0 | 18.4 | 18.2 | 13.0 |
| for scale: 2x faster at all N (existence.md §6) | 1.28 | 1.22 | 1.19 | 1.14 | 1.10 | 2 | 2 | 2 | 2 | 2 | 2 |
| for scale: t ~ N^3 beyond N = 4,000 (§6) | 1.04 | 1.11 | 1.21 | 1.25 | 1.27 | | | | | | |

**T11. Sensitivity** (one setting changed at a time)

| variant | E_new / E_fdb77fc at 1 / 100 / 1,000 CPU-yr | C_fdb77fc / C_new at E = 0.1 / 0.3 / 1 / 2 | CPU-yr for E = 1 (F = 1): fdb77fc -> new |
|---|---|---|---|
| central (stream 10%, > 256-label ratio flat) | 1.11 / 1.19 / 1.24 | 1.21 / 1.49 / 2.34 / 5.6 | 650 -> 280 |
| calibration stream 0% | 1.13 / 1.22 / 1.26 | 1.21 / 1.58 / 2.59 / 6.2 | 650 -> 250 |
| stream 5% | 1.12 / 1.21 / 1.25 | 1.21 / 1.53 / 2.47 / 5.9 | 650 -> 260 |
| stream 20% | 1.10 / 1.16 / 1.21 | 1.21 / 1.40 / 2.10 / 5.0 | 650 -> 310 |
| 3.1-5.5k knot 0.80 instead of 0.663 | 1.08 / 1.18 / 1.23 | 1.12 / 1.39 / 2.31 / 5.6 | 650 -> 280 |
| <= 256-label ratio flat beyond 16.9k | 1.11 / 1.19 / 1.23 | 1.21 / 1.49 / 2.33 / 5.5 | 650 -> 280 |
| > 256-label ratio keeps the N slope (-0.51) | 1.11 / 1.19 / 1.23 | 1.21 / 1.49 / 2.33 / 5.9 | 650 -> 280 |
| above 32k, the new code's own slope 3.53 instead of the matrix path's 5.2 | 1.11 / 1.19 / 1.24 | 1.21 / 1.49 / 2.34 / 6.0 | 650 -> 280 |
| the verifier's time fits as the ratio: 0.93 (N/4000)^-0.96 | 1.05 / 1.22 / 1.24 | 1.06 / 1.33 / 2.60 / 5.1 | 650 -> 250 |
| time level x0.58, both versions (T4) | 1.11 / 1.20 / 1.24 | 1.21 / 1.49 / 2.34 / 5.6 | 380 -> 160 |
| analytic: review factors at their low corner | 1.10 / 1.16 / 1.20 | 1.27 / 1.66 / 4.71 / - (E = 2 not reached) | 22,000 -> 4,600 |
| analytic: review factors at their high corner | 1.12 / 1.22 / 1.28 | 1.16 / 1.39 / 1.94 / 2.9 | 120 -> 62 |
| analytic: extra assignments m = 2.0 | 1.11 / 1.20 / 1.24 | 1.26 / 1.61 / 3.19 / 8.4 | 2,300 -> 730 |

* Under the same analytic variants, plain only gives x1.06-1.08 at 1,000 CPU-years and 1.2-2.1x
  CPU at E = 1.
* The stream variants scale R above 5k by -10% to +12.5%. Their effect is a guide to the fit's
  scatter, which is ±10-15% on a knot: ±0.03 in the E ratios and ±10-15% in the CPU for E = 1.

**T12. Cost per expected magic square by N band** (base assignments x 2.7, as existence.md §6)

| N band | < 2k | 2-4k | 4-8k | 8-16k | 16-32k | 32-64k |
|---|---|---|---|---|---|---|
| E available | 0.03 | 0.14 | 0.39 | 0.70 | 0.85 | 0.71 |
| CPU-yr per expected square: fdb77fc | 160 | 1,400 | 14,000 | 1.9x10^5 | 4.5x10^6 | 2.5x10^8 |
| integrated, best mode | 160 | 1,200 | 8,400 | 7.7x10^4 | 5.8x10^5 | 2.6x10^7 (1.1x10^7 at its own slope) |
| fdb77fc / integrated | 1.0 | 1.2 | 1.7 | 2.5 | 7.9 | 9.6 (22) |

### 7.4 Caveats

* **d-first writes no semi-magic squares**, except through the calibration stream (10% of the
  CPU of the d-first sums is assumed; T11 has 0-20%).
  * Scheduler v2 does not launch d-first units yet. So the "best mode" rows are what the ideal
    frontier collects once it runs `msearch --diag-first --calib-r1-stride k` at N >= 5,000.
  * Until then, the "plain only" rows apply.
  * Engine 3's time law has no term for the d-first time, nor for the 0.26x above 256 labels.
* **The ratios rest on 13 sums, 6 of them of 13 7 4 3 1 1.**
  * The 3.1-5.5k knot is a single sum. Its 137 labels let nearly every r1 drop a word, which
    gives 0.663.
  * At 0.80 (the geometric mean of new / base at 137-252 labels), the gain at 1 CPU-year falls
    from x1.11 to x1.08, and the CPU factor at E = 0.1 from 1.21 to 1.12. From 100 CPU-years
    on, the E ratios change by <= 0.01.
* **Sums with more than 256 labels decide the large-E end.**
  * Their ratio rests on three sums at N = 21-32k and is held flat beyond them. The new code
    therefore keeps the existence model's slope of 5.2 above 2x10^4, which cx/profile traced
    to fdb77fc's matrix path. The new code itself measured N^3.5.
  * Letting the ratio fall further changes only E(10^4) and the CPU at F = 0.5: 9,800-10,000
    instead of 10,500 CPU-years.
  * Up to E ~ 1.35 (1,000 CPU-years) these sums hold <= 13% of the frontier's E.
* **The time model** is the existence model (§3, §4). Its level cancels in the ratios (x0.58
  gives the same factors, T11).
* **The ratios are for this machine's native build**: AVX-512 with VPOPCNTDQ and GFNI. Other
  targets were measured only on the bench sums (N < 2.3k), where cascadelake runs at
  0.97-1.03x fdb77fc.

### 7.5 Files (`integ/forecast/` in `research/retrospective/scripts.tar.xz`)

* `forecast_new.py` imports `retro_lib` and runs `frontier_rev.py` as `retro.py` does. It holds
  the 13 measured sums, the R_new construction and the variants, and writes `results_new.json`
  (central, the variants, R_new at selected N, and where the frontier is).
* The analytic variants are `python3 forecast_new.py tag=_anlow analytic_rev_args='[...]'`
  (with sens.py's arguments). They write `results_new_{anlow,anhigh,m20}.json`.
* A run takes ~45 s, most of it loading the lens cache that gives the > 256-label share.

## Generated tables

Produced by `retro.py`; every setting is in `results.json`. E is the expected number of magic squares found, C the CPU in single-core CPU-years. Unless marked, E is the existence study's analytic model as revised (`synth/rev/frontier_rev.py`, VAR {'SQ12': 1.23, 'PAIR': 0.87, 'NTR': 1.35, 'K7': 0.85, 'X01': 0.92, 'C53': 2.0, 'MULT': 2.7}): ideal order of all sums with S <= 2 S_min, all exponent assignments by E_all(C) = 2.7 E(C/2.7) plus the P*p sums. The time per sum is R_v(N) x the existence time model (scheduler model to N = 4,000, 51.3 s (N/4000)^3.92 above, slope 5.2 above 2x10^4).

### T1. Key table

| version | what it adds | speed-up vs legacy at N = 2k / 10k / 30k extrapolated (>256 labels, N = 31.7k) | E(1 CPU-yr) (x legacy) | E(10) (x legacy) | E(100) (x legacy) | CPU-yr for E = 1 at F = 2 / **1** / 0.5 | marginal CPU-yr per magic square at E = 1 |
|---|---|---|---|---|---|---|---|
| legacy | code of the first search (`arrangement_6`) | 1.00 / 1.00 / 1.00 (1.00) | 0.051 (x1.00) | 0.13 (x1.00) | 0.28 (x1.00) | 710 / **12,000** / 760,000 | 53,000 |
| df352df | new arrange.c: forward checking, MRV, label bitsets, AVX-512 | 3.21 / 2.58 / 2.25 (2.00) | 0.080 (x1.57) | 0.19 (x1.48) | 0.38 (x1.36) | 260 / **4,700** / 350,000 | 21,000 |
| 8d07935 | + byte counters, permute filters, vectorized child choice | 5.48 / 3.77 / 3.05 (2.69) | 0.097 (x1.90) | 0.22 (x1.73) | 0.43 (x1.53) | 170 / **3,200** / 260,000 | 15,000 |
| a72cef3 | + support filter | 6.74 / 6.37 / 7.51 (5.07) | 0.11 (x2.12) | 0.24 (x1.93) | 0.48 (x1.71) | 110 / **1,800** / 130,000 | 8,000 |
| 7cd9bc5 | + tie-break by smallest label* | 6.80 / 6.44 / 7.59 (5.12) | 0.11 (x2.13) | 0.24 (x1.93) | 0.48 (x1.72) | 110 / **1,800** / 130,000 | 7,900 |
| 2d2bc6d | + candidate lists carry their bitsets | 12.6 / 12.2 / 11.7 (5.20) | 0.14 (x2.71) | 0.30 (x2.40) | 0.58 (x2.04) | 58 / **1,000** / 110,000 | 4,800 |
| d023d32 | + degeneracy label order* | 12.7 / 12.4 / 11.9 (5.39) | 0.14 (x2.72) | 0.30 (x2.41) | 0.58 (x2.05) | 57 / **1,000** / 110,000 | 4,700 |
| f61e719 | + cross support, label order, msearch fixes | 15.6 / 17.9 / 19.1 (10.6) | 0.15 (x2.99) | 0.33 (x2.64) | 0.63 (x2.24) | 41 / **690** / 59,000 | 3,100 |
| fdb77fc | + round-2 micro-optimizations (current) | 17.0 / 19.0 / 20.6 (10.4) | 0.16 (x3.06) | 0.34 (x2.68) | 0.64 (x2.28) | 39 / **650** / 59,000 | 2,900 |

* Speed-ups: measured time ratios (part 1) at sums with <= 256 distinct numbers; beyond the last knot (N ~ 18k; the largest such sum measured is N = 23k) they follow each version's power-fit slope, so the 30k column is an extrapolation. In parentheses: the one measured sum with > 256 labels (N = 31.7k), where every version uses the N x N matrices (legacy and df352df-a72cef3 always do; the carried bitsets of 2d2bc6d and later need <= 256 labels); legacy's value there rests on 2 sampled first rows (±50%; 8.3-9.0 when chained through df352df and a72cef3), and 8d07935's is interpolated.
* E has the model factor F (x/÷2 at 68%): E(C) at F is F x the central value. The "x legacy" ratios do not depend on F. CPU-yr for E = 1 at F is the central curve's C(E = 1/F).
* *7cd9bc5 and d023d32 were measured on bench sums only (N <= 2.3k); above that they are a72cef3 and 2d2bc6d times their bench ratio.

### T2. F-independent: E gain at a fixed budget, and CPU saved at a fixed E

Left: E_v(C) / E_legacy(C) (a constant factor on E cancels). Right: C_legacy(E) / C_v(E), the factor by which the CPU needed to reach a given E shrank. It depends on the magic-density model only through the N of the sums the frontier uses at that E.

| version | E_v / E_legacy at C = 1 | 10 | 100 | 1,000 | 10^4 CPU-yr | C_legacy / C_v at E = 0.03 | 0.1 | 0.3 | 1 | 2 |
|---|---|---|---|---|---|---|---|---|---|---|
| legacy | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| df352df | 1.57 | 1.48 | 1.36 | 1.28 | 1.21 | 3.14 | 3.01 | 2.82 | 2.51 | 2.16 |
| 8d07935 | 1.90 | 1.73 | 1.53 | 1.41 | 1.29 | 5.30 | 4.99 | 4.42 | 3.65 | 2.98 |
| a72cef3 | 2.12 | 1.93 | 1.71 | 1.59 | 1.45 | 6.47 | 6.63 | 6.33 | 6.48 | 5.69 |
| 7cd9bc5 | 2.13 | 1.93 | 1.72 | 1.60 | 1.45 | 6.53 | 6.69 | 6.40 | 6.55 | 5.75 |
| 2d2bc6d | 2.71 | 2.40 | 2.04 | 1.81 | 1.55 | 12.1 | 12.6 | 12.5 | 11.4 | 6.91 |
| d023d32 | 2.72 | 2.41 | 2.05 | 1.82 | 1.55 | 12.3 | 12.8 | 12.7 | 11.6 | 7.12 |
| f61e719 | 2.99 | 2.64 | 2.24 | 1.98 | 1.67 | 15.1 | 16.3 | 17.0 | 17.3 | 12.9 |
| fdb77fc | 3.06 | 2.68 | 2.28 | 2.01 | 1.67 | 16.4 | 17.4 | 18.0 | 18.2 | 13.0 |

Where the frontier is (E-weighted N of the sums taken, base assignments): current code 3,000 at 1 CPU-yr, 4,100 at 10 CPU-yr, 5,600 at 100 CPU-yr, 7,500 at 1000 CPU-yr; legacy 2,000 at 1, 2,800 at 10, 3,800 at 100, 5,100 at 1000.

### T3. What each step bought

| step | CPU factor at fixed E = 0.03 / 0.1 / 0.3 / 1 | E gain at C = 1 / 10 / 100 / 1,000 CPU-yr | lens check: CPU factor at E = 0.03 / 0.1 / 0.3 (mono) |
|---|---|---|---|
| legacy -> df352df | 3.14 / 3.01 / 2.82 / 2.51 | 1.57 / 1.48 / 1.36 / 1.28 | 2.81 / 2.61 / 2.36 |
| df352df -> 8d07935 | 1.69 / 1.66 / 1.56 / 1.45 | 1.21 / 1.17 / 1.12 / 1.10 | 1.56 / 1.48 / 1.42 |
| 8d07935 -> a72cef3 | 1.22 / 1.33 / 1.43 / 1.78 | 1.11 / 1.11 / 1.12 / 1.13 | 1.45 / 1.67 / 1.94 |
| a72cef3 -> 2d2bc6d | 1.87 / 1.91 / 1.98 / 1.75 | 1.28 / 1.25 / 1.19 / 1.14 | 1.97 / 1.90 / 1.67 |
| 2d2bc6d -> f61e719 | 1.24 / 1.29 / 1.36 / 1.52 | 1.11 / 1.10 / 1.10 / 1.09 | 1.37 / 1.45 / 1.54 |
| f61e719 -> fdb77fc | 1.08 / 1.06 / 1.06 / 1.06 | 1.02 / 1.02 / 1.01 / 1.01 | 1.06 / 1.06 / 1.06 |
| **legacy -> fdb77fc** | 16.4 / 17.4 / 18.0 / 18.2 | 3.06 / 2.68 / 2.28 / 2.01 | 18.1 / 18.8 / 17.8 |

### T4. The level of the current code's time model

Primary: the existence study's model (msearch wall time, fitted before round 2). Alternative: x0.58, part 1's median ratio of the current code's measured CPU time to that model on its 11 sums at N = 3.8k-31.7k (range x0.38-1.03; own fit 26 s (N/4000)^4.15). The level applies to every version alike (the ratios R_v stay), so it is the same as a budget x1.72.

| version | time level | E(1 CPU-yr) | E(10) | E(100) | E(1,000) | CPU-yr for E = 1 at F = 2 / 1 / 0.5 |
|---|---|---|---|---|---|---|
| legacy | existence model (used) | 0.051 | 0.13 | 0.28 | 0.55 | 710 / 12,000 / 760,000 |
| legacy | x0.58 (part 1) | 0.064 | 0.15 | 0.33 | 0.63 | 410 / 6,900 / 440,000 |
| fdb77fc | existence model (used) | 0.16 | 0.34 | 0.64 | 1.09 | 39 / 650 / 59,000 |
| fdb77fc | x0.58 (part 1) | 0.19 | 0.40 | 0.74 | 1.21 | 23 / 380 / 34,000 |

### T5. Search strategy at fixed code (analytic model as the truth)

First search: its own 1,136 P, each swept from S_min to its max_S; 1.13 CPU-yr with the legacy code (actual). Its E: 0.0057 from the analytic model on its sums; from its data 0.0037 (model-free lower bound, SP counts) to 0.0041 (squares x P(magic | square)). For another code version its CPU is the actual per-P CPU times that version's time-weighted ratio to legacy over the P's sums. The scheduler's decisions use only its own fitted model (`model_6.json`, 10,075 P of `gen_candidates`, sums to 3.3 S_min, taken per P in S order in blocks of non-increasing predicted density); they are scored with the analytic E and the true time.

**Code: legacy** (the first search's sums cost 1.1 CPU-yr with it)

| strategy | CPU-yr to reach E = 0.0057 | E after 1.1 CPU-yr | E(1 CPU-yr) | E(10) | E(100) | E(10^4) | E available |
|---|---|---|---|---|---|---|---|
| first search (hand-picked P, each swept from S_min) | 1.1 | 0.0057 | - | - | - | - | - |
| scheduler as fitted (its model and pool) | 0.074 | 0.019 | 0.018 | 0.044 | 0.084 | 0.20 | 0.24 |
| same pool, true E/t, S in order per P | 0.033 | 0.024 | 0.023 | 0.051 | 0.096 | 0.20 | 0.24 |
| ideal frontier, analytic universe (sorted exponents; x1) | 0.026 | 0.030 | 0.028 | 0.067 | 0.14 | 0.44 | 1.38 |
| ideal frontier, all exponent assignments (x2.7 + P*p) | 0.01 | 0.054 | 0.051 | 0.13 | 0.28 | 0.96 | 3.92 |

**Code: fdb77fc** (the first search's sums cost 0.065 CPU-yr with it)

| strategy | CPU-yr to reach E = 0.0057 | E after 0.065 CPU-yr | E(1 CPU-yr) | E(10) | E(100) | E(10^4) | E available |
|---|---|---|---|---|---|---|---|
| first search (hand-picked P, each swept from S_min) | 0.065 | 0.0057 | - | - | - | - | - |
| scheduler as fitted (its model and pool) | 0.0046 | 0.019 | 0.053 | 0.097 | 0.15 | 0.23 | 0.24 |
| same pool, true E/t, S in order per P | 0.0021 | 0.024 | 0.061 | 0.11 | 0.17 | 0.23 | 0.24 |
| ideal frontier, analytic universe (sorted exponents; x1) | 0.0017 | 0.029 | 0.082 | 0.17 | 0.30 | 0.67 | 1.38 |
| ideal frontier, all exponent assignments (x2.7 + P*p) | 0.00074 | 0.053 | 0.16 | 0.34 | 0.64 | 1.61 | 3.92 |

**Reaching the first search's E (0.0057), code x strategy:**

| code, strategy | CPU-yr | factor |
|---|---|---|
| legacy code, first search's strategy | 1.1 |  |
| current code, same sums | 0.065 | x17.6 |
| + scheduler as fitted | 0.0046 | x13.9 |
| + ideal order, analytic universe | 0.0017 | x2.7 |
| + all exponent assignments | 0.00074 | x2.3 |
| total | | x1532 |

### T6. Secondary check: the model lens

Model lens (existence/model, H_corr mid; P with non-increasing exponents "mono", all orders by E_all(C) = 4.5 E_mono(C/4.5), optimistic). "central": the lens's realistic all-orders cost of E = 1 for the current code (3,000 CPU-yr) scaled by each version's CPU ratio.

| version | E(1 CPU-yr) mono / all | E(10) | E(100) | E_v / E_legacy at 1 / 100 CPU-yr (mono) | CPU-yr for E = 1: mono / all / central | C_legacy / C_v at E = 0.03 / 0.1 / 0.3 / 1 (mono) |
|---|---|---|---|---|---|---|
| legacy | 0.0080 / 0.017 | 0.022 / 0.052 | 0.060 / 0.14 | 1.00 / 1.00 | 140,000 / 12,000 / 44,000 | 1.00 / 1.00 / 1.00 / 1.00 |
| df352df | 0.013 / 0.030 | 0.035 / 0.083 | 0.089 / 0.22 | 1.64 / 1.49 | 67,000 / 4,700 / 20,000 | 2.81 / 2.61 / 2.36 / 2.06 |
| 8d07935 | 0.016 / 0.038 | 0.042 / 0.10 | 0.10 / 0.26 | 2.01 / 1.74 | 49,000 / 3,300 / 14,000 | 4.37 / 3.86 / 3.36 / 2.80 |
| a72cef3 | 0.019 / 0.043 | 0.049 / 0.12 | 0.13 / 0.31 | 2.31 / 2.16 | 26,000 / 1,700 / 7,400 | 6.33 / 6.44 / 6.50 / 5.31 |
| 7cd9bc5 | 0.019 / 0.043 | 0.049 / 0.12 | 0.13 / 0.31 | 2.32 / 2.17 | 26,000 / 1,700 / 7,300 | 6.39 / 6.51 / 6.57 / 5.36 |
| 2d2bc6d | 0.025 / 0.058 | 0.065 / 0.16 | 0.16 / 0.41 | 3.09 / 2.76 | 24,000 / 980 / 5,300 | 12.5 / 12.2 / 10.9 / 5.80 |
| d023d32 | 0.025 / 0.058 | 0.065 / 0.16 | 0.17 / 0.41 | 3.10 / 2.78 | 23,000 / 970 / 5,200 | 12.7 / 12.4 / 11.1 / 6.01 |
| f61e719 | 0.028 / 0.065 | 0.075 / 0.18 | 0.19 / 0.47 | 3.53 / 3.22 | 12,000 / 650 / 3,100 | 17.1 / 17.8 / 16.8 / 11.5 |
| fdb77fc | 0.029 / 0.067 | 0.077 / 0.18 | 0.20 / 0.48 | 3.61 / 3.30 | 12,000 / 610 / 3,000 | 18.1 / 18.8 / 17.8 / 11.4 |

Strategy with the lens as the truth (first search: lens E of its sums 0.0025):

| code | strategy | CPU-yr to reach the first search's E | E(1 CPU-yr) | E(10) | E(100) | E available |
|---|---|---|---|---|---|---|
| legacy | first search | 1.1 | - | - | - | - |
| legacy | scheduler as fitted | 0.12 | 0.0068 | 0.017 | 0.042 | 0.62 |
| legacy | same pool, true E/t | 0.09 | 0.0082 | 0.022 | 0.054 | 0.62 |
| legacy | ideal, monotone P | 0.099 | 0.0080 | 0.022 | 0.060 | 208 |
| legacy | ideal, all orders (x4.5) | 0.029 | 0.017 | 0.052 | 0.14 | 935 |
| fdb77fc | first search | 0.065 | - | - | - | - |
| fdb77fc | scheduler as fitted | 0.0075 | 0.021 | 0.052 | 0.13 | 0.62 |
| fdb77fc | same pool, true E/t | 0.0056 | 0.028 | 0.068 | 0.15 | 0.62 |
| fdb77fc | ideal, monotone P | 0.006 | 0.029 | 0.077 | 0.20 | 208 |
| fdb77fc | ideal, all orders (x4.5) | 0.0022 | 0.067 | 0.18 | 0.48 | 935 |

### T7. Checks

* First search, time model: R_legacy(N) x t_current(N) over its sums (to 3.3 S_min) predicts 2.25 CPU-yr against 1.13 actual. On the 479 P whose sweeps stayed below 3.3 S_min: 1.36 vs 0.94 CPU-yr, per-P median ratio 1.18, sd(ln) 0.52.
* First search, squares: the model lens predicts 434,817 semi-magic squares in its sums against 758,949 found.
* Analytic E of the first search's sums: 0.0057 (0.0057 in the P that gave squares); lens: 0.0025; data: 0.0037-0.0041. Scheduler pool: analytic E 0.239 vs lens 0.621 (sums to 3.3 S_min); 0 pool P failed to evaluate, 0.00% of the pool's analytic E is extrapolated beyond the profile grid.
* The current code's analytic curve reproduces existence.md section 6 (0.156 / 0.338 / 0.641 / 1.09 / 1.61): 0.156 / 0.338 / 0.641 / 1.09 / 1.61.
* E-weighted share of the lens sums with more than 256 labels (label model x 1.13), used for the analytic segments: N 3239: 0.00, N 5134: 0.00, N 8137: 0.00, N 12896: 0.04, N 20439: 0.55, N 32393: 1.00, N 51340: 1.00.


### T8. Sensitivity (analytic model; one setting changed at a time)

| variant | E(1 CPU-yr) legacy / current (ratio) | E(100) legacy / current (ratio) | CPU-yr for E = 1: legacy / current | C_legacy / C_current at E = 0.1 / 1 | first search -> scheduler, CPU factor at the first search's E (legacy code) | ideal frontier / scheduler, E(1 CPU-yr), current code |
|---|---|---|---|---|---|---|
| base (as in T1-T5) | 0.051 / 0.16 (3.06) | 0.28 / 0.64 (2.28) | 12,000 / 650 | 17.4 / 18.2 | 15.2 | 3.0 |
| labels x1.0 (label model as is) | 0.051 / 0.16 (3.06) | 0.28 / 0.64 (2.28) | 12,000 / 650 | 17.4 / 19.0 | - | - |
| no > 256-label regime (all sums on the <= 256-label ratios) | 0.051 / 0.16 (3.06) | 0.28 / 0.64 (2.28) | 12,000 / 650 | 17.4 / 19.1 | - | - |
| legacy > 256-label ratio x0.5 (5.2) | 0.051 / 0.16 (3.06) | 0.28 / 0.64 (2.28) | 11,000 / 650 | 17.4 / 17.7 | - | - |
| legacy > 256-label ratio x1.5 (15.6) | 0.051 / 0.16 (3.06) | 0.28 / 0.64 (2.28) | 12,000 / 650 | 17.4 / 18.7 | - | - |
| ratios flat beyond the last measured N (~18k) | 0.051 / 0.16 (3.06) | 0.28 / 0.64 (2.28) | 12,000 / 650 | 17.4 / 18.2 | - | - |
| time: no steepening above 2x10^4 (slope 3.92 throughout) | 0.051 / 0.16 (3.06) | 0.28 / 0.64 (2.28) | 12,000 / 650 | 17.4 / 18.2 | 15.2 | 3.0 |
| time: part 1's fit for the current code, 26 s (N/4000)^4.15 | 0.067 / 0.20 (2.92) | 0.34 / 0.73 (2.14) | 7,600 / 420 | 17.4 / 18.2 | 30.1 | 3.1 |
| analytic: review factors at their low corner | 0.040 / 0.11 (2.75) | 0.19 / 0.41 (2.15) | 310,000 / 22,000 | 17.5 / 14.0 | 18.9 | 3.0 |
| analytic: review factors at their high corner | 0.064 / 0.21 (3.31) | 0.40 / 0.95 (2.40) | 2,200 / 120 | 17.3 / 18.6 | 12.6 | 3.0 |
| analytic: extra assignments m = 2.0 instead of 2.7 | 0.042 / 0.13 (3.04) | 0.23 / 0.51 (2.23) | 39,000 / 2,300 | 17.5 / 16.5 | 15.2 | 2.5 |
