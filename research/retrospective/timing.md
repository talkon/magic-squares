# Retro forecasts, part 1: time per sum vs N for each version of the search

Written 2026-10-08. This part measures how long each version of the arrangement search takes per (P, S) sum, as a function of the number of vectors N, on one common set of sums. It also measures the ratio of each version to the current code (fdb77fc). The forecasts themselves (E(C) per code version) are a later part. They should use the ratios below, applied to the existence study's time model, which is maintained elsewhere.

Everything is single-core CPU time (thread CPU clock) on this machine: an Intel Xeon (KVM, 4 vCPU, AVX-512 with VBMI/GFNI/VPOPCNTDQ/BITALG), with each version built `-O3 -march=native -flto`. Data and scripts are in `measure/` of `research/retrospective/scripts.tar.xz` (section 5). The tables generated from the data follow section 5.

## 1. Bottom line

* **Speed-up of the current code (fdb77fc) over the code of the first search (legacy `arrangement_6`)**, for the time per sum:
  * ~11x at N < 1,000;
  * ~15x at N = 1-2k;
  * **17-20x at N = 2k-25k**, the range that matters for the forecasts. The fit is R = 17.2 (N/4000)^0.10: 17 at N = 4k, 19 at 10k, 20 at 20k.
  * ~8-10x at the one sum with more than 256 distinct numbers (13 7 4 3 1 1, S = 2650, N = 31,743, 279 labels). There, every version uses the N x N matrix path, so the current code loses its carried bitsets. This value rests on 2 sampled first rows for legacy (section 4).
* **The speed-ups compound over the six steps** (time ratio of consecutive versions, geometric mean per N range):

  | step | N 1-2k | N 2-3k | N 3-5.5k | N 5.5-10k | N 10-25k |
  |---|---:|---:|---:|---:|---:|
  | legacy -> df352df (forward checking, MRV, label bitsets, AVX-512 filters) | 3.2 | 3.2 | 2.9 | 2.7 | 2.4 |
  | df352df -> 8d07935 (byte counters, permute filters, vectorized child selection) | 1.7 | 1.7 | 1.6 | 1.5 | 1.4 |
  | 8d07935 -> a72cef3 (support filter) | 1.2 | 1.3 | 1.4 | 1.5 | 2.1 |
  | a72cef3 -> 2d2bc6d (candidate lists carry their bitsets; with 7cd9bc5 tie-break, d023d32 label order in between on the bench sums) | 1.8 | 1.9 | 1.9 | 2.1 | 1.65 |
  | 2d2bc6d -> f61e719 (cross support, label order, msearch fixes) | 1.2 | 1.3 | 1.3 | 1.4 | 1.5 |
  | f61e719 -> fdb77fc (round-2 micro-optimizations) | 1.09 | 1.09 | 1.04 | 1.06 | 1.07 |
  | **legacy -> fdb77fc** | **15** | **18.5** | **17** | **19** | **20** |

  * The early steps lose some of their effect at large N: legacy -> df352df falls from 3.2x to 2.4x, and df352df -> 8d07935 from 1.7x to 1.4x.
  * The pruning steps gain: the support filter goes from 1.2x to 2.1x, and cross support from 1.2x to 1.5x.
  * Overall, the ratio to legacy is nearly flat in N above 2k.
* **Ratio of each version to the current code** (fit over N = 3k-23k, 10 sums, R = c (N/4000)^d; sums with <= 256 labels):

  | version | c | d | R at N = 4k / 10k / 20k | resid sd (ln) |
  |---|---:|---:|---|---:|
  | legacy | 17.2 | +0.10 | 17 / 19 / 20 | 0.10 |
  | df352df | 6.0 | +0.22 | 6.0 / 7.3 / 8.5 | 0.06 |
  | 8d07935 | 3.7 | +0.32 | 3.7 / 5.0 / 6.2 | 0.09 |
  | a72cef3 | 2.8 | +0.01 | 2.8 / 2.9 / 2.9 | 0.11 |
  | 2d2bc6d | 1.36 | +0.14 | 1.36 / 1.55 / 1.7 | 0.07 |
  | f61e719 | 1.04 | +0.02 | 1.04 / 1.06 / 1.08 | 0.02 |

  * Per sum, the ratios scatter only by 2-11% (sd in ln) around these fits.
  * Time per sum vs N alone scatters by a factor ~1.4 (sd 0.3-0.45 in ln) between P.
  * So **t_v(N) = R_v(N) x t_current(N)** is the robust way to build a version's time model, with t_current from the existence study's model.
  * `ratio_model.json` has the per-bin knots and these fits.
* **The current code against the existence time model.** The model 51.3 s (N/4000)^3.92 over-predicts the current code on these 11 sums by x1.0-2.6 (median ~1.7; table "Current code vs the existence study's time model" below).
  * The model's level is the scheduler's t(4000), fitted to msearch wall times.
  * The existence study's own free fit on its 17 sums was 28 s (N/4000)^4.33.
  * Here the current code gives t = 26 s (N/4000)^4.15 (N = 3.8k-23k, 10 sums, sd 0.34) and f61e719 gives 27 s (N/4000)^4.17.
  * Part 2 should keep the existence model for the current code, or rescale it consistently. The ratios do not depend on that choice.
* **Cost per semi-magic square** (CPU s per square; all squares in a range's sums / total time):

  | N range | squares | legacy | df352df | 8d07935 | a72cef3 | 2d2bc6d | f61e719 | fdb77fc |
  |---|---:|---:|---:|---:|---:|---:|---:|---:|
  | 1-2k (11 bench sums) | 22 | 3.5 | 1.07 | 0.63 | 0.53 | 0.29 | 0.24 | 0.22 |
  | 2-3k (5 bench sums) | 45 | 5.2 | 1.6 | 0.91 | 0.71 | 0.37 | 0.28 | 0.26 |
  | 3.8k-5k (3 sums) | 72 | 25 | 8.5 | 5.3 | 3.8 | 2.0 | 1.5 | 1.4 |
  | 6k-9k (4 sums) | 735 | 36 | 13 | 9.3 | 6.1 | 2.9 | 2.0 | 1.9 |
  | 15k-23k (3 sums) | ~6,300 (sampled) | 236 | 94 | 69 | 29 | 20 | 11 | 10.6 |

  * The bench sums are chosen for time, not for squares, so their per-square costs are only indicative.
  * The large sums are the existence study's.
  * Squares in the 15k-23k sums are pooled sampled estimates, ±20-30%. They are common to all versions, so the ratios between versions are unaffected.
* **Correctness.** Every full search found the expected squares, with the right count and order-independent hash:
  * 22 bench instances (6x6) x 9 versions;
  * 14 7 5 3 S = 1656 for 3 versions (12 squares, as in the existence run).
  * Versions never disagreed.
  * Sampled runs are consistent with the known counts within their sampling error. For example, 12 6 3 2 1 1 S = 988 (245 squares): fdb77fc 256 and 296, a72cef3 240 and 192, from 1/8 and 1/16 samples.
* **Cost of this measurement.** ~2.6 CPU-hours in all:
  * 7,335 CPU-s of kept runs;
  * 1,832 CPU-s of legacy runs discarded because of a patch bug (section 3);
  * ~250 CPU-s of checks.

## 2. Versions and sums

| version | what it adds |
|---|---|
| legacy | `src/c/arrangement.c` + `search.c` (unchanged since 2023-05): the code of the first search (758,949 squares in ~1 CPU-year). Reads vector files written by `bin/enumerate`; N x N intersection bitsets, 16-wide gathers. |
| df352df | new `arrange.c`: forward checking, MRV, label bitsets, AVX-512 filters (start of the second day) |
| 8d07935 | + byte counters, permute filters, vectorized child selection |
| a72cef3 | + support filter |
| 7cd9bc5 | + tie-break by the smallest label (bench sums only) |
| 2d2bc6d | + candidate lists carry their bitsets (up to 256 labels) |
| d023d32 | + degeneracy label order (bench sums only) |
| f61e719 | + cross support (and label order, msearch overhead fixes) |
| fdb77fc | + round-2 micro-optimizations (current) |

* **Builds.** Each version is a `git archive` of its commit (the repo is unchanged), patched (section 3) and built with `gcc 13.3 -O3 -DNDEBUG -march=native -flto` like `build.sh`.
* **Driver.** `bench` (from df352df on) runs the new versions, one instance per process, with `--update`, which prints count and hash. Legacy runs as `arrangement_6 --file <bin/enumerate output> --sum S`.
* **Bench sums.** All 6x6 instances of `bench/quick.txt`, `full.txt` and `prod.txt`: 22 sums, N = 451-2,994, every version, full searches. The two 5x5 instances of `full.txt` are left out.
* **Large sums.** 11 of the existence study's sums (`existence/model/sample/bigN_3000.txt`, `verify-forecast/region-pool`), N = 3,786-31,743:
  * 14 7 5 3 S = 1656;
  * 13 7 4 2 0 1 S = 1166;
  * 14 7 4 4 1 0 0 1 S = 3125 and 3231;
  * 12 8 4 2 1 S = 1302;
  * 12 6 3 2 1 1 S = 988 and 1080;
  * 13 7 4 3 1 1 S = 2200, 2400 and 2650;
  * 12 9 6 2 1 1 S = 3500.

  The 7 main versions ran on all of these, except 8d07935 at S = 2650. 7cd9bc5 and d023d32 ran on the bench sums only, since their effect there is small: 7cd9bc5 ~1.00x a72cef3, d023d32 ~0.96x 2d2bc6d.
* **Timing.** Thread CPU time, so time spent descheduled on the shared machine (load 2-4) does not count.
  * Each version ran the same sum back to back, with the order rotated between instances.
  * Bench sums: min of 3 repeats for the new versions, 1 for legacy (2 for the verification sum).
  * Time per sum = setup (relabelling, tables) + search. Enumeration is not included; it is the same for all versions and ~0.3% of a sum.
  * Repeat-to-repeat noise is a few % (e.g. fdb77fc on the N = 2,994 sum: 5.80-5.98 s).
* **Calibration against earlier measurements.**
  * Legacy takes 3.51 s on `quick.txt` here, against 3.56 s in the README's earlier measurement. Its node counts are identical (31.7M).
  * df352df takes 1.14 s here vs 1.11 s.
  * fdb77fc takes 0.268 / 3.37 / 13.2 s on quick / full (6x6) / prod, against the README's 0.27 / 3.29 / 13.4 s.

## 3. Method: r1 sampling, and checks

* **Sampling patch.** The search is a loop over the first row r1, one independent subproblem per r1. Each square is found exactly once, under its lowest-index vector.
  * The patch runs only r1 = o, o + k, o + 2k, ... (env `SAMPLE_STRIDE`, `SAMPLE_OFFSET`). It is the existence study's patch, ported to each version's `SEARCH_ROOT` loop by `patch_new.py`. Every version has the same one-line root loop.
  * For legacy, `patch_legacy.py` reimplements the root iteration of `search_aux` in `arrangement.c`, leaving the recursion untouched.
  * Squares, nodes and rdtsc cycles are recorded per r1 in memory and written after the search.
  * The estimate of the sum is k x (sampled time), plus setup. This is unbiased for a random offset (Horvitz-Thompson).
  * Standard errors come from successive differences of the per-r1 times. This is a systematic sample, and the cost per r1 falls steeply with r1.
* **Patch overhead.**
  * Patched vs unpatched builds (stride 1): fdb77fc 16.80 vs 16.93 s and 17.06 vs 17.20 s on all bench sums; df352df 20.4 vs 20.6 s and 20.2 vs 20.5 s.
  * The first legacy patch put the sampling test and an rdtsc inside `search_aux`'s loop, at every node. That made legacy 1.75x slower (8.6 vs 4.8 s on 13 6 3 2 S = 699). It was found by timing against the unpatched binary. All of its runs (1,832 CPU-s) were discarded and repeated with the fixed patch, which matches the unpatched program (4.78 s search vs 4.9 s total with setup and I/O, same 46,266,360 nodes).
* **Nested samples.** For each large sum:
  * Strides are powers of two, about (rough ratio x current time) / (60-300 s), so each run took 0.5-5 min.
  * One random offset per sum, taken modulo each version's stride.
  * Every version's r1 sample then contains the samples of the versions with larger strides.
  * The ratio between two versions can therefore also be computed on the common r1 ("paired"). Where that has the smaller standard error, it is used: ratio x current time is the version's best estimate (`seconds_best`).
  * Paired and direct ratios agree within their errors. The pairing gains less than it could because the label order changes at d023d32, so r1 positions only roughly correspond between the frequency-ordered versions (legacy-2d2bc6d) and the degeneracy-ordered ones (f61e719, fdb77fc).
* **Verification against full runs** (every version, on the largest bench sum 12 6 3 2 1 0 1 S = 900, N = 2,994; table "Sampling check" below).
  * Two stride-16 samples per version (18 runs) estimate the full search time within 0-9% for 14 runs, and 10-18% for 4. The mean is +5%.
  * Three repeats of the fdb77fc o5 sample gave +1-3% instead of +12%, so the large deviations were mostly load noise on single runs.
  * Per r1, a sampled run's times are 1.00-1.14x (mean 1.05) the same r1's times within the full run. That is noise plus at most a few % of lost cache warmth, the same for all versions. Absolute sampled times may therefore be ~3% high; ratios are unaffected.
  * Node counts per r1 are identical, except in f61e719 and fdb77fc, where the adaptive cross-support switch learns its kill rate from the sampled r1 only. There 120-130 of 187 r1 have the same count, and the time is unaffected.
  * From the full per-r1 times, the exact sampling CV of the estimator is 1-4% at k = 4-16 and 18-20% at k = 64 (47 r1). Large sums used >= 100 r1 wherever affordable.
  * Independent second offsets agree within their errors: all 7 versions on 12 6 3 2 1 1 S = 988, and legacy on S = 2200 and 3500 (table below the sampling check).
* **Precision.** Relative standard errors of the large-sum estimates:
  * 1-3% at N <= 9k for the new versions, 5-10% for legacy (16% at S = 1080, 35 r1);
  * 5-14% at N = 15-17k;
  * 8-30% at S = 2400 (N = 23k; legacy has only 5 r1);
  * 35-50% at S = 2650 (8 r1; legacy 2 r1, where the se is not meaningful). The r1 cost falls ~10^5-fold across r1 there, so only the ratios are usable, and they are rough.

## 4. Notes per range of N

* **N < 3,000 (bench, full searches).** The ratio to the current code rises with N for every old version. Legacy goes from 6.8x at N = 451 to 21x at N = 2,994, and df352df from 3.2x to 6.7x.
* **N = 3.8k-9k.**
  * Legacy is 17-20x the current code, df352df 6-7.4x, 8d07935 3.6-5.0x, a72cef3 2.5-3.2x, 2d2bc6d 1.35-1.6x and f61e719 1.02-1.08x.
* **N = 15k-23k** (13 7 4 3 1 1 S = 2200 and 2400, 12 9 6 2 1 1 S = 3500).
  * Legacy is 17-24x, df352df 7-9x, 8d07935 4.7-7x, a72cef3 2.4-3.0x, 2d2bc6d 1.5-2.0x and f61e719 1.06-1.07x (the best per-sum ratios).
  * Seconds per sum for the current code: 9,900 (N = 15.2k), 6,200 (17.1k, 256 labels) and 50,500 (23.0k).
* **N = 31.7k, 279 labels (> 256, matrix path in every version).**
  * Ratios to the current code: df352df 5.2, a72cef3 2.05, 2d2bc6d 2.0 and f61e719 0.98, paired on the same 8 r1.
  * 2d2bc6d ~ a72cef3 here, since its carried bitsets need <= 256 labels.
  * Legacy direct is 10.4, from 2 r1 at stride 16,384. Chained through the frequency-ordered versions on the same 2 r1, it is 8.3-9.0: legacy / df352df = 1.60 and legacy / a72cef3 = 4.39 on those r1, times 5.21 and 2.05. These rest on 2 r1, so ±50%.
  * Above ~256 labels, the old versions lose less ground: df352df's bitsets are 8 words wide, while the current code also falls back to matrices.
  * For the forecasts this matters only for sums with > 256 labels, which here start around N ~ 25-30k.
* **Extrapolation beyond N = 23k with <= 256 labels.** None of the sums measures this. The ratio fits are nearly flat for legacy and a72cef3, and grow slowly for df352df and 8d07935, so extrapolating them to N ~ 30k changes little: legacy 21, df352df 9.3, 8d07935 7.0, 2d2bc6d 1.8.

## 5. Files

`measure/` in `research/retrospective/scripts.tar.xz` (the patched sources and binaries are not included; `patch_new.py` and `patch_legacy.py` recreate them from `git archive` of each commit):
* `times.json`: one record per (version, sum). Fields:
  * `version`, `P`, `S`, `N` (vectors enumerated), `N_searched`, `labels`;
  * `seconds` (own estimate) and `seconds_best` (paired-ratio estimate where better), with `se` / `rel_se`;
  * `full`, `sampled_stride`, `r1_sampled`;
  * `squares` (count, or sampled estimate), `squares_in_sample`, `squares_best`, `squares_source`;
  * `nodes`, `setup`, `hash` / `check` (full runs), `parts` (each sampled run: stride, offset, r1, time, squares, log);
  * `seconds_per_square(_best)`.

  Also `ratios_to_current` (direct, paired, best per sum), `fits`, `ratio_fits` and `bins`.
* `ratio_model.json`: R_v(N), per-bin knots (geometric-mean N, geometric-mean ratio) and power fits, plus the > 256-label point.
* `runs.jsonl`: every run (legacy runs before the time in `legacy_fix_time.txt` are excluded by the scripts). `logs/*.log`: per-r1 records (r1, squares, nodes, CPU s).
* `patch_new.py`, `patch_legacy.py`: the patches. `src/<commit>/`: the patched sources. `bin/`: the builds (`bench_<commit>`, `arrangement_6_legacy`; `*_orig*` unpatched, for the overhead checks).
* `run.py`: one run, appending to `runs.jsonl`. `bench_round.sh`: all bench sums for given versions. `plan.py`: nested-stride plans for the large sums. `plans/`: the plans as run.
* `inst/bench6.txt`, `inst/big.txt`: the sums. `big.txt` also carries the existence runs' counts, times and sampled square estimates.
* **Re-running.** `python3 analyze.py && python3 verify_sampling.py && python3 report.py` rebuilds `times.json`, `verify_sampling.json`, `ratio_model.json` and this file's tables from `runs.jsonl`. `measure_head.md` is the text above.
  * To add a sum: `python3 plan.py --sums "S ..." --target 60 > plan.sh`, add its line to `inst/big.txt`, run the plan and re-run the three scripts.
  * To change the time model the ratios are applied to: nothing here depends on it.

## Tables (generated by report.py)

### Benchmark files (6x6 instances; full searches; CPU s, min of 1-3 repeats)

| version | quick (7 sums) | full (7 sums) | prod (9 sums) | nodes, prod | ratio to current, prod |
|---|---:|---:|---:|---:|---:|
| legacy | 3.51 | 61.6 | 248 | 2.38e+03M | 18.8 |
| df352df | 1.14 | 18.3 | 77.5 | 541M | 5.9 |
| 8d07935 | 0.698 | 10.6 | 44.4 | 439M | 3.4 |
| a72cef3 | 0.642 | 8.20 | 35.2 | 133M | 2.7 |
| 7cd9bc5 | 0.638 | 7.97 | 35.1 | 132M | 2.7 |
| 2d2bc6d | 0.344 | 4.33 | 18.7 | 132M | 1.4 |
| d023d32 | 0.327 | 4.18 | 18.2 | 129M | 1.4 |
| f61e719 | 0.298 | 3.70 | 14.3 | 50.4M | 1.1 |
| fdb77fc | 0.268 | 3.37 | 13.2 | 50.4M | 1.0 |

### Time per sum (CPU s; `~` = r1-sampled estimate, with its relative standard error)

| P | S | N | squares | legacy | df352df | 8d07935 | a72cef3 | 7cd9bc5 | 2d2bc6d | d023d32 | f61e719 | fdb77fc |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 4 3 2 | 327 | 451 | 1 | 0.0085 | 0.0040 | 0.0030 | 0.0034 | 0.0035 | 0.0017 | 0.0017 | 0.0014 | 0.0012 |
| 11 5 3 2 | 401 | 714 | 1 | 0.121 | 0.0390 | 0.0234 | 0.0241 | 0.0243 | 0.0137 | 0.0127 | 0.0113 | 0.0104 |
| 10 6 3 1 0 1 | 391 | 760 | 2 | 0.140 | 0.0461 | 0.0282 | 0.0286 | 0.0283 | 0.0153 | 0.0148 | 0.0135 | 0.0120 |
| 13 6 3 2 | 506 | 831 | 0 | 0.280 | 0.0866 | 0.0503 | 0.0517 | 0.0499 | 0.0277 | 0.0269 | 0.0245 | 0.0221 |
| 11 4 3 1 1 | 470 | 957 | 1 | 0.120 | 0.0449 | 0.0274 | 0.0310 | 0.0309 | 0.0161 | 0.0157 | 0.0143 | 0.0128 |
| 13 6 3 2 | 517 | 976 | 2 | 0.627 | 0.180 | 0.108 | 0.107 | 0.104 | 0.0654 | 0.0615 | 0.0500 | 0.0464 |
| 13 6 3 2 | 561 | 1195 | 2 | 1.58 | 0.505 | 0.297 | 0.275 | 0.270 | 0.147 | 0.141 | 0.129 | 0.117 |
| 13 6 3 2 | 593 | 1299 | 2 | 2.60 | 0.775 | 0.452 | 0.395 | 0.381 | 0.216 | 0.204 | 0.182 | 0.169 |
| 14 7 5 3 | 1361 | 1491 | 1 | 6.50 | 1.61 | 0.907 | 0.767 | 0.779 | 0.445 | 0.418 | 0.362 | 0.323 |
| 13 6 3 2 | 699 | 1637 | 2 | 4.86 | 1.54 | 0.907 | 0.742 | 0.714 | 0.392 | 0.378 | 0.346 | 0.319 |
| 10 6 3 1 0 1 | 648 | 1682 | 1 | 2.21 | 0.741 | 0.458 | 0.396 | 0.397 | 0.204 | 0.194 | 0.184 | 0.164 |
| 14 6 4 2 1 | 920 | 1838 | 3 | 10.4 | 3.20 | 1.85 | 1.65 | 1.61 | 0.902 | 0.839 | 0.680 | 0.629 |
| 13 8 3 2 1 | 922 | 1912 | 2 | 11.5 | 3.75 | 2.29 | 1.88 | 1.85 | 1.01 | 0.992 | 0.805 | 0.753 |
| 11 6 3 1 1 | 700 | 1921 | 1 | 4.42 | 1.39 | 0.838 | 0.676 | 0.674 | 0.347 | 0.328 | 0.347 | 0.299 |
| 10 4 2 2 1 1 | 460 | 1932 | 1 | 11.5 | 3.32 | 1.89 | 1.47 | 1.44 | 0.780 | 0.756 | 0.682 | 0.623 |
| 11 6 3 1 1 | 719 | 1969 | 0 | 4.41 | 1.59 | 0.930 | 0.760 | 0.750 | 0.458 | 0.438 | 0.436 | 0.407 |
| 10 6 4 2 1 1 | 838 | 1986 | 7 | 16.9 | 5.23 | 2.98 | 2.62 | 2.57 | 1.52 | 1.41 | 1.08 | 1.04 |
| 13 6 3 2 1 1 | 899 | 2063 | 18 | 19.8 | 6.42 | 3.87 | 3.23 | 3.25 | 1.74 | 1.74 | 1.34 | 1.27 |
| 14 7 4 2 1 | 1085 | 2082 | 3 | 18.7 | 6.52 | 4.01 | 3.31 | 3.31 | 1.75 | 1.72 | 1.37 | 1.23 |
| 12 6 3 2 1 0 1 | 863 | 2191 | 2 | 29.6 | 8.86 | 5.11 | 4.01 | 3.98 | 2.13 | 2.09 | 1.62 | 1.50 |
| 14 7 5 3 | 1460 | 2237 | 4 | 41.5 | 11.7 | 6.71 | 4.96 | 4.78 | 2.57 | 2.50 | 2.08 | 1.90 |
| 12 6 3 2 1 0 1 | 900 | 2994 | 18 | 123 | 38.6 | 21.4 | 16.3 | 16.3 | 8.47 | 8.28 | 6.34 | 5.80 |
| 14 7 5 3 | 1656 | 3786 | 12 | ~463 (6%) | ~153 (3%) | ~92.0 (2%) | ~64.1 (2%) |  | 34.7 |  | 26.7 | 25.6 |
| 13 7 4 2 0 1 | 1166 | 4090 | 29 | ~559 (10%) | ~192 (3%) | ~120 (2%) | ~90.1 (2%) |  | ~44.1 (1%) |  | ~33.3 (1%) | ~32.7 (1%) |
| 14 7 4 4 1 0 0 1 | 3125 | 4956 | 31 | ~772 (8%) | ~267 (2%) | ~173 (3%) | ~121 (1%) |  | ~62.9 (1%) |  | ~48.3 (1%) | ~45.0 (1%) |
| 12 8 4 2 1 | 1302 | 5994 | 87 | ~2954 (9%) | ~1165 (6%) | ~725 (4%) | ~501 (4%) |  | ~246 (1%) |  | ~169 (1%) | ~164 (1%) |
| 12 6 3 2 1 1 | 988 | 6671 | 245 | ~4998 (8%) | ~1867 (4%) | ~1158 (5%) | ~786 (3%) |  | ~382 (1%) |  | ~277 (1%) | ~262 (1%) |
| 14 7 4 4 1 0 0 1 | 3231 | 7901 | 128 | ~6204 (10%) | ~2299 (5%) | ~1526 (6%) | ~991 (3%) |  | ~475 (2%) |  | ~335 (2%) | ~312 (2%) |
| 12 6 3 2 1 1 | 1080 | 8891 | 275 | ~12358 (16%) | ~4428 (7%) | ~3285 (7%) | ~2076 (8%) |  | ~1042 (5%) |  | ~707 (5%) | ~653 (5%) |
| 13 7 4 3 1 1 | 2200 | 15199 | ~2500 | ~1.98e+05 (14%) | ~70608 (12%) | ~46086 (11%) | ~22759 (9%) |  | ~14518 (7%) |  | ~10601 (6%) | ~9896 (6%) |
| 12 9 6 2 1 1 | 3500 | 17142 | ~310 | ~1.1e+05 (10%) | ~46701 (13%) | ~33768 (9%) | ~17012 (10%) |  | ~9640 (5%) |  | ~6632 (5%) | ~6235 (5%) |
| 13 7 4 3 1 1 | 2400 | 22992 | ~3504 | ~1.21e+06 (28%) | ~4.67e+05 (10%) | ~3.52e+05 (9%) | ~1.61e+05 (9%) |  | ~99507 (8%) |  | ~54083 (11%) | ~50527 (11%) |
| 13 7 4 3 1 1 | 2650 | 31743 | ~2967 | ~1.63e+06 (10%) | ~8.18e+05 (35%) |  | ~3.22e+05 (38%) |  | ~3.14e+05 (39%) |  | ~1.54e+05 (43%) | ~1.57e+05 (43%) |

### Sampling of the large sums (stride k / r1 sampled, per version)

| P | S | N | legacy | df352df | 8d07935 | a72cef3 | 2d2bc6d | f61e719 | fdb77fc |
|---|---:|---:|---|---|---|---|---|---|---|
| 14 7 5 3 | 1656 | 3786 | 16/237 | 8/473 | 4/947 | 4/947 | full | full | full |
| 13 7 4 2 0 1 | 1166 | 4090 | 32/128 | 8/511 | 4/1023 | 4/1023 | 2/2045 | 2/2045 | 2/2045 |
| 14 7 4 4 1 0 0 1 | 3125 | 4956 | 32/155 | 8/620 | 8/620 | 4/1239 | 2/2478 | 2/2478 | 2/2478 |
| 12 8 4 2 1 | 1302 | 5994 | 64/94 | 32/188 | 16/375 | 16/375 | 4/1498 | 4/1498 | 4/1498 |
| 12 6 3 2 1 1 | 988 | 6671 | 64/105, 64/104 | 32/209, 32/208 | 32/209, 32/208 | 16/417, 16/417 | 8/834, 8/834 | 8/834, 8/834 | 8/834, 8/834 |
| 14 7 4 4 1 0 0 1 | 3231 | 7901 | 64/124 | 32/247 | 32/247 | 16/494 | 8/988 | 8/988 | 8/988 |
| 12 6 3 2 1 1 | 1080 | 8891 | 256/35 | 128/70 | 64/139 | 64/139 | 32/278 | 32/278 | 32/278 |
| 13 7 4 3 1 1 | 2200 | 15199 | 1024/15, 1024/14 | 512/29 | 512/29 | 256/59 | 128/118 | 128/118 | 128/118 |
| 12 9 6 2 1 1 | 3500 | 17142 | 512/34, 512/33 | 512/34 | 256/67 | 256/67 | 64/268 | 64/268 | 64/268 |
| 13 7 4 3 1 1 | 2400 | 22992 | 4096/5 | 2048/11 | 2048/11 | 1024/22 | 512/44 | 512/44 | 512/44 |
| 13 7 4 3 1 1 | 2650 | 31743 | 16384/2 | 4096/8 |  | 4096/8 | 4096/8 | 4096/8 | 4096/8 |

### Ratio to the current code (fdb77fc) per sum: direct [paired on common r1]

Direct = ratio of the two estimates; paired = ratio of the summed per-r1 times over the r1 both versions sampled (nested samples); `*` = the one used (smaller standard error). Full searches: N < 3,000, and 14 7 5 3 S = 1656 for 2d2bc6d, f61e719, fdb77fc.

| P | S | N | legacy | df352df | 8d07935 | a72cef3 | 7cd9bc5 | 2d2bc6d | d023d32 | f61e719 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 4 3 2 | 327 | 451 | 6.81 | 3.22 | 2.38 | 2.76 | 2.82 | 1.35 | 1.33 | 1.11 |
| 11 5 3 2 | 401 | 714 | 11.66 | 3.77 | 2.26 | 2.33 | 2.35 | 1.32 | 1.22 | 1.09 |
| 10 6 3 1 0 1 | 391 | 760 | 11.72 | 3.86 | 2.36 | 2.39 | 2.37 | 1.28 | 1.24 | 1.13 |
| 13 6 3 2 | 506 | 831 | 12.70 | 3.92 | 2.28 | 2.34 | 2.26 | 1.26 | 1.22 | 1.11 |
| 11 4 3 1 1 | 470 | 957 | 9.37 | 3.51 | 2.14 | 2.42 | 2.41 | 1.26 | 1.23 | 1.11 |
| 13 6 3 2 | 517 | 976 | 13.53 | 3.89 | 2.33 | 2.30 | 2.25 | 1.41 | 1.33 | 1.08 |
| 13 6 3 2 | 561 | 1195 | 13.45 | 4.31 | 2.53 | 2.35 | 2.30 | 1.25 | 1.21 | 1.10 |
| 13 6 3 2 | 593 | 1299 | 15.40 | 4.60 | 2.68 | 2.34 | 2.26 | 1.28 | 1.21 | 1.08 |
| 14 7 5 3 | 1361 | 1491 | 20.14 | 4.99 | 2.81 | 2.37 | 2.41 | 1.38 | 1.30 | 1.12 |
| 13 6 3 2 | 699 | 1637 | 15.26 | 4.82 | 2.85 | 2.33 | 2.24 | 1.23 | 1.19 | 1.09 |
| 10 6 3 1 0 1 | 648 | 1682 | 13.50 | 4.53 | 2.80 | 2.42 | 2.43 | 1.25 | 1.19 | 1.12 |
| 14 6 4 2 1 | 920 | 1838 | 16.48 | 5.08 | 2.94 | 2.62 | 2.56 | 1.43 | 1.33 | 1.08 |
| 13 8 3 2 1 | 922 | 1912 | 15.28 | 4.98 | 3.04 | 2.50 | 2.45 | 1.34 | 1.32 | 1.07 |
| 11 6 3 1 1 | 700 | 1921 | 14.80 | 4.67 | 2.81 | 2.27 | 2.26 | 1.16 | 1.10 | 1.16 |
| 10 4 2 2 1 1 | 460 | 1932 | 18.53 | 5.34 | 3.03 | 2.36 | 2.32 | 1.25 | 1.21 | 1.10 |
| 11 6 3 1 1 | 719 | 1969 | 10.85 | 3.90 | 2.28 | 1.87 | 1.84 | 1.13 | 1.08 | 1.07 |
| 10 6 4 2 1 1 | 838 | 1986 | 16.17 | 5.00 | 2.85 | 2.51 | 2.46 | 1.45 | 1.35 | 1.03 |
| 13 6 3 2 1 1 | 899 | 2063 | 15.52 | 5.04 | 3.03 | 2.53 | 2.55 | 1.36 | 1.36 | 1.05 |
| 14 7 4 2 1 | 1085 | 2082 | 15.27 | 5.32 | 3.27 | 2.70 | 2.70 | 1.43 | 1.40 | 1.12 |
| 12 6 3 2 1 0 1 | 863 | 2191 | 19.82 | 5.92 | 3.41 | 2.68 | 2.66 | 1.42 | 1.40 | 1.08 |
| 14 7 5 3 | 1460 | 2237 | 21.85 | 6.17 | 3.53 | 2.61 | 2.52 | 1.35 | 1.31 | 1.10 |
| 12 6 3 2 1 0 1 | 900 | 2994 | 21.15 | 6.65 | 3.70 | 2.81 | 2.81 | 1.46 | 1.43 | 1.09 |
| 14 7 5 3 | 1656 | 3786 | 18.12 [18.35]* | 5.98 [6.12]* | 3.60 [3.60] | 2.51 [2.50] |  | 1.36 |  | 1.04 |
| 13 7 4 2 0 1 | 1166 | 4090 | 17.10 [16.06]* | 5.89 [5.84] | 3.66 [3.61]* | 2.76 [2.72]* |  | 1.35 [1.35]* |  | 1.02 [1.02]* |
| 14 7 4 4 1 0 0 1 | 3125 | 4956 | 17.15 [17.37]* | 5.93 [5.91]* | 3.85 [3.84]* | 2.69 [2.68]* |  | 1.40 [1.40]* |  | 1.07 [1.07]* |
| 12 8 4 2 1 | 1302 | 5994 | 17.96 [17.64]* | 7.09 [6.70]* | 4.41 [4.42]* | 3.05 [3.05]* |  | 1.50 [1.50]* |  | 1.03 [1.03]* |
| 12 6 3 2 1 1 | 988 | 6671 | 19.10 [18.14]* | 7.13 [6.95]* | 4.43 [4.31]* | 3.00 [3.00]* |  | 1.46 [1.46]* |  | 1.06 [1.06]* |
| 14 7 4 4 1 0 0 1 | 3231 | 7901 | 19.86 [20.43]* | 7.36 [7.43]* | 4.88 [4.93]* | 3.17 [3.19]* |  | 1.52 [1.52]* |  | 1.07 [1.07]* |
| 12 6 3 2 1 1 | 1080 | 8891 | 18.92 [21.71] | 6.78 [8.44] | 5.03 [5.25]* | 3.18 [3.32]* |  | 1.60 [1.60]* |  | 1.08 [1.08]* |
| 13 7 4 3 1 1 | 2200 | 15199 | 20.01 [16.91]* | 7.13 [7.14]* | 4.66 [4.66]* | 2.30 [2.44]* |  | 1.47 [1.47]* |  | 1.07 [1.07]* |
| 12 9 6 2 1 1 | 3500 | 17142 | 17.60 [18.53]* | 7.49 [8.40]* | 5.42 [5.90]* | 2.73 [2.97]* |  | 1.55 [1.55]* |  | 1.06 [1.06]* |
| 13 7 4 3 1 1 | 2400 | 22992 | 23.86 [23.25] | 9.24 [6.97] | 6.97 [5.25] | 3.18 [2.78]* |  | 1.97 [1.97]* |  | 1.07 [1.07]* |
| 13 7 4 3 1 1 | 2650 | 31743 | 10.42 | 5.21 [5.21]* |  | 2.05 [2.05]* |  | 2.00 [2.00]* |  | 0.98 [0.98]* |

### By N: geometric-mean ratio to current, and CPU s per semi-magic square

Sums measured by all main versions; CPU s per square = total time of the bin's sums / their squares (full counts, or pooled sampled estimates for N > 10k).

Cells: geometric-mean ratio (range over the bin's sums) / CPU s per square.

| N range (sums, squares) | legacy | df352df | 8d07935 | a72cef3 | 7cd9bc5 | 2d2bc6d | d023d32 | f61e719 | fdb77fc |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 400-1000 (6, 7) | 10.70 (6.8-13.5) / 0.185 | 3.69 (3.2-3.9) / 0.0573 | 2.29 (2.1-2.4) / 0.0343 | 2.42 (2.3-2.8) / 0.0351 | 2.40 (2.2-2.8) / 0.0344 | 1.31 (1.3-1.4) / 0.0200 | 1.26 (1.2-1.3) / 0.0190 | 1.11 (1.1-1.1) / 0.0164 | 1.00 (1.0-1.0) / 0.0150 |
| 1000-2000 (11, 22) | 15.25 (10.8-20.1) / 3.50 | 4.73 (3.9-5.3) / 1.07 | 2.78 (2.3-3.0) / 0.627 | 2.35 (1.9-2.6) / 0.529 | 2.31 (1.8-2.6) / 0.520 | 1.28 (1.1-1.5) / 0.292 | 1.22 (1.1-1.3) / 0.277 | 1.09 (1.0-1.2) / 0.238 | 1.00 (1.0-1.0) / 0.220 |
| 2000-3100 (5, 45) | 18.51 (15.3-21.9) / 5.16 | 5.79 (5.0-6.6) / 1.60 | 3.38 (3.0-3.7) / 0.914 | 2.66 (2.5-2.8) / 0.706 | 2.65 (2.5-2.8) / 0.703 | 1.40 (1.4-1.5) / 0.370 | 1.38 (1.3-1.4) / 0.363 | 1.09 (1.0-1.1) / 0.283 | 1.00 (1.0-1.0) / 0.260 |
| 3100-5500 (3, 72) | 17.23 (16.1-18.3) / 24.7 | 5.97 (5.9-6.1) / 8.54 | 3.68 (3.6-3.8) / 5.32 | 2.64 (2.5-2.7) / 3.80 |  | 1.37 (1.3-1.4) / 1.97 |  | 1.04 (1.0-1.1) / 1.50 | 1.00 (1.0-1.0) / 1.43 |
| 5500-10000 (4, 735) | 18.75 (17.6-20.4) / 35.9 | 6.96 (6.7-7.4) / 13.2 | 4.71 (4.3-5.2) / 9.28 | 3.14 (3.0-3.3) / 6.05 |  | 1.52 (1.5-1.6) / 2.92 |  | 1.06 (1.0-1.1) / 2.02 | 1.00 (1.0-1.0) / 1.89 |
| 10000-25000 (3, 6314) | 19.55 (16.9-23.9) / 236 | 8.21 (7.1-9.2) / 93.5 | 5.76 (4.7-7.0) / 68.9 | 2.72 (2.4-3.0) / 29.0 |  | 1.65 (1.5-2.0) / 19.6 |  | 1.07 (1.1-1.1) / 11.3 | 1.00 (1.0-1.0) / 10.6 |

### Fits t(N) = a (N/4000)^b (CPU s per sum; least squares in ln t; resid sd in ln)

| version | N>=1000 | N<3000 | N>=3000 | N 3000-10000 | N>=10000 |
|---|---|---|---|---|---|
| legacy | a 337, b 4.57 (sd 0.49, n 26) | a 356, b 4.69 (sd 0.56, n 16) | a 452, b 4.25 (sd 0.35, n 10) | a 493, b 3.99 (sd 0.24, n 7) | a 84.2, b 5.37 (sd 0.75, n 3) |
| df352df | a 116, b 4.69 (sd 0.45, n 26) | a 116, b 4.76 (sd 0.50, n 16) | a 157, b 4.36 (sd 0.34, n 10) | a 173, b 4.07 (sd 0.27, n 7) | a 53.2, b 5.10 (sd 0.67, n 3) |
| 8d07935 | a 72.5, b 4.77 (sd 0.44, n 26) | a 66.7, b 4.73 (sd 0.49, n 16) | a 97.5, b 4.46 (sd 0.32, n 10) | a 104, b 4.29 (sd 0.25, n 7) | a 22.7, b 5.43 (sd 0.65, n 3) |
| a72cef3 | a 51.2, b 4.55 (sd 0.45, n 26) | a 47.1, b 4.52 (sd 0.49, n 16) | a 74.3, b 4.16 (sd 0.31, n 10) | a 75.2, b 4.14 (sd 0.26, n 7) | a 30.1, b 4.75 (sd 0.62, n 3) |
| 7cd9bc5 | a 47.2, b 4.54 (sd 0.49, n 16) | a 47.2, b 4.54 (sd 0.49, n 16) |  |  |  |
| 2d2bc6d | a 27.5, b 4.56 (sd 0.45, n 26) | a 24.7, b 4.48 (sd 0.49, n 16) | a 35.9, b 4.29 (sd 0.36, n 10) | a 39.2, b 4.02 (sd 0.25, n 7) | a 8.32, b 5.27 (sd 0.77, n 3) |
| d023d32 | a 24.4, b 4.52 (sd 0.50, n 16) | a 24.4, b 4.52 (sd 0.50, n 16) |  |  |  |
| f61e719 | a 20.8, b 4.45 (sd 0.39, n 26) | a 18.0, b 4.31 (sd 0.42, n 16) | a 27.4, b 4.17 (sd 0.34, n 10) | a 29.9, b 3.88 (sd 0.23, n 7) | a 16.3, b 4.54 (sd 0.75, n 3) |
| fdb77fc | a 19.4, b 4.46 (sd 0.41, n 26) | a 16.7, b 4.33 (sd 0.43, n 16) | a 26.3, b 4.15 (sd 0.34, n 10) | a 28.9, b 3.84 (sd 0.25, n 7) | a 15.3, b 4.54 (sd 0.74, n 3) |

Sums with N < 25,000 (<= 256 labels); per-version times use `seconds_best`. Time per sum depends on P as well as N (resid sd 0.3-0.5 in ln), so these fits depend on the mix of P; the N >= 10,000 fits (3 sums of 2 P) are not meaningful. 7cd9bc5 and d023d32 were measured on the bench sums only.

Ratio fits R(N) = c (N/4000)^d (best ratio to current per sum, N >= 1000 / N >= 3000, N < 25,000):

| version | N >= 1000 | N >= 3000 |
|---|---|---|
| legacy | c 17.41, d +0.10 (sd 0.15, n 26) | c 17.19, d +0.10 (sd 0.10, n 10) |
| df352df | c 5.98, d +0.23 (sd 0.09, n 26) | c 5.99, d +0.22 (sd 0.06, n 10) |
| 8d07935 | c 3.74, d +0.31 (sd 0.09, n 26) | c 3.71, d +0.32 (sd 0.09, n 10) |
| a72cef3 | c 2.64, d +0.09 (sd 0.10, n 26) | c 2.83, d +0.01 (sd 0.11, n 10) |
| 7cd9bc5 | c 2.83, d +0.21 (sd 0.09, n 16) |  |
| 2d2bc6d | c 1.42, d +0.10 (sd 0.07, n 26) | c 1.36, d +0.14 (sd 0.07, n 10) |
| d023d32 | c 1.46, d +0.19 (sd 0.08, n 16) |  |
| f61e719 | c 1.08, d -0.01 (sd 0.03, n 26) | c 1.04, d +0.02 (sd 0.02, n 10) |

### Current code vs the existence study's time model and runs (N >= 3,000)

Model: 51.3 s (N/4000)^3.92 ("a = 3.92"; the steep variant has slope 5.2 above N = 2x10^4). Existence runs: msearch / sampled msearch with the f61e719 search code, wall time on the loaded machine.

| P | S | N | labels | fdb77fc (CPU s) | f61e719 (CPU s) | existence run (f61e719 code, wall s) | model a = 3.92 | model steep | fdb77fc / model |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 14 7 5 3 | 1656 | 3786 | 130 | 25.6 | 26.7 | 30.0 | 41.4 | 41.4 | 0.62 |
| 13 7 4 2 0 1 | 1166 | 4090 | 137 | 32.7 | 33.3 | 42.0 | 56.0 | 56.0 | 0.58 |
| 14 7 4 4 1 0 0 1 | 3125 | 4956 | 156 | 45.0 | 48.3 | 50.0 | 119 | 119 | 0.38 |
| 12 8 4 2 1 | 1302 | 5994 | 155 | 164 | 169 | 214 | 250 | 250 | 0.66 |
| 12 6 3 2 1 1 | 988 | 6671 | 159 | 262 | 277 | 291 | 381 | 381 | 0.69 |
| 14 7 4 4 1 0 0 1 | 3231 | 7901 | 183 | 312 | 335 | 348 | 740 | 740 | 0.42 |
| 12 6 3 2 1 1 | 1080 | 8891 | 179 | 653 | 707 | 1022 | 1175 | 1175 | 0.56 |
| 13 7 4 3 1 1 | 2200 | 15199 | 211 | 9896 | 10601 | 12389 | 9611 | 9611 | 1.03 |
| 12 9 6 2 1 1 | 3500 | 17142 | 256 | 6235 | 6632 | 7206 | 15402 | 15402 | 0.40 |
| 13 7 4 3 1 1 | 2400 | 22992 | 245 | 50527 | 54083 | 46660 | 48688 | 58200 | 0.87 |
| 13 7 4 3 1 1 | 2650 | 31743 | 279 | 1.57e+05 | 1.54e+05 | 3.75e+05 | 1.72e+05 | 3.11e+05 | 0.50 |

### Sampling check (P = 12 6 3 2 1 0 1, S = 900, N = 2994): full search vs stride-16 samples

| version | full (CPU s) | sampled estimate / full (se) | per-r1 time, sampled run / full run | nodes per r1 equal | exact CV of the estimator, k = 4 / 16 / 64 / 256 |
|---|---:|---|---|---|---|
| legacy | 123 | o5: 1.09 (0.07); o13: 1.00 (0.07) | 1.02; 1.00 | 187/187; 187/187 | 0.01 / 0.03 / 0.20 / 0.30 |
| df352df | 38.6 | o5: 1.03 (0.07); o13: 1.01 (0.07) | 1.00; 1.03 | 187/187; 187/187 | 0.01 / 0.04 / 0.18 / 0.26 |
| 8d07935 | 21.4 | o5: 1.08 (0.07); o13: 1.04 (0.07) | 1.05; 1.06 | 187/187; 187/187 | 0.01 / 0.03 / 0.18 / 0.26 |
| a72cef3 | 16.3 | o5: 1.18 (0.09); o13: 1.00 (0.07) | 1.13; 1.02 | 187/187; 187/187 | 0.01 / 0.04 / 0.19 / 0.27 |
| 7cd9bc5 | 16.3 | o5: 1.05 (0.07); o13: 1.10 (0.08) | 1.02; 1.11 | 187/187; 187/187 | 0.01 / 0.03 / 0.19 / 0.26 |
| 2d2bc6d | 8.46 | o5: 1.08 (0.07); o13: 1.05 (0.07) | 1.04; 1.06 | 187/187; 187/187 | 0.01 / 0.03 / 0.19 / 0.26 |
| d023d32 | 8.27 | o5: 1.01 (0.06); o13: 1.07 (0.06) | 1.01; 1.06 | 187/187; 187/187 | 0.01 / 0.02 / 0.20 / 0.26 |
| f61e719 | 6.34 | o5: 1.01 (0.06); o13: 1.02 (0.06); o5: 1.02 (0.06); o5: 1.00 (0.06); o5: 0.99 (0.06) | 1.00; 1.02; 1.00; 0.98; 0.98 | 130/187; 120/187; 130/187; 130/187; 130/187 | 0.01 / 0.02 / 0.20 / 0.26 |
| fdb77fc | 5.80 | o5: 1.12 (0.07); o13: 1.14 (0.07); o5: 1.01 (0.07); o5: 1.03 (0.06); o5: 1.03 (0.06) | 1.12; 1.14; 1.02; 1.04; 1.04 | 130/187; 120/187; 130/187; 130/187; 130/187 | 0.01 / 0.03 / 0.20 / 0.26 |

Sums sampled at two or more offsets by the same version (estimate of the search, CPU s +- se):

| version | P S | estimates |
|---|---|---|
| legacy | 12 6 3 2 1 1 988 | k 64, o 11, m 105: 5373 +- 581; k 64, o 17, m 104: 4618 +- 600 |
| df352df | 12 6 3 2 1 1 988 | k 32, o 11, m 209: 1924 +- 113; k 32, o 17, m 208: 1809 +- 119 |
| 8d07935 | 12 6 3 2 1 1 988 | k 32, o 11, m 209: 1204 +- 75.5; k 32, o 17, m 208: 1112 +- 76.4 |
| a72cef3 | 12 6 3 2 1 1 988 | k 16, o 1, m 417: 784 +- 28.4; k 16, o 11, m 417: 787 +- 28.5 |
| 2d2bc6d | 12 6 3 2 1 1 988 | k 8, o 1, m 834: 376 +- 6.97; k 8, o 3, m 834: 388 +- 6.88 |
| f61e719 | 12 6 3 2 1 1 988 | k 8, o 1, m 834: 278 +- 5.11; k 8, o 3, m 834: 276 +- 4.87 |
| fdb77fc | 12 6 3 2 1 1 988 | k 8, o 1, m 834: 262 +- 4.88; k 8, o 3, m 834: 262 +- 4.67 |
| legacy | 13 7 4 3 1 1 2200 | k 1024, o 487, m 15: 2.26e+05 +- 32026; k 1024, o 871, m 14: 1.68e+05 +- 47310 |
| legacy | 12 9 6 2 1 1 3500 | k 512, o 133, m 34: 1.04e+05 +- 16895; k 512, o 325, m 33: 1.16e+05 +- 13381 |

