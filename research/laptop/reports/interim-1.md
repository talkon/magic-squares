# Interim stage-1 report 1 (digest seq 1: 7,143 of 21,298 units, 16.3 reference CPU-hours)

Made from research/laptop/digest.json.xz (commit f9681b7) with scripts/decide.py, without and with
the scheduler state the plan was made from (a copy of the calibration search's profiled state,
1,324 files; research/stage1/PREREGISTERED.txt). An interim look: the plan is not finished
(5 of its 134 band units are), so the band's f_rho rests mostly on the prior.

Summary:
- Squares, S and P traversals: 0.970 / 0.965 / 0.942 of the frozen predictions (7,143 units).
- (square, SP traversal) pairs outside N' 3-6k: 344 / 361.7, f_rho 0.950 [0.868, 1.037] (flat prior).
- In N' 3-6k, counting the band sums of units that start below 3k: 10 pairs / 12.92 predicted
  (flat prior f_rho 0.744 [0.399, 1.249]; with the shipped prior 0.947 [0.678, 1.279]).
- Predictive P(>=1 magic square): 11.0% / 16.2% / 23.5% at 1 / 3 / 10 CPU-years
  (before stage 1: 13.8% / 20.1% / 28.7%); the x band of E at 1 CPU-year 4.6 (before 5.8).
- No magic square; best score 10 of 14. M1 Pro (6 performance cores, portable kernels):
  3.74x the plan's pred_time per core overall, 3.97x in the last 30 minutes.

## Without the state (the plan's frozen predictions only)

```
# stage-1 decision report: plan research/stage1/plan-stage1.jsonl.xz, units digests research/laptop/digest.json.xz (no state)
records: plan research/stage1/plan-stage1.jsonl.xz: 7143 unit files, 7143 finished (the fit uses these)

frozen predictions of the 7143 units (at launch / in the plan, f_rho = 1):
  squares    276067 observed /  284682.17 predicted = 0.970 [0.967, 0.973]
  S          282124 observed /  292266.22 predicted = 0.965 [0.962, 0.968]
  P           70844 observed /   75221.28 predicted = 0.942 [0.936, 0.948]
  pairs         354 observed /     392.08 predicted = 0.903 [0.825, 0.986]
  CPU        61.023 h measured / 16.303 h predicted (reference); SP traversals before the class factors predicted 564.58
  best square: best_score 10 (14 = magic)

SP coupling (pairs in N' 3000-6000: the finished units whose sums all lie in it):
  f_rho = 0.949 median, 0.965 mean (90% 0.635-1.351; ln sd 0.230, of f_rho^2 0.459) from 1.0 pairs / 3.45 predicted at f_rho = 1 over 4 P, dispersion 1.18, prior shipped
  the band's pairs alone (flat prior): f_rho = 0.187 median, 0.290 mean (90% 0.009-0.922; ln sd 1.453, of f_rho^2 2.907) from 1.0 pairs / 3.45 predicted at f_rho = 1 over 4 P, dispersion 1.18, prior flat; x band of E at 1 CPU-year 10923.99
  check, outside the band (flat prior, not used): f_rho = 0.935 median, 0.936 mean (90% 0.852-1.024; ln sd 0.056, of f_rho^2 0.112) from 318.0 pairs / 339.61 predicted at f_rho = 1 over 5173 P, dispersion 1.00, prior flat
  data: in the band 1 pairs / 3.45 from 5 units; outside 318 / 339.61; 54 units across an edge of the band left out (35 pairs)

predictive E (E_ship from research/stage1/ecurve.json: forecast --shipped (engine 4, the measured truth, seed 1) on the calibration search's profiled state, f_rho = 1: 0.1-10 CPU-years on 10% of the candidates (0.0427 / 0.0718 / 0.1202 / 0.1843 / 0.2827); 30 and 100 on 1% (0.4240 / 0.6009) x 0.944, the 10% / 1% ratio at 10 CPU-years (the 1% sample runs 6-11% high at 0.1-10); x f_rho^2 x lognormal(class factors, pair 0.2, selection 0.12)):
  CPU-years  E_ship E (f^2 med)      5%  median     95% x band  P(>=1)   flat prior   before stage 1 (shipped prior): P(>=1), x band
          1   0.120       0.108   0.044   0.107   0.245   5.60   11.2%        2.1%    13.8%, x5.77
          3   0.184       0.166   0.067   0.165   0.374   5.56   16.5%        3.1%    20.1%, x5.73
         10   0.283       0.254   0.103   0.252   0.572   5.53   23.9%        4.5%    28.7%, x5.69
  CPU-years for P(>=1) = 25%: 12.2, 50%: 282 (before stage 1: 6.34, 131; the curve is extrapolated beyond 100 CPU-years at its last slope)

decision table (research/stage1/decision-table.json; thresholds set by the user before stage 1; this tool only reports):
  stop after stage 1                         if P(>=1 at 10 CPU-years) < X: P = 23.9%: X not set
  continue to 1 CPU-year                     if P(>=1 at 3 CPU-years) >= X: P = 16.5%: X not set
  continue to 3 CPU-years                    if P(>=1 at 3 CPU-years) >= X: P = 16.5%: X not set
  plan for 10 CPU-years                      if P(>=1 at 10 CPU-years) >= X: P = 23.9%: X not set
```

## With the state the plan was made from

```
# stage-1 decision report: STATE (a copy of the plan's state; units: digests research/laptop/digest.json.xz)
stage 1: off
records: plan research/stage1/plan-stage1.jsonl.xz: 7143 unit files
each plain (P, S) sum once: 814055 sums kept, 0 repeats dropped (0 squares, 0 pairs); 0 d-first sums with a plain record left out

checks (obs / pred, 90% Poisson; pred at the state's learned class factors, and at the shipped ones in brackets):
  <1k    squares 4514/5998 = 0.753 [0.734, 0.771] (0.945); S trav 6045/6172 = 0.979 [0.959, 1.000] (1.081); P trav 1410/1533 = 0.920 [0.880, 0.961] (0.907)
  1-3k   squares 259076/266486 = 0.972 [0.969, 0.975] (1.003); S trav 264763/266620 = 0.993 [0.990, 0.996] (1.042); P trav 66661/67592 = 0.986 [0.980, 0.993] (0.992)
  3-6k   squares 12477/10557 = 1.182 [1.165, 1.199] (1.078); S trav 11316/11976 = 0.945 [0.930, 0.960] (0.966); P trav 2773/2926 = 0.948 [0.918, 0.978] (0.918)

frozen predictions of the 7143 units (at launch / in the plan, f_rho = 1):
  squares    276067 observed /  284682.17 predicted = 0.970 [0.967, 0.973]
  S          282124 observed /  292266.22 predicted = 0.965 [0.962, 0.968]
  P           70844 observed /   75221.28 predicted = 0.942 [0.936, 0.948]
  pairs         354 observed /     392.08 predicted = 0.903 [0.825, 0.986]
  CPU        61.023 h measured / 16.303 h predicted (reference); SP traversals before the class factors predicted 564.58
  best square: best_score 10 (14 = magic)

SP coupling (pairs in N' 3000-6000, cells ['3-6k']):
  f_rho = 0.947 median, 0.959 mean (90% 0.678-1.279; ln sd 0.193, of f_rho^2 0.386) from 10.0 pairs / 12.92 predicted at f_rho = 1 over 56 P, dispersion 1.15, prior shipped
  the band's pairs alone (flat prior): f_rho = 0.744 median, 0.774 mean (90% 0.399-1.249; ln sd 0.349, of f_rho^2 0.697) from 10.0 pairs / 12.92 predicted at f_rho = 1 over 56 P, dispersion 1.15, prior flat; x band of E at 1 CPU-year 11.43
  check, outside the band (flat prior, not used): f_rho = 0.950 median, 0.951 mean (90% 0.868-1.037; ln sd 0.054, of f_rho^2 0.108) from 344.0 pairs / 361.69 predicted at f_rho = 1 over 5174 P, dispersion 1.00, prior flat
  data: plain squares 354 pairs / 374.61; 0 d-first sums searched in full in one record (--dfirst-sums on): 0.0 est. pairs / 0.00
  pairs per N' band (obs / model at f_rho = 1): <1k 5/11.97, 1-3k 339/349.72, 3-6k 10/12.92

predictive E (E_ship from research/stage1/ecurve.json: forecast --shipped (engine 4, the measured truth, seed 1) on the calibration search's profiled state, f_rho = 1: 0.1-10 CPU-years on 10% of the candidates (0.0427 / 0.0718 / 0.1202 / 0.1843 / 0.2827); 30 and 100 on 1% (0.4240 / 0.6009) x 0.944, the 10% / 1% ratio at 10 CPU-years (the 1% sample runs 6-11% high at 0.1-10); x f_rho^2 x lognormal(class factors, pair 0.2, selection 0.12)):
  CPU-years  E_ship E (f^2 med)      5%  median     95% x band  P(>=1)   flat prior   before stage 1 (shipped prior): P(>=1), x band
          1   0.120       0.108   0.049   0.107   0.222   4.57   11.0%        7.7%    13.8%, x5.77
          3   0.184       0.165   0.075   0.164   0.339   4.53   16.2%       11.4%    20.1%, x5.73
         10   0.283       0.253   0.115   0.251   0.519   4.50   23.5%       16.7%    28.7%, x5.69
  CPU-years for P(>=1) = 25%: 12.8, 50%: 287 (before stage 1: 6.34, 131; the curve is extrapolated beyond 100 CPU-years at its last slope)

decision table (research/stage1/decision-table.json; thresholds set by the user before stage 1; this tool only reports):
  stop after stage 1                         if P(>=1 at 10 CPU-years) < X: P = 23.5%: X not set
  continue to 1 CPU-year                     if P(>=1 at 3 CPU-years) >= X: P = 16.2%: X not set
  continue to 3 CPU-years                    if P(>=1 at 3 CPU-years) >= X: P = 16.2%: X not set
  plan for 10 CPU-years                      if P(>=1 at 10 CPU-years) >= X: P = 23.5%: X not set
```
