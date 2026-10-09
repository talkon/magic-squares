"""Per-(P, S) predictions of scheduler v2's model at its shipped calibration
(fresh state): the quantities the calibration plan pre-registers.

Formulas (scripts/scheduler.py AnalyticScorer.eval_modes, Calibration.rates;
scripts/amodel.py profile):

  cell c                   = cell_index(N', k, ratio bin, x = ln(S / smin_approx))
  squares  E_sq            = e^{lEs} g_sq(c) SQ12^[N' >= 12k]            (eval_modes "sq")
  P(magic | square)        = min(e^{lPm} m(c), 1e-6),  e^{lPm} = 5400 KAPPA p_pair,
                             m(c) = PAIR (r_S(c) r_P(c) / (0.86 0.64))^2
  magic                    = E_sq P(magic | square)                       (eval_modes "m")
  per-traversal rates      q_S = r_S(c) e^{lpS},  q_P = r_P(c) e^{lpP}     (Summary._ingest: the
                             S / P traversal GLMs are counts against 720 e^{lpS}, 720 e^{lpP})
  per-traversal SP rate    q_SP = k_SP (r_S(c) / 0.86) (r_P(c) / 0.64) e^{lpSP},
                             k_SP = sqrt(KAPPA / 1.1) = 0.7198 (analytic.md 4.3:
                             KAPPA = k_SP^2 x 1.1; the class factors scale it as they scale
                             P(magic) through m(c), which is quadratic in r_S r_P)
  so that P(magic | square) = 5400 q_SP^2 x pairf,  pairf = 1.1 PAIR p_pair / p_SP^2
                             (the pair factor: 1.1 x 0.87 x the lattice / correlation term)
  (square, SP diagonal) pairs of the sum = E_sq x 720 x q_SP   (= sum of sp_count over the
                             squares = the d-first "pairs" of a complete sum)
  SP+S pairs (unordered {SP, S-only} partner pairs) = pairs x 15 x (q_S - q_SP)
  SP+P pairs                                        = pairs x 15 x (q_P - q_SP)
  (i.i.d. classes over the partner graph at the cell rates, as calibrate.py's
  regression predictor; each d-first (square, d) pair has 15 partner diagonals)
  CPU:  plain law tp (TIME_PRIOR + engine-3 shift), d-first law td (DFIRST_TIME_PRIOR),
        ratio law r(N') = exp(-0.028 - 0.566 ln(N'/4000)): plain by ratio = td / r,
        d-first by ratio = tp r; anchored (unit_charge --truth anchored): d-first sum td,
        plain sum tp at N' <= 3k, td / r at >= 5k, geometric blend between.
"""
import math
import os
import sys

import numpy as np

ROOT = "/home/user/magic-squares"
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import scheduler as S  # noqa: E402
import amodel as am  # noqa: E402

K_SP = math.sqrt(am.KAPPA / 1.1)
T = 720
NPART = 15


def predict(sch, a, Ss):
    """dict of arrays over the sums Ss of candidate a"""
    sc = sch.scorer
    Ss = np.asarray(Ss, float)
    sq, m, tp, td, tcal, dm, lNp, lL, cell = sc.eval_modes(a, Ss)
    P = sch.cands.P(a)
    A, v = sch.store.get(P)
    A = np.asarray(A, float)
    Sg = am.grid_sums(P, 6)[v]
    f = {name: np.interp(Ss, Sg, A[j, v]) for j, name in enumerate(am.FIELDS)}
    lg, lS, lP, lm = sch.calib.rates()
    rS, rP = np.exp(lS[cell]), np.exp(lP[cell])
    qS = rS * np.exp(f["lpS"])
    qP = rP * np.exp(f["lpP"])
    qSP = K_SP * (rS / S.TRAV_BASE[0]) * (rP / S.TRAV_BASE[1]) * np.exp(f["lpSP"])
    pm = np.where(sq > 0, m / np.maximum(sq, 1e-300), 0.0)
    pairs = sq * T * qSP
    lr = am.dfirst_log_ratio(lNp, sc.lr0)
    r = np.exp(lr)
    w = np.clip((lNp - math.log(3000)) / (math.log(5000) - math.log(3000)), 0.0, 1.0)
    tpa = np.exp((1 - w) * np.log(tp) + w * (np.log(td) - lr))
    labels = am.labels_obs(lNp, lL, int(sch.cands.k[a]))
    return dict(S=Ss, sq=sq, magic=m, pmagic=pm, tp=tp, td=td, tcal=tcal, dfirst_auto=dm,
                lNp=lNp, Np=np.exp(lNp), lL=lL, labels=labels, cell=cell, qS=qS, qP=qP, qSP=qSP,
                pairs=pairs, sps=pairs * NPART * np.maximum(qS - qSP, 0),
                spp=pairs * NPART * np.maximum(qP - qSP, 0),
                pairf=np.where(qSP > 0, pm / (5400 * qSP ** 2), 0.0),
                r=r, tp_ratio=td / r, td_ratio=tp * r, tp_anch=tpa,
                x=np.log(Ss / sch.cands.smin[a]), lpSP=f["lpSP"], lpS=f["lpS"], lpP=f["lpP"])
