#!/usr/bin/env python3
"""Propagate the calib-target measurements into scheduler v2's forecast.

1. The class-factor GLMs (squares, S and P traversals) refit on the new sums
   with the shipped priors as the prior (Calibration.fit: what the scheduler
   would learn from these records on a fresh state, with the scheduler's own
   quasi-Poisson dispersion and stream weights).
2. The greedy forecast (as `forecast --shipped --truth anchored`, the
   calibration replaced, 10% of the candidates, seed 1, budget scaled) run
   to 10 CPU-years for: shipped, updated; and the updated calibration with
   the time truth corrected by the measured CPU (variant "time").
3. Per-sum E by cell at 1 and 10 CPU-years, for the bands from draws.

usage: forecast_update.py VARIANT   (shipped | updated | time | glm)
writes analysis/fc_<VARIANT>.npz and prints E at the marks
"""
import argparse
import json
import math
import os
import sys
import time

import numpy as np

CT = "/tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/calib-target"
WT = "/home/user/magic-squares/.claude/worktrees/calib-target"
sys.path.insert(0, os.path.join(WT, "scripts"))
import scheduler as S  # noqa: E402
import amodel as am  # noqa: E402

PRED = json.load(open(os.path.join(CT, "predictions.json")))
AN = json.load(open(os.path.join(CT, "analysis", "results.json")))
OBS = {o["id"]: o for o in AN["per_run"]}


def load(path):
    out = {}
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            out.setdefault(r["type"], []).append(r)
    return out


def cell_table(phi_between_A=None):
    """[NCELLS, 9] (sums, O_sq, E_sq, O_S, E_S, O_P, E_P, O_SP, E_SP) and ee
    [NCELLS, 3], expectations before the class factors (shipped prior), as
    Summary._ingest builds them"""
    lg0, lS0, lP0, _ = S.Calibration().rates()
    T = np.zeros((S.NCELLS, 9))
    ee = np.zeros((S.NCELLS, 3))
    for p in PRED["runs"]:
        c = int(p["cell"])
        R = load(p["out"])
        if p["mode"] == "plain":
            sq = R.get("square", [])
            n = len(sq)
            w = 1.0
            o_sq, e_sq = n, p["squares_sum"]
            nrec = n
            oS, oP, oSP = (sum(q[k] for q in sq) for k in ("s_count", "p_count", "sp_count"))
        else:
            cs = R["csum"][0]
            k = cs["r1_stride"]
            nrec = cs["squares"]
            phi_s = (cs["se_squares"] ** 2 / (k * cs["est_squares"])) if nrec >= 5 else S.CALIB_PHI
            w = min(1.0, S.PHI_SUM / (phi_s + S.PHI_SUM / k))
            o_sq, e_sq = nrec, p["squares_sum"] / k
            oS, oP, oSP = cs["s_trav"], cs["p_trav"], cs["sp_pairs"]
        eb_sq = w * e_sq / math.exp(lg0[c])
        eb_S = w * nrec * p["S_trav_per_square"] / math.exp(lS0[c])
        eb_P = w * nrec * p["P_trav_per_square"] / math.exp(lP0[c])
        T[c] += [1, w * o_sq, eb_sq, w * oS, eb_S, w * oP, eb_P, w * oSP, w * nrec * p["SP_trav_per_square"]]
        ee[c] += [eb_sq ** 2, eb_S ** 2, eb_P ** 2]
    return T, ee


def fit_updated():
    T, ee = cell_table()
    cal = S.Calibration.fit(T, ee)
    return cal, T, ee


def glm_report():
    cal0 = S.Calibration()
    cal, T, ee = fit_updated()
    out = {"effects": list(S.EFFECTS)}
    for key in ("sq", "S", "P"):
        out[key] = {"prior": [float(v) for v in cal0.beta[key]], "prior_sd": [float(math.sqrt(v)) for v in np.diag(cal0.cov[key])],
                    "post": [float(v) for v in cal.beta[key]], "post_sd": [float(math.sqrt(v)) for v in np.diag(cal.cov[key])]}
    lg0, lS0, lP0, lm0 = cal0.rates()
    lg, lS, lP, lm = cal.rates()
    # factor per N' band and class, averaged over the cells with the E shares of the plan (filled in later)
    out["cell_factor_sq"] = (lg - lg0).tolist()
    out["cell_factor_m"] = (lm - lm0).tolist()
    out["cell_sd_post_m"] = []
    rng = np.random.default_rng(3)
    d = cal.draws(rng, 400)
    out["cell_sd_post_m"] = np.std([x[1] for x in d], axis=0).tolist()
    out["cell_sd_post_sq"] = np.std([x[0] for x in d], axis=0).tolist()
    d0 = cal0.draws(rng, 400)
    out["cell_sd_prior_m"] = np.std([x[1] for x in d0], axis=0).tolist()
    out["cell_sd_prior_sq"] = np.std([x[0] for x in d0], axis=0).tolist()
    out["table_rows_with_data"] = int((T[:, 0] > 0).sum())
    out["calib_post"] = cal.to_json()
    with open(os.path.join(CT, "analysis", "glm_update.json"), "w") as f:
        json.dump(out, f)
    print("effects:", S.EFFECTS)
    for key in ("sq", "S", "P"):
        print(key)
        for e, a, sa, b, sb in zip(S.EFFECTS, out[key]["prior"], out[key]["prior_sd"], out[key]["post"], out[key]["post_sd"]):
            print(f"  {e:12s} prior {a:+.3f} ({sa:.3f})  post {b:+.3f} ({sb:.3f})")
    return out


def make_args(state, hours, sample, seed):
    return argparse.Namespace(
        state=state, vec_size=6, cmd="forecast", model="analytic", pool="wide", pool_ratio=1.2,
        pool_s0_max=6000, pool_primes=10, tau_min=1000, tau_max=12000, explore=None,
        profile_workers=2, no_legacy=False, workers=3, unit_time=120, unit_drop=S.UNIT_DROP,
        node_limit=20_000_000_000, active=300, only=None, dfirst="auto",
        dfirst_min_n=S.DFIRST_MIN_NP, calib_frac=S.CALIB_FRAC, hours=hours, sample=sample,
        seed=seed, draws=0, shipped=True, truth="anchored")


# measured CPU / law, by N' band (analysis/results.json "time"): plain sums
# against the plain law, d-first sums against the d-first law (variant "time")
TIME_PLAIN = {0: 1.0, 1: 1.0, 2: None, 3: None, 4: None, 5: None}
TIME_DF = {}


def time_factors():
    tm = AN["time"]
    fp = {2: tm["plain 3-6k"]["tp"]["sum_ratio"], 3: tm["plain 6-12k"]["tp"]["sum_ratio"]}
    # plain at >= 12k: the streams' est_time against the plain law
    fp[4] = tm["dfirst 12-24k: plain whole sum (stream est_time)"]["tp"]["sum_ratio"]
    fp[5] = tm["dfirst 24-45k: plain whole sum (stream est_time)"]["tp"]["sum_ratio"]
    fd = {3: tm["dfirst 6-12k: d-first whole sum"]["td"]["sum_ratio"],
          4: tm["dfirst 12-24k: d-first whole sum"]["td"]["sum_ratio"],
          5: tm["dfirst 24-45k: d-first whole sum"]["td"]["sum_ratio"]}
    # d-first at 3-6k: not measured; the ratio law held (0.95), so the plain factor x 1
    fd[2] = fp[2]
    return fp, fd


def corrected_charge(sc, u, fp, fd):
    """the measured truth: plain sums at the plain law x the measured factor of
    their N' band, d-first sums at the d-first law x its factor, streams at
    the corrected plain cost / k"""
    if u.dlo is not None:
        Ss = np.array([float(u.lo)])
    else:
        Ss = np.arange(u.lo, u.hi + 1, dtype=float)
    sq, m, tp, td, tc, dm, lNp, lL, cell = sc.eval_modes(u.a, Ss)
    nb = np.searchsorted(np.log(S.NBAND_EDGES), lNp, side="right")
    cp = np.array([fp.get(int(b), 1.0) for b in nb])
    cd = np.array([fd.get(int(b), 1.0) for b in nb])
    tpc = tp * cp
    tdc = td * cd
    Sg = sc.prof(u.a)[0]
    below = Ss < (Sg[0] - 1e-9) if len(Sg) else np.ones(len(Ss), bool)
    tpc[below] = tdc[below] = am.SUM_OVERHEAD
    if u.mode == "dfirst":
        c = float((u.frac * tdc + (tpc / u.calib if u.calib else 0.0)).sum())
    else:
        c = float(tpc.sum())
    return c + am.UNIT_OVERHEAD


def run(variant, hours_tot=10 * 8766.0, frac=0.1, seed=1):
    args = make_args(os.path.join(CT, "state"), hours_tot, frac, seed)
    sch = S.SchedulerV2(args, quiet=True)
    cal = S.Calibration() if variant == "shipped" else fit_updated()[0]
    sch.calib = cal
    sch.tm = am.time_prior(S.ENGINE)
    sch.tmd = am.time_prior(S.ENGINE, mode="dfirst")
    sch.lr0 = am.dfirst_ratio_level(None)
    sc = sch.scorer
    sc.set_calibration(cal)
    sc.tm, sc.tmd, sc.lr0 = sch.tm, sch.tmd, sch.lr0
    sc.lnfsq[:] = 0.0
    sc.lnFm[:] = 0.0
    rng = np.random.default_rng(seed)
    drop = rng.random(len(sch.cands)) >= frac
    sc.lnfsq[drop] = -np.inf
    budget = hours_tot * frac
    yr = 8766.0
    fp, fd = time_factors()
    hours = magic = 0.0
    acc = {}
    marks = [0.1 * yr, yr, 3 * yr, 10 * yr]
    E_at = {}
    t0 = time.time()
    for u in sch.simulate(max_hours=None, detail=True):
        if hours >= budget:
            break
        cu = corrected_charge(sc, u, fp, fd) if variant == "time" else S.unit_charge(sc, u, "anchored")
        h_sc = hours / frac
        if u.dlo is not None:
            sums, mv = [u.lo], [u.magic]
        else:
            sums = list(range(u.lo, u.hi + 1))
            mv = list(u.mvec) if u.mvec is not None and len(u.mvec) == len(sums) else [u.magic / len(sums)] * len(sums)
        for s_, m_ in zip(sums, mv):
            e = acc.get((u.a, s_))
            if e is None:
                acc[(u.a, s_)] = [float(m_), h_sc]
            else:
                e[0] += float(m_)
        hours += cu / 3600
        magic += u.magic
        while marks and hours / frac >= marks[0]:
            E_at[marks[0] / yr] = magic / frac
            print(f"  {variant}: {marks[0] / yr:g} CPU-years: E = {magic / frac:.4g} ({time.time() - t0:.0f} s)",
                  flush=True)
            marks.pop(0)
    keys = np.array(list(acc.keys()), np.int64)
    vals = np.array(list(acc.values()), float)
    a = keys[:, 0]
    Ss = keys[:, 1].astype(float)
    cell = np.zeros(len(a), np.int64)
    order = np.argsort(a, kind="stable")
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and a[order[j]] == a[order[i]]:
            j += 1
        idx = order[i:j]
        cell[idx] = sc.eval_modes(int(a[idx[0]]), Ss[idx])[8]
        i = j
    np.savez_compressed(os.path.join(CT, "analysis", f"fc_{variant}.npz"), magic=vals[:, 0], h=vals[:, 1],
                        cell=cell, frac=frac, E_at=np.array(sorted(E_at.items())))
    print(f"{variant}: done, {len(acc)} sums, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    v = sys.argv[1]
    if v == "glm":
        glm_report()
    else:
        run(v)
