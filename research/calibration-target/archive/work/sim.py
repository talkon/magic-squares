#!/usr/bin/env python3
"""Step 1 of the calibration plan: where E sits for 1-10 CPU-years.

Runs scheduler v2's greedy (shipped calibration, --dfirst auto, each unit
charged under the anchored truth, as `forecast --shipped --truth anchored`)
on a uniform sample of the candidates with the budget scaled, and writes
per planned sum (P, S): its predicted magic squares (summed over the
units of the sum), the charged CPU-hours (scaled) at which it was first
planned, and its mode. Output: work/plan_sums_<tag>.npz and a composition
table on stdout.

usage: sim.py HOURS SAMPLE SEED TAG
"""
import argparse
import collections
import math
import os
import sys
import time

import numpy as np

ROOT = "/home/user/magic-squares"
sys.path.insert(0, os.path.join(ROOT, "scripts"))
import scheduler as S  # noqa: E402
import amodel as am  # noqa: E402

CT = "/tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/calib-target"


def make_args(state, hours=8766.0, sample=1.0, seed=1):
    return argparse.Namespace(
        state=state, vec_size=6, cmd="forecast", model="analytic", pool="wide", pool_ratio=1.2,
        pool_s0_max=6000, pool_primes=10, tau_min=1000, tau_max=12000, explore=None,
        profile_workers=2, no_legacy=False, workers=3, unit_time=120, unit_drop=S.UNIT_DROP,
        node_limit=20_000_000_000, active=300, only=None, dfirst="auto",
        dfirst_min_n=S.DFIRST_MIN_NP, calib_frac=S.CALIB_FRAC, hours=hours, sample=sample,
        seed=seed, draws=0, shipped=True, truth="anchored")


def shipped(sch):
    """as v2_forecast --shipped"""
    sch.calib = S.Calibration()
    sch.tm = am.time_prior(S.ENGINE)
    sch.tmd = am.time_prior(S.ENGINE, mode="dfirst")
    sch.lr0 = am.dfirst_ratio_level(None)
    sc = sch.scorer
    sc.set_calibration(sch.calib)
    sc.tm, sc.tmd, sc.lr0 = sch.tm, sch.tmd, sch.lr0
    sc.lnfsq[:] = 0.0
    sc.lnFm[:] = 0.0


def main():
    hours_tot, frac, seed, tag = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    args = make_args(os.path.join(CT, "state"), hours_tot, frac, seed)
    sch = S.SchedulerV2(args, quiet=True)
    shipped(sch)
    sc = sch.scorer
    rng = np.random.default_rng(seed)
    drop = rng.random(len(sch.cands)) >= frac
    sc.lnfsq[drop] = -np.inf
    budget = hours_tot * frac
    yr = 8766.0
    hours = magic = 0.0
    acc = {}          # (a, S) -> [magic, first hours (scaled), mode dfirst?, cpu charged]
    t0 = time.time()
    n = 0
    marks = [0.1 * yr, yr, 3 * yr, 10 * yr, 30 * yr, 100 * yr]
    E_at = {}
    for u in sch.simulate(max_hours=None, detail=True):
        if hours >= budget:
            break
        cu = S.unit_charge(sc, u, "anchored")
        h_sc = hours / frac
        if u.dlo is not None:
            sums = [u.lo]
            mv = [u.magic]
            cpu = [cu]
        else:
            sums = list(range(u.lo, u.hi + 1))
            mv = list(u.mvec) if u.mvec is not None and len(u.mvec) == len(sums) else \
                [u.magic / len(sums)] * len(sums)
            tot = sum(mv) or 1.0
            cpu = [cu * m / tot for m in mv]
        for s_, m_, c_ in zip(sums, mv, cpu):
            e = acc.get((u.a, s_))
            if e is None:
                acc[(u.a, s_)] = [float(m_), h_sc, u.mode == "dfirst", float(c_)]
            else:
                e[0] += float(m_)
                e[3] += float(c_)
        hours += cu / 3600
        magic += u.magic
        n += 1
        while marks and hours / frac >= marks[0]:
            E_at[marks[0]] = magic / frac
            print(f"  {marks[0] / yr:g} CPU-years: E = {magic / frac:.4g}, units {n / frac:.0f}, "
                  f"{time.time() - t0:.0f} s", flush=True)
            marks.pop(0)
    E_at[hours / frac] = magic / frac
    keys = np.array(list(acc.keys()), np.int64)
    vals = np.array(list(acc.values()), float)
    a = keys[:, 0]
    Ss = keys[:, 1].astype(float)
    # per-sum details at the shipped calibration
    lNp = np.zeros(len(a))
    cell = np.zeros(len(a), np.int64)
    order = np.argsort(a, kind="stable")
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and a[order[j]] == a[order[i]]:
            j += 1
        idx = order[i:j]
        out = sc.eval_modes(int(a[idx[0]]), Ss[idx])
        lNp[idx] = out[6]
        cell[idx] = out[8]
        i = j
    P = np.array([sch.cands.exps[x] for x in a], np.uint8)
    np.savez_compressed(os.path.join(CT, "work", f"plan_sums_{tag}.npz"), a=a, S=Ss, magic=vals[:, 0],
                        h=vals[:, 1], dfirst=vals[:, 2].astype(bool), cpu=vals[:, 3], lNp=lNp,
                        cell=cell, P=P, ratio=sch.cands.ratio[a], k=sch.cands.k[a],
                        tau=sch.cands.tau[a], S0=sch.cands.S0[a], smin=sch.cands.smin[a],
                        frac=frac, E_at=np.array(sorted(E_at.items())))
    print(f"{n} units, {len(acc)} sums, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
