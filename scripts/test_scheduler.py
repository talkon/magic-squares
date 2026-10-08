#!/usr/bin/env python3
"""
Tests of scheduler.py and of scheduler v2 (amodel.py, pool.py): the
regression fit on synthetic data, coverage tracking, the analytic model's
reference numbers and guards, the exact pool generator, the profile store,
the lazy planner against brute force, the learning (GLM MAP, per-P factors,
time refit), the incremental summary, and the commands end to end on a few
small values of P in temporary state directories, with both models. Needs
numpy and a build (bin/msearch; bin/enumerate for --model regression).

usage: python3 scripts/test_scheduler.py
"""
import gzip
import json
import math
import os
import random
import re
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scheduler  # noqa: E402

SMALL = "10 4 3 2;11 4 3 1 1;11 5 3 2"
REGRESSION = ("--model", "regression", "--pool", "classic")


def run(*args):
    cmd = [sys.executable, os.path.join(HERE, "scheduler.py"), *args]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        print(out.stdout, out.stderr)
        raise SystemExit(f"failed: {' '.join(cmd)}")
    return out.stdout



def test_fit_poisson():
    import numpy as np
    rng = np.random.default_rng(0)
    n = 5000
    X = np.column_stack([np.ones(n), rng.normal(5, 1, n), rng.uniform(0, 1, n)])
    beta = np.array([-4.0, 0.7, -1.5])
    expo = rng.uniform(0.5, 2, n)
    y = rng.poisson(expo * np.exp(X @ beta))
    b = scheduler.fit_poisson(X.tolist(), y.tolist(), expo.tolist(), ridge=1e-6)
    assert np.allclose(b, beta, atol=0.25), b
    groups = np.repeat(np.arange(n // 10), 10)
    yg = np.bincount(groups, weights=y)
    b = scheduler.fit_poisson(X.tolist(), yg.tolist(), expo.tolist(), groups.tolist(),
                              ridge=1e-6)
    assert np.allclose(b, beta, atol=0.5), b
    print("fit_poisson ok")


def test_fit_model():
    """fit_model on synthetic results drawn from the default model"""
    import random
    rng = random.Random(1)
    default = scheduler.Model(dict(scheduler.DEFAULT_MODEL))
    with tempfile.TemporaryDirectory() as state:
        pinfo = scheduler.PInfo(os.path.join(state, "pinfo.json"), 6)
        results = scheduler.Results(6)
        for P, smin in (((13, 6, 3, 2, 1, 1), 791), ((12, 7, 4, 2, 1), 723),
                        ((13, 7, 5, 2), 700)):
            pinfo.data[pinfo.key(P)] = {"smin": smin, "counts": {}, "counted_to": smin,
                                        "window": True}
            for i, N in enumerate(range(400, 2400, 25)):
                S = smin + 2 * i
                mu = default.squares(P, S, N, smin)
                k = sum(1 for _ in range(200) if rng.random() < mu / 200)
                results.sums[(P, S)] = {"type": "sum", "S": S, "nvecs_raw": N, "squares": k,
                                        "time": default.sum_time(N), "setup_time": 0.0,
                                        "truncated": 0}
                for j in range(k):
                    results.squares[f"{P}{S}{j}"] = {
                        "P": P, "S": S, "s_count": rng.randint(0, 3),
                        "p_count": rng.randint(0, 1), "sp_count": 0}
        model, diag = scheduler.fit_model(results, pinfo)
        assert diag["own_squares"] > 20, diag
        print(scheduler.describe_fit(diag))
        # the fitted model reproduces the data it was drawn from, on average
        P, smin = (12, 7, 4, 2, 1), 723
        a = sum(model.squares(P, smin + 40, N, smin) for N in (800, 1200, 1600))
        b = sum(default.squares(P, smin + 40, N, smin) for N in (800, 1200, 1600))
        assert 0.3 < a / b < 3, (a, b)
        sc = scheduler.Scorer(model, results, pinfo)
        assert sc.score(P, smin + 40, 1500)[0] > 0
    scheduler.CLAMP.update(scheduler.DEFAULT_MODEL["clamp"])
    print("fit_model ok")


def test_coverage():
    """searched sums are tracked as ranges: a unit that stopped early (or a
    --sums run) leaves a gap that is searched next, even when a later unit
    of the same P finished"""
    class FakePInfo:
        def smin(self, P):
            return 171
    P = (10, 4, 3, 2)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "u.jsonl")
        with open(path, "w") as f:
            for r in ({"type": "sum", "S": 171}, {"type": "sum", "S": 172},
                      {"type": "done", "min_sum": 171, "last_sum": 295, "complete": 0},
                      {"type": "sum", "S": 350},
                      {"type": "done", "min_sum": 350, "last_sum": 399, "complete": 1},
                      {"type": "sum", "S": 450}):
                f.write(json.dumps({"n": 6, "P": list(P), **r}) + "\n")
        res = scheduler.Results(6)
        res.add_file(path)
        assert res.frontier(P, FakePInfo()) == 296, res.covered
        assert res.next_covered(P, 296) == 350
        res.mark_covered(P, 296, 349)
        assert res.frontier(P, FakePInfo()) == 400
        assert res.next_covered(P, 400) == 450
        res.mark_covered(P, 400, 449)
        assert res.frontier(P, FakePInfo()) == 451
        assert res.next_covered(P, 451) is None
    print("coverage ok")




def test_commands():
    """--model regression (v1) end to end"""
    with tempfile.TemporaryDirectory() as state:
        # use the built-in default model
        model = scheduler.Model.load(os.path.join(state, "none.json"))
        assert model is not None, "no default model"
        out = run("--state", state, "plan", *REGRESSION, "--workers", "2", "--only", SMALL)
        assert "10 4 3 2" in out, out
        lines = run("--state", state, "emit", *REGRESSION, "--workers", "2", "--only", SMALL,
                    "--units", "4", "--unit-time", "1").strip().split("\n")
        assert len(lines) == 4, lines
        assert all(l.startswith("--vec-size 6 ") for l in lines), lines
        out = run("--state", state, "run", *REGRESSION, "--workers", "2", "--only", SMALL,
                  "--unit-time", "2", "--hours", "0.004")
        units = os.listdir(os.path.join(state, "units"))
        assert units, out
        sums = 0
        for name in units:
            with open(os.path.join(state, "units", name)) as f:
                sums += sum(json.loads(l)["type"] == "sum" for l in f)
        assert sums > 0
        out = run("--state", state, "report", *REGRESSION)
        assert "semi-magic squares" in out, out
        # resuming continues after the sums already searched
        before = scheduler.Results(6)
        for name in units:
            before.add_file(os.path.join(state, "units", name))
        out = run("--state", state, "plan", *REGRESSION, "--workers", "2", "--only", SMALL)
        for line in out.split("\n")[1:]:
            if line.startswith("10 4 3 2 "):
                lo = int(line.split()[4].split("-")[0])
                done = max(hi for lo_, hi in before.covered[(10, 4, 3, 2)])
                assert lo > done, line
        out = run("--state", state, "forecast", *REGRESSION, "--only", SMALL, "--hours", "1")
        assert "magic squares" in out, out
        print(f"commands (regression) ok ({len(units)} units, {sums} sums)")


# --------------------------------------------------------------------------
# scheduler v2


def test_amodel():
    """T1: the analytic model's reference numbers (analytic.md), vectorised
    vs per point, the divisor cap, the guards"""
    import numpy as np
    import amodel as am
    ref = (((10, 4, 3, 2), 327, 446, 0.03, 2.35e-4, 5.03e6),
           ((11, 4, 3, 1, 1), 289, 450, 0.03, 9.18e-4, 3.23e6),
           ((13, 7, 4, 3, 2, 1), 2820, 7540, 0.03, 414, 6.09e8),
           ((13, 7, 4, 3, 1, 1), 2200, None, 0.0, None, 3.46e9))
    for P, S, N, tolN, Es, inv_pm in ref:
        r = am.profile(P, [S])
        assert r["valid"][0], P
        if N:
            assert abs(math.exp(r["lN"][0]) / N - 1) < tolN, (P, math.exp(r["lN"][0]))
            assert abs(math.exp(r["lEs"][0]) / Es - 1) < 0.05, (P, math.exp(r["lEs"][0]))
        assert abs(math.exp(-r["lPm"][0]) / inv_pm - 1) < 0.05, (P, math.exp(-r["lPm"][0]))
    P = (13, 7, 4, 3, 0, 0, 1, 1)
    S = am.grid_sums(P)
    cap = am.CAP_FACTOR * S.max() / 6
    vec = am.profile(P, S, cap=cap)
    for j in (0, 7, 23):
        one = am.profile(P, S[j:j + 1], cap=cap)
        for key in am.FIELDS:
            assert abs(one[key][0] - vec[key][j]) < 1e-9, (key, j)
    wide = am.profile(P, S, cap=80 * S.max() / 6)
    v = vec["valid"] & wide["valid"]
    assert v.all() and np.abs(vec["lEs"] - wide["lEs"])[v].max() < 0.05
    # ill-conditioned points: below S_min, and near u = 0 with few vectors
    P = (3, 3, 2, 2, 1, 0, 1, 1)
    S = am.grid_sums(P)
    r = am.profile(P, S)
    assert not r["valid"][S < am.smin_approx(P)].any(), r["valid"]
    assert r["valid"][-1]
    P = (30, 1, 2, 1, 1, 0, 1, 1, 0, 0, 1)
    r = am.profile(P, am.grid_sums(P))
    assert not r["valid"][:3].any() and r["valid"][-1], r["valid"]
    for key in am.FIELDS:
        assert np.isfinite(r[key][r["valid"]]).all()
    assert abs(am.ratio((13, 7, 4, 3, 0, 0, 1, 1)) - 1.145) < 0.001
    assert am.ratio((13, 7, 4, 3, 1, 1)) == 1.0
    print("amodel ok")


def test_pool():
    """T2: the pool generator against brute force; filters; PRIMES"""
    import pool
    import amodel as am
    p = pool.gen_pool(8, 900, 300, 300, 3, 8, 1.3, 99, 40)
    got = {tuple(int(x) for x in r) for r in p["exps"]}
    assert len(p["exps"]) == len(got) == 50723, len(p["exps"])
    assert got == pool.brute_force(8, 900, 300, 300, 3, 8, 1.3)
    q = pool.gen_pool(10, 2500, 1500, 1000, ratio_max=1.2)
    for row, S0, r, t, k in zip(q["exps"][::97], q["S0"][::97], q["ratio"][::97], q["tau"][::97],
                                q["k"][::97]):
        P = scheduler.norm_p(int(x) for x in row)
        assert 1500 <= am.s0(P) <= 2500 and abs(am.s0(P) / S0 - 1) < 1e-5
        assert am.tau(P) == t >= 1000 and 4 <= k == sum(1 for a in P if a) <= 10
        assert abs(am.ratio(P) - r) < 1e-5 and r <= 1.2 + 1e-6
        assert max(P) <= 31 and sum(1 for a in P if a == 1) <= 6
    rows = {tuple(int(x) for x in r) for r in q["exps"]}
    assert (13, 7, 4, 3, 0, 0, 1, 1, 0, 0) in rows and (7, 13, 4, 3, 1, 1, 0, 0, 0, 0) not in rows
    assert len(rows) == len(q["exps"])
    # the classic pool (gen_candidates) is unchanged
    c = scheduler.gen_candidates(6, 800, 12000)
    assert len(c) == 10075 and c[0] == (8, 3, 2, 1, 0, 0, 1, 1) and c[-1] == (18, 9, 6, 4)
    # PRIMES matches the C build
    assert len(scheduler.PRIMES) == 10 and scheduler.p_value((0,) * 9 + (1,)) == 29
    hdr = os.path.join(os.path.dirname(HERE), "src", "c", "enumerate.h")
    if os.path.exists(hdr):
        m = re.search(r"#define ENUM_MAX_PRIMES (\d+)", open(hdr).read())
        assert int(m.group(1)) == len(scheduler.PRIMES)
    assert scheduler.PRIMES == am.PRIMES[:10] == pool.PRIMES[:10]
    print("pool ok")


def small_pool(n=60, seed=3):
    import pool
    p = pool.gen_pool(10, 1400, 700, 1000, ratio_max=1.2)
    rng = random.Random(seed)
    rows = sorted(rng.sample(range(len(p["exps"])), n))
    return [scheduler.norm_p(int(x) for x in p["exps"][r]) for r in rows]


def test_profile_store():
    """T3: compute, reload, float16 round trip, on demand, version bump"""
    import numpy as np
    import amodel as am
    Ps = small_pool(30)
    pool = scheduler.Pool(scheduler.Pool.arrays_for(Ps), "test", {"n": 30})
    with tempfile.TemporaryDirectory() as state:
        st = scheduler.ProfileStore(state, pool, 6)
        st.fill(range(20), workers=1, quiet=True)
        assert st.done[:20].all() and not st.done[20:].any()
        a5 = np.array(st.arr[5], np.float32)
        st2 = scheduler.ProfileStore(state, pool, 6)
        assert (np.asarray(st2.arr[5], np.float32) == a5).all() and st2.done[:20].all()
        ref = am.profile(Ps[5], am.grid_sums(Ps[5]))
        v = ref["valid"]
        for j, key in enumerate(am.FIELDS):
            assert np.abs(a5[j][v] - ref[key][v]).max() < 0.02, key
        A, vv = st2.get(Ps[25])     # not computed yet: on demand
        assert st2.done[25] and (vv == am.profile(Ps[25], am.grid_sums(Ps[25]))["valid"]).all()
        A, vv = st2.get((12, 6, 3, 2, 1, 1))   # outside the pool: extra
        assert vv.any() and "12_6_3_2_1_1" in st2.extra
        st2.save_extra()
        old = am.AMODEL_VERSION
        am.AMODEL_VERSION = old + 1000
        try:
            st3 = scheduler.ProfileStore(state, pool, 6)
            assert not np.asarray(st3.done).any() and not st3.extra
        finally:
            am.AMODEL_VERSION = old
    print("profile store ok")


def make_v2(Ps, state, f0=None, cover=None, unit_time=120.0, drop=0.5):
    """a candidate set, scorer and planner for the P in Ps (extras store)"""
    import numpy as np
    store = scheduler.ProfileStore(state, None, 6)
    store.fill_extra(Ps)
    arr = scheduler.Pool.arrays_for(Ps)
    cands = scheduler.Cands(arr["exps"], [-1] * len(Ps), arr["ratio"])
    sc = scheduler.AnalyticScorer(cands, store, scheduler.Calibration(), _am_tm())
    if f0 is None:
        f0 = np.ceil(cands.S0 - 1e-9).astype(np.int64)
    plan = scheduler.Planner(cands, sc, f0, cover or {}, unit_time, drop)
    return cands, sc, plan


def _am_tm():
    import amodel as am
    return am.TimeModel()


def test_planner():
    """T4: UB >= exact, lazy heap == brute force, stopping rules, frontier
    without enumerate"""
    import numpy as np
    Ps = small_pool(200, seed=5)
    with tempfile.TemporaryDirectory() as state:
        cands, sc, plan = make_v2(Ps, state)
        ub = plan.upper_bounds(np.arange(len(Ps)))
        for a in range(len(Ps)):
            u = plan.unit(a)
            assert u is None or u.score <= ub[a] * (1 + 1e-9), (Ps[a], u.score, ub[a])
        # lazy greedy == brute-force argmax for the first 50 units
        plan.build()
        lazy = []
        for _ in range(50):
            u = plan.pop()
            lazy.append((u.a, u.lo, u.hi))
            plan.advance(u.a, u.hi)
            plan.push(u.a)
        _, _, bp = make_v2(Ps, state)
        brute = []
        for _ in range(50):
            us = [bp.unit(a) for a in range(len(Ps))]
            best = max((u for u in us if u is not None), key=lambda u: (u.score, -u.a))
            brute.append((best.a, best.lo, best.hi))
            bp.advance(best.a, best.hi)
        assert lazy == brute, (lazy[:5], brute[:5])
    # stopping rules on one P
    P = (12, 6, 3, 2, 1, 1)
    with tempfile.TemporaryDirectory() as state:
        cands, sc, plan = make_v2([P], state, unit_time=1e9, drop=0.5)
        f = int(plan.frontier[0])
        assert f == math.ceil(scheduler._am().s0(P))
        u = plan.unit(0)
        S = np.arange(u.lo, u.hi + 2)
        sq, m, t, _, _, _ = sc.eval_sums(0, S)
        dens = m / t
        assert dens[-1] < 0.5 * dens[:-1].max() and (dens[1:-1] >= 0.5 * np.maximum.accumulate(
            dens[:-2])).all()
        plan.unit_time, plan.drop = 30.0, 0.0
        u = plan.unit(0)
        _, _, t, _, _, _ = sc.eval_sums(0, np.arange(u.lo, u.hi + 2))
        assert t[:-1].sum() <= 30 < t.sum() and abs(u.time - t[:-1].sum() - 0.05) < 1e-9
        z = int(math.ceil(sc.prof(0)[0][0]))   # first sum with predicted squares
        assert z > f + 2
        plan.set_cover(0, [[f, z + 20], [z + 26, z + 30]])
        assert plan.frontier[0] == z + 21 and plan.unit(0).hi == z + 25
        # a gap where the model predicts nothing (below S_min) is skipped
        plan.set_cover(0, [[f + 1, z + 20]])
        assert plan.frontier[0] == z + 21
        plan.set_cover(0, [[z + 5, z + 20]])
        assert plan.frontier[0] == f and plan.unit(0).hi <= z + 4
    print("planner ok")


def test_fit_poisson_map():
    """T5: recovers synthetic effects within 2 se; no data -> prior means"""
    import numpy as np
    rng = np.random.default_rng(2)
    Z = scheduler.cell_design()
    beta = np.array([0.3, 0.2, -0.4, 0.1, 0.0, 0.0, -0.2, -0.3, 0.0, 0.1, 0.2, -0.15, 0.05])
    E = rng.uniform(50, 400, len(Z))
    y = rng.poisson(E * np.exp(Z @ beta))
    mean, sd = np.zeros(13), np.full(13, 10.0)
    b, C = scheduler.fit_poisson_map(Z, y, np.log(E), mean, sd)
    se = np.sqrt(np.diag(C))
    used = np.abs(Z).sum(0) > 0
    assert (np.abs(b - beta)[used] < 2.5 * se[used]).all(), (b, beta, se)
    for key, (m, s) in scheduler.GLM_PRIORS.items():
        b, C = scheduler.fit_poisson_map(Z[:0], [], [], m, s)
        assert np.allclose(b, m) and np.allclose(np.diag(C), np.array(s) ** 2)
    # the shipped calibration = the GLMs with no data: review factors on P(magic)
    cal = scheduler.Calibration.fit(np.zeros((scheduler.NCELLS, 9)))
    lg, lS, lP, lm = cal.rates()
    lg0, _, _, lm0 = scheduler.Calibration().rates()
    assert np.allclose(lg, 0) and np.allclose(lm, lm0)
    import amodel as am
    c_low = scheduler.cell_index(math.log(2000), 6, 0, 0.2)
    c_hi = scheduler.cell_index(math.log(4000), 7, 0, 0.05)
    assert abs(math.exp(lm[c_low]) - am.PAIR) < 1e-9
    assert abs(math.exp(lm[c_hi]) / (am.PAIR * am.NTR * am.K7 * am.X01) - 1) < 0.01
    print("fit_poisson_map ok")


def test_per_p():
    """T6: per-P gamma factors; no exploration by default"""
    assert scheduler.EXPLORE_V2 == 0
    f_sq, f_S, f_P = scheduler.per_p_factors(0, 0, 0, 0, 0, 0)
    assert f_sq == f_S == f_P == 1
    f_sq, f_S, f_P = scheduler.per_p_factors(60, 20, 30, 30, 0, 30)
    assert abs(f_sq - 100 / 60) < 1e-12 and abs(f_S - 1) < 1e-12 and abs(f_P - 0.5) < 1e-12
    hi = scheduler.per_p_factors(60, 20, 30, 30, 0, 30, explore=1.0)
    assert hi[0] > f_sq and hi[2] > f_P
    with tempfile.TemporaryDirectory() as state:
        P = (12, 6, 3, 2, 1, 1)
        cands, sc, plan = make_v2([P], state)
        c = int(scheduler.cell_index(math.log(2000), 6, 0, 0.2))   # class rates 0.86, 0.64
        summ = SimpleNamespace(perP={"12_6_3_2_1_1": {"o": [10, 30, 0, 0], "cells": {
            c: [5, 2, 6.0, 30, 30.0, 0, 30.0, 0, 0.1]}}})
        sc.set_factors(summ)
        assert abs(math.exp(sc.lnfsq[0]) - 50 / 46) < 1e-9
        f_S, f_P = 60 / (30 + 30 * 0.86), 30 / (30 + 30 * 0.64)
        Fm = f_S ** 2 * f_P ** 2 * (1 + 1 / 60) * (1 + 1 / 30) / (1 + 1 / 30) ** 2
        assert abs(math.exp(sc.lnFm[0]) - Fm) < 1e-9
    print("per-P factors ok")


def test_time_model():
    """T7: an engine at 0.5x the law is learned from 50 sums; the hinge
    stays at its prior without sums above N' = 8k"""
    import numpy as np
    import amodel as am
    rng = np.random.default_rng(4)
    lNp = np.log(rng.uniform(1000, 7000, 50))
    lL = np.log(rng.uniform(80, 200, 50))
    prior = am.TimeModel()
    y = prior.log_time(lNp, lL) + math.log(0.5) + rng.normal(0, 0.17, 50)
    st = am.TimeModel.stats(am.time_features(lNp, lL), y)
    models = scheduler.fit_time_models({"3:plain": st})
    tm = scheduler.current_time_model(models, engine=3)
    ratio = np.exp(tm.log_time(lNp, lL) - prior.log_time(lNp, lL)).mean()
    assert abs(ratio - 0.5) < 0.05, ratio
    assert tm.th[3] == prior.th[3]
    # no data for the newest engine: the previous engine's posterior
    assert np.allclose(scheduler.current_time_model(models, engine=4).th, tm.th)
    assert np.allclose(scheduler.current_time_model({}).th, prior.th)
    print("time model ok")


def write_unit(path, P, sums, squares=(), done=None, partial=False):
    recs = []
    for S, sq in sums:
        recs.append({"type": "sum", "n": 6, "P": list(P), "S": S, "nvecs": 2000,
                     "nvecs_raw": 2000 + S, "labels": 120, "nodes": 10 ** 6, "squares": sq,
                     "time": 0.004, "setup_time": 0.001, "enum_time": 0.001, "truncated": 0,
                     "engine": 2})
    for S, s, p, sp, best in squares:
        recs.append({"type": "square", "n": 6, "P": list(P), "S": S, "s_count": s,
                     "p_count": p, "sp_count": sp, "best_score": best, "hash": f"{S}{s}{p}",
                     "grid": [[1] * 6] * 6})
    if done:
        recs.append({"type": "done", "n": 6, "P": list(P), "min_sum": done[0],
                     "last_sum": done[1], "complete": 1, "time": 3.0})
    text = "".join(json.dumps(r) + "\n" for r in recs)
    if partial:
        text = text[:-20]
    with open(path, "w") as f:
        f.write(text)
    return text


def summary_state(s):
    return json.dumps({"cover": s.cover, "perP": s.perP, "time": s.time, "tband": s.tband,
                       "nbias": s.nbias, "totals": s.totals, "types": s.types,
                       "notable": s.notable}, sort_keys=True, default=str)


def test_summary():
    """T8: incremental ingest with a partial last line, idempotent, equal to
    a full parse, .gz files"""
    P = (12, 6, 3, 2, 1, 1)
    sums = [(S, S % 3) for S in range(880, 900)]
    sqs = [(885, 2, 1, 0, 3), (890, 1, 1, 1, 7)]
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        store = scheduler.ProfileStore(state, None, 6)
        path = os.path.join(units, "u1.jsonl")
        full = write_unit(path, P, sums, sqs, done=(880, 899))
        with open(path, "w") as f:
            f.write(full[:len(full) // 2])       # a running unit: partial last line
        inc = scheduler.Summary(state, 6)
        inc.update(units, store)
        n1 = inc.totals["sums"]
        assert 0 < n1 < 20
        with open(path, "w") as f:
            f.write(full)
        inc.update(units, store)
        snap = summary_state(inc)
        inc.update(units, store)          # idempotent
        assert summary_state(inc) == snap
        ref = scheduler.Summary(state, 6)
        ref.update(units, store)
        assert summary_state(ref) == snap
        assert inc.totals["sums"] == 20 and inc.totals["squares"] == 2
        assert inc.cover["12_6_3_2_1_1"] == [[880, 899]]
        assert len(inc.notable) == 1 and inc.types == {3: 1, 7: 1}
        assert sum(v[1] for v in inc.perP["12_6_3_2_1_1"]["cells"].values()) == sum(
            sq for _, sq in sums)
        inc.save()
        again = scheduler.Summary.load(state, 6)
        assert summary_state(again) == snap
        # the same data compressed (compact) reads the same
        with tempfile.TemporaryDirectory() as state2:
            u2 = os.path.join(state2, "units")
            os.makedirs(u2)
            with gzip.open(os.path.join(u2, "u1.jsonl.gz"), "wt") as f:
                f.write(full)
            gz = scheduler.Summary(state2, 6)
            gz.update(u2, store)
            assert summary_state(gz) == snap
        res = scheduler.Results(6)
        res.add_file(path)
        assert len(res.sums) == 20
    print("summary ok")


def test_no_enumerate():
    """A9: the analytic paths never call bin/enumerate (or any subprocess
    other than msearch)"""
    import numpy as np  # noqa: F401
    src = open(os.path.join(HERE, "scheduler.py")).read()
    v2 = src[src.index("# scheduler v2 (--model analytic"):src.index("\ndef main():")]
    for bad in ("PInfo", "pinfo.smin", '"enumerate"', ".counts(", "Results("):
        assert bad not in v2, f"the v2 code uses {bad}"
    saved_run, saved_bin = subprocess.run, scheduler.binary

    def forbidden(*a, **k):
        raise AssertionError(f"subprocess.run called: {a}")

    def binary(name):
        assert name == "msearch", name
        return saved_bin(name)

    subprocess.run = forbidden
    scheduler.binary = binary
    try:
        with tempfile.TemporaryDirectory() as state:
            os.makedirs(os.path.join(state, "units"))
            write_unit(os.path.join(state, "units", "u.jsonl"), (10, 4, 3, 2),
                       [(S, 0) for S in range(171, 200)], done=(171, 199))
            args = SimpleNamespace(state=state, vec_size=6, model="analytic", pool="wide",
                                   pool_ratio=1.2, pool_s0_max=6000, pool_primes=10,
                                   tau_min=1000, tau_max=12000, explore=None, profile_workers=1,
                                   no_legacy=False, workers=1, unit_time=60.0, unit_drop=0.5,
                                   node_limit=0, active=300, only=SMALL, top=5, units=3,
                                   hours=0.5, sample=1.0, seed=1, draws=5, shipped=False)
            import contextlib
            import io
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                scheduler.v2_plan(args)
                scheduler.v2_emit(args)
                scheduler.v2_forecast(args)
                scheduler.v2_report(args)
                scheduler.v2_fit(args)
                scheduler.v2_profile(args)
            out = buf.getvalue()
            assert "10 4 3 2" in out and "--node-limit" in out, out
            plan_line = [l for l in out.split("\n") if l.startswith("10 4 3 2 ")][0]
            assert int(plan_line.split()[4].split("-")[0]) == 200, plan_line
            stats = os.path.join(state, "stats.txt")
            with open(stats, "w") as f:
                f.write("Statistics for each P, sorted by factorization\n header\n"
                        "  10 4 3 2 508032000 171 200 250 9 1.0 3 2 1 0 0 0 0 0\n")
            with contextlib.redirect_stdout(buf):
                scheduler.cmd_import_legacy(SimpleNamespace(stats=stats, vec_size=6, state=state,
                                                            workers=1, model="analytic"))
            assert "imported 1 values" in buf.getvalue(), buf.getvalue()
    finally:
        subprocess.run = saved_run
        scheduler.binary = saved_bin
    print("no enumerate ok")


def test_commands_v2():
    """T9: the analytic scheduler end to end (defaults, --only SMALL)"""
    with tempfile.TemporaryDirectory() as state:
        out = run("--state", state, "plan", "--only", SMALL)
        assert "10 4 3 2" in out, out
        lines = run("--state", state, "emit", "--only", SMALL, "--units", "4",
                    "--unit-time", "1").strip().split("\n")
        assert len(lines) == 4 and all(l.startswith("--vec-size 6 ") for l in lines), lines
        assert all("--node-limit" in l and "--time-limit" in l for l in lines)
        out = run("--state", state, "run", "--workers", "2", "--only", SMALL,
                  "--unit-time", "2", "--hours", "0.004")
        units = os.listdir(os.path.join(state, "units"))
        assert units, out
        launched = open(os.path.join(state, "launched_6.jsonl")).read().strip().split("\n")
        assert len(launched) == len(units), (launched, units)
        first = [json.loads(l) for l in launched]
        lo = {tuple(r["P"]): r["lo"] for r in reversed(first)}
        assert lo[(10, 4, 3, 2)] == 170, lo          # ceil(S0), no S_min needed
        out = run("--state", state, "report", "--only", SMALL)
        assert "semi-magic squares" in out and "calibration" in out, out
        summ = scheduler.Summary.load(state, 6)
        done = max(hi for lo_, hi in summ.cover["10_4_3_2"])
        out = run("--state", state, "plan", "--only", SMALL)
        for line in out.split("\n")[1:]:
            if line.startswith("10 4 3 2 "):
                assert int(line.split()[4].split("-")[0]) > done, line
        out = run("--state", state, "forecast", "--only", SMALL, "--hours", "1")
        assert "magic squares" in out, out
        out = run("--state", state, "profile", "--only", SMALL)
        assert "extra P" in out, out
        out = run("--state", state, "pool", "--stats")
        assert "377908 values of P" in out and "S0" in out, out
        out = run("--state", state, "compact", "--only", SMALL)
        assert "compressed" in out
        out = run("--state", state, "report", "--only", SMALL)
        assert "semi-magic squares" in out
        print(f"commands (analytic) ok ({len(units)} units)")


if __name__ == "__main__":
    t0 = time.time()
    test_fit_poisson()
    test_fit_model()
    test_coverage()
    test_amodel()
    test_pool()
    test_profile_store()
    test_planner()
    test_fit_poisson_map()
    test_per_p()
    test_time_model()
    test_summary()
    test_no_enumerate()
    test_commands_v2()
    test_commands()
    print(f"all ok ({time.time() - t0:.0f} s)")
