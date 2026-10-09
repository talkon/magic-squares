#!/usr/bin/env python3
"""
Tests of scheduler.py and of scheduler v2 (amodel.py, pool.py): the
regression fit on synthetic data, coverage tracking, the analytic model's
reference numbers and guards, the exact pool generator, the profile store,
the lazy planner against brute force, the learning (GLM MAP, per-P factors,
time refit), the incremental summary, d-first units (mode choice, the
d-first time law, d-range units merged and resumed, the calibration
streams in the cells), and the commands end to end on a few
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


def test_dfirst_records():
    """msearch --diag-first output: the d-first sums have no "sum" record but
    are not empty sums, so they stay out of the fits (sums, squares), and
    only a sum searched in full counts as covered, also inside a "done"
    range (a --d-stride sample or a --d-range part does not)"""
    class FakePInfo:
        def smin(self, P):
            return 171
    P = (10, 4, 3, 2)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "u.jsonl")
        with open(path, "w") as f:
            for r in ({"type": "sum", "S": 171, "squares": 0, "time": 1.0, "setup_time": 0.0},
                      {"type": "dsum", "mode": "dfirst", "S": 172, "complete": 1,
                       "time": 2.0},
                      {"type": "dsquare", "S": 172, "hash": "ab", "magic": 1,
                       "partner": 1},
                      {"type": "dsum", "mode": "dfirst", "S": 174, "complete": 0,
                       "time": 1.0},
                      # an older d-first record without "complete": partial
                      {"type": "dsum", "S": 176, "time": 1.0},
                      {"type": "done", "min_sum": 171, "last_sum": 180, "complete": 1,
                       "mode": "dfirst", "diag_first_min_n": 0}):
                f.write(json.dumps({"n": 6, "P": list(P), **r}) + "\n")
        res = scheduler.Results(6)
        res.add_file(path)
        assert list(res.sums) == [(P, 171)], res.sums
        assert not res.squares and list(res.dmagic) == ["ab"]
        assert sorted(res.dsums) == [(P, 172), (P, 174), (P, 176)]
        assert res.frontier(P, FakePInfo()) == 174, res.covered
        assert res.next_covered(P, 174) == 175
        res.mark_covered(P, 174, 174)
        assert res.frontier(P, FakePInfo()) == 176
        # the v1 report lists the magic squares found d-first
        import contextlib
        import io
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            scheduler.report(res, FakePInfo(), None)
        assert "1 magic squares" in buf.getvalue() and "MAGIC (d-first) P=10 4 3 2" in buf.getvalue(), \
            buf.getvalue()
    assert scheduler.split_range(1, 9, [3, 4, 9, 12]) == [(1, 2), (5, 8)]
    # end to end: complete d-first sums are covered, a d-sampled one is not
    msearch = os.path.join(os.path.dirname(HERE), "bin", "msearch")
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "u.jsonl")
        base = [msearch, "--diag-first", "--diag-first-min-n", "0", "--out", path]
        subprocess.run(base + ["--min-sum", "327", "--max-sum", "328", "10", "4", "3", "2"],
                       check=True, capture_output=True)
        subprocess.run(base + ["--d-stride", "3", "--min-sum", "329", "--max-sum", "330",
                               "10", "4", "3", "2"], check=True, capture_output=True)
        res = scheduler.Results(6)
        res.add_file(path)
        assert not res.sums and sorted(S for _, S in res.dsums) == [327, 328, 329, 330]
        assert [r[0]["complete"] for _, r in sorted(res.dsums.items())] == [1, 1, 0, 0]

        class PI:
            def smin(self, P):
                return 327
        assert res.frontier(P, PI()) == 329, res.covered
        assert res.next_covered(P, 329) is None, res.covered
    # --d-range units: a plain sum (below --diag-first-min-n) is searched
    # by the unit whose range starts at 0 only, the others skip it
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, "a.jsonl"), os.path.join(d, "b.jsonl")
        base = [msearch, "--diag-first", "--diag-first-min-n", "455", "--min-sum", "327",
                "--max-sum", "328"]
        subprocess.run(base + ["--d-range", "0:200", "--out", a, "10", "4", "3", "2"],
                       check=True, capture_output=True)
        subprocess.run(base + ["--d-range", "200:", "--out", b, "10", "4", "3", "2"],
                       check=True, capture_output=True)
        types = []
        for path in (a, b):
            with open(path) as f:
                types.append([(json.loads(l)["type"], json.loads(l)["S"])
                              for l in f if json.loads(l)["type"] in ("sum", "skip", "dsum")])
        # S = 327 has 451 vectors after reduction (plain), 328 has 460
        assert types[0] == [("sum", 327), ("dsum", 328)], types
        assert types[1] == [("skip", 327), ("dsum", 328)], types
        res = scheduler.Results(6)
        res.add_file(b)

        class PI:
            def smin(self, P):
                return 327
        assert res.frontier(P, PI()) == 327, res.covered  # skipped, partial
        res.add_file(a)
        # 327 searched by unit a; 328's two d-range parts are not merged
        # (yet), so it is not covered
        assert res.frontier(P, PI()) == 328, res.covered
    print("d-first records ok")


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
        old = am.PROFILE_VERSION
        am.PROFILE_VERSION = old + 1000
        try:
            st3 = scheduler.ProfileStore(state, pool, 6)
            assert not np.asarray(st3.done).any() and not st3.extra
        finally:
            am.PROFILE_VERSION = old
        # a changed constant without a version bump also invalidates
        old = am.KAPPA
        am.KAPPA = old * 1.01
        try:
            st4 = scheduler.ProfileStore(state, pool, 6)
            assert not np.asarray(st4.done).any()
        finally:
            am.KAPPA = old
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
    # stopping rules on one P (a coverage gap of negligible density is
    # skipped, so this test uses explicit gaps)
    P = (12, 6, 3, 2, 1, 1)
    with tempfile.TemporaryDirectory() as state:
        cands, sc, plan = make_v2([P], state, unit_time=1e9, drop=0.5)
        f = int(plan.frontier[0])
        assert f == math.ceil(scheduler._am().s0(P))
        # eval_sums (written out for speed) == grid_density (TimeModel.time)
        # at the grid points, with band offsets and the label-word step set
        sc.tm = scheduler._am().TimeModel(th=[3.8, 5.4, -3.4, 0.9, 0.4, 0.1, -0.2, 0.3, 0.05])
        Sg = sc.prof(0)[0]
        _, m, t, _, _, _ = sc.eval_sums(0, Sg)
        g = sc.grid_density([0])[0]
        assert np.allclose(np.log(m / t), g[np.isfinite(g)], atol=1e-9)
        sc.tm = _am_tm()
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
        # a gap of negligible predicted density (right at S_min) is skipped
        plan.set_cover(0, [[z + 5, z + 20]])
        neg = plan._negligible(0, f, z + 4)
        assert plan.frontier[0] == (z + 21 if neg else f)
        # a gap with real density is not
        zz = int(cands.S0[0] * 1.1)
        assert not plan._negligible(0, zz + 1, zz + 9)
        plan.set_cover(0, [[f, zz], [zz + 10, zz + 20]])
        assert plan.frontier[0] == zz + 1 and plan.unit(0).hi <= zz + 9
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
    # the shipped calibration = the GLMs with no data: the prior means
    cal = scheduler.Calibration.fit(np.zeros((scheduler.NCELLS, 9)))
    lg, lS, lP, lm = cal.rates()
    lg0, _, _, lm0 = scheduler.Calibration().rates()
    assert np.allclose(lg, lg0) and np.allclose(lm, lm0)
    import amodel as am
    c_low = scheduler.cell_index(math.log(2000), 6, 0, 0.2)
    c_hi = scheduler.cell_index(math.log(4000), 7, 0, 0.05)
    c_sm = scheduler.cell_index(math.log(800), 6, 0, 0.05)
    assert abs(math.exp(lm[c_low]) - am.PAIR) < 1e-9 and abs(lg[c_low]) < 1e-12
    # N' >= 3k: e^{2 (0.04 + 0.06)}; k >= 7: e^{-0.16}; x < 0.1: e^{-0.08}
    assert abs(math.exp(lm[c_hi]) / (am.PAIR * math.exp(0.2) * am.K7 * am.X01) - 1) < 0.01
    assert abs(lg[c_sm] - (-0.25 - 0.1)) < 1e-12
    # quasi-Poisson: the same data with dispersion 4 gives a 2x wider se
    _, C1 = scheduler.fit_poisson_map(Z, y, np.log(E), mean, sd)
    _, C4 = scheduler.fit_poisson_map(Z, y, np.log(E), mean, sd, phi=np.full(len(y), 4.0))
    r = np.sqrt(np.diag(C4) / np.diag(C1))[used]
    assert (np.abs(r - 2) < 0.05).all(), r
    print("fit_poisson_map ok")


def test_per_p():
    """T6: per-P gamma factors; no exploration by default"""
    assert scheduler.EXPLORE_V2 == 0
    f_sq, f_S, f_P = scheduler.per_p_factors(0, 0, 0, 0, 0, 0)
    assert f_sq == f_S == f_P == 1
    f_sq, f_S, f_P = scheduler.per_p_factors(60, 20, 30, 30, 0, 30)
    A, phi = scheduler.A_SQ, scheduler.PHI_SUM
    assert abs(f_sq - (A + 60 / phi) / (A + 20 / phi)) < 1e-12
    assert abs(f_S - 1) < 1e-12 and abs(f_P - 0.5) < 1e-12
    # one sum at 0.55 of its 670 predicted squares: a strong, not total, pull
    f1 = scheduler.per_p_factors(366, 670, 0, 0, 0, 0)[0]
    assert 0.55 < f1 < 0.65, f1
    hi = scheduler.per_p_factors(60, 20, 30, 30, 0, 30, explore=1.0)
    assert hi[0] > f_sq and hi[2] > f_P
    with tempfile.TemporaryDirectory() as state:
        P = (12, 6, 3, 2, 1, 1)
        cands, sc, plan = make_v2([P], state)
        c = int(scheduler.cell_index(math.log(2000), 6, 0, 0.2))   # class rates 0.86, 0.64
        # o (all squares of the P) is not used: the counts come from the
        # cells, where the expectations are
        summ = SimpleNamespace(perP={"12_6_3_2_1_1": {"o": [10, 30, 0, 0], "cells": {
            c: [5, 2, 6.0, 30, 30.0, 0, 30.0, 0, 0.1]}}})
        sc.set_factors(summ)
        assert abs(math.exp(sc.lnfsq[0]) - (A + 2 / phi) / (A + 6 / phi)) < 1e-9
        f_S, f_P = 60 / (30 + 30 * 0.86), 30 / (30 + 30 * 0.64)
        Fm = f_S ** 2 * f_P ** 2 * (1 + 1 / 60) * (1 + 1 / 30) / (1 + 1 / 30) ** 2
        assert abs(math.exp(sc.lnFm[0]) - Fm) < 1e-9
    print("per-P factors ok")


def test_time_model():
    """T7: an engine at 0.5x the law is learned from 50 sums; the shape
    (slopes, hinge) stays at its prior; a refit on many sums at N' < 3k
    with a different shape moves predictions at N' >= 3k by < 10%"""
    import numpy as np
    import amodel as am
    rng = np.random.default_rng(4)
    lNp = np.log(rng.uniform(1000, 7000, 50))
    lL = np.log(rng.uniform(80, 200, 50))
    prior = am.TimeModel()
    # engine 3 starts from the shipped (engine 2) prior with its shift
    prior3 = am.time_prior(3)
    dth = np.zeros(am.NT)
    for f, d in am.ENGINE_TIME_SHIFT.get(3, {}).items():
        dth[am.TIME_FEATURES.index(f)] = d
    assert prior3.engine == 3 and np.allclose(prior3.th - prior.th, dth)
    y = prior3.log_time(lNp, lL, 6) + math.log(0.5) + rng.normal(0, 0.17, 50)
    st = am.TimeModel.stats(am.time_features(lNp, lL, 6), y)
    models = scheduler.fit_time_models({"3:plain": st})
    tm = scheduler.current_time_model(models, engine=3)
    ratio = np.exp(tm.log_time(lNp, lL, 6) - prior3.log_time(lNp, lL, 6)).mean()
    assert abs(ratio - 0.5) < 0.05, ratio
    assert np.abs(tm.th[1:4] - prior.th[1:4]).max() < 1e-3
    # 3000 cheap sums at N' 300-3k: +0.3 on ln t, and a label slope of -3.5
    # instead of the prior's (what drifted the v2 refit by 1.4-2x)
    n = 3000
    lN2 = np.log(rng.uniform(300, 3000, n))
    lL2 = np.log(rng.uniform(70, 160, n))
    F2 = am.time_features(lN2, lL2, 6)
    y2 = (prior.log_time(lN2, lL2, 6) + 0.3 + (-3.5 - prior.th[2]) * (lL2 - math.log(110))
          + rng.normal(0, 0.25, n))
    tm2 = scheduler.fit_time_models({"2:plain": am.TimeModel.stats(F2, y2)})["2:plain"]
    lN3 = np.log(rng.uniform(3000, 20000, 200))
    lL3 = np.log(rng.uniform(120, 260, 200))
    shift = np.exp(tm2.log_time(lN3, lL3, 6) - prior.log_time(lN3, lL3, 6))
    assert np.abs(shift - 1).max() < 0.10, (shift.min(), shift.max())
    # ... while the band it saw is learned
    small = np.exp(tm2.log_time(lN2, lL2, 6) - prior.log_time(lN2, lL2, 6)).mean()
    assert 1.2 < small < 1.45, small
    # no data for the newest engine: the previous engine's posterior (with
    # the newer engines' shifts)
    assert np.allclose(scheduler.current_time_model(models, engine=4).th, tm.th)
    m2 = scheduler.fit_time_models({"2:plain": am.TimeModel.stats(F2, y2)})
    assert np.allclose(scheduler.current_time_model(m2, engine=3).th, tm2.th + dth)
    assert np.allclose(scheduler.current_time_model(m2, engine=2).th, tm2.th)
    # engine 3 data after engine 2's: fitted from engine 2's shifted posterior
    m23 = scheduler.fit_time_models({"2:plain": am.TimeModel.stats(F2, y2), "3:plain": st})
    tm23 = scheduler.current_time_model(m23)
    assert tm23.engine == scheduler.ENGINE == 3 and tm23.n == 50
    assert np.abs(tm23.th[1:4] - tm2.th[1:4]).max() < 1e-3
    assert np.allclose(scheduler.current_time_model({}).th, prior3.th)
    assert np.allclose(scheduler.current_time_model({}, engine=2).th, prior.th)
    print("time model ok")


def write_unit(path, P, sums, squares=(), done=None, partial=False):
    # as msearch: a sum's squares precede its sum record; squares at a sum
    # without a record (a killed run) come last
    recs = []
    sq_recs = [{"type": "square", "n": 6, "P": list(P), "S": S, "s_count": s,
                "p_count": p, "sp_count": sp, "best_score": best, "hash": f"{S}{s}{p}",
                "grid": [[1] * 6] * 6} for S, s, p, sp, best in squares]
    for S, sq in sums:
        recs += [q for q in sq_recs if q["S"] == S]
        recs.append({"type": "sum", "n": 6, "P": list(P), "S": S, "nvecs": 2000,
                     "nvecs_raw": 2000 + S, "labels": 120, "nodes": 10 ** 6, "squares": sq,
                     "time": 0.004, "setup_time": 0.001, "enum_time": 0.001, "truncated": 0,
                     "engine": 2})
    recs += [q for q in sq_recs if q["S"] not in {S for S, _ in sums}]
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
        # a unit killed in the middle of S = 902 (its square written, no sum
        # record) and its rerun from 902: the square counts once
        write_unit(os.path.join(units, "u2.jsonl"), P, [(900, 0), (901, 1)], [(902, 0, 1, 0, 3)])
        inc.update(units, store)
        assert inc.totals["squares"] == 2
        write_unit(os.path.join(units, "u3.jsonl"), P, [(902, 1), (903, 0)], [(902, 0, 1, 0, 3)],
                   done=(902, 903))
        inc.update(units, store)
        assert inc.totals["squares"] == 3 and inc.types == {3: 2, 7: 1}
        # sampled records (research builds) are not ingested
        with open(os.path.join(units, "u4.jsonl"), "w") as f:
            f.write(json.dumps({"type": "sum", "n": 6, "P": list(P), "S": 950, "stride": 8,
                                "nvecs": 9, "nvecs_raw": 9, "labels": 9, "nodes": 1,
                                "squares": 5, "time": 1.0, "setup_time": 0, "truncated": 0})
                    + "\n")
        n0 = inc.totals["sums"]
        inc.update(units, store)
        assert inc.totals["sums"] == n0 and inc.totals["sampled"] == 1
    print("summary ok")


def test_dfirst_summary():
    """v2: msearch --diag-first output in the Summary. A d-first sum has no
    "sum" record and no semi-magic squares, so it stays out of the cells
    (squares and traversal fits) and of the time law; a sum searched in
    full is covered, a part of one (a --d-range unit, a --d-stride sample),
    a plain sum that a --d-range unit skipped and an r1-sampled plain sum
    are cut out of the file's "done" range, also when the records arrive in
    separate incremental reads; a magic square found d-first (twice, once
    per diagonal) is notable once"""
    P = (12, 6, 3, 2, 1, 1)
    key = "12_6_3_2_1_1"
    sums = [(S, S % 3) for S in range(880, 886)]
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        store = scheduler.ProfileStore(state, None, 6)
        # the reference: the plain sums alone
        ref_path = os.path.join(state, "ref.jsonl")
        write_unit(ref_path, P, sums, done=(880, 885))
        ref = scheduler.Summary(state, 6)
        ref.update_file(ref_path, store)
        # the same plain sums with d-first records and a done record over
        # 880..895 (mode dfirst)
        text = write_unit(os.path.join(state, "plain.jsonl"), P, sums)
        grid = [[1] * 6] * 6
        dsq = {"type": "dsquare", "n": 6, "P": list(P), "S": 886, "set_count": 2,
               "s_count": 9, "p_count": 9, "sp_count": 2, "best_score": 14, "magic": 1,
               "partner": 1, "hash": "00000000000000ab", "grid": grid}
        extra = [
            dict(dsq, d=3),
            dict(dsq, d=7),  # the same square on its other diagonal
            {"type": "dsum", "mode": "dfirst", "n": 6, "P": list(P), "S": 886, "nvecs": 6000,
             "nvecs_raw": 6100, "labels": 160, "nodes": 10 ** 7, "pairs": 2, "time": 50.0,
             "cpu": 51.0, "truncated": 0, "complete": 1, "engine": 2},
            {"type": "dsum", "mode": "dfirst", "n": 6, "P": list(P), "S": 887, "nvecs": 6000,
             "nvecs_raw": 6100, "labels": 160, "nodes": 10 ** 6, "pairs": 0, "time": 5.0,
             "cpu": 5.5, "truncated": 0, "complete": 0, "d_stride": 3, "engine": 2},
            # 886 searched again: counted once among the sums, its CPU twice
            {"type": "dsum", "mode": "dfirst", "n": 6, "P": list(P), "S": 886, "nvecs": 6000,
             "nvecs_raw": 6100, "labels": 160, "nodes": 10 ** 7, "pairs": 0, "time": 40.0,
             "cpu": 40.0, "truncated": 0, "complete": 1, "engine": 2},
            {"type": "skip", "mode": "dfirst", "n": 6, "P": list(P), "S": 889, "nvecs": 900,
             "nvecs_raw": 950, "d_lo": 200},
            {"type": "csum", "mode": "sampled", "n": 6, "P": list(P), "S": 891, "r1_stride": 4,
             "squares": 3, "time": 1.0, "cpu": 1.0, "truncated": 0, "engine": 2},
            {"type": "done", "n": 6, "P": list(P), "min_sum": 880, "last_sum": 895,
             "complete": 1, "time": 60.0, "mode": "dfirst", "diag_first_min_n": 5000},
        ]
        full = text + "".join(json.dumps(r) + "\n" for r in extra)
        path = os.path.join(units, "u.jsonl")
        # incremental: the d-first records of 887 / 889 / 891 in one read,
        # the done record in the next
        cut = full.index('"type": "done"')
        cut = full.rindex("\n", 0, cut) + 1
        with open(path, "w") as f:
            f.write(full[:cut])
        inc = scheduler.Summary(state, 6)
        inc.update(units, store)
        assert inc.cover[key] == [[880, 886]], inc.cover
        with open(path, "w") as f:
            f.write(full)
        inc.update(units, store)
        whole = scheduler.Summary(state, 6)
        whole.update(units, store)
        assert summary_state(whole) == summary_state(inc)
        for s in (inc, whole):
            assert s.cover[key] == [[880, 886], [888, 888], [890, 890], [892, 895]], s.cover
            # the fits see the plain sums only
            assert s.time == ref.time and s.tband == ref.tband and s.nbias == ref.nbias
            assert s.perP[key]["cells"] == ref.perP[key]["cells"]
            assert s.perP[key]["cpu"] == ref.perP[key]["cpu"]
            assert s.totals["sums"] == len(sums) and s.totals["squares"] == 0
            assert s.totals["cpu"] == ref.totals["cpu"]
            T = s.totals
            assert (T["dfirst_sums"], T["dfirst_partial"], T["dfirst_magic"]) == (1, 1, 1), T
            assert abs(T["dfirst_cpu"] - 96.5) < 1e-9 and T["sampled"] == 1
            assert T["sampled_cpu"] == 1.0
            assert [q["hash"] for q in s.notable] == ["00000000000000ab"]
            assert s.notable[0]["dfirst"] == 1 and not s.types
        # the holes survive a save and load between the reads
        inc.save()
        again = scheduler.Summary.load(state, 6)
        assert again.files["u.jsonl"]["holes"] == {key: [887, 889, 891]}
        assert summary_state(again) == summary_state(inc)
        # an older summary (without the d-first fields) is rebuilt
        with open(os.path.join(state, "summary_6.json")) as f:
            d = json.load(f)
        d["tag"].pop("summary")
        with open(os.path.join(state, "summary_6.json"), "w") as f:
            json.dump(d, f)
        assert not scheduler.Summary.load(state, 6).files
        # the report says what d-first found
        sch = SimpleNamespace(summary=inc, calib=scheduler.Calibration(), dir=state, n=6,
                              time_models={}, tm=scheduler._am().TimeModel(engine=2))
        import io
        buf = io.StringIO()
        scheduler.report_v2(sch, out=buf)
        assert "d-first (in the d-first time law only): 1 sums searched in full, 1 parts" \
            in buf.getvalue()
        assert "(d-first)" in buf.getvalue()
    # a notable square of a sum searched twice (two units) is notable once;
    # the units never inherit bench's SAMPLE_* variables
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        store = scheduler.ProfileStore(state, None, 6)
        for name in "ab":
            write_unit(os.path.join(units, name + ".jsonl"), P, [(880, 1)],
                       squares=[(880, 3, 3, 1, 7)], done=(880, 880))
        s = scheduler.Summary(state, 6)
        s.update(units, store)
        assert s.totals["squares"] == 2 and len(s.notable) == 1, s.notable
    os.environ["SAMPLE_STRIDE"] = "3"
    try:
        env = scheduler.msearch_env()
        assert "SAMPLE_STRIDE" not in env and env.get("PATH") == os.environ.get("PATH")
    finally:
        del os.environ["SAMPLE_STRIDE"]
    # end to end with msearch: complete d-first sums covered, a d-sampled
    # one and an r1-sampled plain run not; the time law never sees d-first
    msearch = os.path.join(os.path.dirname(HERE), "bin", "msearch")
    P = (10, 4, 3, 2)
    key = "10_4_3_2"
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        store = scheduler.ProfileStore(state, None, 6)
        base = [msearch, "--diag-first", "--diag-first-min-n", "455"]
        # after reduction S = 327 has 451 vectors (plain), 328 460 and 329
        # 480 (d-first), 330 438 (plain)
        subprocess.run(base + ["--min-sum", "327", "--max-sum", "328", "--out",
                               os.path.join(units, "a.jsonl"), "10", "4", "3", "2"],
                       check=True, capture_output=True)
        subprocess.run(base + ["--d-stride", "3", "--min-sum", "329", "--max-sum", "330",
                               "--out", os.path.join(units, "b.jsonl"), "10", "4", "3", "2"],
                       check=True, capture_output=True)
        subprocess.run([msearch, "--r1-stride", "4", "--min-sum", "331", "--max-sum", "332",
                        "--out", os.path.join(units, "c.jsonl"), "10", "4", "3", "2"],
                       check=True, capture_output=True)
        # --pretest-min alone stays a plain run (no "mode" in its done record)
        subprocess.run([msearch, "--pretest-min", "5", "--min-sum", "333", "--max-sum", "333",
                        "--out", os.path.join(units, "d.jsonl"), "10", "4", "3", "2"],
                       check=True, capture_output=True)
        recs = {}
        for name in "abcd":
            with open(os.path.join(units, name + ".jsonl")) as f:
                recs[name] = [json.loads(line) for line in f]
        ds = [r for r in recs["a"] + recs["b"] if r["type"] == "dsum"]
        assert [(r["S"], r["complete"]) for r in ds] == [(328, 1), (329, 0)], ds
        assert all(r["cpu"] >= r["time"] > 0 for r in ds), ds
        assert [r["type"] for r in recs["c"] if r["type"] != "csquare"] == ["csum", "csum", "done"]
        assert recs["c"][-1]["mode"] == "sampled" and recs["c"][-1]["r1_sample"] == 1
        assert [r["type"] for r in recs["d"]] == ["sum", "done"] and "mode" not in recs["d"][-1]
        s = scheduler.Summary(state, 6)
        s.update(units, store)
        assert s.cover[key] == [[327, 328], [330, 330], [333, 333]], s.cover
        assert s.totals["sums"] == 3 and s.totals["dfirst_sums"] == 1
        assert s.totals["dfirst_partial"] == 1
        assert sum(st["n"] for k, st in s.time.items() if k.endswith(":plain")) <= 3
        assert set(s.time) <= {"3:plain", "3:dfirst"}, s.time.keys()
        # v1 agrees on the coverage
        res = scheduler.Results(6)
        for name in "abcd":
            res.add_file(os.path.join(units, name + ".jsonl"))

        class PI:
            def smin(self, P):
                return 327
        assert res.frontier(P, PI()) == 329, res.covered
        assert sorted(res.covered[P]) == [(327, 327), (327, 328), (328, 328), (330, 330),
                                          (330, 330), (333, 333), (333, 333)], res.covered
    print("d-first summary ok")


def _dchunk(P, S, lo, hi, raw=None, stride=1, truncated=0):
    r = {"type": "dchunk", "n": 6, "P": list(P), "S": S, "d_lo": lo, "d_hi": hi,
         "d_stride": stride, "d_offset": 0, "nd": (hi - lo + stride - 1) // stride, "nodes": 10,
         "pairs": 0, "partners": 0, "time": 0.01 * (hi - lo), "truncated": truncated}
    if raw is not None:
        r["nvecs_raw"] = raw
    return r


def _dsum(P, S, lo, hi, raw, cpu=10.0, time_=9.9, index_time=0.1, complete=0, stride=1,
          truncated=0):
    nd = (hi - lo + stride - 1) // stride
    return {"type": "dsum", "mode": "dfirst", "n": 6, "P": list(P), "S": S, "nvecs": raw - 1,
            "nvecs_raw": raw, "labels": 150, "d_lo": lo, "d_hi": hi, "d_stride": stride,
            "d_offset": 0, "nd": nd, "nodes": 10 ** 6, "pairs": 0, "partners": 0, "time": time_,
            "index_time": index_time, "cpu": cpu, "truncated": truncated, "complete": complete,
            "engine": 3}


def _done(P, lo, hi):
    return {"type": "done", "n": 6, "P": list(P), "min_sum": lo, "last_sum": hi, "complete": 1,
            "time": 1.0, "mode": "dfirst", "diag_first_min_n": 0}


def _write(path, recs, tail=""):
    with open(path, "w") as f:
        f.write("".join(json.dumps(r) + "\n" for r in recs) + tail)


def test_dfirst_merge():
    """v2: the parts of a d-first sum (--d-range units: dchunk records with
    d_stride 1) merge into dcov, robust to duplicates, overlaps, killed units
    (chunks without a dsum, a partial last line), sampled parts (ignored)
    and chunks without nvecs_raw; the sum is covered once they cover every d
    (counted once), also across incremental reads and a save / load; the
    d-first time law learns the whole-sum CPU of a part (the cost profile
    along d), and the plain law never sees it"""
    saved = scheduler.DFIRST_TIME_MIN_N
    scheduler.DFIRST_TIME_MIN_N = 1000   # (these sums have 1000 d)
    try:
        _test_dfirst_merge()
    finally:
        scheduler.DFIRST_TIME_MIN_N = saved
    # with the default, a sum of 1000 d teaches the law nothing
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        P = (12, 6, 3, 2, 1, 1)
        _write(os.path.join(units, "a.jsonl"), [_dchunk(P, 890, 0, 512, 1000),
                                                _dsum(P, 890, 0, 512, 1000)])
        s = scheduler.Summary(state, 6)
        s.update(units, scheduler.ProfileStore(state, None, 6))
        assert not s.time and s.dcov
    print("d-first merge ok")


def _test_dfirst_merge():
    import amodel as am
    P = (12, 6, 3, 2, 1, 1)
    key = "12_6_3_2_1_1"
    S, raw = 890, 1000
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        store = scheduler.ProfileStore(state, None, 6)
        # unit a: d 0..512 in two chunks; b: a duplicate of the first chunk;
        # c: an overlap 400..700; d: killed after its chunk 700..800 (no
        # dsum, a partial line); e: a d-stride 3 sample of 800..1000
        _write(os.path.join(units, "a.jsonl"),
               [_dchunk(P, S, 0, 256, raw), _dchunk(P, S, 256, 512, raw),
                _dsum(P, S, 0, 512, raw), _done(P, S, S)])
        _write(os.path.join(units, "b.jsonl"),
               [_dchunk(P, S, 0, 256, raw), _dsum(P, S, 0, 256, raw), _done(P, S, S)])
        _write(os.path.join(units, "c.jsonl"),
               [_dchunk(P, S, 400, 656, raw), _dchunk(P, S, 656, 700, raw),
                _dsum(P, S, 400, 700, raw), _done(P, S, S)])
        _write(os.path.join(units, "d.jsonl"), [_dchunk(P, S, 700, 800, raw)],
               tail=json.dumps(_dchunk(P, S, 800, 900, raw))[:40])
        _write(os.path.join(units, "e.jsonl"),
               [_dchunk(P, S, 800, 1000, raw, stride=3), _dsum(P, S, 800, 1000, raw, stride=3),
                _done(P, S, S)])
        s = scheduler.Summary(state, 6)
        s.update(units, store)
        assert not s.cover.get(key), s.cover
        assert s.dcov[key] == {str(S): {"nd": raw, "iv": [[0, 800]]}}, s.dcov
        nd, iv, frac = s.dfirst_part(key, S)
        assert nd == raw and abs(frac - float(am.dfirst_cost_cum(0.8))) < 1e-12
        assert s.totals["dfirst_sums"] == 0 and s.totals["dfirst_partial"] == 4
        # the d-first time law: one row per part with d_stride 1 (a, b, c),
        # its CPU scaled to the whole sum by the cost profile
        T = s.time["3:dfirst"]
        assert T["n"] == 3 and "3:plain" not in s.time and not s.perP[key]["cells"], s.time.keys()
        ys = []
        for lo, hi in ((0, 512), (0, 256), (400, 700)):
            f = am.dfirst_cost_frac(lo, hi, raw)
            ys.append(math.log(0.1 + 0.1 + 9.8 / f))
        assert abs(T["yy"] - sum(y * y for y in ys)) < 1e-9, (T["yy"], ys)
        # the killed unit's run continues from 800 (a rerun of 800..900
        # whose first record is a duplicate), and covers the sum
        _write(os.path.join(units, "f.jsonl"),
               [_dchunk(P, S, 800, 900, raw), _dchunk(P, S, 900, 1000, raw),
                _dsum(P, S, 800, 1000, raw), _done(P, S, S)])
        s.update(units, store)
        assert s.cover[key] == [[S, S]] and key not in s.dcov, (s.cover, s.dcov)
        assert s.totals["dfirst_sums"] == 1 and s.perP[key]["dfirst_S"] == [S]
        # a complete dsum of the same sum later: still one sum
        _write(os.path.join(units, "g.jsonl"),
               [_dsum(P, S, 0, raw, raw, complete=1), _done(P, S, S)])
        s.update(units, store)
        assert s.totals["dfirst_sums"] == 1 and s.cover[key] == [[S, S]]
        # chunks without nvecs_raw (an older msearch): the number of d comes
        # from a dsum; then the sum completes from chunks alone
        S2 = 891
        _write(os.path.join(units, "h.jsonl"), [_dchunk(P, S2, 0, 500)])
        s.update(units, store)
        assert s.dcov[key][str(S2)] == {"nd": None, "iv": [[0, 500]]}
        _write(os.path.join(units, "i.jsonl"), [_dsum(P, S2, 0, 500, 600)])
        s.update(units, store)
        assert s.dcov[key][str(S2)] == {"nd": 600, "iv": [[0, 500]]}
        _write(os.path.join(units, "j.jsonl"), [_dchunk(P, S2, 500, 600)])
        s.update(units, store)
        assert s.cover[key] == [[S, S2]] and key not in s.dcov
        # the same from scratch, and after a save / load
        whole = scheduler.Summary(state, 6)
        whole.update(units, store)
        assert summary_state(whole) == summary_state(s)
        assert whole.dcov == s.dcov
        s.save()
        again = scheduler.Summary.load(state, 6)
        assert summary_state(again) == summary_state(s) and again.dcov == s.dcov
    # a part read before and after a save keeps its intervals; the plain
    # search of a sum searched in part removes it from dcov
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        store = scheduler.ProfileStore(state, None, 6)
        _write(os.path.join(units, "a.jsonl"), [_dchunk(P, S, 0, 300, raw)])
        s = scheduler.Summary(state, 6)
        s.update(units, store)
        s.save()
        s = scheduler.Summary.load(state, 6)
        _write(os.path.join(units, "b.jsonl"), [_dchunk(P, S, 250, 600, raw)])
        s.update(units, store)
        assert s.dcov[key][str(S)]["iv"] == [[0, 600]]
        write_unit(os.path.join(units, "c.jsonl"), P, [(S, 2)], done=(S, S))
        s.update(units, store)
        assert key not in s.dcov and s.cover[key] == [[S, S]]


def test_dfirst_plan():
    """v2: the mode of each sum (the d-first law plus the calibration
    stream against the plain law), the d-first prior and its online level,
    units of d sized to --unit-time by the cost profile, resumed from the
    first d not searched, the calibration stream on the first unit of a sum
    only, the parts adding up to the sum (E and CPU), and the upper bounds
    of the lazy planner with d-first"""
    import numpy as np
    import amodel as am
    # the d-first prior, and a level learned online (the slope held)
    tmd = scheduler.current_time_model({}, mode="dfirst")
    assert tmd.mode == "dfirst" and tmd.engine == 3
    assert np.allclose(tmd.th, am.DFIRST_TIME_PRIOR["th"]) and tmd.sd == am.DFIRST_TIME_PRIOR["sd"]
    rng = np.random.default_rng(7)
    lNp = np.log(rng.uniform(5000, 30000, 40))
    lL = np.log(rng.uniform(150, 280, 40))
    y = tmd.log_time(lNp, lL, 6) + math.log(0.6) + 0.5 * (lNp - math.log(4000)) \
        + rng.normal(0, 0.2, 40)
    y -= 0.5 * np.mean(lNp - math.log(4000))   # the same mean level, another slope
    st = am.TimeModel.stats(am.time_features(lNp, lL, 6), y)
    models = scheduler.fit_time_models({"3:dfirst": st, "3:plain": st})
    td2 = scheduler.current_time_model(models, mode="dfirst")
    assert td2.mode == "dfirst" and td2.n == 40
    assert abs(td2.th[1] - tmd.th[1]) < 1e-3 and np.allclose(td2.th[2:], tmd.th[2:], atol=1e-3)
    assert abs(td2.th[0] - tmd.th[0] - math.log(0.6)) < 0.1, td2.th
    # the plain law of the same data is fitted from the plain prior
    assert np.allclose(models["3:plain"].th[1:4], am.time_prior(3).th[1:4], atol=1e-3)
    assert scheduler.current_time_model(models, engine=4, mode="dfirst").th[0] == td2.th[0]
    P = (13, 7, 4, 3, 1, 1)
    with tempfile.TemporaryDirectory() as state:
        cands, sc, plan = make_v2([P], state, unit_time=600.0)
        sc.tm = am.time_prior(scheduler.ENGINE)
        tm, tmd = sc.tm, sc.tmd
        # the mode choice, from the two laws
        S = np.arange(1850.0, 2700.0, 10.0)
        sq, m, tp, td, tc, dm, lNp, lL, cell = sc.eval_modes(0, S)
        assert np.allclose(tp, tm.time(lNp, lL, 6), rtol=1e-6)
        assert np.allclose(td, tmd.time(lNp, lL, 6), rtol=1e-6)
        assert np.allclose(tc, np.minimum(tp, scheduler.CALIB_FRAC * td))
        assert (dm == ((td + tc < tp) & (lNp >= math.log(scheduler.DFIRST_MIN_NP)))).all()
        # (measured: plain at N 4.1k, S = 1900; d-first at 15.2k, S = 2200)
        assert not dm[S == 1900][0] and dm[S == 2200][0] and dm[S >= 2200].all()
        t_eff = sc.eval_sums(0, S)[2]
        assert np.allclose(t_eff, np.where(dm, td + tc, tp))
        # grid_density agrees with eval_sums at the grid points
        Sg = sc.prof(0)[0]
        g = sc.grid_density([0])[0]
        g = g[np.isfinite(g)]
        sqg, mg, tg = sc.eval_sums(0, Sg)[:3]
        assert np.allclose(g, np.log(mg / tg), rtol=1e-3, atol=1e-3)
        sc.policy = "off"
        assert not sc.eval_modes(0, S)[5].any() and np.allclose(sc.eval_sums(0, S)[2], tp)
        sc.policy, sc.ldmin = "on", 0.0
        assert sc.eval_modes(0, S)[5].all()
        sc.policy, sc.ldmin = "auto", math.log(scheduler.DFIRST_MIN_NP)
        # a sum of ~11k CPU-s at S = 2400 in units of d (600 s each)
        f0 = int(plan.f0[0])
        plan.set_cover(0, [[f0, 2399]])
        assert plan.frontier[0] == 2400
        sq1, m1, tp1, td1, tc1 = [float(v[0]) for v in sc.eval_modes(0, [2400.0])[:5]]
        units = []
        while plan.frontier[0] == 2400:
            u = plan.unit(0)
            units.append(u)
            plan.advance_unit(u)
            assert len(units) < 100
        n_exp = td1 / 600
        assert abs(len(units) - n_exp) <= 1.5, (len(units), n_exp)
        assert all(u.mode == "dfirst" and u.lo == u.hi == 2400 for u in units)
        nd = units[0].nd
        assert abs(nd - math.exp(sc.eval_modes(0, [2400.0])[6][0])) <= 1
        assert units[0].dlo == 0 and units[0].calib > 1 and all(u.calib == 0 for u in units[1:])
        assert all(a.dhi == b.dlo for a, b in zip(units, units[1:]))
        # (the last unit runs to the end of the sum: nd is only predicted)
        assert units[-1].dhi is None and all(u.dhi is not None for u in units[:-1])
        assert abs(sum(u.frac for u in units) - 1) < 1e-9
        assert abs(sum(u.magic for u in units) - m1) < 1e-9 * m1
        k = units[0].calib
        assert abs(k - tp1 / (scheduler.CALIB_FRAC * td1)) <= 0.5
        assert abs(sum(u.time for u in units) - (td1 + tp1 / k + len(units) * am.UNIT_OVERHEAD)) \
            < 1e-6 * td1
        for u in units[:-2]:
            assert abs(u.frac * td1 - 600) < 0.01 * 600 + td1 / nd, (u.frac * td1, u)
        # the first part has the density of the whole sum with its stream
        # (spread over the parts), the others that of the d loop alone
        assert abs(units[0].score / (m1 / (td1 + tc1)) - 1) < 0.01
        assert all(abs(u.score / (m1 / td1) - 1) < 0.01 for u in units[1:-1]), \
            [u.score * td1 / m1 for u in units]
        assert plan.frontier[0] > 2400 and plan.cover[0] == [[f0, 2400]], plan.cover[0]
        assert 2400 not in plan.dpart.get(0, {})
        # resume from the records of a sum searched in part (nd known)
        plan.set_cover(0, [[f0, 2399]], {2400: (23000, [[0, 5000], [9000, 9500]])}, {2400})
        u = plan.unit(0)
        assert u.dlo == 5000 and 5000 < u.dhi <= 9000 and u.calib == 0 and u.nd == 23000
        assert scheduler.dfirst_args(u) == ["--diag-first", "--diag-first-min-n", "0", "--d-range",
                                            f"5000:{u.dhi}"]
        plan.set_cover(0, [[f0, 2399]], {2400: (23000, [[0, 8990], [9000, 9500]])}, set())
        u = plan.unit(0)
        assert (u.dlo, u.dhi) == (8990, 9000) and u.calib > 0
        plan.advance_unit(u)
        u = plan.unit(0)
        assert u.dlo == 9500 and u.calib == 0
        # (--dfirst off: the sum is searched plain, whole)
        sc.policy = "off"
        u = plan.unit(0)
        assert u.mode == "plain" and u.lo == 2400 and u.dlo is None
        assert scheduler.dfirst_args(u) == []
        sc.policy = "auto"
        # whole d-first sums in one unit (S = 2000: ~240 s each)
        plan.set_cover(0, [[f0, 1999]], {}, set())
        u = plan.unit(0)
        assert u.mode == "dfirst" and u.dlo is None and u.hi > u.lo and u.calib > 0
        assert "--d-range" not in scheduler.dfirst_args(u)
    # the measured sums (ideas.md, "Measurements on the integrated binary"):
    # with the shipped engine-3 laws, the cheaper mode measured
    measured = [((10, 6, 4, 2, 1, 1), 838, False), ((12, 6, 3, 2, 1, 0, 1), 900, False),
                ((12, 6, 3, 2, 1, 1), 988, True), ((12, 6, 3, 2, 1, 1), 1200, True),
                ((11, 6, 4, 3, 2, 1), 2174, True), ((14, 7, 4, 4, 1, 0, 0, 1), 3648, True),
                ((9, 6, 4, 3, 1, 1, 1, 1), 2700, True), ((13, 7, 4, 3, 1, 1), 2650, True)]
    with tempfile.TemporaryDirectory() as state:
        cands, sc, plan = make_v2([P for P, _, _ in measured], state)
        sc.tm = am.time_prior(scheduler.ENGINE)
        for a, (P, S, want) in enumerate(measured):
            assert bool(sc.eval_modes(a, [float(S)])[5][0]) == want, (P, S, want)
    # upper bounds hold with d-first (and with --dfirst on)
    Ps = small_pool(60, seed=11)
    with tempfile.TemporaryDirectory() as state:
        for policy in ("auto", "on"):
            cands, sc, plan = make_v2(Ps, state, unit_time=30.0)
            sc.policy = policy
            if policy == "on":
                sc.ldmin = 0.0
            ub = plan.upper_bounds(np.arange(len(Ps)))
            nd = 0
            for a in range(len(Ps)):
                u = plan.unit(a)
                assert u is None or u.score <= ub[a] * (1 + 1e-9), (Ps[a], u, ub[a])
                nd += u is not None and u.mode == "dfirst"
            assert policy == "auto" or nd > 0
    print("d-first plan ok")


def test_calib_cells():
    """v2: the calibration stream of a d-first sum (csum mode calib and its
    csquare records) enters the cells: squares with expectation E_sq / k
    and the traversals of its squares, weighted by its sampling dispersion;
    a stream repeated, or of a sum searched plain before, is not counted
    again; a plain search after it replaces it; r1-sampled research records
    stay out; the csquares wait for their csum across reads"""
    P = (12, 6, 3, 2, 1, 1)
    key = "12_6_3_2_1_1"
    S = 890
    grid = [[1] * 6] * 6

    def csq(i, best=3):
        return {"type": "csquare", "n": 6, "P": list(P), "weight": 10, "S": S, "s_count": 2 + i,
                "p_count": 1, "sp_count": i % 2, "best_score": best, "hash": f"c{i}",
                "grid": grid}

    def csum(k=10, o=3, se=20.0, mode="calib", S_=S):
        r = {"type": "csum", "mode": mode, "n": 6, "P": list(P), "S": S_, "nvecs": 900,
             "nvecs_raw": 910, "labels": 120, "squares": o, "est_squares": k * o,
             "se_squares": se, "cpu": 0.5, "time": 0.5, "truncated": 0, "engine": 3}
        r["r1_stride"] = k
        return r

    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        store = scheduler.ProfileStore(state, None, 6)
        # the reference: the plain sum S and its 3 squares
        ref = scheduler.Summary(state, 6)
        refp = os.path.join(state, "ref.jsonl")
        write_unit(refp, P, [(S, 3)], [(S, 2, 1, 0, 3), (S, 3, 1, 1, 3), (S, 4, 1, 0, 3)],
                   done=(S, S))
        ref.update_file(refp, store)
        (c, row), = ref.perP[key]["cells"].items()
        E = row[2]
        # the stream: 3 sampled squares at k = 10 (one notable), its csum
        # in the next read
        path = os.path.join(units, "u.jsonl")
        recs = [csq(0), csq(1, best=7), csq(2)]
        _write(path, recs + [_dsum(P, S, 0, 910, 910, complete=1)])
        s = scheduler.Summary(state, 6)
        s.update(units, store)
        assert not s.perP[key]["cells"] and s.totals["calib_sums"] == 0
        _write(path, recs + [_dsum(P, S, 0, 910, 910, complete=1), csum(), _done(P, S, S)])
        s.update(units, store)
        w = min(1.0, scheduler.PHI_SUM / (scheduler.CALIB_PHI + scheduler.PHI_SUM / 10))
        cells = s.perP[key]["cells"]
        assert list(cells) == [c], cells
        got = cells[c]
        assert got[0] == 0 and abs(got[1] - 3 * w) < 1e-9 and abs(got[2] - w * E / 10) < 1e-9
        # the traversals of its squares, with the same expectations per
        # square as a plain search's
        for j in range(3, 9):
            assert abs(got[j] - w * row[j]) < 1e-9 * max(1, row[j]), (j, got, row)
        T = s.totals
        assert (T["calib_sums"], T["calib_squares"], T["calib_cpu"]) == (1, 3, 0.5)
        assert abs(T["calib_est"] - 30) < 1e-9 and abs(T["calib_pred"] - E) < 1e-9
        assert T["squares"] == 0 and T["sampled_cpu"] == 0
        assert [q["hash"] for q in s.notable] == ["c1"] and s.notable[0]["calib"] == 1
        assert s.perP[key]["calS"][str(S)][0] == c
        # the same from scratch
        whole = scheduler.Summary(state, 6)
        whole.update(units, store)
        assert summary_state(whole) == summary_state(s)
        # a second stream of the same sum: not counted again
        _write(os.path.join(units, "v.jsonl"), [csq(0), csum(), _done(P, S, S)])
        s.update(units, store)
        assert s.perP[key]["cells"][c] == got and s.totals["calib_dup"] == 1
        # the plain search of the sum replaces the stream
        write_unit(os.path.join(units, "w.jsonl"), P, [(S, 3)],
                   [(S, 2, 1, 0, 3), (S, 3, 1, 1, 3), (S, 4, 1, 0, 3)], done=(S, S))
        s.update(units, store)
        for j in range(9):
            assert abs(s.perP[key]["cells"][c][j] - row[j]) < 1e-9, (s.perP[key]["cells"], row)
        assert s.totals["calib_sums"] == 0
        # a stream of a sum searched plain before: not counted
        _write(os.path.join(units, "x.jsonl"), [csq(0), csum(), _done(P, S, S)])
        s.update(units, store)
        assert s.totals["calib_dup"] == 2 and abs(s.perP[key]["cells"][c][1] - row[1]) < 1e-9
        # the weight from the stream's dispersion where it sampled >= 5
        # squares: phi = se^2 / (k est)
        S3 = 892
        _write(os.path.join(units, "y.jsonl"),
               [dict(csq(i), S=S3) for i in range(6)]
               + [csum(k=4, o=6, se=24.0, S_=S3), _done(P, S3, S3)])
        s.update(units, store)
        phi = 24.0 ** 2 / (4 * 24)
        w3 = min(1.0, scheduler.PHI_SUM / (phi + scheduler.PHI_SUM / 4))
        row3 = s.perP[key]["calS"][str(S3)]
        assert abs(row3[1] - 6 * w3) < 1e-9 and w3 < 1
        # r1-sampled research output: neither cells nor calibration
        n0 = dict(s.totals)
        _write(os.path.join(units, "z.jsonl"),
               [dict(csq(0), S=S3 + 1), csum(mode="sampled", S_=S3 + 1),
                {"type": "done", "n": 6, "P": list(P), "min_sum": S3 + 1, "last_sum": S3 + 1,
                 "complete": 1, "time": 1.0, "mode": "sampled", "r1_sample": 1}])
        s.update(units, store)
        assert s.totals["calib_sums"] == n0["calib_sums"] and s.totals["sampled_cpu"] == 0.5
        assert str(S3 + 1) not in s.perP[key]["calS"]
        # the class factors learn from the streams: squares at 2x the model
        # in calibration streams alone raise the squares intercept
        import numpy as np
        tab = np.zeros((scheduler.NCELLS, 9))
        tab[c, 1], tab[c, 2] = 2 * 400 * w, 400 * w
        cal = scheduler.Calibration.fit(tab)
        assert cal.beta["sq"][0] > 0.3, cal.beta["sq"]
    print("calibration streams ok")


def test_dfirst_e2e():
    """v2 end to end with msearch: --dfirst on with a low threshold and tiny
    units, so that d-first sums are split into units of d, each run resumed
    from the records of the previous ones and merged into coverage, the
    first with the calibration stream; report, emit and forecast with
    d-first"""
    P = (10, 4, 3, 2)
    key = "10_4_3_2"
    with tempfile.TemporaryDirectory() as state:
        units = os.path.join(state, "units")
        os.makedirs(units)
        # sums up to 269 searched before (a plain unit)
        write_unit(os.path.join(units, "seed.jsonl"), P, [(S, 0) for S in range(170, 270)],
                   done=(170, 269))
        opts = ("--only", "10 4 3 2", "--dfirst", "on", "--dfirst-min-n", "0",
                "--unit-time", "0.004")
        out = run("--state", state, "emit", *opts, "--units", "6")
        assert "--diag-first --diag-first-min-n 0 --d-range 0:" in out, out
        assert "--calib-r1-stride" in out and "d0-" in out, out
        run("--state", state, "run", *opts, "--workers", "1", "--hours", "0.0008")
        launched = [json.loads(l) for l in open(os.path.join(state, "launched_6.jsonl"))]
        assert launched and all(r["mode"] == "dfirst" for r in launched), launched[:3]
        by_S = {}
        for r in launched:
            by_S.setdefault(r["lo"], []).append(r)
        split = {S: rs for S, rs in by_S.items() if len(rs) >= 2}
        assert split, by_S
        s = scheduler.Summary.load(state, 6)
        for S, rs in split.items():
            assert all(r["lo"] == r["hi"] == S and r["dlo"] is not None for r in rs)
            assert rs[0]["dlo"] == 0 and rs[0]["calib"] > 0
            assert all(r["calib"] == 0 for r in rs[1:])
            # each unit starts where the records of the previous ones stop
            assert all(a["dhi"] == b["dlo"] for a, b in zip(rs, rs[1:])), rs
            if rs[-1]["dhi"] is None or rs[-1]["dhi"] == rs[-1]["nd"]:
                assert any(lo <= S <= hi for lo, hi in s.cover[key]), (S, s.cover)
                assert str(S) not in s.dcov.get(key, {})
            # the second unit knows the number of d from the first's records
            assert rs[1]["dhi"] is not None and rs[1]["nd"] == s.dcov.get(key, {}).get(
                str(S), {}).get("nd", rs[1]["nd"])
        T = s.totals
        assert T["dfirst_sums"] >= 1 and T["calib_sums"] + T["calib_dup"] >= 1, T
        # (sums of < DFIRST_TIME_MIN_N vectors do not teach the d-first law)
        assert "3:dfirst" not in s.time, s.time.keys()
        out = run("--state", state, "report", *opts)
        assert "d-first (in the d-first time law only)" in out and "calibration streams" in out
        assert "d-first law" in out, out
        out = run("--state", state, "forecast", *opts, "--hours", "0.01", "--draws", "2")
        assert "with and without d-first" in out and "mode     d-first" in out, out
    print("d-first end to end ok")


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

    def binary(name, check=True):
        assert name == "msearch", name
        return saved_bin(name, check)

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
            # run: only msearch is started
            run_args = SimpleNamespace(**vars(args))
            run_args.hours, run_args.unit_time, run_args.refit_every = 0.0003, 1.0, 50
            run_args.announce_score = 7
            with contextlib.redirect_stdout(io.StringIO()):
                scheduler.SchedulerV2(run_args).run()
            assert os.path.exists(os.path.join(state, "launched_6.jsonl"))
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
    test_dfirst_records()
    test_amodel()
    test_pool()
    test_profile_store()
    test_planner()
    test_fit_poisson_map()
    test_per_p()
    test_time_model()
    test_summary()
    test_dfirst_summary()
    test_dfirst_merge()
    test_dfirst_plan()
    test_calib_cells()
    test_dfirst_e2e()
    test_no_enumerate()
    test_commands_v2()
    test_commands()
    print(f"all ok ({time.time() - t0:.0f} s)")
