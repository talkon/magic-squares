#!/usr/bin/env python3
"""
Smoke test of scheduler.py: fitting on synthetic data, and the plan / emit /
run / report commands end to end on a few small values of P, in a temporary
state directory. Needs numpy and a build (bin/msearch, bin/enumerate).

usage: python3 scripts/test_scheduler.py
"""
import json
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scheduler  # noqa: E402

SMALL = "10 4 3 2;11 4 3 1 1;11 5 3 2"


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
    with tempfile.TemporaryDirectory() as state:
        # use the built-in default model
        model = scheduler.Model.load(os.path.join(state, "none.json"))
        assert model is not None, "no default model"
        out = run("--state", state, "plan", "--workers", "2", "--only", SMALL)
        assert "10 4 3 2" in out, out
        lines = run("--state", state, "emit", "--workers", "2", "--only", SMALL,
                    "--units", "4", "--unit-time", "1").strip().split("\n")
        assert len(lines) == 4, lines
        assert all(l.startswith("--vec-size 6 ") for l in lines), lines
        out = run("--state", state, "run", "--workers", "2", "--only", SMALL,
                  "--unit-time", "2", "--hours", "0.004")
        units = os.listdir(os.path.join(state, "units"))
        assert units, out
        sums = 0
        for name in units:
            with open(os.path.join(state, "units", name)) as f:
                sums += sum(json.loads(l)["type"] == "sum" for l in f)
        assert sums > 0
        out = run("--state", state, "report")
        assert "semi-magic squares" in out, out
        # resuming continues after the sums already searched
        before = scheduler.Results(6)
        for name in units:
            before.add_file(os.path.join(state, "units", name))
        out = run("--state", state, "plan", "--workers", "2", "--only", SMALL)
        for line in out.split("\n")[1:]:
            if line.startswith("10 4 3 2 "):
                lo = int(line.split()[4].split("-")[0])
                done = max(hi for lo_, hi in before.covered[(10, 4, 3, 2)])
                assert lo > done, line
        print(f"commands ok ({len(units)} units, {sums} sums)")


if __name__ == "__main__":
    test_fit_poisson()
    test_fit_model()
    test_coverage()
    test_commands()
