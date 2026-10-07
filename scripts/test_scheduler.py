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
                assert lo > before.done_to.get((10, 4, 3, 2), 0), line
        print(f"commands ok ({len(units)} units, {sums} sums)")


if __name__ == "__main__":
    test_fit_poisson()
    test_commands()
