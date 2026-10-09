#!/usr/bin/env python3
"""
Smoke test of calibrate.py and square_heuristic.py on a few known squares:
exact counts against the msearch records (square_diag_stats), the exact
moments behind the heuristic, the mod-3 congruences of partner pairs, the
heuristic's predictions (regression values), the sub-events (exact counts by
brute force; their full sets reproduce the rungs), the conditional mode's
exclusive classes, the exact entries of the partner-graph table used for
best-pair probabilities, and calibrate.py end to end on a temporary msearch
output directory with a small first-search summary (searched-range split,
pooled top rung). Needs numpy and bin/enumerate (for S_min).

usage: python3 scripts/test_calibrate.py
"""
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import calibrate as C  # noqa: E402
import square_heuristic as SH  # noqa: E402

# msearch records (grids with rows/cols reordered as msearch writes them)
SQUARES = [
    # SP+P (best 10): the square shown in README.md
    {"type": "square", "n": 6, "P": [16, 5, 4, 2], "S": 849, "s_count": 1, "p_count": 4, "sp_count": 1,
     "best_score": 10, "hash": "cbff3a5d48130c53",
     "grid": [[324, 120, 7, 128, 200, 70], [25, 64, 216, 420, 40, 84], [300, 147, 50, 48, 16, 288],
              [112, 8, 180, 175, 54, 320], [32, 360, 140, 18, 224, 75], [56, 150, 256, 60, 315, 12]]},
    # SP+S (best 9)
    {"type": "square", "n": 6, "P": [12, 6, 3, 2, 1, 0, 1], "S": 836, "s_count": 2, "p_count": 2,
     "sp_count": 1, "best_score": 9, "hash": "930eb5533faf48a3",
     "grid": [[250, 102, 54, 192, 154, 84], [99, 112, 150, 280, 51, 144], [96, 90, 66, 105, 224, 255],
              [56, 252, 160, 68, 135, 165], [119, 220, 168, 81, 200, 48], [216, 60, 238, 110, 72, 140]]},
    # P+P (best 6)
    {"type": "square", "n": 6, "P": [12, 6, 3, 2, 1], "S": 629, "s_count": 2, "p_count": 2, "sp_count": 0,
     "best_score": 6, "hash": "18123acc1c1481d3",
     "grid": [[28, 250, 99, 24, 144, 84], [44, 64, 98, 225, 18, 180], [72, 21, 120, 210, 30, 176],
              [315, 36, 60, 66, 112, 40], [90, 132, 240, 48, 105, 14], [80, 126, 12, 56, 220, 135]]},
    # no S- or P-traversal (best 0)
    {"type": "square", "n": 6, "P": [10, 4, 2, 2, 1, 1], "S": 392, "s_count": 0, "p_count": 0,
     "sp_count": 0, "best_score": 0, "hash": "b4eaaed9d2bffaa3",
     "grid": [[147, 32, 120, 22, 26, 45], [104, 30, 88, 105, 9, 56], [55, 126, 8, 42, 96, 65],
              [20, 21, 78, 160, 77, 36], [48, 143, 63, 24, 100, 14], [18, 40, 35, 39, 84, 176]]},
]

# square_heuristic.predict (moment mode) on these squares, October 2026: a
# regression check (update when the model changes on purpose)
REFERENCE = [
    {"S": 0.9031204664, "P": 0.5107866939, "SP": 0.001285049658, "S+S": 0.008563639261,
     "P+P": 0.00326507153, "SP+SP": 2.900848135e-08},
    {"S": 1.547943272, "P": 0.1805822749, "SP": 0.001470994018, "S+S": 0.02469714342,
     "P+P": 0.000570525219, "SP+SP": 0.0},
    {"S": 1.374122344, "P": 0.3817565452, "SP": 0.001696665376, "S+S": 0.01980801679,
     "P+P": 0.001549858658, "SP+SP": 2.172653213e-08},
    {"S": 2.295481008, "P": 0.184635199, "SP": 0.001329486723, "S+S": 0.05100874476,
     "P+P": 0.0003474114202, "SP+SP": 0.0},
]


def test_observe():
    """both exact counters reproduce the records (square_diag_stats)"""
    recs = [dict(r, P=tuple(r["P"])) for r in SQUARES]
    ev, best, bad = C.observe(recs)
    assert not bad, bad
    for i, r in enumerate(SQUARES):
        o = SH.observe(r["grid"], r["S"], r["P"])
        assert (o["S"], o["P"], o["SP"], o["best"]) == (r["s_count"], r["p_count"], r["sp_count"],
                                                        r["best_score"]), (r["hash"], o)
        assert (ev["S"][i], ev["P"][i], ev["SP"][i], best[i]) == (r["s_count"], r["p_count"], r["sp_count"],
                                                                  r["best_score"])
        assert all(ev[k][i] == o[k] for k in SH.PAIR_KEYS)
        assert sum(o[k] for k in SH.PAIR_KEYS) == 5400
        # a grid with two entries swapped is no longer semi-magic
        g = [row[:] for row in r["grid"]]
        g[0][0], g[1][0] = g[1][0], g[0][0]     # rows 0 and 1 no longer sum to S
        assert C.observe([dict(r, P=tuple(r["P"]), grid=g)])[2], "broken grid not detected"
    print("observe ok")


def test_moments():
    """E[Z] = (S, v(P)) and Hoeffding's covariance, exactly, over the 720
    traversals; the covariance of partner diagonals; the pair congruences
    (the mod-3 span of Z1 + Z2 has dimension <= 4)"""
    for r in SQUARES:
        m = SH.moments(r["grid"], r["S"], r["P"])
        W = SH.traversal_values(m["Z"]).astype(float)
        assert np.allclose(W.mean(0), m["t"], rtol=0, atol=1e-9)
        X = W - m["t"]
        assert np.allclose(X.T @ X / len(W), m["Sigma"], rtol=1e-12, atol=1e-9)
        for ti in range(SH.NTAU):
            X2 = X[SH.PIDX[:, ti]]
            assert np.allclose(X.T @ X2 / len(W), m["Ctau"][ti], rtol=1e-12, atol=1e-9)
            V = (W + W[SH.PIDX[:, ti]]).astype(np.int64)
            assert rank_mod3(V - V[0]) <= 4
    print("moments ok")


def rank_mod3(A):
    A = np.asarray(A, dtype=np.int64) % 3
    rank, col = 0, 0
    rows, cols = A.shape
    while rank < rows and col < cols:
        piv = np.nonzero(A[rank:, col])[0]
        if not len(piv):
            col += 1
            continue
        p = rank + piv[0]
        A[[rank, p]] = A[[p, rank]]
        A[rank] = (A[rank] * A[rank, col]) % 3      # 1/x = x mod 3
        A[rank + 1:] = (A[rank + 1:] - np.outer(A[rank + 1:, col], A[rank])) % 3
        rank += 1
        col += 1
    return rank


def test_heuristic():
    for r, ref in zip(SQUARES, REFERENCE):
        e = SH.predict(r["grid"], r["S"], r["P"])
        assert abs(sum(e[k] for k in SH.PAIR_KEYS) - 5400) < 1e-6
        assert abs(sum(e[f"best={k}"] for k in SH.SCORES) - 1) < 1e-9
        assert abs(e["S_only"] + e["P_only"] + e["SP"] + e["none"] - 720) < 1e-6
        for k, v in ref.items():
            assert abs(e[k] - v) <= 1e-6 * abs(v) + 1e-15, (r["hash"], k, e[k], v)
        c = SH.predict(r["grid"], r["S"], r["P"], mode="conditional")
        assert c["S"] == r["s_count"] and c["P"] == r["p_count"]
    print("heuristic ok")


def test_sub_events():
    """the sub-events: deterministic balanced windows; exact counts by brute
    force; the full sets reproduce the rungs (exact and predicted)"""
    assert SH.sub_sets(5)[:4] == [(1,), (2,), (3,), (4,)] and (1, 2, 3, 4) in SH.sub_sets(5)
    for k in (4, 5, 7):
        wins = [w for w in SH.sub_sets(k + 1) if 0 not in w and len(w) == 2]
        assert len(wins) == k and all(sum(i in w for w in wins) == 2 for i in range(1, k + 1))
    for r in SQUARES:
        comp = SH.components(r["grid"], r["S"], r["P"], sub_events=True)
        e = SH.combine(comp)
        o = SH.observe(r["grid"], r["S"], r["P"])
        k1 = comp["k1"]
        m = SH.moments(r["grid"], r["S"], r["P"])
        W = SH.traversal_values(m["Z"])
        t = np.round(m["t"]).astype(np.int64)
        full = {(0,): ("S", ["S+S", "SP+S", "SP+SP"]), tuple(range(1, k1)): ("P", ["P+P", "SP+P", "SP+SP"]),
                tuple(range(k1)): ("SP", ["SP+SP"])}
        seen = 0
        for row in comp["sub"]:
            hit = np.all(W[:, row["I"]] == t[row["I"]], axis=1)
            assert row["n1"] == hit.sum()
            assert row["n2"] == sum(hit[i] and hit[j] for i in range(SH.NT) for j in SH.PIDX[i] if i < j)
            assert len(row["e1"]) == len(row["e2"]) == len(SH.SUB_VARIANTS)
            # congruence classes of the 15 taus: a forbidden tau has no pair (exact)
            assert sum(row["ntc"]) == SH.NTAU and sum(row["n2c"]) == row["n2"] and row["n2c"][0] == 0
            assert row["e2c"][0] == 0.0 and abs(sum(row["e2c"]) - row["e2"][2]) <= 1e-9 * row["e2"][2] + 1e-15
            if tuple(row["I"]) in full:
                one, two = full[tuple(row["I"])]
                assert row["n1"] == o[one] and row["n2"] == sum(o[x] for x in two)
                assert abs(row["e1"][2] - e[one]) <= 1e-9 * e[one]
                assert abs(row["e2"][2] - sum(e[x] for x in two)) <= 1e-9 * sum(e[x] for x in two) + 1e-15
                seen += 1
        assert seen == 3
    # conditional mode: an SP diagonal's S_only / P_only partners are exclusive
    nS, nP, rho = 9, 4, 2.5
    sp = nS * nP / SH.NT * rho
    trav = {1: nS - sp, 2: nP - sp, 3: sp, 0: SH.NT - nS - nP + sp}
    incl = rho * nS * (nS - 1) * nP / (SH.NT * (SH.NT - 1) * SH.NT)
    assert abs(SH.cond_base(3, 1, trav, nS, nP, rho) + SH.cond_base(3, 3, trav, nS, nP, rho) - incl) < 1e-15
    # the pair residual r(d): measured where there are events, log-linear beyond
    summ = {"rows": [{"kind": "Y", "d": d, "O2": o, "o_nk": r} for d, o, r in
                     ((1, 10 ** 6, 0.99), (2, 10 ** 5, 0.97), (3, 10 ** 4, 0.94), (4, 800, 0.88), (5, 30, 0.80),
                      (6, 0, 0.0))] + [{"kind": "XY", "d": 1, "O2": 100, "o_nk": 1.2}]}
    rk, slope, _ = C.pair_residual_by_k(summ)
    assert rk[4] == (0.88, True) and not rk[6][1] and rk[7][0] < rk[6][0] < 0.80 and slope < 0
    print("sub-events ok")


def test_partner_graph(tmp):
    """exact entries of the W table and the i.i.d. / null event formulas"""
    C.W_TAB, se = C.build_w_table(tmp, n_small=20000, n_large=5000)
    frac = 15 / 719
    assert abs(C.W_TAB[0, 2] - frac) < 1e-12 and abs(C.W_TAB[1, 1] - frac) < 1e-12
    assert np.all(np.diff(C.W_TAB[:, 1:], axis=0) >= -1e-9)          # increasing in x
    ge = C.best_ge_from_counts(np.array([1, 0]), np.array([0, 0]), np.array([1, 2]))
    assert abs(ge[C.LEVELS.index(9), 0] - frac) < 1e-12     # one S, one SP: partners w.p. 15/719
    assert abs(ge[C.LEVELS.index(14), 1] - frac) < 1e-12    # two SP traversals
    nv = C.null_events(np.array([3.0]), np.array([2.0]), np.array([1.0]))
    assert abs(sum(nv[k][0] for k in C.PAIR_KEYS) - 5400) < 1e-9
    assert abs(nv["S+S"][0] - frac * 3) < 1e-12 and abs(nv["SP+P"][0] - frac * 2) < 1e-12
    q = np.array([0.002, 0.0005]), np.array([0.0007, 0.0002]), np.array([2e-5, 1e-6])
    ev = C.iid_events(*q)
    assert np.allclose(sum(ev[k] for k in C.PAIR_KEYS), 5400)
    assert np.allclose(ev["best>=7"], 1 - (1 - q[2]) ** 720, rtol=1e-9)
    ge = np.array([ev[f"best>={k}"] for k in C.LEVELS])
    assert np.all(np.diff(ge, axis=0) <= 1e-12)
    print("partner graph ok")


def test_cli(tmp):
    units = os.path.join(tmp, "units")
    os.makedirs(units)
    with open(os.path.join(units, "test.jsonl"), "w") as f:
        for r in SQUARES:
            f.write(json.dumps({"type": "sum", "n": 6, "P": r["P"], "S": r["S"], "nvecs_raw": 1500,
                                "squares": 1}) + "\n")
            f.write(json.dumps(r) + "\n")
        f.write(json.dumps(SQUARES[0]) + "\n")                     # a duplicate
    # a first-search summary: SQUARES[0]'s P searched up to S = 900 (>= its S),
    # SQUARES[2]'s P only up to 600 (< its S); the other P not searched
    stats = os.path.join(tmp, "stats_short.txt")
    with open(stats, "w") as f:
        f.write("Statistics for each P, sorted by factorization of P:\n"
                "  P   P_val num_vecs max max_S count time sols #S #P #SP hSS hSP hPP hM best\n")
        for r, maxs in ((SQUARES[0], 900), (SQUARES[2], 600)):
            f.write(f"  {' '.join(map(str, r['P']))}   {C.sch.p_value(tuple(r['P']))} 1000 300 {maxs} 100 1.0 "
                    f"20 30 8 1 2 6 1 0 7\n")
    out = os.path.join(tmp, "out")
    cmd = [sys.executable, os.path.join(HERE, "calibrate.py"), "--state", tmp, "--out", out, "--legacy",
           stats, "--boot", "20", "--jobs", "1", "--w-samples", "20000", "--heuristic-cond", "--eb"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout, res.stderr)
        raise SystemExit("calibrate.py failed")
    with open(os.path.join(out, "results.json")) as f:
        d = json.load(f)
    assert d["data"]["squares"] == len(SQUARES) and d["data"]["duplicates"] == 1, d["data"]
    assert d["data"]["mismatch_vs_records"] == {"s_count": 0, "p_count": 0, "sp_count": 0, "best_score": 0}
    ov = d["overall"]
    assert ov["S"]["obs"] == sum(r["s_count"] for r in SQUARES)
    assert ov["SP"]["obs"] == sum(r["sp_count"] for r in SQUARES)
    assert ov["best>=9"]["obs"] == 2 and ov["best>=10"]["obs"] == 1
    for p in ("regression", "null", "heuristic", "heuristic-cond", "perP", "eb", "heur-mixed", "heur-marg",
              "heur-gauss", "heur-nolat"):
        assert p in ov["SP+SP"]["pred"], p
    assert abs(ov["SP"]["pred"]["heuristic"]["E"] - sum(r["SP"] for r in REFERENCE)) < 1e-9
    # first search: SQUARES[0] inside its searched range, the rest new
    assert d["bin_sizes"]["first search"] == {"searched": 1, "new": 3}, d["bin_sizes"]["first search"]
    assert d["legacy"]["squares_in_searched_range"] == 1
    tp = d["legacy"]["top_rung_pooled"]
    assert tp["obs_legacy"] == 10 and tp["obs_new"] == 1      # the SP+S square is new, the SP+P one searched
    # sub-events: the full single-diagonal sets are the S rung (X alone) etc.
    se = {(r["kind"], r["d"]): r for r in d["sub_events"]["all"]["rows"]}
    assert se[("XY", 1)]["O1"] == ov["S"]["obs"]
    assert abs(se[("XY", 1)]["E1_pair"] - ov["S"]["pred"]["heuristic"]["E"]) < 1e-6
    assert set(d["bin_sizes"]["P source"]) == {"P of " + s for s in d["data"]["per_source"]}
    with open(os.path.join(out, "ladder.md")) as f:
        assert "| SP traversals |" in f.read()
    # the second run uses the cache
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0 and "squares to compute" not in res.stderr, res.stderr
    print("calibrate.py end to end ok")


def test_cli_weighted(tmp):
    """csquare records (a calibration stream) count with their stride as weight;
    dsquare records (d-first pairs) stay out of the ladder and are listed"""
    units = os.path.join(tmp, "units")
    os.makedirs(units)
    with open(os.path.join(units, "plain.jsonl"), "w") as f:
        for r in SQUARES[1:]:
            f.write(json.dumps({"type": "sum", "n": 6, "P": r["P"], "S": r["S"], "nvecs_raw": 1500,
                                "squares": 1}) + "\n")
            f.write(json.dumps(r) + "\n")
    r0 = SQUARES[0]

    def sp_dvec(r):
        """the numbers of the square's (only) SP traversal"""
        G = np.array(r["grid"]).reshape(-1)
        Pint = 1
        for p_, e in zip(C.PRIMES, r["P"]):
            Pint *= int(p_) ** int(e)
        for cells in C.CELLIDX:
            v = [int(x) for x in G[cells]]
            if sum(v) == r["S"] and np.prod(np.array(v, dtype=object)) == Pint:
                return v
        raise AssertionError("no SP traversal")

    r1 = SQUARES[1]
    non_sp = [int(x) for x in np.array(r1["grid"]).reshape(-1)[C.CELLIDX[0]]]
    with open(os.path.join(units, "dfirst.jsonl"), "w") as f:
        # (a) the SP+S square, at its SP diagonal: own level SP+S
        f.write(json.dumps(dict(r1, type="dsquare", d=7, dvec=sp_dvec(r1), set_count=1, magic=0, partner=0)) + "\n")
        # (b) the same square at a d that is not an SP diagonal of it: judged on its own d
        f.write(json.dumps(dict(r1, type="dsquare", d=9, dvec=non_sp, set_count=1, magic=0, partner=0)) + "\n")
        f.write(json.dumps({"type": "dsum", "mode": "dfirst", "n": 6, "P": r1["P"], "S": r1["S"],
                            "nvecs_raw": 1500, "d_stride": 4, "pairs": 2}) + "\n")
        # (c, d) one "magic" square (flagged) found at both of its SP diagonals: one square, not two;
        # its run was killed before the dsum, so the d stride (2) comes from the dchunk record
        for d_ in (3, 11):
            f.write(json.dumps(dict(r0, type="dsquare", d=d_, dvec=sp_dvec(r0), set_count=1, magic=1, partner=1,
                                    hash="ffffffffffffffff")) + "\n")
        f.write(json.dumps({"type": "dchunk", "n": 6, "P": r0["P"], "S": r0["S"], "d_lo": 0, "d_hi": 256,
                            "d_stride": 2, "d_offset": 1, "nvecs_raw": 1600}) + "\n")
        f.write(json.dumps(dict(r0, type="csquare", weight=5)) + "\n")
        f.write(json.dumps({"type": "csum", "mode": "calib", "n": 6, "P": r0["P"], "S": r0["S"],
                            "nvecs_raw": 1600, "r1_stride": 5, "squares": 1, "est_squares": 5}) + "\n")
    out = os.path.join(tmp, "out")
    base = [sys.executable, os.path.join(HERE, "calibrate.py"), "--state", tmp, "--legacy", "none",
            "--boot", "20", "--jobs", "1", "--w-samples", "20000", "--no-heuristic"]
    for extra, w0 in (([], 5), (["--unweighted"], 1)):
        res = subprocess.run(base + ["--out", out] + extra, capture_output=True, text=True)
        if res.returncode != 0:
            print(res.stdout, res.stderr)
            raise SystemExit("calibrate.py failed")
        with open(os.path.join(out, "results.json")) as f:
            d = json.load(f)
        assert d["data"]["squares"] == len(SQUARES) and d["data"]["csquares"] == 1, d["data"]
        ov = d["overall"]
        assert ov["S"]["obs"] == w0 * r0["s_count"] + sum(r["s_count"] for r in SQUARES[1:])
        assert ov["P"]["obs"] == w0 * r0["p_count"] + sum(r["p_count"] for r in SQUARES[1:])
        assert ov["best>=10"]["obs"] == w0           # the SP+P square is the csquare
        E = ov["S"]["pred"]["null"]["E"]
        assert abs(E - ov["S"]["obs"]) < 1e-9        # the null reproduces the (weighted) S rung
        assert ("event_scale" in ov["S"]["pred"]["regression"]) == (w0 != 1)
        ds = d["dsquares"]
        assert ds["pairs"] == 4 and ds["est_pairs"] == 4 + 4 + 2 + 2 and ds["squares"] == 2, ds
        assert ds["by_best_score"] == {"SP+S": 2, "SP+P": 2}, ds          # whole-square best pair
        assert ds["d_not_sp"] == 1 and ds["by_own_d"] == {"SP": 0, "SP+S": 1, "SP+P": 2, "SP+SP": 0}, ds
        nS, nP = ds["partner_pairs"]["SP+S"], ds["partner_pairs"]["SP+P"]
        assert nS >= 1 and nP >= 2 and nP % 2 == 0, ds
        assert ds["est_partner_pairs"]["SP+S"] == 4 * nS and ds["est_partner_pairs"]["SP+P"] == 2 * nP, ds
        assert ds["magic_squares"] == 1 and ds["est_magic_squares"] == 2, ds  # deduped by hash; (2 + 2) / 2
        with open(os.path.join(out, "per_square.jsonl")) as f:
            ws = {json.loads(line)["hash"]: json.loads(line)["weight"] for line in f}
        assert ws[r0["hash"]] == w0 and all(ws[r["hash"]] == 1 for r in SQUARES[1:])
    print("calibrate.py csquare weights and dsquare records ok")


def test_cli_star(tmp):
    """dsquare records of msearch's star cover (--dfirst-star K): a pair on a
    star d (its dvec has x*) weighs K in the estimates, one of a K = -1 run
    (no star d searched) none"""
    units = os.path.join(tmp, "units")
    os.makedirs(units)
    with open(os.path.join(units, "plain.jsonl"), "w") as f:
        for r in SQUARES[1:]:
            f.write(json.dumps({"type": "sum", "n": 6, "P": r["P"], "S": r["S"], "nvecs_raw": 1500,
                                "squares": 1}) + "\n")
            f.write(json.dumps(r) + "\n")
    r0, r1 = SQUARES[0], SQUARES[1]
    G = np.array(r1["grid"]).reshape(-1)
    Pint = 1
    for p_, e in zip(C.PRIMES, r1["P"]):
        Pint *= int(p_) ** int(e)
    dv = next([int(x) for x in G[c]] for c in C.CELLIDX
              if sum(int(x) for x in G[c]) == r1["S"] and np.prod(np.array(G[c], dtype=object)) == Pint)
    with open(os.path.join(units, "star.jsonl"), "w") as f:
        # x* on the pair's d, K = 4, stride 1: weight 4; then the dsum
        f.write(json.dumps(dict(r1, type="dsquare", d=7, dvec=dv, set_count=1, magic=0, partner=0)) + "\n")
        f.write(json.dumps({"type": "dsum", "mode": "dfirst", "n": 6, "P": r1["P"], "S": r1["S"],
                            "nvecs_raw": 1500, "d_stride": 1, "pairs": 1, "star_x": dv[2],
                            "star_k": 4, "star_only": 0}) + "\n")
    with open(os.path.join(units, "nostar.jsonl"), "w") as f:
        # K = -1 (killed before its dsum: the star from the dchunk): no estimate
        f.write(json.dumps(dict(r0, type="dsquare", d=3, dvec=dv, set_count=1, magic=1, partner=1,
                                hash="ffffffffffffffff")) + "\n")
        f.write(json.dumps({"type": "dchunk", "n": 6, "P": r0["P"], "S": r0["S"], "d_lo": 0, "d_hi": 256,
                            "d_stride": 1, "nvecs_raw": 1600, "star_x": 5, "star_k": -1,
                            "star_only": 0}) + "\n")
    out = os.path.join(tmp, "out")
    res = subprocess.run([sys.executable, os.path.join(HERE, "calibrate.py"), "--state", tmp, "--legacy",
                          "none", "--boot", "20", "--jobs", "1", "--w-samples", "20000", "--no-heuristic",
                          "--out", out], capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout, res.stderr)
        raise SystemExit("calibrate.py failed")
    with open(os.path.join(out, "results.json")) as f:
        ds = json.load(f)["dsquares"]
    assert ds["pairs"] == 2 and ds["est_pairs"] == 4, ds
    assert ds["star_records"] == 1 and ds["records_not_estimated"] == 1, ds
    assert ds["magic_squares"] == 1 and ds["est_magic_squares"] == 0, ds
    print("calibrate.py star-cover dsquare weights ok")


if __name__ == "__main__":
    test_observe()
    test_moments()
    test_heuristic()
    test_sub_events()
    with tempfile.TemporaryDirectory() as tmp:
        test_partner_graph(tmp)
    with tempfile.TemporaryDirectory() as tmp:
        test_cli(tmp)
    with tempfile.TemporaryDirectory() as tmp:
        test_cli_weighted(tmp)
    with tempfile.TemporaryDirectory() as tmp:
        test_cli_star(tmp)
    print("all ok")
