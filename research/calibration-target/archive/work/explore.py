#!/usr/bin/env python3
"""Economics per N' band: E-weighted squares, SP pairs, CPU per mode, SP pairs
per CPU-hour, for (a) the planned sums of the 10 CPU-year plan and (b) every
pool sum (S <= 2 S0) weighted by predicted magic per sum (the existence
study's weighting of where N_magic sits). Writes work/exist_grid.npz."""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sim  # noqa: E402
import cm  # noqa: E402
S = cm.S
am = cm.am
W = os.path.dirname(os.path.abspath(__file__))
EDGES = (3000, 6000, 12000, 24000, 45000)
BN = ("3-6k", "6-12k", "12-24k", "24-45k")


def load():
    args = sim.make_args(os.path.join(sim.CT, "state"))
    sch = S.SchedulerV2(args, quiet=True)
    sim.shipped(sch)
    return sch


def exist_grid(sch):
    """per candidate and grid point in N' 3k-45k: magic per sum x dS (the
    integral weight), lN', cell"""
    sc = sch.scorer
    c = sch.cands
    rows_all = np.arange(len(c))
    du = np.gradient(am.GRID_U)
    out = []
    for i in range(0, len(c), 20000):
        rows = rows_all[i:i + 20000]
        A = np.zeros((len(rows), len(am.FIELDS), len(am.GRID_U)), np.float32)
        V = np.zeros((len(rows), len(am.GRID_U)), bool)
        src = c.src[rows]
        inp = src >= 0
        idx = np.nonzero(inp)[0]
        A[idx] = sch.store.arr[src[idx]]
        V[idx] = sch.store.valid[src[idx]].astype(bool)
        for j in np.nonzero(~inp)[0]:
            A[j], V[j] = sch.store.get(c.P(rows[j]))
        lN, lEs, lPm = A[:, 0], A[:, 1], A[:, 2]
        lNp = lN + am.n_bias(lN)
        Sx = c.S0[rows, None] * (1 + am.GRID_U[None, :])
        x = np.log(Sx / c.smin[rows, None])
        cell = S.cell_index(lNp, c.k[rows, None], c.rbin[rows, None], x)
        lem = (lEs + sc.ln_gsq[cell] + math.log(am.SQ12) * (lNp >= S.LN12K)
               + np.minimum(lPm + sc.ln_m[cell], am.LOG_PM_CAP))
        wgt = np.where(V & (lNp >= math.log(3000)) & (lNp < math.log(45000)),
                       np.exp(lem) * c.S0[rows, None] * du[None, :], 0.0)
        r_, g_ = np.nonzero(wgt > 0)
        out.append((rows[r_], g_, wgt[r_, g_], lNp[r_, g_], cell[r_, g_]))
    a = np.concatenate([o[0] for o in out])
    g = np.concatenate([o[1] for o in out])
    w = np.concatenate([o[2] for o in out])
    l = np.concatenate([o[3] for o in out])
    ce = np.concatenate([o[4] for o in out])
    np.savez_compressed(os.path.join(W, "exist_grid.npz"), a=a, g=g, w=w, lNp=l, cell=ce)
    return a, g, w, l, ce


def summarize(sch, a, Ss, w, label):
    nb = np.searchsorted(np.log(EDGES), np.log(np.maximum(1, 1)), side="right")
    rows = []
    for ai, Si in zip(a, Ss):
        p = cm.predict(sch, int(ai), [float(Si)])
        rows.append({k: float(v[0]) for k, v in p.items() if np.ndim(v)})
    keys = ("sq", "magic", "pairs", "td", "tp", "tp_anch", "tp_ratio", "Np", "qSP", "pmagic",
            "sps", "spp", "pairf", "labels")
    arr = {k: np.array([r[k] for r in rows]) for k in keys}
    b = np.searchsorted(EDGES, arr["Np"], side="right") - 1
    print(f"== {label}")
    for i, name in enumerate(BN):
        s = b == i
        if s.sum() < 3:
            print(f"  {name}: {s.sum()} draws")
            continue
        ww = w[s] / w[s].sum()
        def wm(k):
            return float((ww * arr[k][s]).sum())
        # per CPU-hour of d-first: sum pairs / sum td over draws in proportion to magic
        print(f"  {name}: n={s.sum()} N' {wm('Np'):.0f}  squares {wm('sq'):.4g}  P(magic|sq) "
              f"{wm('pmagic'):.3g}  pairs {wm('pairs'):.3g}  SP+S {wm('sps'):.3g} SP+P "
              f"{wm('spp'):.3g}  td {wm('td'):.4g} s  tp {wm('tp'):.4g}  tp_anch {wm('tp_anch'):.4g}"
              f"  pairs/td-hour {3600 * wm('pairs') / wm('td'):.3g}  median pairs/td-hour "
              f"{np.median(3600 * arr['pairs'][s] / arr['td'][s]):.3g}  squares/tp_anch-hour "
              f"{3600 * wm('sq') / wm('tp_anch'):.3g} pairf {wm('pairf'):.3g} labels {wm('labels'):.0f}")
    return arr, b


def main():
    sch = load()
    rng = np.random.default_rng(7)
    d = np.load(os.path.join(W, "plan_sums_10y_s01.npz"))
    m = d["magic"]
    nb = np.searchsorted(EDGES, np.exp(d["lNp"]), side="right") - 1
    pick = []
    for i in range(4):
        idx = np.nonzero(nb == i)[0]
        if len(idx) == 0:
            continue
        p = m[idx] / m[idx].sum()
        pick.append(rng.choice(idx, size=min(300, len(idx)), p=p))
    pick = np.concatenate(pick)
    summarize(sch, d["a"][pick], d["S"][pick], np.ones(len(pick)), "10 CPU-year plan, draws by magic")
    if os.path.exists(os.path.join(W, "exist_grid.npz")):
        e = np.load(os.path.join(W, "exist_grid.npz"))
        a, g, w, l, ce = e["a"], e["g"], e["w"], e["lNp"], e["cell"]
    else:
        a, g, w, l, ce = exist_grid(sch)
    b = np.searchsorted(np.log(EDGES), l, side="right") - 1
    tot = w.sum()
    print("existence weighting (magic per sum, S <= 2 S0, N' 3-45k): E =", tot,
          {BN[i]: round(float(w[b == i].sum() / tot), 3) for i in range(4)})
    pick = []
    for i in range(4):
        idx = np.nonzero(b == i)[0]
        pick.append(rng.choice(idx, size=200, p=w[idx] / w[idx].sum()))
    pick = np.concatenate(pick)
    Ss = np.round(sch.cands.S0[a[pick]] * (1 + am.GRID_U[g[pick]]))
    summarize(sch, a[pick], Ss, np.ones(len(pick)), "existence weighting, draws by magic per sum")


if __name__ == "__main__":
    main()
