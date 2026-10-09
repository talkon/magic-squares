#!/usr/bin/env python3
"""E at 1 and 10 CPU-years with bands, shipped vs updated calibration.

Band as `forecast` builds it: draws of the class factors (Laplace posterior;
the shipped prior for "shipped") per cell, times lognormal(0, 0.2) for the
pair factor and lognormal(0, 0.2) for SP+SP; the plan is not re-optimised
per draw. "updated" also replaces the SP+SP lognormal by the posterior of
the SP coupling f_rho^2 (prior: the same lognormal, sd 0.2 on ln f_rho^2,
updated with the 7 SP pairs), and multiplies the point by its posterior mean.
Also the existence-style level factor F budget.
"""
import json
import math
import os
import sys

import numpy as np

CT = "/tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/calib-target"
WT = "/home/user/magic-squares/.claude/worktrees/calib-target"
sys.path.insert(0, os.path.join(WT, "scripts"))
sys.path.insert(0, os.path.join(CT, "analysis"))
import scheduler as S  # noqa: E402
import forecast_update as FU  # noqa: E402

YR = 8766.0
ND = 4000
AN = json.load(open(os.path.join(CT, "analysis", "results.json")))


def ecell(variant, mark_years):
    z = np.load(os.path.join(CT, "analysis", f"fc_{variant}.npz"))
    frac = float(z["frac"])
    sel = z["h"] < mark_years * YR
    E = np.bincount(z["cell"][sel], weights=z["magic"][sel], minlength=S.NCELLS) / frac
    Eat = dict((float(a), float(b)) for a, b in z["E_at"])
    return E, Eat


def f_rho_posterior(O, E, prior_sd, n=200000, seed=5, extra=None):
    """posterior of ln f (f = f_rho: O ~ Poisson(E f)), normal prior N(0,
    prior_sd^2), by importance sampling; extra = (mean, sd): a further
    normal likelihood on ln f (the heuristic route); returns draws of ln f"""
    rng = np.random.default_rng(seed)
    x = rng.normal(0, prior_sd, n)
    lw = O * x - E * np.exp(x)
    if extra is not None:
        lw = lw - 0.5 * ((x - extra[0]) / extra[1]) ** 2
    w = np.exp(lw - lw.max())
    w /= w.sum()
    idx = rng.choice(n, size=n, p=w)
    return x[idx]


def band(Ecell, cal, rng, sp_draws=None, nd=ND):
    lg0, _, _, lm0 = cal.rates()
    vals = []
    for lg, lm in cal.draws(rng, nd):
        f = np.exp(lg - lg0 + lm - lm0)
        sp = math.exp(rng.normal(0, 0.2)) if sp_draws is None else float(np.exp(rng.choice(sp_draws)))
        vals.append(float((Ecell * f).sum()) * math.exp(rng.normal(0, 0.2)) * sp)
    return np.array(vals)


def main():
    cal0 = S.Calibration()
    cal1 = FU.fit_updated()[0]
    sp = AN["sp"]
    out = {}
    # SP coupling: f_rho = f_SP / (f_S f_P); O = 7 against E = sum E_SP_cond x f_S f_P
    fS = AN["bands"]["all"]["S_per_square"]["ratio"]
    fP = AN["bands"]["all"]["P_per_square"]["ratio"]
    O = sp["posterior_flat"]["O"]
    E = sp["posterior_flat"]["E"] * fS * fP
    out["f_rho_mle"] = O / E
    post = {}
    for psd in (0.1, 0.145, 0.2, 0.3, 10.0):       # sd on ln f_rho; 10 = flat
        d = f_rho_posterior(O, E, psd)
        post[str(psd)] = {"median_f": float(np.exp(np.median(d))), "sd_ln_f": float(d.std()),
                          "mean_f2": float(np.mean(np.exp(2 * d))),
                          "q05_f": float(np.exp(np.quantile(d, 0.05))), "q95_f": float(np.exp(np.quantile(d, 0.95)))}
    out["f_rho_posterior"] = post
    # main: hierarchical prior N(0, 0.29^2) (kappa calibration 0.145 + extrapolation to the
    # target region 0.25), the direct count and the heuristic route ln f = ln 1.35 +- 0.25
    HEUR = (math.log(1.35), 0.25)
    hier = {}
    for name, psd, ex in (("count only, prior 0.29", 0.29, None), ("count + heuristic, prior 0.29", 0.29, HEUR),
                          ("count + heuristic, prior 0.145", 0.145, HEUR), ("count + heuristic, flat", 10.0, HEUR)):
        d = f_rho_posterior(O, E, psd, extra=ex)
        hier[name] = {"median_f": float(np.exp(np.median(d))), "mean_ln": float(d.mean()), "sd_ln_f": float(d.std()),
                      "median_f2": float(np.exp(2 * np.median(d))), "mean_f2": float(np.mean(np.exp(2 * d))),
                      "q16_f": float(np.exp(np.quantile(d, 0.16))), "q84_f": float(np.exp(np.quantile(d, 0.84))),
                      "q05_f": float(np.exp(np.quantile(d, 0.05))), "q95_f": float(np.exp(np.quantile(d, 0.95)))}
    out["f_rho_hier"] = hier
    rng = np.random.default_rng(11)
    main_ln = f_rho_posterior(O, E, 0.29, extra=HEUR)
    f2_point = float(np.exp(2 * np.median(main_ln)))
    sp_draws = 2 * (main_ln - np.median(main_ln))     # spread of ln f_rho^2 about its point
    for mark in (1, 10):
        E0, Eat0 = ecell("shipped", mark)
        E1, Eat1 = ecell("updated", mark)
        res = {"E_shipped_point": Eat0.get(float(mark)), "E_updated_point": Eat1.get(float(mark))}
        try:
            E2, Eat2 = ecell("time", mark)
            res["E_updated_time_point"] = Eat2.get(float(mark))
        except FileNotFoundError:
            E2 = None
        b0 = band(E0, cal0, rng)
        b1 = band(E1, cal1, rng)
        b1s = band(E1, cal1, rng, sp_draws=sp_draws) * f2_point
        res["E_updated_sp_point"] = res["E_updated_point"] * f2_point
        res["f_rho2_point"] = f2_point
        for name, b in (("shipped", b0), ("updated", b1), ("updated_with_sp_posterior", b1s)):
            q = np.percentile(b, [5, 16, 50, 84, 95])
            res[name] = {"q05": q[0], "q16": q[1], "median": q[2], "q84": q[3], "q95": q[4], "mean": float(b.mean()),
                         "sd_ln": float(np.log(b).std()), "ratio_90": float(q[4] / q[0])}
        # share of E by N' band and the per-band factors
        nb = np.arange(S.NCELLS) // 36
        names = ["<1k", "1-3k", "3-6k", "6-12k", "12-24k", ">=24k"]
        res["E_by_band"] = {names[i]: {"shipped": float(E0[nb == i].sum()), "updated": float(E1[nb == i].sum()),
                                       "time": (float(E2[nb == i].sum()) if E2 is not None else None)}
                            for i in range(6)}
        kb = (np.arange(S.NCELLS) // 9) % 4
        res["E_by_k"] = {["k<=5", "k=6", "k=7", "k>=8"][i]: {"shipped": float(E0[kb == i].sum()),
                                                            "updated": float(E1[kb == i].sum())} for i in range(4)}
        xb = np.arange(S.NCELLS) % 3
        res["E_by_x"] = {["x<0.1", "0.1-0.25", ">=0.25"][i]: {"shipped": float(E0[xb == i].sum()),
                                                             "updated": float(E1[xb == i].sum())} for i in range(3)}
        rb = (np.arange(S.NCELLS) // 3) % 3
        res["E_by_ratio"] = {["sorted", "r<=1.1", "r>1.1"][i]: {"shipped": float(E0[rb == i].sum()),
                                                               "updated": float(E1[rb == i].sum())} for i in range(3)}
        out[f"{mark}y"] = res
    with open(os.path.join(CT, "analysis", "bands.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
