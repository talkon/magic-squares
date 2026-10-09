#!/usr/bin/env python3
"""E at 1 and 10 CPU-years with bands, shipped vs updated calibration, the
SP-coupling posterior, and the existence study's level factor F, before and
after, on the same components.

SP coupling f_rho = f_SP / (f_S f_P): the SP rate per square beyond what
the S and P rates (already refit by the GLMs) imply. Its likelihood is the
count of (square, SP diagonal) pairs: O = 7 against E = E_SP(per square,
3-12k) x f_S x f_P. P(magic | square) scales with f_rho^2. The main
posterior (pre-specified: the count alone) uses the prior N(0, 0.29^2) on
ln f_rho: the kappa calibration (0.145, existence.md) and 0.25 for carrying
it from near S_min to the target region (an assumption). Variants are
sensitivity analyses. The "+ heuristic route" variant of the first analysis
multiplied the count by a likelihood built from the j = 4-5 sub-events,
which contain the same SP traversals (32 of 35 events at j = 5, 40 of 109
at j = 4); it double-counts and is listed only to show the size of that
effect. subevents_nonsp.py gives the sub-event coupling without them.

Band as `forecast` builds it: draws of the class factors (Laplace posterior;
the shipped prior for "shipped") per cell, times lognormal(0, 0.2) for the
pair factor and, for the shipped band, lognormal(0, 0.2) for SP+SP. The
updated band replaces that SP term by the f_rho^2 posterior and adds a
selection term lognormal(0, SEL_SD) in place of the old point discount
x0.8 (analyze.py "selection": no winner's curse detectable within band;
top density quartile at 3-6k 0.96 [0.82, 1.11]). The plan is not
re-optimised per draw.

usage: CALIB_TARGET_DIR=... python3 bands.py   (after analyze.py and
forecast_update.py shipped / updated / time; writes analysis/bands.json)
"""
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ctpaths  # noqa: E402

CT, OUT = ctpaths.CT, ctpaths.OUT
sys.path.insert(0, ctpaths.SCRIPTS)
import scheduler as S  # noqa: E402
import forecast_update as FU  # noqa: E402

YR = 8766.0
ND = 4000
PRIOR_SD = 0.29            # ln f_rho: kappa calibration 0.145 (+) extrapolation 0.25
SEL_SD = 0.12              # selection term (ln sd) in the updated band
AN = json.load(open(os.path.join(OUT, "results.json")))


def ecell(variant, mark_years):
    z = np.load(os.path.join(OUT, f"fc_{variant}.npz"))
    frac = float(z["frac"])
    sel = z["h"] < mark_years * YR
    E = np.bincount(z["cell"][sel], weights=z["magic"][sel], minlength=S.NCELLS) / frac
    Eat = dict((float(a), float(b)) for a, b in z["E_at"])
    return E, Eat


def f_rho_posterior(O, E, prior_sd, n=200000, seed=5, extra=None):
    """posterior draws of ln f (O ~ Poisson(E f), prior N(0, prior_sd^2)) by
    importance sampling; extra = (mean, sd): a further normal likelihood on ln f"""
    rng = np.random.default_rng(seed)
    x = rng.normal(0, prior_sd, n)
    lw = O * x - E * np.exp(x)
    if extra is not None:
        lw = lw - 0.5 * ((x - extra[0]) / extra[1]) ** 2
    w = np.exp(lw - lw.max())
    w /= w.sum()
    return x[rng.choice(n, size=n, p=w)]


def describe(d):
    return {"median_f_rho": float(np.exp(np.median(d))), "sd_ln_f_rho": float(d.std()),
            "median_f_rho2": float(np.exp(2 * np.median(d))), "mean_f_rho2": float(np.mean(np.exp(2 * d))),
            "sd_ln_f_rho2": float(2 * d.std()),
            "q16_f_rho": float(np.exp(np.quantile(d, 0.16))), "q84_f_rho": float(np.exp(np.quantile(d, 0.84))),
            "q05_f_rho": float(np.exp(np.quantile(d, 0.05))), "q95_f_rho": float(np.exp(np.quantile(d, 0.95)))}


def band(Ecell, cal, rng, sp_draws=None, sel_sd=0.0, nd=ND):
    lg0, _, _, lm0 = cal.rates()
    vals = []
    for lg, lm in cal.draws(rng, nd):
        f = np.exp(lg - lg0 + lm - lm0)
        sp = math.exp(rng.normal(0, 0.2)) if sp_draws is None else float(np.exp(rng.choice(sp_draws)))
        sel = math.exp(rng.normal(0, sel_sd)) if sel_sd else 1.0
        vals.append(float((Ecell * f).sum()) * math.exp(rng.normal(0, 0.2)) * sp * sel)
    return np.array(vals)


def qd(b):
    q = np.percentile(b, [5, 16, 50, 84, 95])
    return {"q05": q[0], "q16": q[1], "median": q[2], "q84": q[3], "q95": q[4], "mean": float(b.mean()),
            "sd_ln": float(np.log(b).std()), "ratio_90": float(q[4] / q[0])}


def class_sd(Ecell, cal, nd=3000, seed=11):
    rng = np.random.default_rng(seed)
    lg0, _, _, lm0 = cal.rates()
    v = [float((Ecell * np.exp(lg - lg0 + lm - lm0)).sum()) for lg, lm in cal.draws(rng, nd)]
    return float(np.log(v).std())


def main():
    cal0 = S.Calibration()
    cal1 = FU.fit_updated()[0]
    out = {"A_SQ": FU.A_SQ_MAIN, "prior_sd_ln_f_rho": PRIOR_SD, "selection_sd": SEL_SD}
    b312 = AN["bands"]["3-12k"]
    fS, fP = b312["S_per_square"]["ratio"], b312["P_per_square"]["ratio"]
    O = AN["sp"]["per_square_jeffreys"]["O"]
    E = AN["sp"]["per_square_jeffreys"]["E"] * fS * fP
    out["sp_count"] = {"O": O, "E_model": AN["sp"]["per_square_jeffreys"]["E"], "f_S": fS, "f_P": fP,
                       "E_given_fS_fP": E, "f_rho_mle": O / E}
    nonsp = json.load(open(os.path.join(OUT, "subevents_nonsp.json")))
    j4 = next(r for r in nonsp["rows"] if r["j"] == 4)
    variants = {
        "main: count only, prior 0.29": f_rho_posterior(O, E, PRIOR_SD),
        "count only, prior 0.145 (kappa only, no extrapolation term)": f_rho_posterior(O, E, 0.145),
        "count only, flat prior": f_rho_posterior(O, E, 10.0),
        "REJECTED (double counts): count + heuristic route N(ln 1.35, 0.25), prior 0.29":
            f_rho_posterior(O, E, PRIOR_SD, extra=(math.log(1.35), 0.25)),
    }
    out["f_rho"] = {k: describe(v) for k, v in variants.items()}
    out["subevent_j4_without_sp"] = {"coupling": j4["coupling_nonsp"], "ci90_approx": j4["ci90_nonsp_approx"],
                                     "O": j4["O_nonsp"], "E": j4["E_nonsp"],
                                     "note": "cross-check of the direction only (not combined)"}
    main_ln = variants["main: count only, prior 0.29"]
    f2_point = float(np.exp(2 * np.median(main_ln)))
    sp_draws = 2 * (main_ln - np.median(main_ln))
    rng = np.random.default_rng(11)
    for mark in (1, 10):
        E0, Eat0 = ecell("shipped", mark)
        E1, Eat1 = ecell("updated", mark)
        res = {"E_shipped_point": Eat0.get(float(mark)), "E_updated_point": Eat1.get(float(mark))}
        try:
            E2, Eat2 = ecell("time", mark)
            res["E_updated_time_point"] = Eat2.get(float(mark))
        except FileNotFoundError:
            E2 = None
        try:
            _, Eat3 = ecell("updated_a15", mark)
            res["E_updated_a15_point"] = Eat3.get(float(mark))
        except FileNotFoundError:
            pass
        res["E_main_point"] = res["E_updated_point"] * f2_point
        res["E_main_time_point"] = (res["E_updated_time_point"] * f2_point) if E2 is not None else None
        res["f_rho2_point"] = f2_point
        for name, v in variants.items():
            res.setdefault("E_point_by_f_rho_variant", {})[name] = res["E_updated_point"] * float(
                np.exp(2 * np.median(v)))
        b0 = band(E0, cal0, rng)
        # the shipped band with an honest SP term: the f_rho^2 prior (ln sd 2 x 0.29) in place of 0.2
        b0x = band(E0, cal0, rng, sp_draws=rng.normal(0, 2 * PRIOR_SD, 200000))
        b1 = band(E1, cal1, rng)
        b1s = band(E1, cal1, rng, sp_draws=sp_draws, sel_sd=SEL_SD) * f2_point
        for name, b in (("shipped", b0), ("shipped_with_sp_prior", b0x), ("updated_class_factors_only", b1),
                        ("main", b1s)):
            res[name] = qd(b)
        res["class_factor_sd_ln"] = {"shipped": class_sd(E0, cal0), "updated": class_sd(E1, cal1)}
        nb = np.arange(S.NCELLS) // 36
        names = ["<1k", "1-3k", "3-6k", "6-12k", "12-24k", ">=24k"]
        res["E_by_band"] = {names[i]: {"shipped": float(E0[nb == i].sum()), "updated": float(E1[nb == i].sum()),
                                       "time": (float(E2[nb == i].sum()) if E2 is not None else None)}
                            for i in range(6)}
        out[f"{mark}y"] = res

    # the existence study's level factor F (ln sd, 68%), component by component, before
    # (existence.md 4.1) and after, on the same basis (no extrapolation term for kappa; existence.md
    # then rounds 0.55 up to 0.7 for what is not modelled: the same factor is applied after)
    s312 = AN["bands"]["3-12k"]
    sq_lo, sq_hi = AN["bands"]["3-6k"]["squares"]["od90"]
    sq_sd = math.log(sq_hi / sq_lo) / (2 * 1.645)
    bS, bP = s312["S_per_square"]["boot90"], s312["P_per_square"]["boot90"]
    sp_sd = 2 * math.sqrt((math.log(bS[1] / bS[0]) / 3.29) ** 2 + (math.log(bP[1] / bP[0]) / 3.29) ** 2)
    kappa_only = out["f_rho"]["count only, prior 0.145 (kappa only, no extrapolation term)"]["sd_ln_f_rho2"]
    with_extrap = out["f_rho"]["main: count only, prior 0.29"]["sd_ln_f_rho2"]
    comp_before = {"squares level (K calibration, x0.85-1.25)": math.log(1.25 / 0.85) / 3.29,
                   "S / P rates (diagonal review factors)": 0.25,
                   "SP coupling (base kappa 0.46-1.2)": math.log(1.2 / 0.46) / 3.29,
                   "pair factor": 0.11, "model form": 0.30, "sampling and coverage": 0.11}
    comp_after = dict(comp_before)
    comp_after["squares level (K calibration, x0.85-1.25)"] = sq_sd
    comp_after["S / P rates (diagonal review factors)"] = sp_sd
    comp_after["SP coupling (base kappa 0.46-1.2)"] = kappa_only

    def tot(c):
        return math.sqrt(sum(v * v for v in c.values()))
    roundup = 0.7 / 0.55
    F = {"components_before": comp_before, "components_after": comp_after,
         "quadrature_before": tot(comp_before), "quadrature_after": tot(comp_after),
         "sigma_before (existence.md, rounded up)": 0.7, "sigma_after (same round-up)": tot(comp_after) * roundup,
         "with the SP extrapolation term on both sides": {
             "before": math.sqrt(tot(comp_before) ** 2 - comp_before["SP coupling (base kappa 0.46-1.2)"] ** 2
                                 + (2 * PRIOR_SD) ** 2),
             "after": math.sqrt(tot(comp_after) ** 2 - kappa_only ** 2 + with_extrap ** 2)},
         "note": "model form kept at 0.30: the squares do not test the magic-rung form (pair factor, "
                 "higher-order coupling); the large-X widening of existence.md (0.75 / 1.0 at X >= 2e4) is "
                 "not reduced by this search (7 sums at 12-45k)"}
    F["x_div_before"] = math.exp(0.7)
    F["x_div_after"] = math.exp(F["sigma_after (same round-up)"])
    out["F"] = F
    with open(os.path.join(OUT, "bands.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
