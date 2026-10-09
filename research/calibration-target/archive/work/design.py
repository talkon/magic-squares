#!/usr/bin/env python3
"""Steps 1-5 of the pre-registered calibration search: draw (P, S) per
stratum, choose modes and strides, write the predictions, plan.json and
run.sh. Nothing is run here.

Strata: N' band (3-6k, 6-12k, 12-24k, 24-45k) x assignment class (sorted,
ratio <= 1.1, ratio > 1.1) x k (<= 5, 6, >= 7).
Population and weights:
  3-6k, 6-12k   the sums of scheduler v2's 10 CPU-year plan (sim.py: shipped
                calibration, --dfirst auto, anchored charge, 10% candidate
                sample seed 1), weight = predicted magic squares of the sum;
  12-24k, 24-45k  the plan has (almost) no E there (1.6% / 0% at 10 CPU-years),
                so every pool sum with S <= 2 S0 (the existence study's
                weighting of where N_magic sits), weight = predicted magic per sum.
Draws: probability proportional to weight within substratum (class x k), the
substrata filled in proportion to their E shares (largest deficit first; a
thin substratum holding >= 1% of the band's E, k >= 7 or ratio > 1.1, gets
one draw first), one sum per P across the whole design, (P, S) already
searched excluded, until the band's CPU budget is used.
"""
import hashlib
import json
import math
import os
import random
import shlex
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cm  # noqa: E402
import explore  # noqa: E402
import calibrate as CAL  # noqa: E402  (scripts/, on the path through cm)

S = cm.S
am = cm.am
CT = os.path.dirname(HERE)
SEED = 20261009
YR = 8766.0
BANDS = ("3-6k", "6-12k", "12-24k", "24-45k")
EDGES = (3000.0, 6000.0, 12000.0, 24000.0, 45000.0)
RB = ("sorted", "r<=1.1", "r>1.1")
KB = ("k<=5", "k=6", "k>=7")
MSEARCH = os.path.join(CT, "bin", "msearch")

# the design's knobs (pre-registered)
TOTAL_BUDGET_H = 7.5         # predicted CPU-hours of the whole design (anchored truth)
COVER_H = {"12-24k": 0.85, "24-45k": 0.60}   # fixed coverage budgets
Q_STREAM = {"3-6k": 40, "6-12k": 40, "12-24k": 25, "24-45k": 8}   # target stream squares
D_CAP = {"3-6k": 1800.0, "6-12k": 1800.0, "12-24k": 300.0, "24-45k": 120.0}  # d-loop CPU / sum
S_CAP = {"3-6k": 1800.0, "6-12k": 1800.0, "12-24k": 650.0, "24-45k": 600.0}  # stream CPU / sum
PLAIN_BELOW = 4000.0         # N' below: plain
RATIO_RULE = 0.95            # plain also where r(N') + 1/k_cal >= this (d-first + stream ~ plain)
MIN_KCAL = 3
TL_FACTOR, TL_ADD = 2.5, 60.0     # msearch --time-limit (wall) = 2.5 pred + 60 (d-first runs)
CPU_FACTOR, CPU_ADD = 3.0, 90.0   # RLIMIT_CPU of the process = 3 pred + 90
HARD_CAP_H = 8.0
K_SUB = {"j3": 50.0, "j4": 10.0}   # sum + j exponents sub-events per SP traversal (new squares,
                                   # research/calibration.md: 149 and 31 vs 3 SP)


def band_of(Np):
    return int(np.searchsorted(EDGES, Np, side="right")) - 1


def kbin(k):
    return 0 if k <= 5 else (1 if k == 6 else 2)


def load_excluded(sch):
    with open(os.path.join(HERE, "searched.json")) as f:
        sr = {k: set(v) for k, v in json.load(f).items()}
    leg = {S.p_str(P, "_"): lg["maxS"] for P, lg in sch.legacy.items()}
    return sr, leg


def is_excluded(key, Sv, sr, leg):
    if key in sr and int(Sv) in sr[key]:
        return True
    if key in leg and int(Sv) <= leg[key]:
        return True
    return False


def run_design(p, band):
    """mode, strides and predicted CPU of one sum (p: cm.predict at that sum, scalars)"""
    Np, sq, td, tpa, r = p["Np"], p["sq"], p["td"], p["tp_anch"], p["r"]
    Q = Q_STREAM[BANDS[band]]
    kc = max(1, int(round(sq / Q)))
    plain = Np < PLAIN_BELOW or kc < MIN_KCAL or r + 1.0 / kc >= RATIO_RULE
    if BANDS[band] in ("12-24k", "24-45k"):
        plain = False
    if plain:
        return dict(mode="plain", kd=1, kc=0, cpu_dloop=0.0, cpu_stream=0.0, cpu=tpa + 0.1)
    # the stream: ~Q squares, but at most S_CAP CPU-seconds (fewer squares then)
    kc = max(kc, 1, int(math.ceil(tpa / S_CAP[BANDS[band]])))
    kd = max(1, int(math.ceil(td / D_CAP[BANDS[band]])))
    cd = td / kd
    cs = tpa / kc
    return dict(mode="dfirst", kd=kd, kc=kc, cpu_dloop=cd, cpu_stream=cs, cpu=cd + cs + 0.2)


def scal(p):
    return {k: (float(v[0]) if np.ndim(v) else v) for k, v in p.items()}


def plan_population(d, band, used_P, sr, leg, sch):
    """plan sums of a band: indices into d, weights, substratum"""
    b = np.searchsorted(np.log(EDGES), d["lNp"], side="right") - 1
    s = (b == band) & (d["h"] < 10 * YR) & (d["magic"] > 0)
    idx = np.nonzero(s)[0]
    keys = {}
    keep = []
    for i, ai, Si in zip(idx, d["a"][idx], d["S"][idx]):
        key = keys.get(ai)
        if key is None:
            key = keys[ai] = sch.cands.key(int(ai))
        if not is_excluded(key, Si, sr, leg):
            keep.append(i)
    idx = np.array(keep, int)
    aa = d["a"][idx]
    kk = sch.cands.k[aa]
    sub = 3 * sch.cands.rbin[aa].astype(int) + np.where(kk <= 5, 0, np.where(kk == 6, 1, 2))
    return idx, d["magic"][idx], sub


def draw_design(sch, rng, prng):
    d = dict(np.load(os.path.join(HERE, "plan_sums_10y_s01.npz")))
    e = dict(np.load(os.path.join(HERE, "exist_grid.npz")))
    sr, leg = load_excluded(sch)
    eb = np.searchsorted(np.log(EDGES), e["lNp"], side="right") - 1
    # --- step 1 numbers: E composition (for plan.json)
    comp = {}
    for lim in (1, 3, 10):
        s = d["h"] < lim * YR
        M = d["magic"][s].sum()
        bb = np.searchsorted(np.log(EDGES), d["lNp"][s], side="right") - 1
        comp[f"{lim}y"] = {"E": float(M / d["frac"]),
                           "N'": {nm: float(d["magic"][s][bb == i].sum() / M)
                                  for i, nm in enumerate(BANDS)},
                           "N'<3k": float(d["magic"][s][bb < 0].sum() / M)}
    exist_share = {BANDS[i]: float(e["w"][eb == i].sum() / e["w"].sum()) for i in range(4)}
    # --- budgets: coverage bands fixed, the rest by Neyman allocation for the
    # E-weighted SP factor (CPU_h proportional to w_h / sqrt(rho_h)), w_h the 10 CPU-year
    # plan's E shares, rho_h the SP pairs per CPU-hour of the band's design (pilot draws)
    pilot = {}
    for band in (0, 1):
        idx, w, sub = plan_population(d, band, set(), sr, leg, sch)
        pk = rng.choice(len(idx), size=150, p=w / w.sum())
        sp = cpu = 0.0
        for j in pk:
            i = idx[j]
            p = scal(cm.predict(sch, int(d["a"][i]), [float(d["S"][i])]))
            rd = run_design(p, band)
            sp += p["pairs"] / rd["kd"]
            cpu += rd["cpu"]
        pilot[BANDS[band]] = 3600 * sp / cpu
    w10 = comp["10y"]["N'"]
    rest = TOTAL_BUDGET_H - sum(COVER_H.values())
    a = {b: w10[b] / math.sqrt(pilot[b]) for b in ("3-6k", "6-12k")}
    budget = {b: rest * a[b] / sum(a.values()) for b in a}
    budget.update(COVER_H)
    # --- draws
    used_P = set()
    runs = []
    for band in range(4):
        bn = BANDS[band]
        if band < 2:
            idx, w, sub = plan_population(d, band, used_P, sr, leg, sch)
            cand_a = d["a"][idx]
            cand_S = d["S"][idx]
        else:
            s = eb == band
            ea, ew = e["a"][s], e["w"][s]
            # per candidate: its weight in the band
            ua, inv = np.unique(ea, return_inverse=True)
            wa = np.bincount(inv, weights=ew)
            idx = np.arange(len(ua))
            w = wa
            sub = np.array([3 * int(sch.cands.rbin[x]) + kbin(int(sch.cands.k[x])) for x in ua])
            cand_a, cand_S = ua, None
        share = np.array([w[sub == c].sum() for c in range(9)]) / w.sum()
        forced = [c for c in range(9) if share[c] >= 0.01 and (c % 3 == 2 or c // 3 == 2)]
        prng.shuffle(forced)
        n_c = np.zeros(9, int)
        spent = 0.0
        tries = 0
        nb = 0
        while spent < budget[bn] * 3600 and tries < 400:
            tries += 1
            if forced:
                c = forced.pop(0)
            else:
                c = int(np.argmax(share * (nb + 1) - n_c))
            pool = np.nonzero(sub == c)[0]
            if len(pool) == 0:
                share[c] = 0
                continue
            j = pool[rng.choice(len(pool), p=w[pool] / w[pool].sum())]
            aa = int(cand_a[j])
            key = sch.cands.key(aa)
            if key in used_P:
                continue
            if cand_S is not None:
                Sv = float(cand_S[j])
                p = scal(cm.predict(sch, aa, [Sv]))
            else:
                # the sums of this P in the band: S drawn by predicted magic per sum
                Sg = sch.cands.S0[aa] * (1 + am.GRID_U)
                lo = int(math.floor(Sg[0]))
                hi = int(math.floor(2 * sch.cands.S0[aa]))
                Ss = np.arange(lo, hi + 1, dtype=float)
                pr = cm.predict(sch, aa, Ss)
                okb = (pr["Np"] >= EDGES[band]) & (pr["Np"] < EDGES[band + 1]) & (pr["magic"] > 0)
                okb &= np.array([not is_excluded(key, x, sr, leg) for x in Ss])
                if not okb.any():
                    continue
                ii = np.nonzero(okb)[0]
                Sv = float(Ss[ii[rng.choice(len(ii), p=pr["magic"][ii] / pr["magic"][ii].sum())]])
                p = scal(cm.predict(sch, aa, [Sv]))
            rd = run_design(p, band)
            if spent + rd["cpu"] > budget[bn] * 3600:
                # too big for what is left of the band: try a smaller draw (a few times)
                if tries > 300:
                    break
                continue
            used_P.add(key)
            n_c[c] += 1
            nb += 1
            spent += rd["cpu"]
            runs.append(dict(band=bn, sub=c, a=aa, key=key, P=list(sch.cands.P(aa)), S=int(Sv),
                             p=p, rd=rd, share_sub=float(share[c])))
        for r in runs:
            if r["band"] == bn:
                r["n_sub"] = int(n_c[r["sub"]])
                r["design_weight"] = r["share_sub"] / max(r["n_sub"], 1)
    return runs, comp, exist_share, pilot, budget


def poisson90(n):
    lo, hi = CAL.poisson_ci(n, 0.90)
    return lo, hi


def precision(E, phi=1.0, nP=1, between=0.0):
    """90% interval on obs/pred if obs came in at its expectation: exact Poisson
    (Garwood) at phi = 1 and no between-P term; else a lognormal approximation
    with variance phi / E + between^2 / nP"""
    if E < 1.0:
        return None     # less than one expected event: no information
    if phi == 1.0 and between == 0.0:
        n = int(round(E))
        lo, hi = poisson90(n)
        return [lo / E if n > 0 else 0.0, hi / E]
    sd = math.sqrt(phi / E + between ** 2 / max(nP, 1))
    return [math.exp(-1.645 * sd), math.exp(1.645 * sd)]


def main():
    rng = np.random.default_rng(SEED)
    prng = random.Random(SEED)
    sch = explore.load()
    t0 = time.time()
    runs, comp, exist_share, pilot, budget = draw_design(sch, rng, prng)
    # order: interleave the strata by priority (band order 6-12k, 3-6k, 12-24k, 24-45k,
    # round robin) so that a budget stop truncates every band alike
    byb = {b: [r for r in runs if r["band"] == b] for b in BANDS}
    for b in byb:
        prng.shuffle(byb[b])
    order = []
    prio = ("6-12k", "3-6k", "12-24k", "24-45k")
    # interleave in proportion to each band's predicted CPU
    tot = {b: sum(r["rd"]["cpu"] for r in byb[b]) for b in BANDS}
    acc = {b: 0.0 for b in BANDS}
    while any(byb.values()):
        cand = [b for b in prio if byb[b]]
        b = min(cand, key=lambda x: (acc[x] / max(tot[x], 1e-9), prio.index(x)))
        r = byb[b].pop(0)
        acc[b] += r["rd"]["cpu"]
        order.append(r)
    preds = []
    lines = []
    for i, r in enumerate(order):
        rid = f"R{i + 1:03d}"
        p, rd = r["p"], r["rd"]
        out = os.path.join(CT, "runs", f"{rid}.jsonl")
        seed = int(rng.integers(1, 2 ** 31 - 1))
        cmd = ["nice", "-n", "19", MSEARCH, "--min-sum", str(r["S"]), "--max-sum", str(r["S"])]
        cpu_limit = int(math.ceil(CPU_FACTOR * rd["cpu"] + CPU_ADD))
        if rd["mode"] == "dfirst":
            nd_pred = max(int(round(p["Np"])), 1)
            d_off = int(rng.integers(0, rd["kd"])) if rd["kd"] > 1 else 0
            d_chunk = max(8, int(math.ceil(nd_pred / 16)))
            tl = int(math.ceil(TL_FACTOR * rd["cpu"] + TL_ADD))
            cmd += ["--diag-first", "--diag-first-min-n", "0", "--d-chunk", str(d_chunk),
                    "--calib-r1-stride", str(rd["kc"]), "--time-limit", str(tl),
                    "--sample-seed", str(seed)]
            if rd["kd"] > 1:
                cmd += ["--d-stride", str(rd["kd"]), "--d-offset", str(d_off)]
        else:
            d_off, d_chunk, tl = None, None, None
        cmd += ["--out", out] + [str(x) for x in r["P"]]
        kd, kc = rd["kd"], rd["kc"]
        sq = p["sq"]
        rec_sq = sq if rd["mode"] == "plain" else sq / kc
        pairs_found = p["pairs"] / kd if rd["mode"] == "dfirst" else p["pairs"]
        stream_sp = rec_sq * cm.T * p["qSP"] if rd["mode"] == "dfirst" else 0.0
        if rd["mode"] == "plain":
            mfound = p["magic"]
        else:
            pd = 1 - (1 - 1 / kd) ** 2 if kd > 1 else 1.0
            mfound = p["magic"] * (1 - (1 - pd) * (1 - 1 / kc))
        frac_d = 1.0 / kd if rd["mode"] == "dfirst" else 0.0
        pred = {
            "id": rid, "out": out, "band": r["band"], "class": RB[r["sub"] // 3], "kbin": KB[r["sub"] % 3],
            "substratum_share_of_band_E": r["share_sub"], "draws_in_substratum": r["n_sub"],
            "design_weight": r["design_weight"],
            "P": r["P"], "P_key": r["key"], "S": r["S"], "k": int(sch.cands.k[r["a"]]),
            "tau": int(sch.cands.tau[r["a"]]), "ratio": float(sch.cands.ratio[r["a"]]),
            "S0": float(sch.cands.S0[r["a"]]), "smin_approx": float(sch.cands.smin[r["a"]]),
            "S_over_smin": float(r["S"] / sch.cands.smin[r["a"]]),
            "x": p["x"], "cell": int(p["cell"]), "cell_parts": list(S.cell_parts(p["cell"])),
            "N_pred": p["Np"], "labels_pred": p["labels"],
            "mode": rd["mode"], "d_stride": kd, "d_offset": d_off, "d_chunk": d_chunk,
            "calib_r1_stride": kc, "sample_seed": seed if rd["mode"] == "dfirst" else None,
            "time_limit": tl, "cpu_limit": cpu_limit,
            "cpu": {
                "sum_plain_law": p["tp"], "sum_dfirst_law": p["td"], "ratio_r": p["r"],
                "sum_plain_by_ratio": p["tp_ratio"], "sum_dfirst_by_ratio": p["td_ratio"],
                "sum_anchored_plain": p["tp_anch"],
                "sum_anchored_this_mode": p["td"] if rd["mode"] == "dfirst" else p["tp_anch"],
                "run_dloop": rd["cpu_dloop"], "run_stream": rd["cpu_stream"],
                "run_pred_anchored": rd["cpu"],
                "run_pred_plain_law": (rd["cpu_dloop"] + p["tp"] / kc + 0.2) if rd["mode"] == "dfirst"
                else p["tp"] + 0.1,
            },
            "squares_sum": sq, "squares_recorded": rec_sq,
            "q_S": p["qS"], "q_P": p["qP"], "q_SP": p["qSP"],
            "S_trav_per_square": cm.T * p["qS"], "P_trav_per_square": cm.T * p["qP"],
            "SP_trav_per_square": cm.T * p["qSP"],
            "S_trav_recorded": rec_sq * cm.T * p["qS"], "P_trav_recorded": rec_sq * cm.T * p["qP"],
            "sp_pairs_sum": p["pairs"], "sp_pairs_found_dfirst": pairs_found if rd["mode"] == "dfirst" else 0.0,
            "sp_pairs_found_plain": p["pairs"] if rd["mode"] == "plain" else 0.0,
            "sp_pairs_found": pairs_found, "sp_in_stream_squares": stream_sp,
            "spS_pairs_sum": p["sps"], "spP_pairs_sum": p["spp"],
            "spS_pairs_found": p["sps"] * (frac_d if rd["mode"] == "dfirst" else 1.0),
            "spP_pairs_found": p["spp"] * (frac_d if rd["mode"] == "dfirst" else 1.0),
            "pmagic_per_square": p["pmagic"], "pair_factor": p["pairf"],
            "magic_sum": p["magic"], "magic_expected_found": mfound,
            "subevents_rough": {"sum+3exp": K_SUB["j3"] * (rec_sq * cm.T * p["qSP"]),
                                "sum+4exp": K_SUB["j4"] * (rec_sq * cm.T * p["qSP"])},
        }
        preds.append(pred)
        lines.append(f"{rid} pred={rd['cpu']:.1f} limit={cpu_limit} :: " + shlex.join(cmd))
    # --- per band and total expectations and precision
    summ = {}
    for b in BANDS + ("all",):
        rs = [q for q in preds if b == "all" or q["band"] == b]
        if not rs:
            continue
        nP = len({q["P_key"] for q in rs})
        Esp = sum(q["sp_pairs_found"] for q in rs)
        Esq = sum(q["squares_recorded"] for q in rs)
        ES = sum(q["S_trav_recorded"] for q in rs)
        EP = sum(q["P_trav_recorded"] for q in rs)
        Esps = sum(q["spS_pairs_found"] for q in rs)
        Espp = sum(q["spP_pairs_found"] for q in rs)
        Esub3 = sum(q["subevents_rough"]["sum+3exp"] for q in rs)
        Esub4 = sum(q["subevents_rough"]["sum+4exp"] for q in rs)
        Esp_sq = sum(q["sp_in_stream_squares"] + q["sp_pairs_found_plain"] for q in rs)
        cpu = sum(q["cpu"]["run_pred_anchored"] for q in rs)
        summ[b] = {
            "runs": len(rs), "distinct_P": nP, "dfirst_runs": sum(q["mode"] == "dfirst" for q in rs),
            "d_sampled_runs": sum(q["d_stride"] > 1 for q in rs),
            "cpu_hours_pred_anchored": cpu / 3600,
            "cpu_hours_pred_plain_law": sum(q["cpu"]["run_pred_plain_law"] for q in rs) / 3600,
            "E_sp_pairs": Esp, "E_squares_recorded": Esq, "E_S_trav": ES, "E_P_trav": EP,
            "E_spS_pairs": Esps, "E_spP_pairs": Espp, "E_magic_found": sum(q["magic_expected_found"] for q in rs),
            "E_sp_in_recorded_squares": Esp_sq,
            "E_subevents_sum+3exp_rough": Esub3, "E_subevents_sum+4exp_rough": Esub4,
            "sp_pairs_per_cpu_hour": Esp / max(cpu / 3600, 1e-9),
            "precision90": {
                "sp_pairs_obs_over_pred_poisson": precision(Esp),
                "sp_pairs_obs_over_pred_with_between_P": precision(Esp, 1.3, nP, 0.36),
                "squares_obs_over_pred": precision(Esq, S.PHI_SUM, nP, 0.26),
                "S_trav_per_square": precision(ES, 1.0, nP, 0.18),
                "P_trav_per_square": precision(EP, 1.0, nP, 0.18),
                "spS_plus_spP_pairs": precision(Esps + Espp),
            },
        }
    rho_best = max(pilot.values())
    n_all_best = TOTAL_BUDGET_H * rho_best
    lo_b, hi_b = poisson90(int(round(n_all_best)))
    verdict = {
        "sp_count_useful": False,
        "why": (f"the design expects {summ['all']['E_sp_pairs']:.1f} (square, SP diagonal) pairs in "
                f"{summ['all']['cpu_hours_pred_anchored']:.2f} predicted CPU-hours: if they came in at "
                f"their expectation, obs/pred would have a 90% interval of "
                f"{summ['all']['precision90']['sp_pairs_obs_over_pred_poisson']} (P(magic) moves with its "
                f"square). Even the whole 7.5 h in the band with the most SP pairs per CPU-hour "
                f"({rho_best:.2f}/h, 3-6k) gives ~{n_all_best:.1f} pairs, interval "
                f"[{lo_b / n_all_best:.2f}, {hi_b / n_all_best:.2f}]. A +-25% (90%) SP rate needs ~45 "
                f"pairs (~{45 / rho_best:.0f} CPU-hours at 3-6k), +-15% ~120 pairs "
                f"(~{120 / rho_best:.0f} CPU-hours); 12-45k gives < 0.05 pairs per CPU-hour"),
        "sp_rungs_useful": False,
        "sp_rungs_why": (f"SP+S and SP+P pairs: {summ['all']['E_spS_pairs']:.3f} and "
                         f"{summ['all']['E_spP_pairs']:.3f} expected in the whole design (each SP "
                         "diagonal has 15 partners, each S-only / P-only with probability ~1e-3 / "
                         "3e-4): no information"),
        "lean_on": [
            f"the squares: ~{summ['all']['E_squares_recorded']:.0f} recorded (plain sums and streams), "
            f"squares obs/pred to {summ['all']['precision90']['squares_obs_over_pred']} (90%), 3-6k "
            f"{summ['3-6k']['precision90']['squares_obs_over_pred']}, 6-12k "
            f"{summ['6-12k']['precision90']['squares_obs_over_pred']}",
            f"S and P traversals per square: ~{summ['all']['E_S_trav']:.0f} and "
            f"~{summ['all']['E_P_trav']:.0f} events, rates to "
            f"{summ['all']['precision90']['S_trav_per_square']} and "
            f"{summ['all']['precision90']['P_trav_per_square']}",
            f"the sub-events on the same squares (calibrate.py): roughly "
            f"{summ['all']['E_subevents_sum+3exp_rough']:.0f} sum + 3 exponent and "
            f"{summ['all']['E_subevents_sum+4exp_rough']:.0f} sum + 4 exponent events (50x / 10x the "
            "expected SP traversals in the recorded squares, the ratio on the 2,280 new squares of "
            "research/calibration.md), which test the sum-product coupling, the part of the SP rate "
            "beyond the S and P rates, to roughly +-12% / +-25%, and the pair factor below the magic "
            "rung (d exponents on both diagonals)",
            "the direct SP count as a check of those (a factor-2 excess, as on the re-found seed P, "
            "would show as ~11 pairs against ~5.5: p ~ 0.02)",
            "the time laws: the run CPU against the plain, d-first, ratio and anchored predictions, "
            "with d-first and stream (plain) estimates at N' 15-38k where the laws rest on few sums",
        ],
        "analysis_cpu": "calibrate.py's heuristic with sub-events takes ~0.25 s per square: ~45 "
                        "CPU-minutes for the expected squares (outside the search budget; run niced, "
                        "--jobs 2)",
        "sp_s_alternative": "SP+S, SP+P alternatives: x1.19 (the S+S residual of calibration.md) and "
                            "x1/0.87 (squares with an SP diagonal are S-richer: per-P i.i.d. rates "
                            "give 0.87 of the own-count null there)",
    }
    plan = {
        "verdict": verdict,
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "question": "have we run a quick search around the ranges of N where we expect the first "
                    "magic square, to calibrate our estimates?",
        "repo": {"path": "/home/user/magic-squares", "commit": "e0ace5436d2fa69f2f33966ee5333a8a5146f9cc",
                 "msearch": MSEARCH, "msearch_sha256": hashlib.sha256(open(MSEARCH, "rb").read()).hexdigest(),
                 "msearch_source": "/home/user/magic-squares/cmake-build-release/src/c/msearch (engine 3)"},
        "state": os.path.join(CT, "state"),
        "seed": SEED,
        "step1_where_E_sits": {
            "method": "work/sim.py: scheduler v2 greedy on a fresh state + legacy import, shipped "
                      "calibration, --dfirst auto, units charged under --truth anchored, 10% "
                      "candidate sample (seed 1), budget scaled",
            "E_and_composition_by_budget": comp,
            "existence_weighting_share_N'_3-45k": exist_share,
        },
        "strata": {"N'_bands": BANDS, "class": RB, "k": KB,
                   "population": {"3-6k": "10 CPU-year plan sums, weight predicted magic",
                                  "6-12k": "10 CPU-year plan sums, weight predicted magic",
                                  "12-24k": "all pool sums S <= 2 S0, weight predicted magic per sum",
                                  "24-45k": "all pool sums S <= 2 S0, weight predicted magic per sum"}},
        "budget_hours": budget, "pilot_sp_pairs_per_cpu_hour": pilot,
        "budget_rule": f"{TOTAL_BUDGET_H} predicted CPU-hours (anchored); 12-24k and 24-45k fixed "
                       "coverage budgets; the rest split between 3-6k and 6-12k by Neyman allocation "
                       "for the E-weighted SP factor: CPU_h ~ w_h / sqrt(rho_h), w_h = the 10 CPU-year "
                       "plan's E shares, rho_h = expected SP pairs per CPU-hour of the band's design",
        "mode_rule": {"plain": f"N' < {PLAIN_BELOW:.0f}, or k_cal = round(squares / Q) < {MIN_KCAL}, or "
                               f"r(N') + 1/k_cal >= {RATIO_RULE} (d-first plus its stream would cost "
                               "about the plain search, which records every square and every SP pair "
                               "too); never in 12-45k",
                      "dfirst": "msearch --diag-first --diag-first-min-n 0, calibration stream "
                                "--calib-r1-stride k_cal = round(predicted squares / Q_band) for "
                                "~Q_band squares; --d-stride k_d = ceil(t_dfirst / D_cap_band) with a "
                                "random --d-offset where the whole d loop exceeds D_cap_band",
                      "Q_band": Q_STREAM, "D_cap_band_cpu_s": D_CAP},
        "limits": {"time_limit": f"d-first runs: --time-limit {TL_FACTOR} x pred + {TL_ADD:.0f} s (wall; "
                                 "stops between chunks of d, the sum is then recorded as not complete "
                                 "and its stream skipped)",
                   "cpu_limit": f"every process: RLIMIT_CPU {CPU_FACTOR} x pred + {CPU_ADD:.0f} s",
                   "hard_cap_hours": HARD_CAP_H,
                   "runner": "runner.py: at most 2 processes, nice 19; before each launch, spent + the "
                             "running processes' CPU limits + this run's limit <= 8 CPU-hours (the "
                             "limit is shrunk to fit, the run skipped if that leaves < 1.25 x its "
                             "prediction); ledger runs/ledger.jsonl; resumable"},
        "exclusions": "(P, S) with any record in the scratchpad JSONL/tar outputs or in "
                      "research/wip/data/all_records.tsv.xz (work/searched.json), and legacy sums "
                      "S <= maxS of the first search (state/legacy.json); one sum per P",
        "analysis_plan": {
            "primary": [
                "SP rate: per band and pooled, R_SP = sum obs / sum pred of the (square, SP diagonal) "
                "pairs found (d-first: dsum 'pairs' (d_stride 1) or 'pairs' against pred/k_d (sampled); "
                "plain: sum of sp_count over the square records), Garwood 90% interval and a P bootstrap; "
                "pooled over bands with the 10 CPU-year E shares (3-6k, 6-12k) as weights",
                "SP per square: pairs / squares, the squares from plain sums and the streams' est_squares "
                "(the stream's own csquares give a second, independent SP-per-square estimate)",
                "squares: obs (plain) or est_squares (stream) / pred per band, quasi-Poisson phi 2.5, P "
                "bootstrap; S and P traversals per square obs / pred (720 q_S, 720 q_P)",
                "SP+S and SP+P pairs: counted from the dsquare / square grids (partners of each SP "
                "diagonal), obs / pred; expected to be ~0 events (stated below)",
                "time: run CPU obs / pred per law (plain law, d-first law, ratio law, anchored); the "
                "d-first law from dsum est_time, the plain law from plain sums and the streams' est_time",
            ],
            "secondary": [
                "scripts/calibrate.py on all square and csquare records: the ladder and the sub-events "
                "(sum + j exponents pinned, d exponents on both diagonals) against the heuristic, by band; "
                "the coupling at j = 3, 4, 5 tests the SP rate indirectly with 10-50x the events",
                "design-weighted (Hansen-Hurwitz, weights design_weight) versions of every ratio",
                "magic squares found (expected ~0): any is reported at once",
            ],
            "do_not_change": "predictions.json is frozen (sha256 in PREREGISTERED.txt); the analysis "
                             "uses it as written, also for runs that are truncated, stopped or skipped "
                             "(those are reported as such, with their partial counts against the "
                             "prediction scaled by the d fraction searched)",
        },
        "summary_expectations": summ,
        "runs": [{k: q[k] for k in ("id", "band", "class", "kbin", "P", "S", "mode", "d_stride", "d_offset",
                                    "calib_r1_stride", "time_limit", "cpu_limit", "out")} for q in preds],
        "commands": lines,
    }
    predictions = {
        "created": plan["created"],
        "frozen_before_any_run": True,
        "model": {"scheduler": "scripts/scheduler.py v2 (analytic), shipped calibration (Calibration() "
                               "priors, time_prior(engine 3) plain and d-first, ratio level -0.028), no "
                               "per-P factors",
                  "amodel_version": am.AMODEL_VERSION, "profile_hash": am.profile_hash(),
                  "model_hash": am.model_hash(), "commit": plan["repo"]["commit"],
                  "formulas": cm.__doc__},
        "summary": summ,
        "runs": preds,
    }
    with open(os.path.join(CT, "predictions.json"), "w") as f:
        json.dump(predictions, f, indent=1, sort_keys=False, default=float)
    with open(os.path.join(CT, "plan.json"), "w") as f:
        json.dump(plan, f, indent=1, default=float)
    with open(os.path.join(HERE, "commands.txt"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(json.dumps({"budget": budget, "pilot": pilot, "summary": summ}, indent=1, default=float))
    print(f"{len(preds)} runs, {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
