#!/usr/bin/env python3
"""Pre-registered analysis of the calib-target runs: observed vs the frozen
predictions.json, per run and pooled per stratum (band, band x class x k),
with Poisson and overdispersed (between-sum) intervals.

Writes analysis/results.json and analysis/per_run.tsv; prints markdown tables.
"""
import collections
import hashlib
import json
import math
import os
import sys

import numpy as np

CT = "/tmp/claude-0/-home-user-magic-squares/f8940ae0-7961-577c-be6c-6db2a2de5f8e/scratchpad/calib-target"
WT = "/home/user/magic-squares/.claude/worktrees/calib-target"
sys.path.insert(0, os.path.join(WT, "scripts"))
import calibrate as CAL  # noqa: E402  (poisson_ci, poisson_tail, gamma_quantile)

SHA = "c1e673b595d7de9452a5c2c01113967c4c6579220c885a3d25e15ed270b06c08"
BANDS = ["3-6k", "6-12k", "12-24k", "24-45k"]
E10 = {"3-6k": 0.54, "6-12k": 0.26}          # share of the 10 CPU-year E (plan: pooling weights)
B = 4000
rng = np.random.default_rng(20261009)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def garwood(k, E, conf=0.90):
    if E <= 0:
        return None
    lo, hi = CAL.poisson_ci(int(round(k)), conf)
    return [lo / E, hi / E]


def lognorm_ci(r, sd):
    return [r * math.exp(-1.645 * sd), r * math.exp(1.645 * sd)]


def fmt(x, d=3):
    if x is None:
        return "-"
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(fmt(v, d) for v in x) + "]"
    if x == 0:
        return "0"
    if abs(x) >= 1000:
        return f"{x:,.0f}"
    return f"{x:.{d}g}"


# ---------------------------------------------------------------------------
pred_path = os.path.join(CT, "predictions.json")
assert sha256(pred_path) == SHA, "predictions.json changed"
PRED = json.load(open(pred_path))
pre = {}
for line in open(os.path.join(CT, "PREREGISTERED.txt")):
    if "sha256" in line and ":" in line:
        name, h = line.split("sha256")[0].strip(), line.split(":")[-1].strip().split()[0]
        pre[name] = h
for name, rel in (("plan.json", "plan.json"), ("run.sh", "run.sh"), ("runner.py", "runner.py"),
                  ("bin/msearch", "bin/msearch")):
    assert sha256(os.path.join(CT, rel)) == pre[name], name

runs = PRED["runs"]


def load(path):
    recs = collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            recs[r["type"]].append(r)
    return recs


rows = []
for p in runs:
    R = load(p["out"])
    o = {"id": p["id"], "band": p["band"], "class": p["class"], "kbin": p["kbin"], "mode": p["mode"],
         "k": p["k"], "tau": p["tau"], "x": p["x"], "S_over_smin": p["S_over_smin"], "Np": p["N_pred"],
         "dw": p["design_weight"], "P_key": p["P_key"], "S": p["S"]}
    assert len(R["done"]) == 1
    if p["mode"] == "plain":
        s = R["sum"][0]
        sq = R["square"]
        assert s["squares"] == len(sq) and not s["truncated"]
        o.update(N=s["nvecs_raw"], labels=s["labels"], cpu=s["cpu"], sq_obs=len(sq), sq_rec=len(sq),
                 sq_var=len(sq), w_rec=1.0,
                 S_obs=sum(q["s_count"] for q in sq), P_obs=sum(q["p_count"] for q in sq),
                 S_raw=sum(q["s_count"] for q in sq), P_raw=sum(q["p_count"] for q in sq),
                 SP_obs=sum(q["sp_count"] for q in sq),
                 spS_sq=sum(1 for q in sq if q["best_score"] == 9),
                 spP_sq=sum(1 for q in sq if q["best_score"] == 10),
                 magic=sum(1 for q in sq if q["best_score"] >= 14),
                 sp_stream=0, cpu_dloop=0.0, cpu_stream=0.0, d_stride=1,
                 t_plain_whole=s["cpu"], t_dfirst_whole=None)
        # expectations
        o.update(E_sq=p["squares_sum"], E_sq_rec=p["squares_recorded"],
                 E_SP_found=p["sp_pairs_found"],
                 E_SP_cond=len(sq) * p["SP_trav_per_square"],
                 E_S_cond=len(sq) * p["S_trav_per_square"], E_P_cond=len(sq) * p["P_trav_per_square"],
                 E_S_rec=p["S_trav_recorded"], E_P_rec=p["P_trav_recorded"])
    else:
        ds, cs = R["dsum"][0], R["csum"][0]
        csq = R["csquare"]
        k = cs["r1_stride"]
        assert cs["squares"] == len(csq) and not ds["truncated"] and not cs["truncated"]
        kd = ds["d_stride"]
        o.update(N=ds["nvecs_raw"], labels=ds["labels"], cpu=ds["cpu"] + cs["cpu"],
                 sq_obs=cs["est_squares"], sq_rec=len(csq), sq_var=cs["se_squares"] ** 2, w_rec=float(k),
                 S_obs=k * cs["s_trav"], P_obs=k * cs["p_trav"], SP_obs=ds["pairs"],
                 S_raw=cs["s_trav"], P_raw=cs["p_trav"],
                 spS_sq=len([q for q in R["dsquare"] if q["best_score"] == 9]),
                 spP_sq=len([q for q in R["dsquare"] if q["best_score"] == 10]),
                 magic=ds["magic_pairs"], sp_stream=cs["sp_pairs"], cpu_dloop=ds["cpu"], cpu_stream=cs["cpu"],
                 d_stride=kd, r1_stride=k, nd=ds["nd"],
                 t_dfirst_whole=ds["cpu"] - ds["time"] + ds["est_time"],
                 t_plain_whole=cs["est_time"] + cs["enum_time"] + cs["reduce_time"],
                 se_tplain=cs["se_time"])
        assert sum(1 for d in range(ds["d_offset"], ds["nvecs_raw"], kd)) == ds["nd"]
        o.update(E_sq=p["squares_sum"], E_sq_rec=p["squares_recorded"],
                 E_SP_found=p["sp_pairs_found_dfirst"], E_SP_stream=p["sp_in_stream_squares"],
                 E_SP_cond=cs["est_squares"] * p["SP_trav_per_square"] / kd,
                 E_S_cond=k * len(csq) * p["S_trav_per_square"], E_P_cond=k * len(csq) * p["P_trav_per_square"],
                 E_S_rec=p["S_trav_recorded"], E_P_rec=p["P_trav_recorded"])
    c = p["cpu"]
    o.update(cpu_pred=c["run_pred_anchored"], cpu_pred_plainlaw=c["run_pred_plain_law"],
             tp=c["sum_plain_law"], td=c["sum_dfirst_law"], r=c["ratio_r"], tp_ratio=c["sum_plain_by_ratio"],
             td_ratio=c["sum_dfirst_by_ratio"], tp_anch=c["sum_anchored_plain"],
             pred_dloop=c["run_dloop"], pred_stream=c["run_stream"],
             E_spS=p["spS_pairs_found"], E_spP=p["spP_pairs_found"], E_magic=p["magic_expected_found"],
             magic_sum=p["magic_sum"], pmagic=p["pmagic_per_square"], qS=p["q_S"], qP=p["q_P"], qSP=p["q_SP"],
             Spsq=p["S_trav_per_square"], Ppsq=p["P_trav_per_square"], SPpsq=p["SP_trav_per_square"],
             dens=p["magic_sum"] / max(c["sum_anchored_this_mode"], 1e-9))
    rows.append(o)

ID = {o["id"]: o for o in rows}


def arr(key, sel):
    return np.array([o[key] if o.get(key) is not None else np.nan for o in sel], float)


def boot_ratio(sel, okey, ekey, strat=True):
    """percentile 90% interval of sum O / sum E over a bootstrap of the sums
    (stratified by band)"""
    O, E = arr(okey, sel), arr(ekey, sel)
    bands = np.array([o["band"] for o in sel])
    idx_by = [np.nonzero(bands == b)[0] for b in BANDS if (bands == b).any()] if strat else [np.arange(len(sel))]
    res = np.zeros(B)
    num = np.zeros(B)
    den = np.zeros(B)
    for ix in idx_by:
        n = len(ix)
        pick = ix[rng.integers(0, n, size=(B, n))]
        num += O[pick].sum(1)
        den += E[pick].sum(1)
    with np.errstate(divide="ignore", invalid="ignore"):
        res = num / den
    res = res[np.isfinite(res)]
    return [float(np.quantile(res, 0.05)), float(np.quantile(res, 0.95))]


def hh(sel, okey, ekey):
    """design-weighted (Hansen-Hurwitz) mean of the per-sum ratios, weights
    design_weight (the sums were drawn in proportion to predicted magic)"""
    w = arr("dw", sel)
    O, E = arr(okey, sel), arr(ekey, sel)
    ok = E > 0
    return float((w[ok] * O[ok] / E[ok]).sum() / w[ok].sum())


def pearson_phi(sel, okey, ekey, vkey=None):
    O, E = arr(okey, sel), arr(ekey, sel)
    R = O.sum() / E.sum()
    mu = R * E
    ok = mu > 0
    if ok.sum() < 3:
        return None
    return float(((O[ok] - mu[ok]) ** 2 / mu[ok]).sum() / (ok.sum() - 1))


def between_sd(sel, okey, ekey):
    """between-sum sd of ln(O/E) beyond Poisson (method of moments on the
    sums with E >= 5)"""
    O, E = arr(okey, sel), arr(ekey, sel)
    ok = (E >= 5) & (O > 0)
    if ok.sum() < 5:
        return None
    lr = np.log(O[ok] / E[ok])
    v = lr.var(ddof=1) - np.mean(1.0 / E[ok])
    return float(math.sqrt(max(v, 0.0)))


def summary(sel, okey, ekey, between=None, phi_pre=1.0, count_key=None):
    O = float(np.nansum(arr(okey, sel)))
    E = float(np.nansum(arr(ekey, sel)))
    out = {"n": len(sel), "O": O, "E": E, "ratio": O / E if E > 0 else None}
    if E <= 0:
        return out
    # Poisson: on the effective count (d-first estimates carry the sampling variance)
    cnt = O if count_key is None else float(np.nansum(arr(count_key, sel)))
    out["poisson90"] = garwood(cnt, E * cnt / O) if (O > 0 and cnt != O) else garwood(O, E)
    out["p_ge"], out["p_le"] = CAL.poisson_tail(int(round(cnt)), E * (cnt / O if O > 0 else 1.0))
    if len(sel) >= 8:
        out["boot90"] = boot_ratio(sel, okey, ekey)
        out["phi"] = pearson_phi(sel, okey, ekey)
        if out["phi"] and O > 0:
            out["quasi90"] = lognorm_ci(O / E, math.sqrt(max(out["phi"], 1.0) / O))
    if between is not None and O > 0:
        sd = math.sqrt(phi_pre / max(O, 1e-9) + between ** 2 / len(sel))
        out["prereg90"] = lognorm_ci(O / E, sd)
    out["hh"] = hh(sel, okey, ekey)
    return out


def effective(sel, okey, vkey):
    """effective Poisson count of a total with per-sum variances (plain: O;
    d-first est_squares: se^2): O^2 / sum var"""
    O = arr(okey, sel)
    V = arr(vkey, sel)
    return float(O.sum() ** 2 / V.sum()) if V.sum() > 0 else 0.0


res = {"hash_ok": True, "bands": {}, "strata": {}, "per_run": rows}
groups = {b: [o for o in rows if o["band"] == b] for b in BANDS}
groups["all"] = rows
groups["3-12k"] = groups["3-6k"] + groups["6-12k"]
groups["12-45k"] = groups["12-24k"] + groups["24-45k"]
groups[">=6k"] = groups["6-12k"] + groups["12-45k"]
groups["6-12k plain"] = [o for o in groups["6-12k"] if o["mode"] == "plain"]
groups["6-12k dfirst"] = [o for o in groups["6-12k"] if o["mode"] == "dfirst"]
groups["3-12k plain"] = [o for o in groups["3-12k"] if o["mode"] == "plain"]

for o in rows:
    o["sq_eff"] = o["sq_obs"] ** 2 / o["sq_var"] if o["sq_var"] > 0 else 0.0
    o["E_sq_eff"] = o["E_sq"] * (o["sq_eff"] / o["sq_obs"]) if o["sq_obs"] > 0 else o["E_sq"] / o["w_rec"]

out = {}
for g, sel in groups.items():
    if not sel:
        continue
    d = {}
    d["squares"] = summary(sel, "sq_obs", "E_sq", between=0.26, phi_pre=2.5, count_key="sq_eff")
    # Poisson interval on the effective count (d-first streams: est^2 / se^2)
    eff = sum(o["sq_eff"] for o in sel)
    Esum = sum(o["E_sq"] for o in sel)
    Osum = sum(o["sq_obs"] for o in sel)
    d["squares"]["eff_count"] = eff
    d["squares"]["poisson90"] = [x for x in garwood(eff, Esum * eff / Osum)] if eff > 0 else None
    d["squares"]["between_sd"] = between_sd([o for o in sel if o["mode"] == "plain"], "sq_obs", "E_sq")
    d["squares_recorded"] = summary(sel, "sq_rec", "E_sq_rec")
    d["S_per_square"] = summary(sel, "S_obs", "E_S_cond", between=0.18, count_key="S_raw")
    d["P_per_square"] = summary(sel, "P_obs", "E_P_cond", between=0.18, count_key="P_raw")
    d["SP_found"] = summary(sel, "SP_obs", "E_SP_found", between=0.36, phi_pre=1.3)
    d["SP_per_square"] = summary(sel, "SP_obs", "E_SP_cond", between=0.36, phi_pre=1.3)
    d["cpu_run_anchored"] = {"O": float(arr("cpu", sel).sum()), "E": float(arr("cpu_pred", sel).sum())}
    d["cpu_run_anchored"]["ratio"] = d["cpu_run_anchored"]["O"] / d["cpu_run_anchored"]["E"]
    d["spS"] = {"O": float(arr("spS_sq", sel).sum()), "E": float(arr("E_spS", sel).sum())}
    d["spP"] = {"O": float(arr("spP_sq", sel).sum()), "E": float(arr("E_spP", sel).sum())}
    d["magic"] = {"O": float(arr("magic", sel).sum()), "E": float(arr("E_magic", sel).sum())}
    d["sp_stream"] = {"O": float(arr("sp_stream", sel).sum()), "E": float(np.nansum(arr("E_SP_stream", sel)))}
    out[g] = d
res["bands"] = out

# strata band x class x k (squares, S, P per square)
strata = collections.defaultdict(list)
for o in rows:
    strata[(o["band"], o["class"], o["kbin"])].append(o)
st = {}
for key, sel in sorted(strata.items(), key=lambda kv: (BANDS.index(kv[0][0]), kv[0][1], kv[0][2])):
    st[" / ".join(key)] = {"n": len(sel), "squares": summary(sel, "sq_obs", "E_sq", between=0.26, phi_pre=2.5, count_key="sq_eff"),
                           "S": summary(sel, "S_obs", "E_S_cond", count_key="S_raw"), "P": summary(sel, "P_obs", "E_P_cond", count_key="P_raw"),
                           "SP": {"O": float(arr("SP_obs", sel).sum()), "E": float(arr("E_SP_found", sel).sum())},
                           "E_share": float(arr("dw", sel).sum())}
res["strata"] = st
# marginal splits within 3-12k: class, k, x, tau
marg = {}
sel312 = groups["3-12k"]
for name, f in (("class", lambda o: o["class"]), ("k", lambda o: o["kbin"]),
                ("S/S_min", lambda o: "<1.1" if o["S_over_smin"] < 1.1 else ("1.1-1.2" if o["S_over_smin"] < 1.2
                                                                            else ("1.2-1.4" if o["S_over_smin"] < 1.4 else ">=1.4"))),
                ("labels", lambda o: "<=128" if o["labels"] <= 128 else ("129-256" if o["labels"] <= 256 else ">256")),
                ("N actual", lambda o: "<4k" if o["N"] < 4000 else ("4-6k" if o["N"] < 6000 else ("6-9k" if o["N"] < 9000 else ">=9k")))):
    dd = collections.defaultdict(list)
    for o in sel312:
        dd[f(o)].append(o)
    marg[name] = {kk: {"n": len(v), "squares": summary(v, "sq_obs", "E_sq", count_key="sq_eff"), "S": summary(v, "S_obs", "E_S_cond", count_key="S_raw"),
                       "P": summary(v, "P_obs", "E_P_cond", count_key="P_raw"),
                       "SP": summary(v, "SP_obs", "E_SP_cond")}
                  for kk, v in sorted(dd.items())}
res["marginals_3_12k"] = marg

# ---------------------------------------------------------------------------
# the SP rate: pooled over bands with the 10 CPU-year E shares, posterior
sp = {}
for g in ("3-6k", "6-12k"):
    sp[g] = {"O": out[g]["SP_per_square"]["O"], "E": out[g]["SP_per_square"]["E"]}
w = E10
num = sum(w[g] * sp[g]["O"] / sp[g]["E"] for g in w) / sum(w.values())
sp["E10_weighted"] = num
# Jeffreys gamma posterior for f_SP (O, E pooled 3-12k + 12-45k): Gamma(O + 1/2, E)
O_all = out["all"]["SP_per_square"]["O"]
E_all = out["all"]["SP_per_square"]["E"]
qs = [0.05, 0.16, 0.5, 0.84, 0.95]
sp["posterior_flat"] = {"O": O_all, "E": E_all,
                        "q": {str(q): CAL.gamma_quantile(O_all + 0.5, q) / E_all for q in qs}}
# found (pre-registered expectation)
sp["found"] = {"O": out["all"]["SP_found"]["O"], "E": out["all"]["SP_found"]["E"]}
res["sp"] = sp

# ---------------------------------------------------------------------------
# time: per law
tm = {}
plain = [o for o in rows if o["mode"] == "plain"]
dfr = [o for o in rows if o["mode"] == "dfirst"]


def gm(v):
    v = np.asarray(v, float)
    return float(np.exp(np.mean(np.log(v))))


def tsum(sel, okey, ekey):
    O, E = arr(okey, sel), arr(ekey, sel)
    lr = np.log(O / E)
    return {"n": len(sel), "sum_ratio": float(O.sum() / E.sum()), "gmean": float(np.exp(lr.mean())),
            "sd_ln": float(lr.std(ddof=1)) if len(sel) > 1 else None}


for b in BANDS + ["3-12k"]:
    sel = [o for o in plain if (o["band"] == b if b != "3-12k" else o["band"] in ("3-6k", "6-12k"))]
    if sel:
        tm[f"plain {b}"] = {law: tsum(sel, "t_plain_whole", law) for law in ("tp", "tp_ratio", "tp_anch")}
# plain by N' sub-bands (actual model N')
for lo, hi in ((3000, 4000), (4000, 5000), (5000, 6000), (6000, 8000), (8000, 12000)):
    sel = [o for o in plain if lo <= o["Np"] < hi]
    if sel:
        tm[f"plain N' {lo // 1000}-{hi // 1000}k"] = {law: tsum(sel, "t_plain_whole", law)
                                                    for law in ("tp", "tp_ratio", "tp_anch")}
for b in ("6-12k", "12-24k", "24-45k", ">=6k"):
    sel = [o for o in dfr if (o["band"] == b if b != ">=6k" else True)]
    if sel:
        tm[f"dfirst {b}: d-first whole sum"] = {law: tsum(sel, "t_dfirst_whole", law) for law in ("td", "td_ratio")}
        tm[f"dfirst {b}: plain whole sum (stream est_time)"] = {law: tsum(sel, "t_plain_whole", law)
                                                               for law in ("tp", "tp_ratio", "tp_anch")}
        for o in sel:
            o["ratio_meas"] = o["t_dfirst_whole"] / o["t_plain_whole"]
        tm[f"dfirst {b}: ratio d-first / plain"] = tsum(sel, "ratio_meas", "r")
        tm[f"dfirst {b}: d loop run"] = tsum(sel, "cpu_dloop", "pred_dloop")
        tm[f"dfirst {b}: stream run"] = tsum(sel, "cpu_stream", "pred_stream")
for b in BANDS + ["all"]:
    sel = rows if b == "all" else [o for o in rows if o["band"] == b]
    tm[f"run CPU {b} (anchored, as planned)"] = tsum(sel, "cpu", "cpu_pred")
res["time"] = tm
res["N_over_Npred"] = {b: float(np.median([o["N"] / o["Np"] for o in groups[b]])) for b in BANDS}
print("N / N' median by band", res["N_over_Npred"])

# ---------------------------------------------------------------------------
# selection effect: rank by predicted magic density (magic per anchored CPU-s)
# and by predicted magic per sum, within 3-12k (plain + d-first), quartiles
sel_out = {}
for rank_key in ("dens", "magic_sum", "pmagic"):
    sel = groups["3-12k"]
    v = arr(rank_key, sel)
    order = np.argsort(-v)
    qn = np.array_split(order, 4)
    rr = []
    for qi, ix in enumerate(qn):
        ss = [sel[i] for i in ix]
        sq = summary(ss, "sq_obs", "E_sq", count_key="sq_eff")
        S_ = summary(ss, "S_obs", "E_S_cond", count_key="S_raw")
        P_ = summary(ss, "P_obs", "E_P_cond", count_key="P_raw")
        SP_ = summary(ss, "SP_obs", "E_SP_cond")
        m = (S_["ratio"] * P_["ratio"]) ** 2
        rr.append({"quartile": qi + 1, "n": len(ss), "range": [float(v[ix].min()), float(v[ix].max())],
                   "squares": sq, "S": S_, "P": P_, "SP": SP_, "rSrP2": m,
                   "E_magic_share": float(arr("magic_sum", ss).sum() / arr("magic_sum", sel).sum())})
    # regression of ln(O/E squares) on ln(rank value) over plain sums with E >= 5
    ss = [o for o in sel if o["E_sq"] >= 5 and o["sq_obs"] > 0]
    x = np.log(arr(rank_key, ss))
    y = np.log(arr("sq_obs", ss) / arr("E_sq", ss))
    wts = 1.0 / (1.0 / arr("E_sq", ss) + 0.26 ** 2)
    X = np.vstack([np.ones_like(x), x - np.average(x, weights=wts)]).T
    Wm = np.diag(wts)
    beta = np.linalg.solve(X.T @ Wm @ X, X.T @ Wm @ y)
    resid = y - X @ beta
    s2 = float((wts * resid ** 2).sum() / (len(y) - 2))
    cov = np.linalg.inv(X.T @ Wm @ X) * s2
    sel_out[rank_key] = {"quartiles": rr, "slope_ln_sq_ratio": float(beta[1]), "slope_se": float(math.sqrt(cov[1, 1])),
                         "n": len(ss)}
# the same for the per-square traversal residuals (S, P): per-sum O/E vs rank
res["selection"] = sel_out

# ---------------------------------------------------------------------------
with open(os.path.join(CT, "analysis", "results.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
cols = ["id", "band", "class", "kbin", "mode", "Np", "N", "labels", "S_over_smin", "sq_obs", "E_sq", "sq_rec",
        "E_sq_rec", "S_obs", "E_S_cond", "P_obs", "E_P_cond", "SP_obs", "E_SP_found", "E_SP_cond", "sp_stream",
        "cpu", "cpu_pred", "t_plain_whole", "tp", "tp_ratio", "tp_anch", "t_dfirst_whole", "td", "r",
        "magic_sum", "dens", "dw"]
with open(os.path.join(CT, "analysis", "per_run.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for o in rows:
        f.write("\t".join(fmt(o.get(c)) if isinstance(o.get(c), float) else str(o.get(c)) for c in cols) + "\n")

# ---------------------------------------------------------------------------
# print
print("## Bands\n")
print("| group | sums | squares O / E = ratio | Poisson90 (eff.) | boot90 (sums) | prereg90 | phi | between sd | HH |")
print("|---|---:|---|---|---|---|---:|---:|---:|")
for g in ["3-6k", "6-12k", "6-12k plain", "6-12k dfirst", "12-24k", "24-45k", "3-12k", "3-12k plain", "12-45k", ">=6k", "all"]:
    s = out[g]["squares"]
    print(f"| {g} | {s['n']} | {fmt(s['O'])} / {fmt(s['E'])} = {fmt(s['ratio'])} | {fmt(s.get('poisson90'))} | "
          f"{fmt(s.get('boot90'))} | {fmt(s.get('prereg90'))} | {fmt(s.get('phi'))} | {fmt(s.get('between_sd'))} | "
          f"{fmt(s.get('hh'))} |")
for q in ("S_per_square", "P_per_square", "SP_found", "SP_per_square"):
    print(f"\n### {q}\n")
    print("| group | O / E = ratio | Poisson90 | boot90 | prereg90 | phi | HH | p(>=O) | p(<=O) |")
    print("|---|---|---|---|---|---:|---:|---:|---:|")
    for g in ["3-6k", "6-12k", "6-12k plain", "6-12k dfirst", "12-24k", "24-45k", "3-12k", "3-12k plain", "12-45k", ">=6k", "all"]:
        s = out[g][q]
        print(f"| {g} | {fmt(s['O'])} / {fmt(s['E'])} = {fmt(s['ratio'])} | {fmt(s.get('poisson90'))} | "
              f"{fmt(s.get('boot90'))} | {fmt(s.get('prereg90'))} | {fmt(s.get('phi'))} | {fmt(s.get('hh'))} | "
              f"{fmt(s.get('p_ge'))} | {fmt(s.get('p_le'))} |")
print("\n### rungs and other counts\n")
for g in ["3-6k", "6-12k", "12-24k", "24-45k", "all"]:
    d = out[g]
    print(f"{g}: SP+S {d['spS']}, SP+P {d['spP']}, magic {d['magic']}, SP in streams {d['sp_stream']}, "
          f"CPU {d['cpu_run_anchored']}")
print("\n## SP\n")
print(json.dumps(sp, indent=1, default=float))
print("\n## Strata\n")
print("| stratum | sums | E share in band | squares O/E [boot or Poisson] | S/sq | P/sq | SP O/E |")
print("|---|---:|---:|---|---|---|---|")
for k, v in st.items():
    s = v["squares"]
    print(f"| {k} | {v['n']} | {fmt(v['E_share'])} | {fmt(s['ratio'])} {fmt(s.get('boot90') or s.get('poisson90'))} | "
          f"{fmt(v['S']['ratio'])} | {fmt(v['P']['ratio'])} | {fmt(v['SP']['O'])}/{fmt(v['SP']['E'])} |")
print("\n## Marginals 3-12k\n")
for name, m in marg.items():
    print(f"\n{name}:")
    for kk, v in m.items():
        print(f"  {kk}: n {v['n']}, squares {fmt(v['squares']['ratio'])} {fmt(v['squares'].get('boot90'))}, "
              f"S {fmt(v['S']['ratio'])} {fmt(v['S'].get('poisson90'))}, P {fmt(v['P']['ratio'])} "
              f"{fmt(v['P'].get('poisson90'))}, SP {fmt(v['SP']['O'])}/{fmt(v['SP']['E'])}")
print("\n## Time\n")
for k, v in tm.items():
    print(k, json.dumps(v, default=lambda x: round(x, 3)))
print("\n## Selection\n")
for k, v in sel_out.items():
    print(f"rank by {k}: slope of ln(O/E squares) on ln(rank) {v['slope_ln_sq_ratio']:.3f} +- {v['slope_se']:.3f} "
          f"(n {v['n']})")
    for q in v["quartiles"]:
        print(f"  Q{q['quartile']} n {q['n']} E-share {q['E_magic_share']:.2f}: squares {fmt(q['squares']['ratio'])} "
              f"{fmt(q['squares'].get('boot90'))}, S {fmt(q['S']['ratio'])}, P {fmt(q['P']['ratio'])} "
              f"{fmt(q['P'].get('poisson90'))}, (rSrP)^2 {fmt(q['rSrP2'])}, SP {fmt(q['SP']['O'])}/{fmt(q['SP']['E'])}")
