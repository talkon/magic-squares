#!/usr/bin/env python3
"""Pre-registered analysis of the calibration-target runs: observed vs the
frozen predictions.json, per run and pooled per stratum (band, band x class
x k), with Poisson, bootstrap-over-sums, overdispersed and design-weighted
(Hansen-Hurwitz) estimates.

Totals that mix plain sums with calibration-stream estimates (est_squares,
stride x traversals) are weighted: a stream's events count with weight k
(its r1 stride) but carry the uncertainty of one event each. Every interval
on such a total uses the effective count O / s against E / s, with the
event scale s = sum_i w_i E_i / sum_i E_i (w_i = 1 for plain sums, k for a
stream: Var(sum w o) = sum w^2 mu = s x E under the model), and the
overdispersion phi is computed on the raw counts. The between-sum term of
the overdispersed intervals uses the E-weighted effective number of sums
(sum E)^2 / sum E^2.

usage: CALIB_TARGET_DIR=... python3 analyze.py     (writes $CALIB_TARGET_DIR/analysis/
       results.json, per_run.tsv; prints the tables; run it with nice -n 19)
"""
import collections
import hashlib
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ctpaths  # noqa: E402

CT, OUT = ctpaths.CT, ctpaths.OUT
sys.path.insert(0, ctpaths.SCRIPTS)
import calibrate as CAL  # noqa: E402  (poisson_ci, poisson_tail, gamma_quantile, observe, d_rooted)

SHA = "c1e673b595d7de9452a5c2c01113967c4c6579220c885a3d25e15ed270b06c08"
BANDS = ["3-6k", "6-12k", "12-24k", "24-45k"]
E10 = {"3-6k": 0.54, "6-12k": 0.26}          # share of the 10 CPU-year E (plan: pooling weights)
B = 4000
BETWEEN_SQ_MEASURED = 0.41                   # between-sum sd of squares obs/pred (ratio scale, plain sums)
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


def norm_cdf(z):
    return 0.5 * math.erfc(-z / math.sqrt(2))


def fmt(x, d=3):
    if x is None:
        return "-"
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(fmt(v, d) for v in x) + "]"
    if isinstance(x, (int, np.integer)):
        return f"{x:,}" if abs(x) >= 1000 else str(x)
    if x == 0:
        return "0"
    if not math.isfinite(x):
        return str(x)
    if abs(x) >= 1000:
        return f"{x:,.0f}"
    return f"{x:.{d}g}"


# ---------------------------------------------------------------------------
# the frozen files
pred_path = os.path.join(CT, "predictions.json")
assert sha256(pred_path) == SHA, "predictions.json changed"
PRED = json.load(open(pred_path))
pre = {}
for line in open(os.path.join(CT, "PREREGISTERED.txt")):
    if "sha256" in line and ":" in line:
        name, h = line.split("sha256")[0].strip(), line.split(":")[-1].strip().split()[0]
        pre[name] = h
checked = []
for name, rel in (("plan.json", "plan.json"), ("run.sh", "run.sh"), ("runner.py", "runner.py"),
                  ("bin/msearch", "bin/msearch")):
    path = os.path.join(CT, rel)
    if not os.path.exists(path) and name == "bin/msearch":
        continue            # the binary is not archived (built from e0ace54)
    assert sha256(path) == pre[name], name
    checked.append(name)
print("verified sha256:", ", ".join(["predictions.json"] + checked))

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
    R = load(ctpaths.run_path(p["out"]))
    o = {"id": p["id"], "band": p["band"], "class": p["class"], "kbin": p["kbin"], "mode": p["mode"],
         "k": p["k"], "tau": p["tau"], "x": p["x"], "S_over_smin": p["S_over_smin"], "Np": p["N_pred"],
         "dw": p["design_weight"], "P_key": p["P_key"], "S": p["S"]}
    assert len(R["done"]) == 1
    if p["mode"] == "plain":
        s = R["sum"][0]
        sq = R["square"]
        assert s["squares"] == len(sq) and not s["truncated"]
        # the rungs above SP, counted from the grids (pre-registered: the partners of each SP
        # diagonal): unordered {SP, S-only} / {SP, P-only} / {SP, SP} partner pairs
        if sq:
            ev, best, bad = CAL.observe(sq)
            assert not bad and int(ev["SP"].sum()) == sum(q["sp_count"] for q in sq)
            spS, spP = (float(ev[k_].sum()) for k_ in ("SP+S", "SP+P"))
            nmagic = int((ev["SP+SP"] > 0).sum())
        else:
            spS = spP = 0.0
            nmagic = 0
        o.update(N=s["nvecs_raw"], labels=s["labels"], cpu=s["cpu"], sq_obs=len(sq), sq_rec=len(sq),
                 w_rec=1.0,
                 S_obs=sum(q["s_count"] for q in sq), P_obs=sum(q["p_count"] for q in sq),
                 S_raw=sum(q["s_count"] for q in sq), P_raw=sum(q["p_count"] for q in sq),
                 SP_obs=sum(q["sp_count"] for q in sq),
                 spS_obs=spS, spP_obs=spP, magic=nmagic,
                 sp_stream=0, cpu_dloop=0.0, cpu_stream=0.0, d_stride=1,
                 t_plain_whole=s["cpu"], t_dfirst_whole=None)
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
        # d-first pairs: each record judged on its own SP diagonal; magic squares by unique hash
        dsq = R.get("dsquare", [])
        dr = CAL.d_rooted(dsq) if dsq else []
        magic_h = {r["hash"] for r, x in zip(dsq, dr) if x["level"] == 14 or r.get("magic") or r.get("partner")}
        o.update(N=ds["nvecs_raw"], labels=ds["labels"], cpu=ds["cpu"] + cs["cpu"],
                 sq_obs=cs["est_squares"], sq_rec=len(csq), w_rec=float(k),
                 S_obs=k * cs["s_trav"], P_obs=k * cs["p_trav"], SP_obs=ds["pairs"],
                 S_raw=cs["s_trav"], P_raw=cs["p_trav"],
                 spS_obs=float(sum(x["nS"] for x in dr)), spP_obs=float(sum(x["nP"] for x in dr)),
                 magic=len(magic_h), sp_stream=cs["sp_pairs"], cpu_dloop=ds["cpu"], cpu_stream=cs["cpu"],
                 d_stride=kd, r1_stride=k, nd=ds["nd"],
                 t_dfirst_whole=ds["cpu"] - ds["time"] + ds["est_time"],
                 t_plain_whole=cs["est_time"] + cs["enum_time"] + cs["reduce_time"],
                 se_tplain=cs["se_time"])
        assert sum(1 for d in range(ds["d_offset"], ds["nvecs_raw"], kd)) == ds["nd"]
        assert ds["pairs"] == len(dsq)
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
    o["w_one"] = 1.0
    rows.append(o)

ID = {o["id"]: o for o in rows}


def arr(key, sel):
    return np.array([o[key] if o.get(key) is not None else np.nan for o in sel], float)


def boot_M(sel, seed=0):
    """(B, n) multiplicities of a bootstrap over the sums of sel, stratified by band"""
    r_ = np.random.default_rng(20261009 + seed)
    bands = np.array([o["band"] for o in sel])
    M = np.zeros((B, len(sel)))
    for b in BANDS:
        ix = np.nonzero(bands == b)[0]
        if not len(ix):
            continue
        pick = ix[r_.integers(0, len(ix), size=(B, len(ix)))]
        for j in range(B):
            M[j] += np.bincount(pick[j], minlength=len(sel))
    return M


_MCACHE = {}


def boot_ratio(sel, okey, ekey):
    """percentile 90% interval of sum O / sum E over a bootstrap of the sums
    (stratified by band)"""
    key = tuple(o["id"] for o in sel)
    if key not in _MCACHE:
        _MCACHE[key] = boot_M(sel)
    M = _MCACHE[key]
    O, E = np.nan_to_num(arr(okey, sel)), np.nan_to_num(arr(ekey, sel))
    with np.errstate(divide="ignore", invalid="ignore"):
        res = (M @ O) / (M @ E)
    res = res[np.isfinite(res)]
    return [float(np.quantile(res, 0.05)), float(np.quantile(res, 0.95))]


def hh(sel, okey, ekey):
    """design-weighted (Hansen-Hurwitz) mean of the per-sum ratios, weights
    design_weight (the sums were drawn in proportion to predicted magic)"""
    w = arr("dw", sel)
    O, E = arr(okey, sel), arr(ekey, sel)
    ok = E > 0
    return float((w[ok] * O[ok] / E[ok]).sum() / w[ok].sum())


def pearson_phi(sel, okey, ekey, wkey):
    """Pearson dispersion on the raw counts (O_i / w_i against R E_i / w_i)"""
    O, E, w = arr(okey, sel), arr(ekey, sel), arr(wkey, sel)
    R = O.sum() / E.sum()
    mu = R * E
    ok = mu > 0
    if ok.sum() < 3:
        return None
    return float(((O[ok] - mu[ok]) ** 2 / (w[ok] * mu[ok])).sum() / (ok.sum() - 1))


def between_sd(sel, okey, ekey):
    """between-sum sd of ln(O/E) beyond Poisson (method of moments on the
    plain sums with E >= 5)"""
    O, E = arr(okey, sel), arr(ekey, sel)
    ok = (E >= 5) & (O > 0)
    if ok.sum() < 5:
        return None
    lr = np.log(O[ok] / E[ok])
    v = lr.var(ddof=1) - np.mean(1.0 / E[ok])
    return float(math.sqrt(max(v, 0.0)))


def summary(sel, okey, ekey, wkey="w_one", between=None, phi_pre=1.0, between_meas=None):
    O = float(np.nansum(arr(okey, sel)))
    Ei = np.nan_to_num(arr(ekey, sel))
    E = float(Ei.sum())
    out = {"n": len(sel), "O": O, "E": E, "ratio": O / E if E > 0 else None}
    if E <= 0:
        return out
    w = arr(wkey, sel)
    s = float((w * Ei).sum() / E)                   # event scale (1 without stream estimates)
    cnt, Es = O / s, E / s
    nsum = E ** 2 / float((Ei ** 2).sum())          # E-weighted effective number of sums
    out.update(event_scale=s, eff_count=cnt, eff_sums=nsum)
    out["poisson90"] = garwood(cnt, Es)
    out["p_ge"], out["p_le"] = CAL.poisson_tail(int(round(cnt)), Es)
    if len(sel) >= 8:
        out["boot90"] = boot_ratio(sel, okey, ekey)
        out["phi"] = pearson_phi(sel, okey, ekey, wkey)
        if out["phi"] and O > 0:
            out["quasi90"] = lognorm_ci(O / E, math.sqrt(max(out["phi"], 1.0) * s / O))
    if between is not None and O > 0:
        out["prereg90"] = lognorm_ci(O / E, math.sqrt(phi_pre * s / O + between ** 2 / nsum))
    if between_meas is not None and O > 0:
        out["od90"] = lognorm_ci(O / E, math.sqrt(s / O + between_meas ** 2 / nsum))
    elif between_meas is not None:
        out["od90"] = [0.0, garwood(0, Es)[1] * math.exp(1.645 * between_meas / math.sqrt(nsum))]
    out["hh"] = hh(sel, okey, ekey)
    return out


def sq_summary(sel):
    d = summary(sel, "sq_obs", "E_sq", wkey="w_rec", between=0.26, phi_pre=2.5, between_meas=BETWEEN_SQ_MEASURED)
    d["between_sd"] = between_sd([o for o in sel if o["mode"] == "plain"], "sq_obs", "E_sq")
    return d


def tr_summary(sel, q):
    return summary(sel, f"{q}_obs", f"E_{q}_cond", wkey="w_rec", between=0.18)


res = {"hash_ok": True, "bands": {}, "strata": {}, "per_run": rows}
groups = {b: [o for o in rows if o["band"] == b] for b in BANDS}
groups["3-12k"] = groups["3-6k"] + groups["6-12k"]
groups["12-45k"] = groups["12-24k"] + groups["24-45k"]
groups["6-12k plain"] = [o for o in groups["6-12k"] if o["mode"] == "plain"]
groups["6-12k dfirst"] = [o for o in groups["6-12k"] if o["mode"] == "dfirst"]
groups["3-12k plain"] = [o for o in groups["3-12k"] if o["mode"] == "plain"]
# "all" pools differently defined populations (plan sums at 3-12k, pool sums S <= 2 S0 above 12k):
# kept for the raw totals only, not an estimate of a defined quantity
groups["all (raw totals only)"] = rows
GROUP_ORDER = ["3-6k", "6-12k", "6-12k plain", "6-12k dfirst", "3-12k", "3-12k plain", "12-24k", "24-45k", "12-45k"]

out = {}
for g, sel in groups.items():
    d = {}
    d["squares"] = sq_summary(sel)
    d["squares_recorded"] = summary(sel, "sq_rec", "E_sq_rec")
    d["S_per_square"] = tr_summary(sel, "S")
    d["P_per_square"] = tr_summary(sel, "P")
    d["SP_found"] = summary(sel, "SP_obs", "E_SP_found", between=0.36, phi_pre=1.3)
    d["SP_per_square"] = summary(sel, "SP_obs", "E_SP_cond", between=0.36, phi_pre=1.3)
    d["cpu_run_anchored"] = {"O": float(arr("cpu", sel).sum()), "E": float(arr("cpu_pred", sel).sum())}
    d["cpu_run_anchored"]["ratio"] = d["cpu_run_anchored"]["O"] / d["cpu_run_anchored"]["E"]
    d["spS"] = {"O": float(arr("spS_obs", sel).sum()), "E": float(arr("E_spS", sel).sum())}
    d["spP"] = {"O": float(arr("spP_obs", sel).sum()), "E": float(arr("E_spP", sel).sum())}
    d["magic"] = {"O": float(arr("magic", sel).sum()), "E": float(arr("E_magic", sel).sum())}
    d["sp_stream"] = {"O": float(arr("sp_stream", sel).sum()), "E": float(np.nansum(arr("E_SP_stream", sel)))}
    out[g] = d
res["bands"] = out

# the S deficit at 3-6k (and 3-12k): Poisson p and a dispersion-adjusted p (phi of the per-sum
# S counts; normal approximation on ln(O/E))
sdef = {}
for g in ("3-6k", "3-12k"):
    s_ = out[g]["S_per_square"]
    phi = max(s_.get("phi") or 1.0, 1.0)
    z = math.log(s_["ratio"]) / math.sqrt(phi * s_["event_scale"] / s_["O"])
    sdef[g] = {"ratio": s_["ratio"], "p_le_poisson": s_["p_le"], "phi": phi, "z": z, "p_le_overdispersed": norm_cdf(z)}
res["S_deficit"] = sdef

# strata band x class x k (squares, S, P per square)
strata = collections.defaultdict(list)
for o in rows:
    strata[(o["band"], o["class"], o["kbin"])].append(o)
st = {}
for key, sel in sorted(strata.items(), key=lambda kv: (BANDS.index(kv[0][0]), kv[0][1], kv[0][2])):
    st[" / ".join(key)] = {"n": len(sel), "squares": sq_summary(sel), "S": tr_summary(sel, "S"),
                           "P": tr_summary(sel, "P"),
                           "SP": {"O": float(arr("SP_obs", sel).sum()), "E": float(arr("E_SP_found", sel).sum())},
                           "E_share": float(arr("dw", sel).sum())}
res["strata"] = st

# marginal splits within 3-12k, each against the overall 3-12k rate: bootstrap (over sums,
# stratified by band) of the subgroup's ratio over the overall ratio, two-sided p
sel312 = groups["3-12k"]
M312 = boot_M(sel312, seed=1)
SPLITS = (("class", lambda o: o["class"]), ("k", lambda o: o["kbin"]),
          ("S/S_min", lambda o: "<1.1" if o["S_over_smin"] < 1.1 else ("1.1-1.2" if o["S_over_smin"] < 1.2
                                                                       else ("1.2-1.4" if o["S_over_smin"] < 1.4
                                                                             else ">=1.4"))),
          ("labels", lambda o: "<=128" if o["labels"] <= 128 else ("129-256" if o["labels"] <= 256 else ">256")),
          ("N actual", lambda o: "<4k" if o["N"] < 4000 else ("4-6k" if o["N"] < 6000 else
                                                               ("6-9k" if o["N"] < 9000 else ">=9k"))))
QTY = {"squares": ("sq_obs", "E_sq"), "S": ("S_obs", "E_S_cond"), "P": ("P_obs", "E_P_cond")}
marg = {}
ntests = 0
pvals = []
for name, f in SPLITS:
    labs = np.array([f(o) for o in sel312])
    marg[name] = {}
    for lab in sorted(set(labs.tolist())):
        m = labs == lab
        v = [o for o, mm in zip(sel312, m) if mm]
        e = {"n": len(v), "squares": sq_summary(v), "S": tr_summary(v, "S"), "P": tr_summary(v, "P"),
             "SP": summary(v, "SP_obs", "E_SP_cond"), "contrast": {}}
        for q, (ok_, ek_) in QTY.items():
            O, E = np.nan_to_num(arr(ok_, sel312)), np.nan_to_num(arr(ek_, sel312))
            with np.errstate(divide="ignore", invalid="ignore"):
                rb = ((M312 @ (O * m)) / (M312 @ (E * m))) / ((M312 @ O) / (M312 @ E))
            rb = rb[np.isfinite(rb)]
            r0 = (O[m].sum() / E[m].sum()) / (O.sum() / E.sum())
            pv = float(min(1.0, 2 * min((rb <= 1).mean(), (rb >= 1).mean())))
            e["contrast"][q] = {"rel": float(r0), "boot90": [float(np.quantile(rb, 0.05)), float(np.quantile(rb, 0.95))],
                                "p_two_sided": pv}
            ntests += 1
            pvals.append(pv)
        marg[name][lab] = e
res["marginals_3_12k"] = marg
res["marginal_tests"] = {"n_tests": ntests, "p_lt_0.05": int(sum(p < 0.05 for p in pvals)),
                         "p_lt_0.10": int(sum(p < 0.10 for p in pvals)),
                         "expected_lt_0.05_by_chance": 0.05 * ntests,
                         "holm_significant_0.05": 0}
ps = sorted(pvals)
for i, pv in enumerate(ps):
    if pv <= 0.05 / (ntests - i):
        res["marginal_tests"]["holm_significant_0.05"] += 1
    else:
        break

# ---------------------------------------------------------------------------
# the SP rate (the posterior of the SP coupling is in bands.py)
sp = {}
for g in ("3-6k", "6-12k"):
    sp[g] = {"O": out[g]["SP_per_square"]["O"], "E": out[g]["SP_per_square"]["E"]}
sp["E10_weighted"] = sum(E10[g] * sp[g]["O"] / sp[g]["E"] for g in E10) / sum(E10.values())
O_all = out["all (raw totals only)"]["SP_per_square"]["O"]
E_all = out["all (raw totals only)"]["SP_per_square"]["E"]
qs = [0.05, 0.16, 0.5, 0.84, 0.95]
sp["per_square_jeffreys"] = {"O": O_all, "E": E_all,
                             "q": {str(q): CAL.gamma_quantile(O_all + 0.5, q) / E_all for q in qs}}
sp["found"] = {"O": out["all (raw totals only)"]["SP_found"]["O"], "E": out["all (raw totals only)"]["SP_found"]["E"]}
Ef = sp["found"]["E"]
sp["found"].update({"p_ge": CAL.poisson_tail(7, Ef)[0], "p_ge_if_half": CAL.poisson_tail(7, Ef / 2)[0],
                    "p_le_if_double": CAL.poisson_tail(7, 2 * Ef)[1]})
sp["plain_only"] = summary([o for o in rows if o["mode"] == "plain"], "SP_obs", "E_SP_cond")
sp["runs_with_sp"] = [o["id"] for o in rows if o["SP_obs"] > 0]
res["sp"] = sp

# ---------------------------------------------------------------------------
# time: per law
tm = {}
plain = [o for o in rows if o["mode"] == "plain"]
dfr = [o for o in rows if o["mode"] == "dfirst"]


def tsum(sel, okey, ekey):
    O, E = arr(okey, sel), arr(ekey, sel)
    lr = np.log(O / E)
    return {"n": len(sel), "sum_ratio": float(O.sum() / E.sum()), "gmean": float(np.exp(lr.mean())),
            "sd_ln": float(lr.std(ddof=1)) if len(sel) > 1 else None}


for b in BANDS + ["3-12k"]:
    sel = [o for o in plain if (o["band"] == b if b != "3-12k" else o["band"] in ("3-6k", "6-12k"))]
    if sel:
        tm[f"plain {b}"] = {law: tsum(sel, "t_plain_whole", law) for law in ("tp", "tp_ratio", "tp_anch")}
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

# ---------------------------------------------------------------------------
# selection effect (winner's curse), within each band, against the stratum-adjusted model:
# each sum's E (squares, S, P) times its band x class x k stratum's O/E, so that what is left
# is the part of the residual that goes with the predicted density beyond what a recalibration
# by stratum removes. Quartiles of predicted magic per anchored CPU-second; the quantity is
# squares x (S P)^2 (P(magic) without the SP coupling), relative to the band.
for o in rows:
    sel = strata[(o["band"], o["class"], o["kbin"])]
    for q, (ok_, ek_) in QTY.items():
        Oq, Eq = np.nansum(arr(ok_, sel)), np.nansum(arr(ek_, sel))
        o[f"{ek_}_adj"] = o[ek_] * (Oq / Eq if Eq > 0 and Oq > 0 else 1.0)


def m_rel(sel, M=None):
    """squares x (S P)^2 against the stratum-adjusted expectations (and its bootstrap)"""
    def tot(key):
        v = np.nan_to_num(arr(key, sel))
        return v.sum() if M is None else M @ v
    r = tot("sq_obs") / tot("E_sq_adj")
    rs = tot("S_obs") / tot("E_S_cond_adj")
    rp = tot("P_obs") / tot("E_P_cond_adj")
    return r * (rs * rp) ** 2


sel_out = {}
for b in ("3-6k", "6-12k", "3-12k"):
    sel = groups[b]
    v = arr("dens", sel)
    if b == "3-12k":       # rank within band, then pool the band quartiles
        qi_of = {}
        for bb in ("3-6k", "6-12k"):
            ss = [o for o in sel if o["band"] == bb]
            vv = arr("dens", ss)
            for qi, ix in enumerate(np.array_split(np.argsort(-vv), 4)):
                for i in ix:
                    qi_of[ss[i]["id"]] = qi
        qn = [[i for i, o in enumerate(sel) if qi_of[o["id"]] == qi] for qi in range(4)]
    else:
        qn = [list(ix) for ix in np.array_split(np.argsort(-v), 4)]
    rr = []
    for qi, ix in enumerate(qn):
        ss = [sel[i] for i in ix]
        M = boot_M(ss, seed=10 + qi)
        mb = m_rel(ss, M)
        mb = mb[np.isfinite(mb)]
        rr.append({"quartile": qi + 1, "n": len(ss), "rel": float(m_rel(ss)),
                   "boot90": [float(np.quantile(mb, 0.05)), float(np.quantile(mb, 0.95))],
                   "sd_ln": float(np.log(mb[mb > 0]).std()),
                   "squares_rel": float(np.nansum(arr("sq_obs", ss)) / np.nansum(arr("E_sq_adj", ss))),
                   "SP": {"O": float(arr("SP_obs", ss).sum()), "E": float(arr("E_SP_cond", ss).sum())},
                   "E_magic_share": float(arr("magic_sum", ss).sum() / arr("magic_sum", sel).sum())})
    # regression of ln(O/E_adj squares) on ln(density) over the plain sums with E >= 5
    ss = [o for o in sel if o["E_sq"] >= 5 and o["sq_obs"] > 0 and o["mode"] == "plain"]
    x = np.log(arr("dens", ss))
    y = np.log(arr("sq_obs", ss) / arr("E_sq_adj", ss))
    wts = 1.0 / (1.0 / arr("E_sq", ss) + BETWEEN_SQ_MEASURED ** 2)
    X = np.vstack([np.ones_like(x), x - np.average(x, weights=wts)]).T
    Wm = np.diag(wts)
    beta = np.linalg.solve(X.T @ Wm @ X, X.T @ Wm @ y)
    resid = y - X @ beta
    s2 = float((wts * resid ** 2).sum() / (len(y) - 2))
    cov = np.linalg.inv(X.T @ Wm @ X) * max(s2, 1.0)
    sel_out[b] = {"quartiles": rr, "slope_ln_sq_ratio": float(beta[1]), "slope_se": float(math.sqrt(cov[1, 1])),
                  "n_slope": len(ss)}
res["selection"] = sel_out

# ---------------------------------------------------------------------------
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "results.json"), "w") as f:
    json.dump(res, f, indent=1, default=float)
cols = ["id", "band", "class", "kbin", "mode", "Np", "N", "labels", "S_over_smin", "sq_obs", "E_sq", "w_rec",
        "sq_rec", "E_sq_rec", "S_obs", "E_S_cond", "P_obs", "E_P_cond", "SP_obs", "E_SP_found", "E_SP_cond",
        "spS_obs", "E_spS", "spP_obs", "E_spP", "magic", "sp_stream", "cpu", "cpu_pred", "t_plain_whole", "tp",
        "tp_ratio", "tp_anch", "t_dfirst_whole", "td", "r", "magic_sum", "dens", "dw"]
with open(os.path.join(OUT, "per_run.tsv"), "w") as f:
    f.write("\t".join(cols) + "\n")
    for o in rows:
        f.write("\t".join(fmt(o.get(c)) if isinstance(o.get(c), float) else str(o.get(c)) for c in cols) + "\n")

# ---------------------------------------------------------------------------
# print
print("N / N' median by band", res["N_over_Npred"])
print("\n## Squares per sum (obs / pred)\n")
print("| group | sums | O / E = ratio | HH | eff. count | Poisson90 (eff.) | boot90 (sums) | od90 (between "
      f"{BETWEEN_SQ_MEASURED}) | prereg90 | phi (raw) | between sd |")
print("|---|---:|---|---:|---:|---|---|---|---|---:|---:|")
for g in GROUP_ORDER:
    s = out[g]["squares"]
    print(f"| {g} | {s['n']} | {fmt(s['O'])} / {fmt(s['E'])} = {fmt(s['ratio'])} | {fmt(s.get('hh'))} | "
          f"{fmt(s.get('eff_count'))} | {fmt(s.get('poisson90'))} | {fmt(s.get('boot90'))} | {fmt(s.get('od90'))} | "
          f"{fmt(s.get('prereg90'))} | {fmt(s.get('phi'))} | {fmt(s.get('between_sd'))} |")
for q in ("S_per_square", "P_per_square", "SP_found", "SP_per_square"):
    print(f"\n### {q}\n")
    print("| group | O / E = ratio | eff. count | Poisson90 | boot90 | quasi90 | prereg90 | phi (raw) | HH | "
          "p(>=O) | p(<=O) |")
    print("|---|---|---:|---|---|---|---|---:|---:|---:|---:|")
    for g in GROUP_ORDER:
        s = out[g][q]
        print(f"| {g} | {fmt(s['O'])} / {fmt(s['E'])} = {fmt(s['ratio'])} | {fmt(s.get('eff_count'))} | "
              f"{fmt(s.get('poisson90'))} | {fmt(s.get('boot90'))} | {fmt(s.get('quasi90'))} | "
              f"{fmt(s.get('prereg90'))} | {fmt(s.get('phi'))} | {fmt(s.get('hh'))} | "
              f"{fmt(s.get('p_ge'))} | {fmt(s.get('p_le'))} |")
print("\nS deficit:", json.dumps(sdef, default=float))
print("\n### rungs above SP (partner pairs from the grids) and other counts\n")
for g in ["3-6k", "6-12k", "12-24k", "24-45k", "all (raw totals only)"]:
    d = out[g]
    print(f"{g}: SP+S {d['spS']}, SP+P {d['spP']}, magic {d['magic']}, SP in streams {d['sp_stream']}, "
          f"CPU {d['cpu_run_anchored']}")
print("\n## SP\n")
print(json.dumps(sp, indent=1, default=float))
print("\n## Strata\n")
print("| stratum | sums | E share in band | squares O/E [boot or Poisson] | S/sq | P/sq | SP O/E |")
print("|---|---:|---:|---|---|---|---|")
for k_, v in st.items():
    s = v["squares"]
    print(f"| {k_} | {v['n']} | {fmt(v['E_share'])} | {fmt(s['ratio'])} {fmt(s.get('boot90') or s.get('poisson90'))} | "
          f"{fmt(v['S']['ratio'])} | {fmt(v['P']['ratio'])} | {fmt(v['SP']['O'])}/{fmt(v['SP']['E'])} |")
print(f"\n## Marginals 3-12k, against the overall 3-12k rate ({ntests} contrasts; "
      f"{res['marginal_tests']['p_lt_0.05']} with p < 0.05, {0.05 * ntests:.1f} expected by chance; "
      f"{res['marginal_tests']['holm_significant_0.05']} after Holm)\n")
print("| split | group | sums | squares (rel, p) | S (rel, p) | P (rel, p) | SP O/E |")
print("|---|---|---:|---|---|---|---|")
for name, m in marg.items():
    for kk, v in m.items():
        c_ = v["contrast"]
        print(f"| {name} | {kk} | {v['n']} | " + " | ".join(
            f"{fmt(v[q]['ratio'])} (rel {fmt(c_[q]['rel'])} {fmt(c_[q]['boot90'])}, p {fmt(c_[q]['p_two_sided'], 2)})"
            for q in ("squares", "S", "P")) + f" | {fmt(v['SP']['O'])}/{fmt(v['SP']['E'])} |")
print("\n## Time\n")
for k_, v in tm.items():
    print(k_, json.dumps(v, default=lambda x: round(x, 3)))
print("\n## Selection (within band, against the stratum-adjusted model; squares x (S P)^2)\n")
for k_, v in sel_out.items():
    print(f"{k_}: slope of ln(O/E_adj squares) on ln(density) {v['slope_ln_sq_ratio']:.3f} +- {v['slope_se']:.3f} "
          f"(n {v['n_slope']})")
    for q in v["quartiles"]:
        print(f"  Q{q['quartile']} n {q['n']} E-share {q['E_magic_share']:.2f}: rel {fmt(q['rel'])} "
              f"{fmt(q['boot90'])} (sd ln {q['sd_ln']:.3f}), squares rel {fmt(q['squares_rel'])}, "
              f"SP {fmt(q['SP']['O'])}/{fmt(q['SP']['E'])}")
