#!/usr/bin/env python3
"""
Stage 1's decision report (research/stage1.md): what the pairs found so far
say about the SP coupling f_rho, and what that means for spending 1, 3 or 10
CPU-years. It only reports; it never decides spend.

usage: decide.py STATE [--select stage1|all] [--plan PLAN --units-dir DIR]
                 [--prior shipped|pre|flat] [--band LO:HI]
                 [--ecurve FILE] [--table FILE] [--json OUT]
       decide.py --simulate PLAN [--reps 1000] [--pairs 45] [--phi 1.41]

STATE is a scheduler state (scripts/scheduler.py --state). The records used:
  --select stage1 (the default when the state has a stage-1 clock): the
      units `scheduler.py run --stage1` launched while stage 1 was on
      (launched_6.jsonl "stage1"), with the predictions frozen at launch;
  --select all: every record of the state (the calibration search's replay);
  --plan PLAN --units-dir DIR: a static stage-1 plan (scripts/laptop/
      make_plan.py --stage1) run by scripts/laptop/run.py into DIR
      (U<i>.jsonl = plan line i), with the plan's frozen predictions.

The fit (from worktree f1/p1-learn-the-sp-coup, scratchpad frho_block.py,
"as written", on the plain records only): per P its (square, SP traversal)
pairs Y_P in the cells of the N' band [LO, HI) (default the stage-1 band
3-6k) against the model's M_P = sum over its squares of k_SP r_S r_P 720
e^{lpSP} (scheduler.sp_kappa; r_S, r_P the state's learned class rates);
quasi-Poisson with a gamma prior: posterior Gamma(Y / phi + a, M / phi + a /
f0), phi = the Pearson dispersion of the per-P totals x G / (G - 1), floored
at 1; priors (f0, ln sd of f_rho): shipped = the calibration search's
posterior (1.06, 0.235: f_rho^2 ln sd 0.47; for records that are not the
calibration search's), pre = its prior (1, 0.29: to replay its records),
flat. Pairs elsewhere (other N' bands) are reported as a check only: f_rho
at small N is what the model's kappa was calibrated on, and carrying it to
N' 3-12k is the extrapolation term this stage measures.

The predictive E at C CPU-years: E_ship(C) (forecast --shipped, f_rho = 1;
--ecurve, research/stage1/ecurve.json) x f_rho^2 (posterior draws) x
lognormal(0, s(C)^2 + 0.2^2 + 0.12^2) (the class factors' posterior s(C):
0.115 at 1 CPU-year, 0.097 at 10, calibration-target.md 5; the pair factor;
selection), applied to all of E, as the calibration search's forecast did;
P(>=1) = mean over draws of 1 - exp(-E).
"""
import argparse
import collections
import glob
import gzip
import json
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

FRHO_PRIORS = {"shipped": (1.06, 0.47 / 2), "pre": (1.0, 0.29), "flat": (1.0, None)}
PAIR_SD, SELECTION_SD = 0.2, 0.12
CLASS_SD = ((1.0, 0.115), (10.0, 0.097))     # (CPU-years, ln sd), calibration-target.md 5
MARKS = (1.0, 3.0, 10.0)


# --------------------------------------------------------------------------
# the fit


def trigamma(x):
    """psi'(x), x > 0 (recurrence to x >= 6, then the asymptotic series)"""
    v = 0.0
    while x < 6:
        v += 1.0 / (x * x)
        x += 1
    x2 = 1.0 / (x * x)
    return v + 1 / x + x2 / 2 + (1 / x) * x2 * (1 / 6 - x2 * (1 / 30 - x2 * (1 / 42 - x2 / 30)))


def gamma_quantile(shape, rate, z):
    """the quantile at standard normal z of Gamma(shape, rate) (Wilson-
    Hilferty; shape >= ~1)"""
    if shape <= 0 or rate <= 0:
        return float("nan")
    c = 1.0 / (9.0 * shape)
    return shape * max(1.0 - c + z * math.sqrt(c), 0.0) ** 3 / rate


def fit_frho(Y, M, prior="shipped", phi_min=1.0):
    """f_rho from per-P pairs Y and model pairs M (sequences; the frho_block
    fit, as written): dict with the gamma posterior and summaries, or None
    (flat prior and no pairs)"""
    f0, sd0 = FRHO_PRIORS[prior]
    a = 0.0 if sd0 is None else 1.0 / sd0 ** 2 + 0.5
    Yt, Mt = float(sum(Y)), float(sum(M))
    if a == 0 and Yt <= 0:
        return None
    fm = Yt / Mt if Yt > 0 and Mt > 0 else f0
    G = sum(1 for m in M if m > 0)
    phi = max(1.0, phi_min)
    if G >= 2 and Yt > 0:
        phi = max(phi, sum((y - fm * m) ** 2 for y, m in zip(Y, M)) / (fm * Mt) * G / (G - 1))
    shape = Yt / phi + a
    rate = Mt / phi + (a / f0 if a else 0.0)
    out = {"prior": prior, "Y": Yt, "M": Mt, "nP": G, "phi": phi, "shape": shape, "rate": rate,
           "mean": shape / rate, "median": gamma_quantile(shape, rate, 0.0),
           "q05": gamma_quantile(shape, rate, -1.645), "q95": gamma_quantile(shape, rate, 1.645),
           "ln_sd": math.sqrt(trigamma(shape)) if shape > 0 else float("inf")}
    out["ln_sd_f2"] = 2 * out["ln_sd"]
    return out


def describe_fit(fr, what=""):
    if fr is None:
        return f"f_rho{what}: no pairs (flat prior: no estimate)"
    return (f"f_rho{what} = {fr['median']:.3f} median, {fr['mean']:.3f} mean (90% "
            f"{fr['q05']:.3f}-{fr['q95']:.3f}; ln sd {fr['ln_sd']:.3f}, of f_rho^2 "
            f"{fr['ln_sd_f2']:.3f}) from {fr['Y']:.1f} pairs / {fr['M']:.2f} predicted at f_rho = 1 "
            f"over {fr['nP']} P, dispersion {fr['phi']:.2f}, prior {fr['prior']}")


# --------------------------------------------------------------------------
# the predictive E


def load_ecurve(path):
    with open(path) as f:
        d = json.load(f)
    pts = sorted((float(c), float(e)) for c, e in d["E"].items())
    return d, pts


def e_at(pts, C):
    """E_ship(C): log-log interpolation in the curve, extrapolated with the
    last segment's slope"""
    if C <= pts[0][0]:
        (c0, e0), (c1, e1) = pts[0], pts[1]
    elif C >= pts[-1][0]:
        (c0, e0), (c1, e1) = pts[-2], pts[-1]
    else:
        i = max(j for j in range(len(pts)) if pts[j][0] <= C)
        (c0, e0), (c1, e1) = pts[i], pts[min(i + 1, len(pts) - 1)]
        if c1 == c0:
            return e0
    b = math.log(e1 / e0) / math.log(c1 / c0)
    return e0 * (C / c0) ** b


def class_sd(C):
    (c0, s0), (c1, s1) = CLASS_SD
    w = min(max(math.log(C / c0) / math.log(c1 / c0), 0.0), 1.0)
    return s0 + w * (s1 - s0)


def predictive(fr, pts, Cs, draws=20000, seed=1):
    """per C: E point (median f_rho^2), E draws' 5/50/95% and mean, P(>=1)"""
    import numpy as np
    rng = np.random.default_rng(seed)
    if fr is None:
        f2 = np.ones(draws)
        f2med = 1.0
    else:
        f2 = rng.gamma(fr["shape"], 1.0 / fr["rate"], draws) ** 2
        f2med = fr["median"] ** 2
    z = rng.standard_normal(draws)
    out = {}
    for C in Cs:
        s = math.sqrt(class_sd(C) ** 2 + PAIR_SD ** 2 + SELECTION_SD ** 2)
        E = e_at(pts, C) * f2 * np.exp(s * z)
        q = np.percentile(E, [5, 50, 95])
        out[C] = {"E_ship": e_at(pts, C), "E_point": e_at(pts, C) * f2med,
                  "q05": float(q[0]), "q50": float(q[1]), "q95": float(q[2]),
                  "mean": float(E.mean()), "band": float(q[2] / q[0]),
                  "P1": float((1 - np.exp(-E)).mean())}
    return out


def cpu_for(fr, pts, p, draws=20000, lo=0.01, hi=1e5):
    """CPU-years at which P(>=1) = p (bisection in ln C; the curve is
    extrapolated beyond its last point)"""
    def P(C):
        return predictive(fr, pts, [C], draws)[C]["P1"]
    if P(hi) < p:
        return float("inf")
    for _ in range(40):
        mid = math.sqrt(lo * hi)
        if P(mid) >= p:
            hi = mid
        else:
            lo = mid
    return hi


# --------------------------------------------------------------------------
# records


def scheduler_args(state):
    """the scheduler's argument namespace (emit's defaults), with the
    state's stage-1 spec (stage1_6.json) when it has one"""
    import scheduler as sch
    got = {}
    emit = sch.v2_emit
    sch.v2_emit = lambda a: got.setdefault("args", a)
    argv = sys.argv
    st1 = []
    try:
        with open(os.path.join(state, "stage1_6.json")) as f:
            st1 = ["--stage1", ":".join(f"{v:g}" for v in json.load(f)["spec"])]
    except (OSError, ValueError, KeyError):
        pass
    sys.argv = ["scheduler.py", "--state", state, "emit", "--units", "1"] + st1
    try:
        sch.main()
    finally:
        sys.argv = argv
        sch.v2_emit = emit
    return got["args"]


def open_any(path):
    if path.endswith(".xz"):
        import lzma
        return lzma.open(path, "rt")
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def file_counts(path):
    """observed records of a unit file: plain squares of completed sums and
    calibration-stream squares, with their S, P and SP traversals; sums;
    CPU (sum records' cpu, else time)"""
    o = collections.Counter()
    pend = []
    try:
        with open_any(path) as f:
            for line in f:
                try:
                    r = json.loads(line)
                except ValueError:
                    continue
                t = r.get("type")
                if t in ("square", "csquare"):
                    pend.append(r)
                elif t in ("sum", "csum"):
                    mine = [q for q in pend if q["S"] == r["S"]]
                    pend = [q for q in pend if q["S"] != r["S"]]
                    if t == "csum" and r.get("mode") != "calib":
                        continue
                    o["sums" if t == "sum" else "csums"] += 1
                    o["cpu"] += float(r.get("cpu", r.get("time", 0.0)))
                    for q in mine:
                        o["squares"] += 1
                        o["S"] += q["s_count"]
                        o["P"] += q["p_count"]
                        o["pairs"] += q["sp_count"]
                        o["best"] = max(o["best"], q.get("best_score", 0))
                elif t == "dsum":
                    o["cpu"] += float(r.get("cpu", 0.0))
                    o["dpairs"] += float(r.get("pairs", 0))
    except OSError:
        pass
    return o


def selected_units(state, sch_obj, select, plan, units_dir):
    """[(path, frozen predictions or {})] of the selected records"""
    if plan:
        lines = [json.loads(x) for x in open_any(plan) if x.strip()]
        out = []
        for path in sorted(glob.glob(os.path.join(units_dir, "U*.jsonl"))):
            i = int(os.path.basename(path)[1:-6])
            out.append((path, lines[i]))
        return out
    units = sch_obj.units
    if select == "all":
        import scheduler as sch
        return [(p, {}) for p in sch.unit_files(units)]
    path = os.path.join(state, f"launched_{sch_obj.n}.jsonl")
    out = []
    have = {}
    for p in glob.glob(os.path.join(units, "*.jsonl*")):
        b = os.path.basename(p)
        have[b[:-3] if b.endswith(".gz") else b] = p
    try:
        with open(path) as f:
            for line in f:
                r = json.loads(line)
                if r.get("stage1") and r["file"] in have:
                    out.append((have[r["file"]], r))
    except OSError:
        pass
    return out


def band_cells(lo, hi):
    """the N' bands (NBAND index) of the cells inside [lo, hi)"""
    import scheduler as sch
    edges = (0,) + tuple(sch.NBAND_EDGES) + (float("inf"),)
    return [b for b in range(len(edges) - 1) if edges[b] >= lo - 1e-9 and edges[b + 1] <= hi + 1e-9]


def dfirst_sums(paths, s):
    """the d-first sums searched in full in one "dsum" record (complete,
    d_stride 1, not truncated, with an estimate of its pairs: est_pairs over
    every d, the star d weighed K x) of the files: [(P key, S, cell, est
    pairs, the model's pairs before the class factors m0 = E_sq 720
    e^{lpSP})], once per sum (as frho_block's _sp_dfirst; sums split into
    units of d are left out)"""
    import numpy as np
    import scheduler as sch
    am = sch._am()
    seen = {}
    for path in paths:
        try:
            with open_any(path) as f:
                for line in f:
                    if '"dsum"' not in line:
                        continue
                    r = json.loads(line)
                    if (r.get("type") != "dsum" or not r.get("complete") or r.get("truncated")
                            or r.get("d_stride", 1) != 1 or r.get("n") != s.n):
                        continue
                    est = sch.dsum_est_pairs(r)
                    if est is None:
                        continue
                    P = sch.norm_p(r["P"])
                    seen.setdefault((P, int(r["S"])), est)
        except (OSError, ValueError):
            continue
    out = []
    for (P, S), est in sorted(seen.items()):
        try:
            A, v = s.store.get(P)
        except Exception:
            continue
        if not np.asarray(v).any():
            continue
        v = np.asarray(v, bool)
        Sg = am.grid_sums(P, s.n)[v]
        if not Sg[0] - 1e-9 <= S <= Sg[-1] + 1e-9:
            continue
        f = {fl: float(np.interp(S, Sg, np.asarray(A[j], float)[v])) for j, fl in enumerate(am.FIELDS)}
        lNp = f["lN"] + float(am.n_bias(f["lN"]))
        k6 = sum(1 for x in P if x)
        rb = int(sch.ratio_bin(np.float32(am.ratio(P, s.n))))
        cell = int(sch.cell_index(lNp, k6, rb, math.log(S / am.smin_approx(P, s.n))))
        esq = math.exp(f["lEs"]) * (am.SQ12 if lNp >= sch.LN12K else 1.0)
        out.append((sch.p_str(P, "_"), S, cell, est, esq * sch.NUM_TRAVERSALS[s.n] * math.exp(f["lpSP"])))
    return out


def pair_data(summary, calib, bands, dsums=()):
    """per P with model pairs: (Y_P, M_P) inside the N' bands and outside:
    its plain squares' pairs (cells; the calibration streams of the d-first
    sums in dsums left out, their pairs are in the d-first estimate) and
    its d-first sums' estimated pairs (dsums) against the model's"""
    import numpy as np
    import scheduler as sch
    lg, lS, lP, _ = calib.rates()
    kq = sch.sp_kappa() * np.exp(lS + lP)
    dby = collections.defaultdict(list)
    for key, S, c, est, m0 in dsums:
        dby[key].append((S, c, est, m0))
    inside, outside = [], []
    per_band = collections.defaultdict(lambda: [0.0, 0.0])
    per_mode = collections.defaultdict(lambda: [0.0, 0.0])
    for key in set(summary.perP) | set(dby):
        ps = summary.perP.get(key, {})
        acc = collections.defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])   # band -> y, m plain, y, m d-first
        for c, v in ps.get("cells", {}).items():
            c = int(c)
            acc[c // 36][0] += float(v[7])
            acc[c // 36][1] += float(v[8]) * float(kq[c])
        fin = {S for S, _, _, _ in dby.get(key, ())}
        for S, row in (ps.get("calS") or {}).items():
            if row and int(S) in fin:
                c = int(row[0])
                acc[c // 36][0] -= float(row[7])
                acc[c // 36][1] -= float(row[8]) * float(kq[c])
        for S, c, est, m0 in dby.get(key, ()):
            acc[c // 36][2] += est
            acc[c // 36][3] += m0 * float(np.exp(lg[c])) * float(kq[c])
        yi = mi = yo = mo = 0.0
        for b, (y, m, yd, md) in acc.items():
            per_band[b][0] += y + yd
            per_band[b][1] += m + md
            per_mode["plain"][0] += y
            per_mode["plain"][1] += m
            per_mode["d-first"][0] += yd
            per_mode["d-first"][1] += md
            if b in bands:
                yi, mi = yi + y + yd, mi + m + md
            else:
                yo, mo = yo + y + yd, mo + m + md
        if mi > 0 or yi > 0:
            inside.append((yi, mi))
        if mo > 0 or yo > 0:
            outside.append((yo, mo))
    return inside, outside, per_band, per_mode


def poisson_interval(o, e, z=1.645):
    """approximate 90% interval of o / e (Wilson-Hilferty on the count)"""
    if e <= 0:
        return float("nan"), float("nan"), float("nan")
    lo = 0.0 if o <= 0 else o * (1 - 1 / (9 * o) - z / (3 * math.sqrt(o))) ** 3
    o1 = o + 1
    hi = o1 * (1 - 1 / (9 * o1) + z / (3 * math.sqrt(o1))) ** 3
    return o / e, lo / e, hi / e


# --------------------------------------------------------------------------
# the report


def report(a):
    import numpy as np
    import scheduler as sch
    args = scheduler_args(a.state)
    s = sch.SchedulerV2(args, quiet=True)
    pr = print
    spec = s.stage1.spec
    select = a.select or ("stage1" if spec else "all")
    lo, hi = (map(float, a.band.split(":")) if a.band else
              ((spec[1], spec[2]) if spec else (sch.STAGE1_DEFAULT[1], sch.STAGE1_DEFAULT[2])))
    bands = band_cells(lo, hi)
    pr(f"# stage-1 decision report: {os.path.abspath(a.state)}")
    pr(s.stage1_status())
    units = selected_units(a.state, s, select, a.plan, a.units_dir)
    pr(f"records: {'plan ' + a.plan if a.plan else '--select ' + select}: {len(units)} unit files")
    if select == "all" and not a.plan:
        summ = s.summary
    else:
        summ = sch.Summary(tempfile.mkdtemp(prefix="decide_"), s.n)
        for path, _ in units:
            summ.update_file(path, s.store)
    calib = s.calib
    # checks: observed / predicted squares, S and P traversals (shipped and
    # learned class factors), per N' band
    tab = summ.cell_table()
    out = {"state": os.path.abspath(a.state), "select": select, "units": len(units),
           "band": [lo, hi]}
    pr(f"\nchecks (obs / pred, 90% Poisson; pred at the state's learned class factors, and at "
       f"the shipped ones in brackets):")
    for name, cal in (("learned", calib), ("shipped", sch.Calibration())):
        lg, lS, lP, _ = cal.rates()
        rows = []
        for b in range(len(sch.NBAND_NAMES)):
            m = (np.arange(sch.NCELLS) // 36) == b
            t = tab[m]
            if t[:, 0].sum() <= 0:
                continue
            rows.append((sch.NBAND_NAMES[b], t[:, 1].sum(), (t[:, 2] * np.exp(lg[m])).sum(),
                         t[:, 3].sum(), (t[:, 4] * np.exp(lS[m])).sum(),
                         t[:, 5].sum(), (t[:, 6] * np.exp(lP[m])).sum()))
        out["checks_" + name] = rows
    for (bn, osq, esq, oS, eS, oP, eP), r2 in zip(out["checks_learned"], out["checks_shipped"]):
        f = lambda o, e, e2: "{:.0f}/{:.0f} = {:.3f} [{:.3f}, {:.3f}] ({:.3f})".format(
            o, e, *poisson_interval(o, e), o / e2 if e2 > 0 else float("nan"))
        pr(f"  {bn:6} squares {f(osq, esq, r2[2])}; S trav {f(oS, eS, r2[4])}; "
           f"P trav {f(oP, eP, r2[6])}")
    # frozen per-unit predictions against the unit files
    if units and any(u for _, u in units):
        O = collections.Counter()
        Pd = collections.Counter()
        for path, u in units:
            fc = file_counts(path)
            best = max(fc.pop("best", 0), O.get("best", 0))
            O.update(fc)
            O["best"] = best
            for k in ("squares", "S", "P", "pairs", "SP"):
                Pd[k] += float(u.get("pred_" + k, 0.0) or 0.0)
            Pd["time"] += float(u.get("pred_time", u.get("time", 0.0)) or 0.0)
        out["frozen"] = {"observed": dict(O), "predicted": dict(Pd)}
        pr(f"\nfrozen predictions of the {len(units)} units (at launch / in the plan, f_rho = 1):")
        for k in ("squares", "S", "P", "pairs"):
            r, l, h = poisson_interval(O[k], Pd[k])
            pr(f"  {k:8} {O[k]:8.0f} observed / {Pd[k]:10.2f} predicted = {r:.3f} [{l:.3f}, {h:.3f}]")
        pr(f"  CPU      {O['cpu'] / 3600:8.3f} h measured / {Pd['time'] / 3600:.3f} h predicted "
           f"(reference); SP traversals before the class factors predicted {Pd['SP']:.2f}")
        if O["best"] >= 9:
            pr(f"  best square: best_score {O['best']} (14 = magic)")
    # the fit
    dsums = dfirst_sums([p for p, _ in units], s) if a.dfirst == "on" else []
    ins, outs, per_band, per_mode = pair_data(summ, calib, bands, dsums)
    fr = fit_frho([y for y, _ in ins], [m for _, m in ins], a.prior, a.phi_min)
    fro = fit_frho([y for y, _ in outs], [m for _, m in outs], "flat")
    out["fit"], out["fit_outside_flat"] = fr, fro
    pr(f"\nSP coupling (plain pairs in N' {lo:g}-{hi:g}, cells {[sch.NBAND_NAMES[b] for b in bands]}):")
    pr("  " + describe_fit(fr))
    pr("  check, outside the band (flat prior, not used): " + describe_fit(fro))
    if fro is not None and fr is not None and fro["phi"] > fr["phi"]:
        # pairs come clustered (by P, and by SP vector within a sum): with
        # few pairs in the band its own dispersion is poorly estimated
        frc = fit_frho([y for y, _ in ins], [m for _, m in ins], a.prior, fro["phi"])
        pc = predictive(frc, load_ecurve(a.ecurve)[1], [1.0], a.draws)[1.0]
        out["fit_phi_outside"] = frc
        pr(f"  sensitivity, at the dispersion outside the band ({fro['phi']:.2f}): "
           + describe_fit(frc) + f"; x band of E at 1 CPU-year {pc['band']:.2f}")
    pr(f"  data: plain squares {per_mode['plain'][0]:.0f} pairs / {per_mode['plain'][1]:.2f}; "
       f"{len(dsums)} d-first sums searched in full in one record (--dfirst-sums "
       f"{a.dfirst}): {per_mode['d-first'][0]:.1f} est. pairs / "
       f"{per_mode['d-first'][1]:.2f}")
    pr("  pairs per N' band (obs / model at f_rho = 1): " + ", ".join(
        f"{sch.NBAND_NAMES[b]} {v[0]:.0f}/{v[1]:.2f}" for b, v in sorted(per_band.items())
        if v[1] > 0 or v[0] > 0))
    # predictive E
    meta, pts = load_ecurve(a.ecurve)
    Cs = sorted(set(MARKS) | {float(x) for x in a.marks.split(",")} if a.marks else set(MARKS))
    pred = predictive(fr, pts, Cs, a.draws)
    pred0 = predictive(fit_frho([], [], "shipped"), pts, Cs, a.draws)
    out["predictive"] = {str(C): v for C, v in pred.items()}
    out["predictive_prior"] = {str(C): v for C, v in pred0.items()}
    pr(f"\npredictive E (E_ship from {os.path.relpath(a.ecurve, ROOT)}: {meta.get('source', '')}; "
       f"x f_rho^2 x lognormal(class factors, pair {PAIR_SD}, selection {SELECTION_SD})):")
    pr(f"  {'CPU-years':>9} {'E_ship':>7} {'E (f^2 med)':>11} {'5%':>7} {'median':>7} {'95%':>7} "
       f"{'x band':>6} {'P(>=1)':>7}   before stage 1 (shipped prior): P(>=1), x band")
    for C in Cs:
        v, v0 = pred[C], pred0[C]
        pr(f"  {C:9g} {v['E_ship']:7.3f} {v['E_point']:11.3f} {v['q05']:7.3f} {v['q50']:7.3f} "
           f"{v['q95']:7.3f} {v['band']:6.2f} {100 * v['P1']:6.1f}%   {100 * v0['P1']:5.1f}%, "
           f"x{v0['band']:.2f}")
    need = {p: cpu_for(fr, pts, p, a.draws // 4) for p in (0.25, 0.5)}
    need0 = {p: cpu_for(fit_frho([], [], "shipped"), pts, p, a.draws // 4) for p in (0.25, 0.5)}
    out["cpu_for"] = {str(p): v for p, v in need.items()}
    ext = max(c for c, _ in pts)
    pr("  CPU-years for P(>=1) = 25%: {:.3g}, 50%: {:.3g} (before stage 1: {:.3g}, {:.3g}; the "
       "curve is extrapolated beyond {:g} CPU-years at its last slope)".format(
           need[0.25], need[0.5], need0[0.25], need0[0.5], ext))
    # the decision table: evaluated, never acted on
    try:
        with open(a.table) as f:
            table = json.load(f)
    except OSError:
        table = None
    if table:
        pr(f"\ndecision table ({os.path.relpath(a.table, ROOT)}; thresholds set by the user "
           f"before stage 1; this tool only reports):")
        rows = []
        for row in table.get("rows", []):
            C = float(row["C"])
            val = predictive(fr, pts, [C], a.draws)[C]["P1"]
            X = row.get("X")
            op = row.get("op", ">=")
            if X is None:
                verdict = "X not set"
            else:
                ok = val >= X if op == ">=" else val < X
                verdict = "condition holds" if ok else "condition does not hold"
            rows.append(dict(row, value=val, verdict=verdict))
            pr(f"  {row['action']:42} if P(>=1 at {C:g} CPU-years) {op} "
               f"{'X' if X is None else f'{X:.0%}'}: P = {100 * val:.1f}%: {verdict}")
        out["table"] = rows
    if a.json:
        with open(a.json, "w") as f:
            json.dump(out, f, indent=1, default=float)


# --------------------------------------------------------------------------
# the simulation gate


def simulate(a):
    """pair counts drawn from the model at a known f_rho, clustered by P as
    the plan's units are (M_P per P from the plan's band units, truncated
    where the band's predicted pairs reach --pairs), with gamma per-P effects
    of the variance that gives the Pearson dispersion --phi (the
    calibration search's 1.41 of S traversals; the pairs' own was 1.00 on 7
    pairs); fitted as decide.py fits; coverage of the 90% interval and bias
    of the posterior median over --reps replicates per f_rho, flat and
    shipped prior, and the predictive band of E at 1 CPU-year"""
    import numpy as np
    lines = [json.loads(x) for x in open_any(a.simulate) if x.strip()]
    MP = collections.defaultdict(float)
    cum = 0.0
    hours = 0.0
    for r in lines:
        hours += r["pred_time"] / 3600
        if not r.get("stage1"):
            continue
        MP[tuple(r["P"])] += r["pred_pairs"]
        cum += r["pred_pairs"]
        if cum >= a.pairs:
            break
    M = np.array([m for m in MP.values() if m > 0])
    k2 = float((M ** 2).sum() / M.sum())
    _, pts = load_ecurve(a.ecurve)
    rng = np.random.default_rng(a.seed)
    print(f"plan {a.simulate}: the band's first {cum:.1f} predicted pairs (f_rho = 1) after "
          f"{hours:.1f} reference CPU-hours, {len(M)} P, sum M_P^2 / sum M_P = {k2:.3f}")
    res = {}
    for phi in a.phi:
        for f in (0.5, 1.0, 2.0):
            # phi = 1 + f tau^2 k2 (Pearson dispersion of the per-P totals)
            tau2 = max(phi - 1.0, 0.0) / (f * k2)
            for prior in ("flat", "shipped"):
                cover = 0
                meds = []
                bands = []
                phis = []
                for _ in range(a.reps):
                    u = rng.gamma(1 / tau2, tau2, len(M)) if tau2 > 0 else 1.0
                    Y = rng.poisson(f * M * u)
                    fr = fit_frho(Y, M, prior)
                    if fr is None:
                        continue
                    cover += fr["q05"] <= f <= fr["q95"]
                    meds.append(fr["median"])
                    phis.append(fr["phi"])
                    if prior == "shipped" and f == 1.0 and len(bands) < 200:
                        bands.append(predictive(fr, pts, [1.0], 4000, seed=len(bands))[1.0]["band"])
                n = len(meds)
                r = {"cover": cover / max(n, 1), "bias": float(np.mean(meds)) / f - 1,
                     "phi_hat": float(np.mean(phis)), "n": n}
                if bands:
                    r["band_1y"] = float(np.median(bands))
                res[(phi, f, prior)] = r
                print(f"phi {phi:.2f} f_rho {f:3.1f} prior {prior:7}: 90% interval covers "
                      f"{100 * r['cover']:.1f}%, mean posterior median / truth - 1 = "
                      f"{100 * r['bias']:+.1f}%, mean phi_hat {r['phi_hat']:.2f} ({n} replicates)"
                      + (f"; median x band of E at 1 CPU-year {r['band_1y']:.2f}" if bands else ""))
    # Bayesian check: f_rho drawn from the shipped prior
    f0, sd0 = FRHO_PRIORS["shipped"]
    a_ = 1 / sd0 ** 2 + 0.5
    for phi in a.phi:
        cover = 0
        for _ in range(a.reps):
            f = rng.gamma(a_, f0 / a_)
            tau2 = max(phi - 1.0, 0.0) / (f * k2)
            u = rng.gamma(1 / tau2, tau2, len(M)) if tau2 > 0 else 1.0
            fr = fit_frho(rng.poisson(f * M * u), M, "shipped")
            cover += fr["q05"] <= f <= fr["q95"]
        print(f"phi {phi:.2f}, f_rho drawn from the shipped prior: the shipped-prior interval "
              f"covers {100 * cover / a.reps:.1f}%")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("state", nargs="?")
    ap.add_argument("--select", choices=("stage1", "all"), default=None)
    ap.add_argument("--plan")
    ap.add_argument("--units-dir")
    ap.add_argument("--prior", choices=tuple(FRHO_PRIORS), default="shipped")
    ap.add_argument("--band", help="LO:HI of the N' band whose pairs are fitted (default the "
                                    "state's stage-1 band, else 3000:6000)")
    ap.add_argument("--ecurve", default=os.path.join(ROOT, "research", "stage1", "ecurve.json"))
    ap.add_argument("--table", default=os.path.join(ROOT, "research", "stage1",
                                                    "decision-table.json"))
    ap.add_argument("--dfirst-sums", dest="dfirst", choices=("on", "off"), default="on",
                    help="also fit the estimated pairs of d-first sums searched in full in one "
                         "record (default on)")
    ap.add_argument("--phi-min", type=float, default=1.0,
                    help="floor of the pairs' dispersion in the fit (default 1: the estimate)")
    ap.add_argument("--marks", help="extra CPU-years for the predictive table, e.g. 0.3,30")
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--json")
    ap.add_argument("--simulate", metavar="PLAN")
    ap.add_argument("--reps", type=int, default=1000)
    ap.add_argument("--pairs", type=float, default=45)
    ap.add_argument("--phi", type=float, nargs="+", default=[1.41])
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.simulate:
        simulate(a)
        return
    if not a.state:
        ap.error("STATE or --simulate PLAN")
    if bool(a.plan) != bool(a.units_dir):
        ap.error("--plan and --units-dir go together")
    report(a)


if __name__ == "__main__":
    main()
