#!/usr/bin/env python3
"""The sum-product coupling of the sub-events (calibrate.py, "X + j exponents
pinned" over "j exponents alone") with the SP traversals taken out.

An SP traversal satisfies every window of the square's k exponents, so it is
counted in the "X + j" sub-event once per window: k times for j < k, once
for j = k. On the 10,677 plain squares at 3-12k the 7 SP traversals make 32
of the 35 events at j = 5 and 40 of the 109 at j = 4, so the sub-event
coupling at j >= 4 is not independent of the SP count. Here they are
removed from the observed side (O' = O - sum m(j, k) n_SP) and their
heuristic expectation from the expected side (E' = E - sum m(j, k) E_SP,
per square, the heuristic's SP rate), on both the X + j and the j-alone
totals. This is a cross-check of the direction of the SP coupling, not a
measurement of it: going from j = 4 to j = k needs a model of how the
coupling grows with j, which is not derived here, and E_SP uses the
"heuristic" variant while the sub-event E uses the "pair" variant
(they differ by ~10-30%), so E' at j = 5 is not reliable (the SP part is
~95% of it).

Inputs: $CALIB_TARGET_DIR/calibrate_unw/{results.json, per_square.jsonl}
(calibrate.py --unweighted on the 307 plain sums at 3-12k). Output:
analysis/subevents_nonsp.json, printed table.
"""
import collections
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ctpaths  # noqa: E402

sys.path.insert(0, ctpaths.SCRIPTS)
import calibrate as CAL  # noqa: E402

D = os.path.join(ctpaths.CT, "calibrate_unw")
res = json.load(open(os.path.join(D, "results.json")))
rows = {(r["kind"], r["d"]): r for r in res["sub_events"]["all"]["rows"]}

nsp, esp = collections.Counter(), collections.Counter()
with open(os.path.join(D, "per_square.jsonl")) as f:
    for line in f:
        r = json.loads(line)
        k = sum(1 for e in r["P"] if e)
        nsp[k] += r["obs"]["SP"]
        esp[k] += r["heuristic"]["SP"] if r["heuristic"].get("SP") is not None else 0.0


def mult(j, k):
    """windows of j of the k exponent coordinates that one SP traversal satisfies"""
    return k if j < k else (1 if j == k else 0)


out = {"n_sp_by_k": dict(nsp), "E_sp_heuristic_by_k": dict(esp), "rows": []}
print("| j | X + j: O / E | SP part O / E | without SP: O' / E' | j alone O / E | coupling | without SP | "
      "Poisson 90% of O' (x1.5 for clustering) |")
print("|---:|---|---|---|---|---:|---:|---|")
for j in range(1, 7):
    xy, y = rows.get(("XY", j + 1)), rows.get(("Y", j))
    if not xy or not y:
        continue
    o_sp = sum(mult(j, k) * n for k, n in nsp.items())
    e_sp = sum(mult(j, k) * e for k, e in esp.items())
    O, E = xy["O1"], xy["E1_pair"]
    YO, YE = y["O1"], y["E1_pair"]
    c = (O / E) / (YO / YE)
    O2, E2 = O - o_sp, E - e_sp
    YO2, YE2 = YO - o_sp, YE - e_sp
    c2 = (O2 / E2) / (YO2 / YE2) if E2 > 0 and O2 >= 0 else float("nan")
    lo, hi = CAL.poisson_ci(int(round(O2)), 0.90) if O2 >= 0 else (float("nan"),) * 2
    # widen the Poisson interval by sqrt(1.5^2)~ for the clustering of a square's windows
    # (calibrate.py's bootstrap over squares is ~1.5x the Poisson width at j = 4)
    if O2 > 0 and E2 > 0:
        sd = 1.5 * (math.log(max(hi, 1e-9)) - math.log(max(lo, 1e-9))) / (2 * 1.645)
        ci = [c2 * math.exp(-1.645 * sd), c2 * math.exp(1.645 * sd)]
    else:
        ci = [float("nan")] * 2
    row = {"j": j, "O": O, "E": E, "O_sp": o_sp, "E_sp": e_sp, "O_nonsp": O2, "E_nonsp": E2,
           "YO": YO, "YE": YE, "coupling": c, "coupling_nonsp": c2, "ci90_nonsp_approx": ci,
           "sp_share_of_E": e_sp / E if E > 0 else None}
    out["rows"].append(row)
    print(f"| {j} | {O:.0f} / {E:.1f} | {o_sp:.0f} / {e_sp:.1f} | {O2:.0f} / {E2:.1f} | {YO:.0f} / {YE:.0f} | "
          f"{c:.2f} | {c2:.2f} | [{ci[0]:.2f}, {ci[1]:.2f}] |")
os.makedirs(ctpaths.OUT, exist_ok=True)
with open(os.path.join(ctpaths.OUT, "subevents_nonsp.json"), "w") as f:
    json.dump(out, f, indent=1)
