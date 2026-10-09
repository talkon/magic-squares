"""calibrate.py plugin: scheduler v2's pre-registered per-traversal rates
(predictions.json, frozen) for the square's (P, S), as an i.i.d. model."""
import json
import os

PRED = os.path.join(os.environ.get("CALIB_TARGET_DIR", os.path.join(os.path.dirname(os.path.abspath(__file__)), "work")),
                    "predictions.json")
_T = None


def _norm(P):
    P = list(P)
    while P and P[-1] == 0:
        P.pop()
    return tuple(P)


def predict(grid, S, P):
    global _T
    if _T is None:
        _T = {(_norm(r["P"]), int(r["S"])): r for r in json.load(open(PRED))["runs"]}
    r = _T[(_norm(P), int(S))]
    return {"q_S_only": r["q_S"] - r["q_SP"], "q_P_only": r["q_P"] - r["q_SP"], "q_SP": r["q_SP"]}
