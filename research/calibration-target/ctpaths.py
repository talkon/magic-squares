"""Paths of the calibration-target analysis scripts.

CALIB_TARGET_DIR: the directory with predictions.json, plan.json,
PREREGISTERED.txt, run.sh, runner.py and runs/ (default: work/ next to this
file, which prepare.sh fills from archive/). The scripts write to
$CALIB_TARGET_DIR/analysis. The repository's scripts/ (scheduler.py,
calibrate.py, amodel.py) are imported from the checkout this file is in.
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SCRIPTS = os.path.join(REPO, "scripts")
CT = os.path.abspath(os.environ.get("CALIB_TARGET_DIR", os.path.join(HERE, "work")))
OUT = os.path.join(CT, "analysis")


def run_path(p):
    """a run's output file: predictions.json names the absolute path it was
    written to; the file is looked up by its name under CT/runs"""
    return os.path.join(CT, "runs", os.path.basename(p))
