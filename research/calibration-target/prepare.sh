#!/bin/sh
# Unpack the archived calibration-target data (archive/) into $CALIB_TARGET_DIR
# (default: work/ next to this file) for analyze.py, forecast_update.py,
# bands.py and subevents_nonsp.py, and check the frozen files against
# PREREGISTERED.txt.
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
W=${CALIB_TARGET_DIR:-$HERE/work}
mkdir -p "$W/work" "$W/analysis"
tar -xJf "$HERE/archive/runs.tar.xz" -C "$W"
for f in predictions.json plan.json run.sh runner.log; do
  xz -dc "$HERE/archive/$f.xz" > "$W/$f"
done
cp "$HERE/archive/PREREGISTERED.txt" "$HERE/archive/runner.py" "$W/"
cp "$HERE/archive/work/"* "$W/work/"
cd "$W"
for f in predictions.json plan.json run.sh runner.py work/sim.py work/cm.py work/explore.py work/design.py work/searched.py; do
  want=$(grep -E "^$f sha256:" PREREGISTERED.txt | sed 's/.*: *//; s/ .*//')
  got=$(sha256sum "$f" | cut -d' ' -f1)
  [ "$want" = "$got" ] || { echo "sha256 mismatch: $f" >&2; exit 1; }
done
echo "unpacked into $W ($(ls runs/R*.jsonl | wc -l) runs); frozen files match PREREGISTERED.txt"
