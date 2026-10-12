#!/bin/sh
# Pack the outputs of scripts/laptop/run.py that were not packed before:
# research/laptop/results-<seq>[-partNN].tar.xz, each part under ~45 MB,
# with the build's record, units.tsv, progress.log and the meta files
# (scrubbed: no host name, OS version, path, time zone or clock time). laptop-runs/collected.json remembers what was packed. It prints
# the files to git add. See scripts/laptop/collect.py.
#
# usage: collect.sh [DIR]   (run.py's --out, default laptop-runs)
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
exec python3 "$ROOT/scripts/laptop/collect.py" "$@"
