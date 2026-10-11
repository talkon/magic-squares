#!/bin/sh
# Pack the outputs of scripts/laptop/run.py for analysis:
# research/laptop/results-<host>-<date>.tar.xz (only finished units), with
# the build's record (build-laptop/build_info.txt: compiler, flags, search
# kernels, bench/quick.txt time; and check_quick.txt) and units.tsv (each
# unit's CPU and wall time). Nothing else needs to be sent. The units enter
# the scheduler's state as counts only (scheduler.py ingest --counts-only).
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OUT=${1:-$ROOT/laptop-runs}
NAME=results-$(hostname -s 2>/dev/null || echo host)-$(date +%Y%m%d-%H%M).tar.xz
mkdir -p "$ROOT/research/laptop"
cd "$OUT"
for f in build_info.txt check_quick.txt; do
  if [ -f "$ROOT/build-laptop/$f" ]; then cp "$ROOT/build-laptop/$f" "build-laptop_$f"; fi
done
ls | grep -E '^U[0-9]+\.jsonl$|^meta_.*\.json$|^progress\.log$|^units\.tsv$|^MAGIC\.txt$|^build-laptop_.*\.txt$' > .collect_list
tar -cJf "$ROOT/research/laptop/$NAME" -T .collect_list
rm .collect_list
echo "wrote research/laptop/$NAME ($(du -h "$ROOT/research/laptop/$NAME" | cut -f1))"
echo "(here the counts go to decide.py and into the scheduler's state; the laptop's times stay out of its time laws)"
echo "to send it back: git add research/laptop/$NAME && git commit -m 'laptop run results' && git push"
