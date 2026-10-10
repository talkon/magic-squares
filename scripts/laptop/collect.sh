#!/bin/sh
# Pack the outputs of scripts/laptop/run.py for analysis:
# research/laptop/results-<host>-<date>.tar.xz (only finished units).
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
OUT=${1:-$ROOT/laptop-runs}
NAME=results-$(hostname -s 2>/dev/null || echo host)-$(date +%Y%m%d-%H%M).tar.xz
mkdir -p "$ROOT/research/laptop"
cd "$OUT"
ls | grep -E '^U[0-9]+\.jsonl$|^meta_.*\.json$|^progress\.log$|^MAGIC\.txt$' > .collect_list
tar -cJf "$ROOT/research/laptop/$NAME" -T .collect_list
rm .collect_list
echo "wrote research/laptop/$NAME ($(du -h "$ROOT/research/laptop/$NAME" | cut -f1))"
echo "to send it back: git add research/laptop/$NAME && git commit -m 'laptop run results' && git push"
