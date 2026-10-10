#!/bin/sh
# Build msearch, bench and fuzz_arrange without CMake (macOS or Linux, x86 or
# ARM) into build-laptop/, then check the build: every bench/quick.txt and
# bench/full.txt instance must give the expected squares and hash, and the
# differential fuzz test must pass. Prints the time of bench/quick.txt.
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
B=$ROOT/build-laptop
mkdir -p "$B"
CC=${CC:-cc}
cd "$ROOT/src/c"
probe() { echo 'int main(void){return 0;}' > "$B/probe.c"; $CC $1 -o "$B/probe" "$B/probe.c" 2>/dev/null; }
ARCH=""
for f in -march=native -mcpu=native -mcpu=apple-m1; do
  if probe "$f"; then ARCH=$f; break; fi
done
LTO=""
if probe -flto; then LTO=-flto; fi
DEFS=""
if [ "$(uname)" = Darwin ]; then DEFS=-D_DARWIN_C_SOURCE; fi
FLAGS="-O3 -DNDEBUG -std=gnu17 $ARCH $LTO $DEFS"
echo "compiler: $($CC --version | head -1)"
echo "flags: $FLAGS"
CORE="enumerate.c arrange.c square.c dfirst.c"
for t in msearch bench fuzz_arrange; do
  $CC $FLAGS -o "$B/$t" $t.c $CORE -lm
done
{ echo "compiler: $($CC --version | head -1)"; echo "flags: $FLAGS"; uname -a;
  git -C "$ROOT" rev-parse HEAD 2>/dev/null || true; } > "$B/build_info.txt"
cd "$ROOT"
echo "checking bench/quick.txt and bench/full.txt (expected squares and hashes)"
"$B/bench" bench/quick.txt > "$B/check_quick.txt"
"$B/bench" bench/full.txt > "$B/check_full.txt"
if grep -v '^#' "$B/check_quick.txt" "$B/check_full.txt" | grep -q 'MISMATCH\|FAIL\| bad'; then
  echo "BUILD CHECK FAILED: see $B/check_*.txt"; exit 1
fi
for f in quick full; do
  if ! tail -1 "$B/check_$f.txt" | grep -q ' ok '; then
    echo "BUILD CHECK FAILED: see $B/check_$f.txt"; exit 1
  fi
done
echo "checking the search against brute force (fuzz_arrange, 60 seeds, ~1 min)"
cd "$B"
./fuzz_arrange 7300000 30 > check_fuzz.txt 2>&1 || { echo "FUZZ FAILED: see $B/check_fuzz.txt"; exit 1; }
./fuzz_arrange --dfirst 7400000 30 >> check_fuzz.txt 2>&1 || { echo "FUZZ FAILED: see $B/check_fuzz.txt"; exit 1; }
if grep -q 'MISMATCH\|fails [1-9]' check_fuzz.txt; then echo "FUZZ FAILED: see $B/check_fuzz.txt"; exit 1; fi
cd "$ROOT"
echo "all checks passed"
tail -1 "$B/check_quick.txt" | awk '{print "bench/quick.txt search time: " $4 " s (fast x86 build: ~0.29 s)"}'
