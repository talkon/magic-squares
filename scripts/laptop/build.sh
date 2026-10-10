#!/bin/sh
# Build msearch, bench and fuzz_arrange without CMake (macOS or Linux, x86 or
# ARM) into build-laptop/, then check the build: every bench/quick.txt and
# bench/full.txt instance must give the expected squares and hash, with the
# node counts of the fast x86 build (the same search, also without
# AVX-512), and the differential fuzz test must pass, with the node count
# of every seed and variant of the fast x86 build. Prints the time of
# bench/quick.txt and the time per core it implies against the fast x86
# build, and adds them, with the search kernels, to
# build-laptop/build_info.txt (run.py copies that file into its meta file,
# collect.sh packs it).
#
# The target is the host's CPU by default: on an M1 with Apple clang,
# `-mcpu=native -flto -D_DARWIN_C_SOURCE` and the portable kernels; on x86,
# -march=native. CC picks the compiler and may name a target itself (e.g.
# CC='gcc -march=x86-64-v3'); ARCH replaces the target flags (ARCH= for
# none). On an x86 machine with AVX-512, ARCH=-march=x86-64-v2 (or -v3)
# builds the portable kernels, to rehearse the laptop's build.
set -e
ROOT=$(cd "$(dirname "$0")/../.." && pwd)
B=$ROOT/build-laptop
mkdir -p "$B"
CC=${CC:-cc}
cd "$ROOT/src/c"
probe() { echo 'int main(void){return 0;}' > "$B/probe.c"; $CC $1 -o "$B/probe" "$B/probe.c" 2>/dev/null; }
# the target: the host's (-mcpu=native on ARM, e.g. Apple clang on an M1;
# -march=native on x86), unless CC already names one (e.g.
# CC='gcc -march=x86-64-v3') or ARCH is set (ARCH= for none)
if [ -n "${ARCH+set}" ]; then
  if [ -n "$ARCH" ]; then
    probe "$ARCH" || { echo "the compiler does not take ARCH=$ARCH"; exit 1; }
  fi
else
  ARCH=""
  case " $CC " in
    *" -march="* | *" -mcpu="*) ;;
    *)
      case $(uname -m) in
        arm64 | aarch64) CANDS="-mcpu=native -mcpu=apple-m1 -march=native" ;;
        *) CANDS="-march=native" ;;
      esac
      for f in $CANDS; do
        if probe "$f"; then ARCH=$f; break; fi
      done ;;
  esac
fi
LTO=""
if probe -flto; then LTO=-flto; fi
DEFS=""
if [ "$(uname)" = Darwin ]; then DEFS=-D_DARWIN_C_SOURCE; fi
FLAGS="-O3 -DNDEBUG -std=gnu17 $ARCH $LTO $DEFS"
echo "compiler: $CC ($($CC --version | head -1))"
echo "flags: $FLAGS"
CORE="enumerate.c arrange.c square.c dfirst.c"
for t in msearch bench fuzz_arrange; do
  $CC $FLAGS -o "$B/$t" $t.c $CORE -lm
done
{ echo "compiler: $CC ($($CC --version | head -1))"; echo "flags: $FLAGS"; uname -a;
  git -C "$ROOT" rev-parse HEAD 2>/dev/null || true; } > "$B/build_info.txt"
cd "$ROOT"
echo "checking bench/quick.txt and bench/full.txt (expected squares and hashes)"
"$B/bench" bench/quick.txt > "$B/check_quick.txt"
"$B/bench" bench/full.txt > "$B/check_full.txt"
if grep -v '^#' "$B/check_quick.txt" "$B/check_full.txt" | grep -Eq 'MISMATCH|FAIL| bad'; then
  echo "BUILD CHECK FAILED: see $B/check_*.txt"; exit 1
fi
for f in quick full; do
  if ! tail -1 "$B/check_$f.txt" | grep -q ' ok '; then
    echo "BUILD CHECK FAILED: see $B/check_$f.txt"; exit 1
  fi
done
# the search is the carried path (with AVX-512BW its AVX-512 kernels,
# otherwise the plain C ones of src/c/arrange_carry.h, e.g. on ARM), which
# visits exactly the nodes of the fast x86 build
if $CC $FLAGS -dM -E -x c /dev/null 2>/dev/null | grep -q __AVX512BW__; then
  KERNELS="search kernels: AVX-512"
else
  KERNELS="search kernels: portable (plain C, no AVX-512BW)"
fi
echo "$KERNELS"
for fn in quick:1770779 full:14958507; do
  f=${fn%%:*}; n=${fn#*:}
  got=$(tail -1 "$B/check_$f.txt" | awk '{print $2}')
  if [ "$got" != "$n" ]; then
    echo "BUILD CHECK FAILED: bench/$f.txt took $got nodes, the fast x86 build $n (not the same search)"; exit 1
  fi
done
echo "checking the search against brute force (fuzz_arrange, 60 seeds, ~1 min)"
cd "$B"
./fuzz_arrange -v 7300000 30 > check_fuzz.txt 2>&1 || { echo "FUZZ FAILED: see $B/check_fuzz.txt"; exit 1; }
./fuzz_arrange -v --dfirst 7400000 30 >> check_fuzz.txt 2>&1 || { echo "FUZZ FAILED: see $B/check_fuzz.txt"; exit 1; }
if grep -Eq 'MISMATCH|fails [1-9]' check_fuzz.txt; then echo "FUZZ FAILED: see $B/check_fuzz.txt"; exit 1; fi
# the node count of every seed and variant, as in the fast x86 build (the
# totals of bench quick/full miss some divergences that keep the squares,
# e.g. every label count off by one; the fuzz seeds do not); the lines of
# the matrix path (more than 512 labels, or forced 16/64 words) left out:
# its counters differ between builds
grep -E '^seed .* (squares|pairs)' check_fuzz.txt | sed -E 's/, [0-9.]+s \(oracle [0-9.]+s\)//' |
  grep -vE 'variant (w16|w64|xcross_w16|w16_nomrv)[:_ ]' |
  awk '{for (i = 1; i < NF; i++) if ($i == "L" && $(i + 1) > 512) next; print}' > check_fuzz_nodes.txt
got=$(cksum < check_fuzz_nodes.txt | awk '{print $1, $2}')
if [ "$got" != "3761910603 92405" ]; then
  echo "BUILD CHECK FAILED: the fuzz_arrange node counts differ from the fast x86 build's (cksum $got): not the same search, see $B/check_fuzz_nodes.txt"; exit 1
fi
cd "$ROOT"
echo "all checks passed"
{ echo "$KERNELS"; echo "all checks passed"; } >> "$B/build_info.txt"
# the speed per core: the best of three runs of bench/quick.txt against the
# fast x86 build's 0.27 s (an estimate of the plan's CPU too: clang without
# AVX-512 took 2.5x on bench/quick.txt and 2.55-3.0x on plain sums)
for i in 1 2; do "$B/bench" bench/quick.txt | tail -1; done > "$B/time_quick.txt"
tail -1 "$B/check_quick.txt" >> "$B/time_quick.txt"
awk 'NR == 1 || $4 < t { t = $4 } END {
  r = t / 0.27
  printf "bench/quick.txt search time: %.3f s (best of 3; fast x86 build: ~0.27 s)\n", t
  printf "time per core: ~%.1fx the fast x86 build; the stage-1 plan (60 CPU-hours of that build) would take ~%.0f CPU-hours here\n", r, 60 * r
}' "$B/time_quick.txt" | tee -a "$B/build_info.txt"
