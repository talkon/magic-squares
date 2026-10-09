# New pipeline: portable C, no AVX-512 required.
add_library(magic_core OBJECT enumerate.c arrange.c square.c dfirst.c)
add_executable(enumerate enumerate_main.c $<TARGET_OBJECTS:magic_core>)
add_executable(msearch msearch.c $<TARGET_OBJECTS:magic_core>)
add_executable(bench bench.c $<TARGET_OBJECTS:magic_core>)
add_executable(dsearch dsearch.c $<TARGET_OBJECTS:magic_core>)
foreach(t enumerate msearch bench dsearch)
  target_link_libraries(${t} m)
endforeach()

# Legacy arrangement program (uses AVX-512 intrinsics, so only built where
# the compiler targets AVX-512).
include(CheckCSourceCompiles)
set(CMAKE_REQUIRED_FLAGS "${CMAKE_C_FLAGS}")
check_c_source_compiles("
#include <immintrin.h>
int main(void) { __m512i x = _mm512_set1_epi32(1); return _mm512_reduce_add_epi32(x) == 16 ? 0 : 1; }
" HAVE_AVX512)
if (HAVE_AVX512)
  add_executable(arrangement_5 arrangement.c bitset.c search.c labelling.c)
  add_executable(arrangement_6 arrangement.c bitset.c search.c labelling.c)
  add_executable(arrangement_7 arrangement.c bitset.c search.c labelling.c)
  target_compile_definitions(arrangement_5 PUBLIC VEC_SIZE=5)
  target_compile_definitions(arrangement_6 PUBLIC VEC_SIZE=6)
  target_compile_definitions(arrangement_7 PUBLIC VEC_SIZE=7)
else()
  message(STATUS "No AVX-512: skipping legacy arrangement_{5,6,7}")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
add_executable(enumeration enumeration.cpp)

# Tests of the new pipeline (fast; run with `ctest -R fast_`).
set(BENCH_DIR ${PROJECT_SOURCE_DIR}/bench)
add_test(NAME fast_bench_quick COMMAND bench ${BENCH_DIR}/quick.txt)
# the same instances, with wider label bitsets than needed (carried up to
# 8 words with AVX-512BW, each r1 at the width it needs unless
# --no-r1-width)
add_test(NAME fast_bench_quick_wide3 COMMAND bench --min-words 3 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_wide5 COMMAND bench --min-words 5 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_wide7 COMMAND bench --min-words 7 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_wide8 COMMAND bench --min-words 8 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_wide6_fixed COMMAND bench --no-r1-width --min-words 6 ${BENCH_DIR}/quick.txt)
# the benchmark built with the N x N intersection matrices for all widths
# (the filter used without AVX-512BW, here still with byte counters), and
# with the fallbacks for AVX-512 CPUs without GFNI / VPOPCNTDQ (Skylake-X,
# Cascade Lake), to test those paths
add_executable(bench_matrices bench.c enumerate.c arrange.c square.c)
target_compile_definitions(bench_matrices PRIVATE CARRY_MAX_W=0)
target_link_libraries(bench_matrices m)
add_test(NAME fast_bench_quick_matrices COMMAND bench_matrices ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_matrices_no_support COMMAND bench_matrices --no-support ${BENCH_DIR}/quick.txt)
# the matrices with 8 words (as for 257-512 labels without AVX-512BW), and
# with gathers as for N > 3072
add_test(NAME fast_bench_quick_matrices_wide8 COMMAND bench_matrices --min-words 8 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_gather COMMAND bench_matrices --min-words 8 --gather ${BENCH_DIR}/quick.txt)
include(CheckCCompilerFlag)
check_c_compiler_flag("-mno-gfni -mno-avx512vpopcntdq" HAVE_NO_GFNI_FLAGS)
if (HAVE_NO_GFNI_FLAGS)
  add_executable(bench_no_gfni bench.c enumerate.c arrange.c square.c)
  target_compile_options(bench_no_gfni PRIVATE -mno-gfni -mno-avx512vpopcntdq)
  target_link_libraries(bench_no_gfni m)
  add_test(NAME fast_bench_quick_no_gfni COMMAND bench_no_gfni ${BENCH_DIR}/quick.txt)
  add_test(NAME fast_bench_quick_no_gfni_wide3 COMMAND bench_no_gfni --min-words 3 ${BENCH_DIR}/quick.txt)
  add_test(NAME fast_bench_quick_no_gfni_no_support COMMAND bench_no_gfni --no-support ${BENCH_DIR}/quick.txt)
  # the fallbacks on the carried path with 7 words, and with the pretest
  # off
  add_test(NAME fast_bench_quick_no_gfni_wide7 COMMAND bench_no_gfni --min-words 7 ${BENCH_DIR}/quick.txt)
  add_test(NAME fast_bench_quick_no_gfni_no_pretest COMMAND bench_no_gfni --pretest-min 0 ${BENCH_DIR}/quick.txt)
endif()
# builds for CPUs without AVX-512BW (bit-sliced counters, 8 slices with
# AVX-512F) and without AVX-512 (plain C filters, 4 slices), which the
# native build never compiles on an AVX-512BW machine
foreach(variant "no_bw:-mno-avx512bw" "portable:-mno-avx512f")
  string(REPLACE ":" ";" parts ${variant})
  list(GET parts 0 vname)
  list(GET parts 1 vflag)
  check_c_compiler_flag(${vflag} HAVE_FLAG_${vname})
  if (HAVE_FLAG_${vname})
    add_executable(bench_${vname} bench.c enumerate.c arrange.c square.c)
    target_compile_options(bench_${vname} PRIVATE ${vflag})
    target_link_libraries(bench_${vname} m)
    add_test(NAME fast_bench_quick_${vname} COMMAND bench_${vname} ${BENCH_DIR}/quick.txt)
    add_test(NAME fast_bench_quick_${vname}_wide8 COMMAND bench_${vname} --min-words 8 ${BENCH_DIR}/quick.txt)
    add_test(NAME fast_bench_quick_${vname}_no_support COMMAND bench_${vname} --no-support ${BENCH_DIR}/quick.txt)
  endif()
endforeach()
# the generic bit-sliced counters (more than 512 labels), on the smaller
# instances
add_test(NAME fast_bench_quick_wide16 COMMAND bench --min-words 16 --only "S=3" ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_no_fc COMMAND bench --no-fc --no-mrv ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_no_support COMMAND bench --no-support ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_no_cross COMMAND bench --no-cross ${BENCH_DIR}/quick.txt)
# cross support in the portable code (not used by default on the matrix
# path with 8 words)
add_test(NAME fast_bench_quick_cross_wide8 COMMAND bench_matrices --cross --min-words 8 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_matrices_no_cross COMMAND bench_matrices --no-cross ${BENCH_DIR}/quick.txt)
# msearch finds Morgenstern's square (P = 2^10 3^4 5^3 7^2, S = 327)
add_test(NAME fast_msearch_327 COMMAND sh -c
  "$<TARGET_FILE:msearch> --sums 327 10 4 3 2 | grep -c '\"type\":\"square\"' | grep -qx 1")
# differential test against a brute-force oracle on random synthetic
# families (n = 3..8, all label widths, every option set), in the native
# build and without AVX-512
add_executable(fuzz_arrange fuzz_arrange.c enumerate.c arrange.c square.c dfirst.c)
target_link_libraries(fuzz_arrange m)
add_test(NAME fast_fuzz_arrange COMMAND fuzz_arrange 1 150
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# 257-500 labels: the carried path with 5-8 words, its r1 at 2-8 words
add_test(NAME fast_fuzz_arrange_wide COMMAND fuzz_arrange --mode 7 1 150
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if (HAVE_FLAG_portable)
  add_executable(fuzz_arrange_portable fuzz_arrange.c enumerate.c arrange.c square.c dfirst.c)
  target_compile_options(fuzz_arrange_portable PRIVATE -mno-avx512f)
  target_link_libraries(fuzz_arrange_portable m)
  add_test(NAME fast_fuzz_arrange_portable COMMAND fuzz_arrange_portable 1000 60
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  # the d-first search with the plain C filters
  add_test(NAME fast_fuzz_dfirst_portable COMMAND fuzz_arrange_portable --dfirst 7000 40
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
if (HAVE_NO_GFNI_FLAGS)
  # the carried path with 5-8 words (and the pretest) with the fallbacks
  # for AVX-512 CPUs without GFNI / VPOPCNTDQ
  add_executable(fuzz_arrange_no_gfni fuzz_arrange.c enumerate.c arrange.c square.c dfirst.c)
  target_compile_options(fuzz_arrange_no_gfni PRIVATE -mno-gfni -mno-avx512vpopcntdq)
  target_link_libraries(fuzz_arrange_no_gfni m)
  add_test(NAME fast_fuzz_arrange_no_gfni_wide COMMAND fuzz_arrange_no_gfni --mode 7 3000 60
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME fast_fuzz_dfirst_no_gfni COMMAND fuzz_arrange_no_gfni --dfirst 3100 40
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
add_executable(fuzz_arrange_matrices fuzz_arrange.c enumerate.c arrange.c square.c dfirst.c)
target_compile_definitions(fuzz_arrange_matrices PRIVATE CARRY_MAX_W=0)
target_link_libraries(fuzz_arrange_matrices m)
add_test(NAME fast_fuzz_arrange_matrices COMMAND fuzz_arrange_matrices 2000 60
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME fast_fuzz_arrange_matrices_wide COMMAND fuzz_arrange_matrices --mode 7 2000 40
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# the diagonal-first search (dfirst.c) against brute force: every (square,
# traversal in a random diagonal set) pair, the partner flag, a planted
# diagonal + anti-diagonal pair, d-stride splits; with the carried lists
# and with the matrices
add_test(NAME fast_fuzz_dfirst COMMAND fuzz_arrange --dfirst 1 150
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME fast_fuzz_dfirst_matrices COMMAND fuzz_arrange_matrices --dfirst 5000 60
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# the class support of the V_d searches through the API with the plain
# root, whose r1 runs can be narrower than the class labels; built with
# ARRANGE_DEBUG (the class blocks are checked, and the words a narrower run
# must not read are poisoned)
add_executable(class_support_test class_support_test.c enumerate.c arrange.c square.c)
target_compile_definitions(class_support_test PRIVATE ARRANGE_DEBUG)
target_link_libraries(class_support_test m)
add_test(NAME fast_class_support_plain_root COMMAND class_support_test 1 400)
# msearch --diag-first finds the one (square, SP traversal) pair of
# P = 2^16 3^5 5^4 7^2, S = 849 (as bin/dsearch does)
add_test(NAME fast_msearch_dfirst_849 COMMAND sh -c
  "$<TARGET_FILE:msearch> --diag-first --diag-first-min-n 0 --sums 849 16 5 4 2 | grep -c '\"type\":\"dsquare\"' | grep -qx 1")
# ... also without the class support of the V_d searches
add_test(NAME fast_msearch_dfirst_849_no_class_support COMMAND sh -c
  "$<TARGET_FILE:msearch> --diag-first --diag-first-min-n 0 --no-class-support --sums 849 16 5 4 2 | grep -c '\"type\":\"dsquare\"' | grep -qx 1")
# msearch ignores bench's SAMPLE_* environment (a leftover export must not
# turn scheduler units into sampled runs), and rejects malformed options
add_test(NAME fast_msearch_no_sample_env COMMAND sh -c
  "SAMPLE_STRIDE=3 SAMPLE_LOG=/dev/null $<TARGET_FILE:msearch> --sums 327 10 4 3 2 | grep -c '\"type\":\"square\"' | grep -qx 1")
add_test(NAME fast_msearch_bad_args COMMAND sh -c
  "m=$<TARGET_FILE:msearch>; for a in '--diag-first --d-range 5:3' '--diag-first --d-range x' '--d-stride 3' '--calib-r1-stride 2' '--r1-strata 1,2,3,4,5,6,7,8,9' '--min-sum 3x' '--diag-first --diag-first-min-n abc' '--diag-first --d-log /nonexistent/x' '--no-class-support'; do $m $a --sums 327 10 4 3 2 > /dev/null 2>&1; [ $? -eq 2 ] || { echo \"accepted: $a\"; exit 1; }; done")
# saturated byte counters at the node below r1 (every unmatched cell in
# >= 255 candidates, so the cell is chosen among the saturated counts; the
# oracle is too slow there, so the variants are compared with each other,
# with MRV off and with bit-sliced counters, and must find the planted
# square), with carried bitsets and with the matrices
add_test(NAME fast_fuzz_arrange_saturated COMMAND fuzz_arrange --mode 6
  --noracle --variants default,nomrv,w3,w4,w8 1 2
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME fast_fuzz_arrange_matrices_saturated COMMAND fuzz_arrange_matrices
  --mode 6 --noracle --variants default,nomrv,w4,w8 3 1
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# n = 8 (SQ_MAX_N) searches run (they used to overflow the per-depth state)
add_test(NAME fast_msearch_n8 COMMAND sh -c
  "$<TARGET_FILE:msearch> --vec-size 8 --min-sum 320 --max-sum 325 12 6 3 2 1 1 | grep -c '\"type\":\"sum\"' | grep -qx 6")
# the new enumerator writes exactly what the old one did
add_test(NAME fast_enumerate_matches_cpp COMMAND sh -c
  "$<TARGET_FILE:enumeration> old_10_4_3_2.txt 10 4 3 2 > /dev/null && $<TARGET_FILE:enumerate> --file new_10_4_3_2.txt 10 4 3 2 2> /dev/null && cmp old_10_4_3_2.txt new_10_4_3_2.txt"
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
# the scheduler (needs python3 with numpy; skipped if python3 is missing)
find_package(Python3 COMPONENTS Interpreter)
if (Python3_FOUND)
  add_test(NAME fast_scheduler COMMAND ${Python3_EXECUTABLE}
    ${PROJECT_SOURCE_DIR}/scripts/test_scheduler.py)
  # the calibration ladder (scripts/calibrate.py, research/calibration.md)
  add_test(NAME fast_calibrate COMMAND ${Python3_EXECUTABLE}
    ${PROJECT_SOURCE_DIR}/scripts/test_calibrate.py)
endif()
