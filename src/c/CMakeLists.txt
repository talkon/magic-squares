# New pipeline: portable C, no AVX-512 required.
add_library(magic_core OBJECT enumerate.c arrange.c square.c)
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
# the same instances, with wider label bitsets than needed
add_test(NAME fast_bench_quick_wide3 COMMAND bench --min-words 3 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_wide8 COMMAND bench --min-words 8 ${BENCH_DIR}/quick.txt)
# the same instances, filtering with the intersection matrices and gathers
# as for more than 256 labels and N > 3072
add_test(NAME fast_bench_quick_gather COMMAND bench --min-words 8 --gather ${BENCH_DIR}/quick.txt)
# the benchmark built with the N x N intersection matrices for all widths
# (the filter used without AVX-512BW, here still with byte counters), and
# with the fallbacks for AVX-512 CPUs without GFNI / VPOPCNTDQ (Skylake-X,
# Cascade Lake), to test those paths
add_executable(bench_matrices bench.c enumerate.c arrange.c square.c)
target_compile_definitions(bench_matrices PRIVATE CARRY_MAX_W=0)
target_link_libraries(bench_matrices m)
add_test(NAME fast_bench_quick_matrices COMMAND bench_matrices ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_matrices_no_support COMMAND bench_matrices --no-support ${BENCH_DIR}/quick.txt)
include(CheckCCompilerFlag)
check_c_compiler_flag("-mno-gfni -mno-avx512vpopcntdq" HAVE_NO_GFNI_FLAGS)
if (HAVE_NO_GFNI_FLAGS)
  add_executable(bench_no_gfni bench.c enumerate.c arrange.c square.c)
  target_compile_options(bench_no_gfni PRIVATE -mno-gfni -mno-avx512vpopcntdq)
  target_link_libraries(bench_no_gfni m)
  add_test(NAME fast_bench_quick_no_gfni COMMAND bench_no_gfni ${BENCH_DIR}/quick.txt)
  add_test(NAME fast_bench_quick_no_gfni_wide3 COMMAND bench_no_gfni --min-words 3 ${BENCH_DIR}/quick.txt)
  add_test(NAME fast_bench_quick_no_gfni_no_support COMMAND bench_no_gfni --no-support ${BENCH_DIR}/quick.txt)
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
# cross support in the portable code (not used by default with 8 words)
add_test(NAME fast_bench_quick_cross_wide8 COMMAND bench --cross --min-words 8 ${BENCH_DIR}/quick.txt)
add_test(NAME fast_bench_quick_matrices_no_cross COMMAND bench_matrices --no-cross ${BENCH_DIR}/quick.txt)
# msearch finds Morgenstern's square (P = 2^10 3^4 5^3 7^2, S = 327)
add_test(NAME fast_msearch_327 COMMAND sh -c
  "$<TARGET_FILE:msearch> --sums 327 10 4 3 2 | grep -c '\"type\":\"square\"' | grep -qx 1")
# differential test against a brute-force oracle on random synthetic
# families (n = 3..8, all label widths, every option set), in the native
# build and without AVX-512
add_executable(fuzz_arrange fuzz_arrange.c enumerate.c arrange.c square.c)
target_link_libraries(fuzz_arrange m)
add_test(NAME fast_fuzz_arrange COMMAND fuzz_arrange 1 150
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
if (HAVE_FLAG_portable)
  add_executable(fuzz_arrange_portable fuzz_arrange.c enumerate.c arrange.c square.c)
  target_compile_options(fuzz_arrange_portable PRIVATE -mno-avx512f)
  target_link_libraries(fuzz_arrange_portable m)
  add_test(NAME fast_fuzz_arrange_portable COMMAND fuzz_arrange_portable 1000 60
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
add_executable(fuzz_arrange_matrices fuzz_arrange.c enumerate.c arrange.c square.c)
target_compile_definitions(fuzz_arrange_matrices PRIVATE CARRY_MAX_W=0)
target_link_libraries(fuzz_arrange_matrices m)
add_test(NAME fast_fuzz_arrange_matrices COMMAND fuzz_arrange_matrices 2000 60
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
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
