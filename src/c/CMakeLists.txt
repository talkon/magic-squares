# New pipeline: portable C, no AVX-512 required.
add_library(magic_core OBJECT enumerate.c arrange.c square.c)
add_executable(enumerate enumerate_main.c $<TARGET_OBJECTS:magic_core>)
add_executable(msearch msearch.c $<TARGET_OBJECTS:magic_core>)
add_executable(bench bench.c $<TARGET_OBJECTS:magic_core>)
foreach(t enumerate msearch bench)
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
