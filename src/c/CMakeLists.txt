add_library(bitset OBJECT bitset.c)
add_executable(arrangement arrangement.c)
target_link_libraries(arrangement bitset)
