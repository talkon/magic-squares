add_library(bitset OBJECT bitset.c search.c labelling.c)
add_executable(arrangement arrangement.c)
target_link_libraries(arrangement bitset)
