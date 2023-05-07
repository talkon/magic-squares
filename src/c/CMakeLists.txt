add_library(bitset STATIC bitset.c search.c labelling.c)
add_executable(arrangement arrangement.c)
target_link_libraries(arrangement bitset)
