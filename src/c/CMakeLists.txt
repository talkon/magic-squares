add_library(bitset STATIC bitset.c search.c labelling.c)
add_executable(arrangement arrangement.c)
target_link_libraries(arrangement bitset)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
add_executable(enumeration enumeration.cpp)
