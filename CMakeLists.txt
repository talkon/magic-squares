cmake_minimum_required (VERSION 3.13)
project (MAGIC)

add_compile_options(-O3 -flto)
add_link_options(-flto)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory (src)