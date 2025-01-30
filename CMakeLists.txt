cmake_minimum_required(VERSION 3.29)
project(dag_partitioning)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_library(SCOTCH_LIBRARY scotch HINTS "${CMAKE_BINARY_DIR}/third_party/scotch/lib" REQUIRED)
find_library(SCOTCHERR_LIBRARY scotcherr HINTS "${CMAKE_BINARY_DIR}/third_party/scotch/lib" REQUIRED)
find_library(METIS_LIBRARY metis HINTS "${CMAKE_BINARY_DIR}/third_party/metis/lib" REQUIRED)
find_library(GKLIB_LIBRARY GKlib HINTS "${CMAKE_BINARY_DIR}/third_party/gklib/lib" REQUIRED)

find_path(SCOTCH_INCLUDE_DIR scotch.h HINTS "${CMAKE_BINARY_DIR}/third_party/scotch/include" REQUIRED)
find_path(METIS_INCLUDE_DIR metis.h HINTS "${CMAKE_BINARY_DIR}/third_party/metis/include" REQUIRED)

include_directories(${SCOTCH_INCLUDE_DIR})
include_directories(${METIS_INCLUDE_DIR})

add_subdirectory(src)