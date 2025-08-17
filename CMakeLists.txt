cmake_minimum_required(VERSION 3.28)

# Compiler configuration (must be before project())
if(NOT CMAKE_CXX_COMPILER)
    find_program(CLANG_CXX_COMPILER clang++)
    if(CLANG_CXX_COMPILER)
        set(CMAKE_CXX_COMPILER ${CLANG_CXX_COMPILER})
    endif()
endif()

if(NOT CMAKE_C_COMPILER)
    find_program(CLANG_C_COMPILER clang)
    if(CLANG_C_COMPILER)
        set(CMAKE_C_COMPILER ${CLANG_C_COMPILER})
    endif()
endif()

project(dag_partitioning VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release build
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Options
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_PIC "Enable position independent code" ON)
option(METIS_USE_I64 "Use 64-bit indices for METIS" ON)
option(METIS_USE_R64 "Use 64-bit reals for METIS" ON)
set(EXTRA_CXX_FLAGS "" CACHE STRING "Additional CXX flags")
set(GKLIB_EXTRA_C_FLAGS "" CACHE STRING "Additional C flags for GKlib")
set(METIS_EXTRA_C_FLAGS "" CACHE STRING "Additional C flags for METIS")

# Compiler launcher (ccache)
if(NOT CMAKE_CXX_COMPILER_LAUNCHER)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE_PROGRAM})
    endif()
endif()

# Position independent code
if(ENABLE_PIC)
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Find required packages
find_package(Threads REQUIRED)

# External dependencies
include(cmake/Dependencies.cmake)

# Add subdirectories
add_subdirectory(src)

# Add executables
add_executable(rand-dag test/rand-dag.cpp)

add_executable(test test/test.cpp)
target_link_libraries(test PRIVATE dag_partitioning)

# Installation
include(cmake/Install.cmake)