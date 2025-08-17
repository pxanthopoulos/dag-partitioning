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
option(GKLIB_GDB "Enable GDB support for GKlib" OFF)
option(GKLIB_ASSERT "Turn GKlib asserts on" OFF)
option(GKLIB_ASSERT2 "Additional GKlib assertions" OFF)
option(GKLIB_DEBUG "Add GKlib debugging support" OFF)
option(GKLIB_GPROF "Add GKlib gprof support" OFF)
option(GKLIB_VALGRIND "Add GKlib valgrind support" OFF)
option(GKLIB_OPENMP "Enable GKlib OpenMP support" OFF)
option(GKLIB_PCRE "Enable GKlib PCRE support" OFF)
option(GKLIB_GKREGEX "Enable GKlib GKREGEX support" OFF)
option(GKLIB_GKRAND "Enable GKlib GKRAND support" OFF)
option(METIS_USE_I64 "Use 64-bit indices for METIS" ON)
option(METIS_USE_R64 "Use 64-bit reals for METIS" ON)
option(METIS_GDB "Enable GDB support for METIS" OFF)
option(METIS_ASSERT "Turn METIS asserts on" OFF)
option(METIS_ASSERT2 "Additional METIS assertions" OFF)
option(METIS_DEBUG "Add METIS debugging support" OFF)
option(METIS_GPROF "Add METIS gprof support" OFF)
option(METIS_VALGRIND "Add METIS valgrind support" OFF)
option(METIS_OPENMP "Enable METIS OpenMP support" OFF)
option(SCOTCH_USE_INTSIZE64 "Use 64-bit integers for SCOTCH" ON)
option(SCOTCH_BUILD_PTSCOTCH "Build PT-SCOTCH (parallel version)" OFF)
option(SCOTCH_BUILD_FORTRAN "Build Fortran interface for SCOTCH" OFF)
option(SCOTCH_BUILD_LIBSCOTCHMETIS "Build libscotchmetis compatibility library" OFF)
option(SCOTCH_METIS_PREFIX "Use METIS prefix for SCOTCH compatibility" OFF)
option(SCOTCH_INSTALL_METIS_HEADERS "Install libScotchMeTiS headers" ON)
option(SCOTCH_THREADS "Use multi-threading in Scotch and PT-Scotch" ON)
option(SCOTCH_MPI_THREAD_MULTIPLE "Use multi-threading with MPI in PT-Scotch" ON)
option(SCOTCH_BUILD_LIBESMUMPS "Build EsMUMPS library" ON)
option(SCOTCH_USE_ZLIB "Use ZLIB compression format if found" ON)
option(SCOTCH_USE_LZMA "Use LZMA compression format if found" ON)
option(SCOTCH_USE_BZ2 "Use BZ2 compression format if found" ON)
option(SCOTCH_ENABLE_TESTS "Enable SCOTCH testing" ON)
option(DAG_PARTITIONING_OPENMP "Enable OpenMP support" OFF)
set(SCOTCH_DETERMINISTIC "FIXED_SEED" CACHE STRING "Level of determinism in Scotch (NONE/FIXED_SEED/FULL)")
set(SCOTCH_IDXSIZE "" CACHE STRING "Index size: empty string stands for C compiler int type size")
set(SCOTCH_NAME_SUFFIX "" CACHE STRING "Suffix to add to every symbol in the Scotch/PT-Scotch libraries")
set(SCOTCH_LIBSCOTCHERR "" CACHE STRING "Which error library is linked by default")
set(SCOTCH_LIBPTSCOTCHERR "" CACHE STRING "Which parallel error library is linked by default")
set(EXTRA_CXX_FLAGS "" CACHE STRING "Additional CXX flags")
set(GKLIB_EXTRA_C_FLAGS "" CACHE STRING "Additional C flags for GKlib")
set(METIS_EXTRA_C_FLAGS "" CACHE STRING "Additional C flags for METIS")
set(SCOTCH_EXTRA_C_FLAGS "" CACHE STRING "Additional C flags for SCOTCH")

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

# Find OpenMP if enabled
if(DAG_PARTITIONING_OPENMP)
    find_package(OpenMP REQUIRED)
else()
    find_package(Threads REQUIRED)
endif()

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