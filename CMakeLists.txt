cmake_minimum_required(VERSION 3.28)
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

# Compiler configuration
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

# Additional CXX flags
if(EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif()

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# External dependencies
include(FetchContent)
include(ProcessorCount)
ProcessorCount(NCPU)

# GKlib dependency
FetchContent_Declare(
    GKlib
    GIT_REPOSITORY https://github.com/pxanthopoulos/GKlib.git
    GIT_TAG        master
)

# Configure GKlib build
set(SHARED $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,1,0>)
set(GKLIB_PREFIX ${CMAKE_BINARY_DIR}/_deps/gklib-build)

FetchContent_GetProperties(GKlib)
if(NOT gklib_POPULATED)
    FetchContent_Populate(GKlib)
    
    # Build GKlib using external command
    execute_process(
        COMMAND make config cc=${CMAKE_C_COMPILER} prefix=${GKLIB_PREFIX} shared=${SHARED}
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        RESULT_VARIABLE GKLIB_CONFIG_RESULT
    )
    
    if(NOT GKLIB_CONFIG_RESULT EQUAL 0)
        message(FATAL_ERROR "GKlib configuration failed")
    endif()
    
    execute_process(
        COMMAND make -j${NCPU}
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        RESULT_VARIABLE GKLIB_BUILD_RESULT
    )
    
    if(NOT GKLIB_BUILD_RESULT EQUAL 0)
        message(FATAL_ERROR "GKlib build failed")
    endif()
    
    execute_process(
        COMMAND make install
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        RESULT_VARIABLE GKLIB_INSTALL_RESULT
    )
    
    if(NOT GKLIB_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR "GKlib install failed")
    endif()
endif()

# METIS dependency
FetchContent_Declare(
    METIS
    GIT_REPOSITORY https://github.com/pxanthopoulos/METIS.git
    GIT_TAG        master
)

# Configure METIS build
set(METIS_I64 $<IF:$<BOOL:${METIS_USE_I64}>,1,0>)
set(METIS_R64 $<IF:$<BOOL:${METIS_USE_R64}>,1,0>)
set(METIS_PREFIX ${CMAKE_BINARY_DIR}/_deps/metis-build)

FetchContent_GetProperties(METIS)
if(NOT metis_POPULATED)
    FetchContent_Populate(METIS)
    
    # Build METIS using external command
    execute_process(
        COMMAND make config shared=${SHARED} cc=${CMAKE_C_COMPILER} prefix=${METIS_PREFIX} i64=${METIS_I64} r64=${METIS_R64} gklib_path=${GKLIB_PREFIX}
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        RESULT_VARIABLE METIS_CONFIG_RESULT
    )
    
    if(NOT METIS_CONFIG_RESULT EQUAL 0)
        message(FATAL_ERROR "METIS configuration failed")
    endif()
    
    execute_process(
        COMMAND make -j${NCPU}
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        RESULT_VARIABLE METIS_BUILD_RESULT
    )
    
    if(NOT METIS_BUILD_RESULT EQUAL 0)
        message(FATAL_ERROR "METIS build failed")
    endif()
    
    execute_process(
        COMMAND make install
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        RESULT_VARIABLE METIS_INSTALL_RESULT
    )
    
    if(NOT METIS_INSTALL_RESULT EQUAL 0)
        message(FATAL_ERROR "METIS install failed")
    endif()
endif()

# Create custom targets for building dependencies
add_custom_target(build_gklib
    DEPENDS ${GKLIB_PREFIX}/lib/libGKlib.a
)

add_custom_target(build_metis
    DEPENDS ${METIS_PREFIX}/lib/libmetis.a
    DEPENDS build_gklib
)

# Create imported targets
add_library(GKlib::GKlib STATIC IMPORTED)
set_target_properties(GKlib::GKlib PROPERTIES
    IMPORTED_LOCATION ${GKLIB_PREFIX}/lib/libGKlib.a
    INTERFACE_INCLUDE_DIRECTORIES ${GKLIB_PREFIX}/include
)
add_dependencies(GKlib::GKlib build_gklib)

add_library(METIS::METIS STATIC IMPORTED)
set_target_properties(METIS::METIS PROPERTIES
    IMPORTED_LOCATION ${METIS_PREFIX}/lib/libmetis.a
    INTERFACE_INCLUDE_DIRECTORIES ${METIS_PREFIX}/include
    INTERFACE_LINK_LIBRARIES GKlib::GKlib
)
add_dependencies(METIS::METIS build_metis)

# Add subdirectories
add_subdirectory(src)

# Add executables
add_executable(rand-dag test/rand-dag.cpp)

add_executable(test test/test.cpp)
target_link_libraries(test PRIVATE dag_partitioning)

# Install configuration
include(GNUInstallDirs)

install(TARGETS dag_partitioning rand-dag test
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)