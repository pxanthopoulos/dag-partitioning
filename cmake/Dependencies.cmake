# External dependencies for dag_partitioning

include(FetchContent)
include(ProcessorCount)

# Detect default parallel job count
ProcessorCount(NCPU)
if(NOT NCPU)
    set(NCPU 1)
endif()

# GKlib dependency
FetchContent_Declare(
    GKlib
    GIT_REPOSITORY https://github.com/pxanthopoulos/GKlib.git
    GIT_TAG        master
)

# Configure GKlib build
set(SHARED $<IF:$<BOOL:${BUILD_SHARED_LIBS}>,1,0>)
set(GKLIB_PREFIX ${CMAKE_BINARY_DIR}/_deps/gklib-build)

# Build GKlib C flags (PIC + custom flags)
set(GKLIB_C_FLAGS "${CMAKE_C_FLAGS}")
if(ENABLE_PIC)
    set(GKLIB_C_FLAGS "${GKLIB_C_FLAGS} -fPIC")
endif()
if(GKLIB_EXTRA_C_FLAGS)
    set(GKLIB_C_FLAGS "${GKLIB_C_FLAGS} ${GKLIB_EXTRA_C_FLAGS}")
endif()
message(STATUS "GKlib C flags: ${GKLIB_C_FLAGS}")

FetchContent_GetProperties(GKlib)
if(NOT gklib_POPULATED)
    FetchContent_Populate(GKlib)
    
    # Configure GKlib during CMake configure phase
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${GKLIB_C_FLAGS}
                make config cc=${CMAKE_C_COMPILER} prefix=${GKLIB_PREFIX} shared=${SHARED}
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        RESULT_VARIABLE GKLIB_CONFIG_RESULT
    )
    
    if(NOT GKLIB_CONFIG_RESULT EQUAL 0)
        message(FATAL_ERROR "GKlib configuration failed")
    endif()
endif()

# Build GKlib during build phase (respects -j flag)
# Emulate GKlib Makefile logic: BUILDDIR = build/$(systype)-$(cputype)
set(GKLIB_BUILDDIR ${gklib_SOURCE_DIR}/build/${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})

# Use $(MAKE) to inherit jobserver from parent make, or use cmake --build with inherited settings
if(CMAKE_GENERATOR MATCHES "Make")
    # For Makefiles, use $(MAKE) to inherit the jobserver
    add_custom_command(
        OUTPUT ${GKLIB_PREFIX}/lib/libGKlib.a
        COMMAND $(MAKE) -C ${GKLIB_BUILDDIR}
        COMMAND $(MAKE) -C ${GKLIB_BUILDDIR} install
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        COMMENT "Building GKlib"
        USES_TERMINAL
    )
else()
    # For other generators (Ninja, etc.), use CMAKE_BUILD_PARALLEL_LEVEL
    add_custom_command(
        OUTPUT ${GKLIB_PREFIX}/lib/libGKlib.a
        COMMAND ${CMAKE_COMMAND} -E env CMAKE_BUILD_PARALLEL_LEVEL=$ENV{CMAKE_BUILD_PARALLEL_LEVEL} 
                ${CMAKE_COMMAND} --build ${GKLIB_BUILDDIR} --parallel
        COMMAND ${CMAKE_COMMAND} --install ${GKLIB_BUILDDIR}
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        COMMENT "Building GKlib"
        USES_TERMINAL
    )
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

# Build METIS C flags (PIC + custom flags)
set(METIS_C_FLAGS "${CMAKE_C_FLAGS}")
if(ENABLE_PIC)
    set(METIS_C_FLAGS "${METIS_C_FLAGS} -fPIC")
endif()
if(METIS_EXTRA_C_FLAGS)
    set(METIS_C_FLAGS "${METIS_C_FLAGS} ${METIS_EXTRA_C_FLAGS}")
endif()
message(STATUS "METIS C flags: ${METIS_C_FLAGS}")

FetchContent_GetProperties(METIS)
if(NOT metis_POPULATED)
    FetchContent_Populate(METIS)
    
    # Configure METIS during CMake configure phase
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${METIS_C_FLAGS}
                make config shared=${SHARED} cc=${CMAKE_C_COMPILER} prefix=${METIS_PREFIX} i64=${METIS_I64} r64=${METIS_R64} gklib_path=${GKLIB_PREFIX}
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        RESULT_VARIABLE METIS_CONFIG_RESULT
    )
    
    if(NOT METIS_CONFIG_RESULT EQUAL 0)
        message(FATAL_ERROR "METIS configuration failed")
    endif()
endif()

# Build METIS during build phase (respects -j flag)
if(CMAKE_GENERATOR MATCHES "Make")
    # For Makefiles, use $(MAKE) to inherit the jobserver
    add_custom_command(
        OUTPUT ${METIS_PREFIX}/lib/libmetis.a
        COMMAND $(MAKE) -C ${metis_SOURCE_DIR}/build
        COMMAND $(MAKE) -C ${metis_SOURCE_DIR}/build install
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        DEPENDS ${GKLIB_PREFIX}/lib/libGKlib.a
        COMMENT "Building METIS"
        USES_TERMINAL
    )
else()
    # For other generators (Ninja, etc.), use CMAKE_BUILD_PARALLEL_LEVEL
    add_custom_command(
        OUTPUT ${METIS_PREFIX}/lib/libmetis.a
        COMMAND ${CMAKE_COMMAND} -E env CMAKE_BUILD_PARALLEL_LEVEL=$ENV{CMAKE_BUILD_PARALLEL_LEVEL}
                ${CMAKE_COMMAND} --build ${metis_SOURCE_DIR}/build --parallel
        COMMAND ${CMAKE_COMMAND} --install ${metis_SOURCE_DIR}/build
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        DEPENDS ${GKLIB_PREFIX}/lib/libGKlib.a
        COMMENT "Building METIS"
        USES_TERMINAL
    )
endif()

# Create directories for imported targets (needed at configure time)
file(MAKE_DIRECTORY ${GKLIB_PREFIX}/include)
file(MAKE_DIRECTORY ${GKLIB_PREFIX}/lib)
file(MAKE_DIRECTORY ${METIS_PREFIX}/include)
file(MAKE_DIRECTORY ${METIS_PREFIX}/lib)

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