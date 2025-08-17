# METIS dependency configuration

# METIS dependency
FetchContent_Declare(
    METIS
    GIT_REPOSITORY https://github.com/pxanthopoulos/METIS.git
    GIT_TAG        master
)

# Configure METIS build
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
    # Convert CMAKE bool to 0/1 for make config
    if(BUILD_SHARED_LIBS)
        set(METIS_SHARED_VALUE 1)
    else()
        set(METIS_SHARED_VALUE 0)
    endif()
    
    if(DAG_PARTITIONING_OPENMP)
        set(METIS_OPENMP_VALUE 1)
    else()
        set(METIS_OPENMP_VALUE 0)
    endif()
    
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${METIS_C_FLAGS}
                make config shared=${METIS_SHARED_VALUE} cc=${CMAKE_C_COMPILER} prefix=${METIS_PREFIX} 
                i64=${METIS_USE_I64} r64=${METIS_USE_R64} 
                gdb=${METIS_GDB} assert=${METIS_ASSERT} assert2=${METIS_ASSERT2}
                debug=${METIS_DEBUG} gprof=${METIS_GPROF} valgrind=${METIS_VALGRIND}
                openmp=${METIS_OPENMP_VALUE} 
                gklib_path=${GKLIB_PREFIX}
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        RESULT_VARIABLE METIS_CONFIG_RESULT
    )
    
    if(NOT METIS_CONFIG_RESULT EQUAL 0)
        message(FATAL_ERROR "METIS configuration failed")
    endif()
endif()

# Set library file based on BUILD_SHARED_LIBS
if(BUILD_SHARED_LIBS)
    set(METIS_LIB_FILE ${METIS_PREFIX}/lib/libmetis${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(METIS_LIB_FILE ${METIS_PREFIX}/lib/libmetis${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

# Build METIS during build phase (respects -j flag)
if(CMAKE_GENERATOR MATCHES "Make")
    # For Makefiles, use $(MAKE) to inherit the jobserver
    add_custom_command(
        OUTPUT ${METIS_LIB_FILE}
        COMMAND $(MAKE) -C ${metis_SOURCE_DIR}/build
        COMMAND $(MAKE) -C ${metis_SOURCE_DIR}/build install
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        DEPENDS ${GKLIB_LIB_FILE}
        COMMENT "Building METIS"
        USES_TERMINAL
    )
else()
    # For other generators (Ninja, etc.), use CMAKE_BUILD_PARALLEL_LEVEL
    add_custom_command(
        OUTPUT ${METIS_LIB_FILE}
        COMMAND ${CMAKE_COMMAND} -E env CMAKE_BUILD_PARALLEL_LEVEL=$ENV{CMAKE_BUILD_PARALLEL_LEVEL}
                ${CMAKE_COMMAND} --build ${metis_SOURCE_DIR}/build --parallel
        COMMAND ${CMAKE_COMMAND} --install ${metis_SOURCE_DIR}/build
        WORKING_DIRECTORY ${metis_SOURCE_DIR}
        DEPENDS ${GKLIB_LIB_FILE}
        COMMENT "Building METIS"
        USES_TERMINAL
    )
endif()

# Create directories for imported targets (needed at configure time)
file(MAKE_DIRECTORY ${METIS_PREFIX}/include)
file(MAKE_DIRECTORY ${METIS_PREFIX}/lib)

# Create custom targets for building dependencies
add_custom_target(build_metis
    DEPENDS ${METIS_LIB_FILE}
    DEPENDS build_gklib
)

# Create imported targets
if(BUILD_SHARED_LIBS)
    add_library(METIS::METIS SHARED IMPORTED)
else()
    add_library(METIS::METIS STATIC IMPORTED)
endif()
set_target_properties(METIS::METIS PROPERTIES
    IMPORTED_LOCATION ${METIS_LIB_FILE}
    INTERFACE_INCLUDE_DIRECTORIES ${METIS_PREFIX}/include
    INTERFACE_LINK_LIBRARIES GKlib::GKlib
)
add_dependencies(METIS::METIS build_metis)