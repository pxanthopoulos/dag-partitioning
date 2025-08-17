# GKlib dependency configuration

# GKlib dependency
FetchContent_Declare(
    GKlib
    GIT_REPOSITORY https://github.com/pxanthopoulos/GKlib.git
    GIT_TAG        master
)

# Configure GKlib build
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
    # Convert CMAKE bool to 0/1 for make config
    if(BUILD_SHARED_LIBS)
        set(GKLIB_SHARED_VALUE 1)
    else()
        set(GKLIB_SHARED_VALUE 0)
    endif()
    
    execute_process(
        COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${GKLIB_C_FLAGS}
                make config cc=${CMAKE_C_COMPILER} prefix=${GKLIB_PREFIX} shared=${GKLIB_SHARED_VALUE}
                gdb=${GKLIB_GDB} assert=${GKLIB_ASSERT} assert2=${GKLIB_ASSERT2}
                debug=${GKLIB_DEBUG} gprof=${GKLIB_GPROF} valgrind=${GKLIB_VALGRIND}
                openmp=${GKLIB_OPENMP} pcre=${GKLIB_PCRE} gkregex=${GKLIB_GKREGEX} gkrand=${GKLIB_GKRAND}
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

# Set library file based on BUILD_SHARED_LIBS
if(BUILD_SHARED_LIBS)
    set(GKLIB_LIB_FILE ${GKLIB_PREFIX}/lib/libGKlib${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(GKLIB_LIB_FILE ${GKLIB_PREFIX}/lib/libGKlib${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

# Use $(MAKE) to inherit jobserver from parent make, or use cmake --build with inherited settings
if(CMAKE_GENERATOR MATCHES "Make")
    # For Makefiles, use $(MAKE) to inherit the jobserver
    add_custom_command(
        OUTPUT ${GKLIB_LIB_FILE}
        COMMAND $(MAKE) -C ${GKLIB_BUILDDIR}
        COMMAND $(MAKE) -C ${GKLIB_BUILDDIR} install
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        COMMENT "Building GKlib"
        USES_TERMINAL
    )
else()
    # For other generators (Ninja, etc.), use CMAKE_BUILD_PARALLEL_LEVEL
    add_custom_command(
        OUTPUT ${GKLIB_LIB_FILE}
        COMMAND ${CMAKE_COMMAND} -E env CMAKE_BUILD_PARALLEL_LEVEL=$ENV{CMAKE_BUILD_PARALLEL_LEVEL} 
                ${CMAKE_COMMAND} --build ${GKLIB_BUILDDIR} --parallel
        COMMAND ${CMAKE_COMMAND} --install ${GKLIB_BUILDDIR}
        WORKING_DIRECTORY ${gklib_SOURCE_DIR}
        COMMENT "Building GKlib"
        USES_TERMINAL
    )
endif()

# Create directories for imported targets (needed at configure time)
file(MAKE_DIRECTORY ${GKLIB_PREFIX}/include)
file(MAKE_DIRECTORY ${GKLIB_PREFIX}/lib)

# Create custom targets for building dependencies
add_custom_target(build_gklib
    DEPENDS ${GKLIB_LIB_FILE}
)

# Create imported targets
if(BUILD_SHARED_LIBS)
    add_library(GKlib::GKlib SHARED IMPORTED)
else()
    add_library(GKlib::GKlib STATIC IMPORTED)
endif()
set_target_properties(GKlib::GKlib PROPERTIES
    IMPORTED_LOCATION ${GKLIB_LIB_FILE}
    INTERFACE_INCLUDE_DIRECTORIES ${GKLIB_PREFIX}/include
)
add_dependencies(GKlib::GKlib build_gklib)