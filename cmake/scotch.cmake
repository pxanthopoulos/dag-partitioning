# SCOTCH dependency configuration

# SCOTCH dependency
FetchContent_Declare(
    SCOTCH
    URL https://gitlab.inria.fr/scotch/scotch/-/archive/master/scotch-master.tar.gz
)

# Configure SCOTCH build
if(SCOTCH_USE_INTSIZE64)
    set(SCOTCH_INTSIZE 64)
else()
    set(SCOTCH_INTSIZE 32)
endif()
set(SCOTCH_PREFIX ${CMAKE_BINARY_DIR}/_deps/scotch-build)

# Build SCOTCH C flags (PIC + custom flags)
set(SCOTCH_C_FLAGS "${CMAKE_C_FLAGS}")
if(ENABLE_PIC)
    set(SCOTCH_C_FLAGS "${SCOTCH_C_FLAGS} -fPIC")
endif()
if(SCOTCH_EXTRA_C_FLAGS)
    set(SCOTCH_C_FLAGS "${SCOTCH_C_FLAGS} ${SCOTCH_EXTRA_C_FLAGS}")
endif()
message(STATUS "SCOTCH C flags: ${SCOTCH_C_FLAGS}")

FetchContent_GetProperties(SCOTCH)
if(NOT scotch_POPULATED)
    FetchContent_Populate(SCOTCH)
    
    # Configure SCOTCH during CMake configure phase
    execute_process(
        COMMAND ${CMAKE_COMMAND} -S ${scotch_SOURCE_DIR} -B ${scotch_SOURCE_DIR}/build
                -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                -DCMAKE_C_FLAGS=${SCOTCH_C_FLAGS}
                -DCMAKE_INSTALL_PREFIX=${SCOTCH_PREFIX}
                -DINTSIZE=${SCOTCH_INTSIZE}
                -DIDXSIZE=${SCOTCH_IDXSIZE}
                -DBUILD_PTSCOTCH=${SCOTCH_BUILD_PTSCOTCH}
                -DBUILD_FORTRAN=${SCOTCH_BUILD_FORTRAN}
                -DBUILD_LIBSCOTCHMETIS=${SCOTCH_BUILD_LIBSCOTCHMETIS}
                -DSCOTCH_METIS_PREFIX=${SCOTCH_METIS_PREFIX}
                -DINSTALL_METIS_HEADERS=${SCOTCH_INSTALL_METIS_HEADERS}
                -DTHREADS=${SCOTCH_THREADS}
                -DMPI_THREAD_MULTIPLE=${SCOTCH_MPI_THREAD_MULTIPLE}
                -DBUILD_LIBESMUMPS=${SCOTCH_BUILD_LIBESMUMPS}
                -DUSE_ZLIB=${SCOTCH_USE_ZLIB}
                -DUSE_LZMA=${SCOTCH_USE_LZMA}
                -DUSE_BZ2=${SCOTCH_USE_BZ2}
                -DENABLE_TESTS=${SCOTCH_ENABLE_TESTS}
                -DSCOTCH_DETERMINISTIC=${SCOTCH_DETERMINISTIC}
                -DSCOTCH_NAME_SUFFIX=${SCOTCH_NAME_SUFFIX}
                -DLIBSCOTCHERR=${SCOTCH_LIBSCOTCHERR}
                -DLIBPTSCOTCHERR=${SCOTCH_LIBPTSCOTCHERR}
                -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
        RESULT_VARIABLE SCOTCH_CONFIG_RESULT
    )
    
    if(NOT SCOTCH_CONFIG_RESULT EQUAL 0)
        message(FATAL_ERROR "SCOTCH configuration failed")
    endif()
endif()

# Set library files based on BUILD_SHARED_LIBS
if(BUILD_SHARED_LIBS)
    set(SCOTCH_LIB_FILE ${SCOTCH_PREFIX}/lib/libscotch${CMAKE_SHARED_LIBRARY_SUFFIX})
    set(SCOTCHERR_LIB_FILE ${SCOTCH_PREFIX}/lib/libscotcherr${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
    set(SCOTCH_LIB_FILE ${SCOTCH_PREFIX}/lib/libscotch${CMAKE_STATIC_LIBRARY_SUFFIX})
    set(SCOTCHERR_LIB_FILE ${SCOTCH_PREFIX}/lib/libscotcherr${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

# Build SCOTCH during build phase (respects -j flag)
if(CMAKE_GENERATOR MATCHES "Make")
    # For Makefiles, use $(MAKE) to inherit the jobserver
    add_custom_command(
        OUTPUT ${SCOTCH_LIB_FILE} ${SCOTCHERR_LIB_FILE}
        COMMAND $(MAKE) -C ${scotch_SOURCE_DIR}/build
        COMMAND $(MAKE) -C ${scotch_SOURCE_DIR}/build install
        WORKING_DIRECTORY ${scotch_SOURCE_DIR}
        COMMENT "Building SCOTCH"
        USES_TERMINAL
    )
else()
    # For other generators (Ninja, etc.), use CMAKE_BUILD_PARALLEL_LEVEL
    add_custom_command(
        OUTPUT ${SCOTCH_LIB_FILE} ${SCOTCHERR_LIB_FILE}
        COMMAND ${CMAKE_COMMAND} --build ${scotch_SOURCE_DIR}/build --parallel
        COMMAND ${CMAKE_COMMAND} --install ${scotch_SOURCE_DIR}/build
        WORKING_DIRECTORY ${scotch_SOURCE_DIR}
        COMMENT "Building SCOTCH"
        USES_TERMINAL
    )
endif()

# Create directories for imported targets (needed at configure time)
file(MAKE_DIRECTORY ${SCOTCH_PREFIX}/include)
file(MAKE_DIRECTORY ${SCOTCH_PREFIX}/lib)

# Create custom targets for building dependencies
add_custom_target(build_scotch
    DEPENDS ${SCOTCH_LIB_FILE} ${SCOTCHERR_LIB_FILE}
)

# Create imported targets
if(BUILD_SHARED_LIBS)
    add_library(SCOTCH::scotcherr SHARED IMPORTED)
    add_library(SCOTCH::SCOTCH SHARED IMPORTED)
else()
    add_library(SCOTCH::scotcherr STATIC IMPORTED)
    add_library(SCOTCH::SCOTCH STATIC IMPORTED)
endif()

set_target_properties(SCOTCH::scotcherr PROPERTIES
    IMPORTED_LOCATION ${SCOTCHERR_LIB_FILE}
)
add_dependencies(SCOTCH::scotcherr build_scotch)

set_target_properties(SCOTCH::SCOTCH PROPERTIES
    IMPORTED_LOCATION ${SCOTCH_LIB_FILE}
    INTERFACE_INCLUDE_DIRECTORIES ${SCOTCH_PREFIX}/include
    INTERFACE_LINK_LIBRARIES SCOTCH::scotcherr
)
add_dependencies(SCOTCH::SCOTCH build_scotch)