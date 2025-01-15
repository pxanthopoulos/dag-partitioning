cmake_minimum_required(VERSION 3.29)
project(dag_partitioning)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package(LLVM REQUIRED CONFIG)
include_directories(${LLVM_INCLUDE_DIRS})
add_definitions(${LLVM_DEFINITIONS})
llvm_map_components_to_libnames(llvm_libs Support Core)

include(ExternalProject)
ExternalProject_Add(scotch
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/third-party/scotch
        PREFIX ${CMAKE_BINARY_DIR}/third-party/scotch
        CONFIGURE_COMMAND
        cmake -G ${CMAKE_GENERATOR}
        -B <BINARY_DIR>
        -S <SOURCE_DIR>
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER}
        -DBUILD_SHARED_LIBS=OFF
        -DBUILD_PTSCOTCH=OFF
        -DBUILD_FORTRAN=OFF
        -DBUILD_LIBSCOTCHMETIS=OFF
        -DSCOTCH_METIS_PREFIX=OFF
        -DINTSIZE=64
        -DIDXSIZE64=ON
        -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        BUILD_COMMAND
        cmake --build <BINARY_DIR>
        INSTALL_COMMAND
        cmake --install <BINARY_DIR>
        BUILD_IN_SOURCE 0
        BUILD_BYPRODUCTS
        <INSTALL_DIR>/lib/libscotch.a
        <INSTALL_DIR>/lib/libscotcherr.a
)
include_directories(${CMAKE_BINARY_DIR}/third-party/scotch/include)

ExternalProject_Add(gklib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/third-party/gklib
        PREFIX ${CMAKE_BINARY_DIR}/third-party/gklib
        CONFIGURE_COMMAND
        cd <SOURCE_DIR> && make config cc=${CMAKE_C_COMPILER}
        prefix=<INSTALL_DIR>
        openmp=1
        BUILD_COMMAND
        cd <SOURCE_DIR> && make
        INSTALL_COMMAND
        cd <SOURCE_DIR> && make install
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS
        <INSTALL_DIR>/lib/libGKlib.a
)

ExternalProject_Add(metis
        DEPENDS gklib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/third-party/metis
        PREFIX ${CMAKE_BINARY_DIR}/third-party/metis
        CONFIGURE_COMMAND
        cd <SOURCE_DIR> && make config cc=${CMAKE_C_COMPILER}
        prefix=<INSTALL_DIR>
        gklib_path=${CMAKE_BINARY_DIR}/third-party/gklib
        i64=1
        r64=1
        gdb=1
        assert=1
        openmp=1
        BUILD_COMMAND
        cd <SOURCE_DIR> && make
        INSTALL_COMMAND
        cd <SOURCE_DIR> && make install
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS
        <INSTALL_DIR>/lib/libmetis.a
)
include_directories(${CMAKE_BINARY_DIR}/third-party/metis/include)

add_subdirectory(src)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE
        partitioner
        ${llvm_libs}
        m
)