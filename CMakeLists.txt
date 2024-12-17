cmake_minimum_required(VERSION 3.29)
project(dag_partitioning)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O3")

find_package(LLVM REQUIRED CONFIG)
include_directories(${LLVM_INCLUDE_DIRS})
add_definitions(${LLVM_DEFINITIONS})
llvm_map_components_to_libnames(llvm_libs Support Core)

add_subdirectory(Graph)
add_subdirectory(Clustering)
add_subdirectory(Bisection)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE Clustering Graph Bisection ${llvm_libs})