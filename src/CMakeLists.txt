add_library(partitioner STATIC)

add_subdirectory(Graph)
add_subdirectory(Clustering)
add_subdirectory(Bisection)
add_subdirectory(Driver)
add_subdirectory(Refinement)

target_link_libraries(partitioner
        PUBLIC
        Clustering
        Graph
        Bisection
        Driver
        Refinement
)

target_include_directories(partitioner
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/Graph/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Clustering/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Bisection/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Driver/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Refinement/include
)