add_library(dag_partitioning STATIC)

add_subdirectory(Graph)
add_subdirectory(Clustering)
add_subdirectory(Bisection)
add_subdirectory(Driver)
add_subdirectory(Refinement)

target_link_libraries(dag_partitioning
        PUBLIC
        Clustering
        Graph
        Bisection
        Driver
        Refinement
)

target_include_directories(dag_partitioning
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/Graph/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Clustering/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Bisection/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Driver/include
        ${CMAKE_CURRENT_SOURCE_DIR}/Refinement/include
)

install(TARGETS dag_partitioning
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        RUNTIME DESTINATION bin
)
