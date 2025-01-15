add_library(Clustering OBJECT
        src/ClusteringForbiddenEdges.cpp
        src/Clustering.cpp
        src/ClusteringCycleDetection.cpp
        src/ClusteringHybrid.cpp
)

target_include_directories(Clustering
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Clustering
        PUBLIC Graph
)