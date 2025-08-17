add_library(Driver OBJECT
        src/MultilevelBisectioner.cpp
        src/RecursivePartitioner.cpp
)

target_include_directories(Driver
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(Driver
        PRIVATE
        Bisection
        Clustering
        Graph
        Refinement
)