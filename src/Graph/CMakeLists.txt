add_library(Graph OBJECT
        src/Graph.cpp
)

target_include_directories(Graph
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)