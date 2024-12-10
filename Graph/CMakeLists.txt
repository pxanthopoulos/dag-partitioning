add_library(Graph
        src/Graph.cpp
)

target_include_directories(Graph
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)