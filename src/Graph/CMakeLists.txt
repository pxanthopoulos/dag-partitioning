add_library(Graph OBJECT
        src/Graph.cpp
)

target_include_directories(Graph
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)