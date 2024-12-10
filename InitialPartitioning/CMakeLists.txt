add_library(InitialPartitioning
        src/GreedyDirectedGraphGrowing.cpp
)

target_include_directories(InitialPartitioning
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(InitialPartitioning
        PUBLIC Graph
)