add_library(Bisection
        src/GreedyDirectedGraphGrowing.cpp
        src/Bisection.cpp
)

target_include_directories(Bisection
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Bisection
        PUBLIC Graph
)