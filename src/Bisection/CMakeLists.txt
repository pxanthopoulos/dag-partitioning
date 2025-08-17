add_library(Bisection OBJECT
        src/GreedyDirectedGraphGrowing.cpp
        src/Bisection.cpp
        src/UndirectedFix.cpp
)

target_include_directories(Bisection
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(Bisection
        PUBLIC
        METIS::METIS
        SCOTCH::SCOTCH
        PRIVATE 
        Graph 
        Refinement
)