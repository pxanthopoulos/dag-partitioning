add_library(Bisection OBJECT
        src/GreedyDirectedGraphGrowing.cpp
        src/Bisection.cpp
        src/UndirectedFix.cpp
)

target_include_directories(Bisection
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        ${SCOTCH_INCLUDE_DIR}
)

target_link_libraries(Bisection
        PUBLIC
        METIS::METIS
        PRIVATE 
        Graph 
        Refinement 
        ${SCOTCH_LIBRARY}
        ${SCOTCHERR_LIBRARY}
)