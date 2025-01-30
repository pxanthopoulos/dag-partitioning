add_library(Bisection OBJECT
        src/GreedyDirectedGraphGrowing.cpp
        src/Bisection.cpp
        src/UndirectedFix.cpp
)

target_include_directories(Bisection
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${SCOTCH_INCLUDE_DIR}
        ${METIS_INCLUDE_DIR}
)

target_link_libraries(Bisection
        PUBLIC Graph Refinement
        ${SCOTCH_LIBRARY}
        ${SCOTCHERR_LIBRARY}
        ${METIS_LIBRARY}
)