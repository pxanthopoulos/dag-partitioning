add_library(Bisection
        src/GreedyDirectedGraphGrowing.cpp
        src/Bisection.cpp
        src/UndirectedFix.cpp
)

target_include_directories(Bisection
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Bisection
        PUBLIC Graph
        ${CMAKE_BINARY_DIR}/third-party/scotch/lib/libscotch.a
        ${CMAKE_BINARY_DIR}/third-party/scotch/lib/libscotcherr.a
        ${CMAKE_BINARY_DIR}/third-party/metis/lib/libmetis.a
        ${CMAKE_BINARY_DIR}/third-party/gklib/lib/libGKlib.a
)
add_dependencies(Bisection scotch metis)
