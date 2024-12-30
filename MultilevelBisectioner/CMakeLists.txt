add_library(MultilevelBisectioner
        src/MultilevelBisectioner.cpp
)

target_include_directories(MultilevelBisectioner
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(MultilevelBisectioner
        PUBLIC Bisection Clustering Graph Refinement
)