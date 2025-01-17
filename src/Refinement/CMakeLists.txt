add_library(Refinement OBJECT
        src/Refinement.cpp
        src/BoundaryFM.cpp
        src/BoundaryKL.cpp
        src/Mixed.cpp
        src/RefinementWrapper.cpp
)

target_include_directories(Refinement
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Refinement
        PUBLIC Graph
)