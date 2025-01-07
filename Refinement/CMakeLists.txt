add_library(Refinement
        src/Refinement.cpp
        src/BoundaryFM.cpp
        src/BoundaryKL.cpp
        src/Mixed.cpp
)

target_include_directories(Refinement
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Refinement
        PUBLIC Graph
)