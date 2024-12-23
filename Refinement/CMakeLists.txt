add_library(Refinement
        src/Refinement.cpp
)

target_include_directories(Refinement
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(Refinement
        PUBLIC Graph
)