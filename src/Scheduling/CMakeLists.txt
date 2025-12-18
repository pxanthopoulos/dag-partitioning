add_library(Scheduling OBJECT
        src/Scheduling.cpp
        src/Utils.cpp
)

target_include_directories(Scheduling
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(Scheduling
        PRIVATE
        Graph
        robin_hood::robin_hood
        ortools::ortools
)