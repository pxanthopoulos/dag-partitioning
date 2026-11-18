add_library(Scheduling OBJECT
        src/BruteForceSolver.cpp
        src/CPSAT.cpp
        src/CPSATMemorySolver.cpp
        src/CPSATTransferMemorySolver.cpp
        src/HyperGraph.cpp
        src/PartitionFootprints.cpp
        src/RPOScheduler.cpp
        src/Scheduler.cpp
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
