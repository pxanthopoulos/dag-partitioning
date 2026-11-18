include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.15.2.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

add_executable(unit-tests
        unit/BruteForceSolverTest.cpp
        unit/CPSATMemorySolverTest.cpp
        unit/CPSATTransferMemorySolverTest.cpp
        unit/DotFileTest.cpp
        unit/PartitionFootprintsTest.cpp
        unit/RPOSchedulerTest.cpp
        unit/SchedulerTest.cpp
        unit/UtilsTest.cpp
)
target_link_libraries(unit-tests PRIVATE dag_partitioning GTest::gtest_main)
gtest_discover_tests(unit-tests)

add_executable(scheduler-fuzz e2e/scheduler-fuzz.cpp)
target_link_libraries(scheduler-fuzz PRIVATE dag_partitioning)
add_test(NAME scheduler-fuzz COMMAND scheduler-fuzz 1 300)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME fuzz
        COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/e2e/fuzz.py
                $<TARGET_FILE_DIR:dag-test> ${CMAKE_CURRENT_BINARY_DIR}/fuzz 20)
add_test(NAME fuzz2
        COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/e2e/fuzz2.py
                $<TARGET_FILE_DIR:dag-test> ${CMAKE_CURRENT_BINARY_DIR}/fuzz2 40)
