add_library(dag_partitioning)

add_subdirectory(Graph)
add_subdirectory(Clustering)
add_subdirectory(Bisection)
add_subdirectory(Driver)
add_subdirectory(Refinement)
add_subdirectory(Scheduling)

# Include all object library sources in the main library
target_sources(dag_partitioning PRIVATE
    $<TARGET_OBJECTS:Graph>
    $<TARGET_OBJECTS:Clustering>
    $<TARGET_OBJECTS:Bisection>
    $<TARGET_OBJECTS:Driver>
    $<TARGET_OBJECTS:Refinement>
    $<TARGET_OBJECTS:Scheduling>
)

# Set include directories for public interface
target_include_directories(dag_partitioning
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Graph/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Driver/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Scheduling/include>
)

# Link only external dependencies and robin_hood for public headers
target_link_libraries(dag_partitioning
        PUBLIC
        robin_hood::robin_hood
        ortools::ortools
        PRIVATE
        METIS::METIS
        GKlib::GKlib
        SCOTCH::SCOTCH
)

# Link OpenMP if enabled
if(DAG_PARTITIONING_OPENMP)
    target_link_libraries(dag_partitioning PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(dag_partitioning PRIVATE DAG_PARTITIONING_OPENMP_ENABLED)
else()
    target_link_libraries(dag_partitioning PRIVATE Threads::Threads)
endif()
