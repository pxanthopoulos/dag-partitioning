add_library(dag_partitioning STATIC)

add_subdirectory(Graph)
add_subdirectory(Clustering)
add_subdirectory(Bisection)
add_subdirectory(Driver)
add_subdirectory(Refinement)

# Include all object library sources in the main library
target_sources(dag_partitioning PRIVATE
    $<TARGET_OBJECTS:Graph>
    $<TARGET_OBJECTS:Clustering>
    $<TARGET_OBJECTS:Bisection>
    $<TARGET_OBJECTS:Driver>
    $<TARGET_OBJECTS:Refinement>
)

# Set include directories for public interface
target_include_directories(dag_partitioning
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Graph/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Driver/include>
)

# Link only external dependencies
target_link_libraries(dag_partitioning
        PRIVATE
        METIS::METIS
        GKlib::GKlib
        SCOTCH::SCOTCH
        Threads::Threads
)
