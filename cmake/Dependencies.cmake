# External dependencies for dag_partitioning

include(FetchContent)

# Include individual dependency configurations
include(${CMAKE_CURRENT_LIST_DIR}/GKlib.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/METIS.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/SCOTCH.cmake)