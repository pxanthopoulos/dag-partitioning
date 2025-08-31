# Installation configuration for dag_partitioning

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Export targets for installation
install(TARGETS dag_partitioning
    EXPORT dag_partitioningTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(TARGETS rand-dag dag-test
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install only public API headers
install(DIRECTORY 
    src/Graph/include/ 
    src/Driver/include/ 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install export targets
install(EXPORT dag_partitioningTargets
    FILE dag_partitioningTargets.cmake
    NAMESPACE dag_partitioning::
    DESTINATION lib/cmake/dag_partitioning
)

# Generate and install package config files
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/dag_partitioningConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/dag_partitioningConfig.cmake
    INSTALL_DESTINATION lib/cmake/dag_partitioning
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/dag_partitioningConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/dag_partitioningConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/dag_partitioningConfigVersion.cmake
    DESTINATION lib/cmake/dag_partitioning
)