/**
 * @file Utils.h
 * @brief Utility functions for memory-aware scheduling
 *
 * Provides helper functions and data structures used in memory-aware scheduling
 * of DAG operations to minimize peak memory usage.
 */
#ifndef DAG_PARTITIONING_UTILS_H
#define DAG_PARTITIONING_UTILS_H

#include <cstdint>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace utils {

/**
 * @brief Tensor representation for memory allocation
 */
struct Tensor {
    uint64_t size;   // Memory size required
    uint64_t birth;  // Local position when tensor is created
    uint64_t death;  // Local position when tensor is last used
    uint64_t offset; // Memory offset assigned by packing
};

/**
 * @brief Performs best-fit bin packing on tensors
 *
 * Assigns memory offsets to tensors while minimizing peak memory usage.
 * Uses a best-fit algorithm that finds the smallest gap for each tensor.
 *
 * @param tensors Vector of tensors to pack (modified in-place)
 * @return Peak memory usage after packing
 */
[[nodiscard]] uint64_t packBestFit(std::vector<Tensor> &tensors);

} // namespace utils

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_UTILS_H
