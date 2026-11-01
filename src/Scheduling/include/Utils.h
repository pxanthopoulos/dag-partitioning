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
#include <functional>
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
 * Sizes are rounded up to rangeAlignment and offsets to offsetAlignment, as
 * IREE does when it packs transient slices.
 *
 * @param tensors Vector of tensors to pack (modified in-place, sizes become
 * the aligned sizes)
 * @param offsetAlignment Alignment of every assigned offset
 * @param rangeAlignment Alignment of every size and of the peak
 * @return Peak memory usage after packing
 */
[[nodiscard]] uint64_t packBestFit(std::vector<Tensor> &tensors,
                                   uint64_t offsetAlignment = 1,
                                   uint64_t rangeAlignment = 1);

/**
 * @brief Packs the tensors of one partition into a slab
 *
 * Gets the partition's tensors sorted by birth, then death. Assigns each an
 * offset (sizes may be rounded up in place) and returns the slab size, which
 * becomes the partition's weight. Lets a compiler give the scheduler the
 * packer it will itself use for the partition.
 */
using SlabPacker = std::function<uint64_t(std::vector<Tensor> &tensors)>;

} // namespace utils

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_UTILS_H
