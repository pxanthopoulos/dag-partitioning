/**
 * @file Scheduling.h
 * @brief TODO
 */

#ifndef DAG_PARTITIONING_SCHEDULING_H
#define DAG_PARTITIONING_SCHEDULING_H

#include <Graph.h>

#include <cstdint>
#include <memory>
#include <vector>

namespace dag_partitioning {

namespace core {
class Graph;
}

namespace scheduling {

/**
 * @brief Bin packing utilities for memory allocation
 */
namespace packing {

/**
 * @brief Tensor representation for memory allocation
 */
struct Tensor {
    uint64_t producer; // Node that produces this tensor
    uint64_t size;     // Memory size required
    uint64_t birth;    // Local position when tensor is created
    uint64_t death;    // Local position when tensor is last used
    uint64_t offset;   // Memory offset assigned by packing
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

} // namespace packing

class Scheduler {
  protected:
    const core::Graph &originalGraph;              // Original graph to schedule
    const std::vector<uint64_t> &partitionMapping; // Node to partition mapping
    std::unique_ptr<core::Graph> coarseGraph;      // Resulting coarse graph

    /**
     * @brief Calculates memory requirements for each partition
     *
     * Uses best-fit memory packing to compute peak memory usage for
     * each partition based on tensor lifetimes.
     *
     * @return Vector of partition weights (peak memory usage)
     */
    [[nodiscard]] std::vector<uint64_t> calculatePartitionWeights() const;

    /**
     * @brief Builds the coarse graph from partitioned original graph
     *
     * Creates a new graph where each node represents a partition.
     * Node weights are computed using memory packing, and edge weights
     * represent inter-partition communication.
     */
    void buildCoarseGraph();

  public:
    /**
     * @brief Constructor for Scheduler
     * @param originalGraph Graph to schedule
     * @param partitionMapping Vector mapping each node to its partition
     */
    Scheduler(const core::Graph &originalGraph,
              std::vector<uint64_t> &partitionMapping);

    /**
     * @brief Virtual destructor for proper cleanup
     */
    virtual ~Scheduler() = default;

    /**
     * @brief Executes the scheduling algorithm
     *
     * Builds the coarse graph with partition weights computed via memory
     * packing and inter-partition edge weights.
     */
    void run();

    /**
     * @brief Gets the coarse graph
     * @return Const reference to the unique_ptr of the coarse graph
     */
    [[nodiscard]] const std::unique_ptr<core::Graph> &getCoarseGraph() const;
};

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_H
