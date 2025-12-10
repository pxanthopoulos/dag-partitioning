/**
 * @file Scheduling.h
 * @brief TODO
 */

#ifndef DAG_PARTITIONING_SCHEDULING_H
#define DAG_PARTITIONING_SCHEDULING_H

#include <Graph.h>

#include "ortools/linear_solver/linear_solver.h"
#include "robin_hood.h"
#include <cstdint>
#include <iostream>
#include <memory>
#include <queue>
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

namespace ilp {

class ILPGraph {
  private:
    uint64_t size;
    std::vector<uint64_t> tensorSizes; // S_i: output tensor size
    std::vector<uint64_t> extraSizes;  // ES_i: extra memory during execution
    std::vector<robin_hood::unordered_set<uint64_t>>
        inputTensors; // IN_i: input tensors

    std::vector<robin_hood::unordered_set<uint64_t>> ancestors;
    std::vector<robin_hood::unordered_set<uint64_t>> descendants;
    std::vector<robin_hood::unordered_set<uint64_t>> tensorUsers;

    void computeTopology();

    void computeAncestors();

    void computeDescendants();

    void print(std::ostream &os) const;

  public:
    ILPGraph(const core::Graph &graph);

    virtual ~ILPGraph() = default;

    uint64_t getSize() const;

    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getAncestors() const;

    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getDescendants() const;

    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getTensorUsers() const;

    const std::vector<uint64_t> &getTensorSizes() const;

    const std::vector<uint64_t> &getExtraSizes() const;

    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getInputTensors() const;

    friend std::ostream &operator<<(std::ostream &os, const ILPGraph &graph);
};
class ILPSolver {
  protected:
    const ILPGraph &graph;
    std::unique_ptr<operations_research::MPSolver> solver;
    bool debug = false;

    // ILP variables
    std::vector<std::vector<operations_research::MPVariable *>> O;
    std::vector<std::vector<operations_research::MPVariable *>> T;
    operations_research::MPVariable *mem = nullptr;

    void computePruningBound(std::vector<uint64_t> &earliestOp,
                             std::vector<uint64_t> &latestOp,
                             std::vector<uint64_t> &earliestTensor,
                             std::vector<uint64_t> &latestTensor) const;

    void createVariables(const std::vector<uint64_t> &earliestOp,
                         const std::vector<uint64_t> &latestOp,
                         const std::vector<uint64_t> &earliestTensor,
                         const std::vector<uint64_t> &latestTensor);

  public:
    ILPSolver(const ILPGraph &graph, bool debug);

    std::vector<uint64_t> solve(uint64_t timeLimitSeconds = 60);
};

} // namespace ilp

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
     *
     * @return Vector of partition numbers ordered by their scheduled execution
     */
    [[nodiscard]] std::vector<uint64_t> run();

    /**
     * @brief Gets the coarse graph
     * @return Const reference to the unique_ptr of the coarse graph
     */
    [[nodiscard]] const std::unique_ptr<core::Graph> &getCoarseGraph() const;
};

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_H
