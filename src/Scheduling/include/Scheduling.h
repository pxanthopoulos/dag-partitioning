/**
 * @file Scheduling.h
 * @brief Memory-aware scheduling for DAGs
 *
 * Performs scheduling of DAG operations to minimize peak memory
 * usage. Uses an ILP-based solver (brute-force solver for small
 * graphs for verification).
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

    /**
     * @brief Computes topological relationships
     */
    void computeTopology();

    /**
     * @brief Computes transitive closure of predecessors
     */
    void computeAncestors();

    /**
     * @brief Computes transitive closure of successors
     */
    void computeDescendants();

    /**
     * @brief Outputs graph structure for debugging
     * @param os Output stream
     */
    void print(std::ostream &os) const;

  public:
    /**
     * @brief Constructs ILP graph from a DAG
     * @param graph Input graph to analyze
     */
    ILPGraph(const core::Graph &graph);

    virtual ~ILPGraph() = default;

    /**
     * @brief Returns number of nodes
     * @return Graph size
     */
    uint64_t getSize() const;

    /**
     * @brief Returns ancestor sets for all nodes
     * @return Vector of ancestor sets
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getAncestors() const;

    /**
     * @brief Returns descendant sets for all nodes
     * @return Vector of descendant sets
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getDescendants() const;

    /**
     * @brief Returns tensor consumer information
     * @return Vector of sets indicating which operations consume each tensor
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getTensorUsers() const;

    /**
     * @brief Returns output tensor sizes
     * @return Vector of tensor sizes
     */
    const std::vector<uint64_t> &getTensorSizes() const;

    /**
     * @brief Returns extra memory required during operation execution
     * @return Vector of workspace memory sizes
     */
    const std::vector<uint64_t> &getExtraSizes() const;

    /**
     * @brief Returns input tensor dependencies
     * @return Vector of sets indicating which tensors each operation needs
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getInputTensors() const;

    /**
     * @brief Stream output operator
     * @param os Output stream
     * @param graph Graph to output
     * @return Output stream
     */
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

    /**
     * @brief Computes earliest/latest valid positions for pruning
     *
     * Uses topological information to determine valid time windows,
     * reducing the number of ILP variables.
     *
     * @param earliestOp Earliest position each operation can run
     * @param latestOp Latest position each operation can run
     * @param earliestTensor Earliest time tensor can exist
     * @param latestTensor Latest time tensor can be needed
     */
    void computePruningBound(std::vector<uint64_t> &earliestOp,
                             std::vector<uint64_t> &latestOp,
                             std::vector<uint64_t> &earliestTensor,
                             std::vector<uint64_t> &latestTensor) const;

    /**
     * @brief Creates ILP variables and constraints
     *
     * Sets up O[i][j] (operation i at step j), T[i][j] (tensor i alive
     * at step j), and memory bound variable. Adds all scheduling and
     * memory constraints.
     *
     * @param earliestOp Earliest position each operation can run
     * @param latestOp Latest position each operation can run
     * @param earliestTensor Earliest time tensor can exist
     * @param latestTensor Latest time tensor is needed
     */
    void createVariables(const std::vector<uint64_t> &earliestOp,
                         const std::vector<uint64_t> &latestOp,
                         const std::vector<uint64_t> &earliestTensor,
                         const std::vector<uint64_t> &latestTensor);

  public:
    /**
     * @brief Constructor
     * @param graph ILP graph representation
     * @param debug Enable debug output
     */
    ILPSolver(const ILPGraph &graph, bool debug);

    virtual ~ILPSolver() = default;

    /**
     * @brief Solves the scheduling problem using ILP
     *
     * @param timeLimitSeconds Time limit for solver in seconds
     * @return Tuple of (solver status, schedule, peak memory)
     */
    [[nodiscard]] std::tuple<operations_research::MPSolver::ResultStatus,
                             std::vector<uint64_t>, uint64_t>
    solve(uint64_t timeLimitSeconds = 60);

    /**
     * @brief Outputs variable values for debugging
     * @param os Output stream
     */
    void printVariables(std::ostream &os) const;
};

} // namespace ilp

namespace bruteforce {

class BruteForceSolver {
  protected:
    const ilp::ILPGraph &graph;
    bool debug = false;

    /**
     * @brief Computes peak memory for a given execution order
     *
     * @param order Execution order of operations
     * @return Peak memory usage
     */
    [[nodiscard]] uint64_t
    computePeakMemory(const std::vector<uint64_t> &order) const;

    /**
     * @brief Checks if a node can be placed given current placement state
     *
     * @param node Node to check
     * @param placed Vector indicating which nodes have been placed
     * @return true if all dependencies are satisfied
     */
    [[nodiscard]] bool canPlace(uint64_t node,
                                const std::vector<bool> &placed) const;

    /**
     * @brief Recursively enumerates all valid topological orders
     *
     * @param current Current partial order being built
     * @param placed Vector tracking which nodes have been placed
     * @param bestPeak Best peak memory found so far
     * @param bestOrder Best order found so far
     * @param count Counter for total orders enumerated
     * @param outputs Storage for debug output strings
     */
    void enumerate(std::vector<uint64_t> &current, std::vector<bool> &placed,
                   uint64_t &bestPeak, std::vector<uint64_t> &bestOrder,
                   uint64_t &count, std::vector<std::string> &outputs) const;

  public:
    /**
     * @brief Constructor
     *
     * @param graph ILP graph representing the scheduling problem
     * @param debug Enable debug output
     */
    BruteForceSolver(const ilp::ILPGraph &graph, bool debug = false);

    virtual ~BruteForceSolver() = default;

    /**
     * @brief Solves the scheduling problem via brute-force enumeration
     *
     * @return Pair of (schedule, peak memory)
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> solve() const;
};

} // namespace bruteforce

class Scheduler {
  protected:
    const core::Graph &originalGraph;              // Original graph to schedule
    const std::vector<uint64_t> &partitionMapping; // Node to partition mapping
    std::unique_ptr<core::Graph> coarseGraph;      // Resulting coarse graph
    bool debug = false;
    bool verify = false;

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
     * @brief Constructor
     * @param originalGraph Graph to schedule
     * @param partitionMapping Vector mapping each node to its partition
     * @param debug Enable debug output
     * @param verify If true, run brute-force solver and verify peak memory
     * matches
     */
    Scheduler(const core::Graph &originalGraph,
              std::vector<uint64_t> &partitionMapping, bool debug = false,
              bool verify = false);

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
     * @return Pair of a vector with partition numbers ordered by their
     * scheduled execution and the peak memory usage of the schedule
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> run();
};

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_H
