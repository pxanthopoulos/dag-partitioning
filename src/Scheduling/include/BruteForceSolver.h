/**
 * @file BruteForceSolver.h
 * @brief Exhaustive search over all orders, to verify the other solvers
 */

#ifndef DAG_PARTITIONING_SCHEDULING_BRUTEFORCESOLVER_H
#define DAG_PARTITIONING_SCHEDULING_BRUTEFORCESOLVER_H

#include "HyperGraph.h"
#include "Objective.h"
#include "PartitionFootprints.h"

#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace bruteforce {

class BruteForceSolver {
  protected:
    const hypergraph::HyperGraph &graph;
    bool debug = false;
    Objective objective = Objective::Memory;
    const PartitionFootprints *footprints = nullptr;
    uint64_t memoryLimit = UINT64_MAX;

    /**
     * @brief Cost of an execution order; UINT64_MAX means not allowed
     */
    using Cost = std::function<uint64_t(const std::vector<uint64_t> &)>;

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
     * @param bestCost Best cost found so far
     * @param bestOrder Best order found so far
     * @param count Counter for total orders enumerated
     * @param outputs Storage for debug output strings
     * @param cost Cost of a complete order
     */
    void enumerate(std::vector<uint64_t> &current, std::vector<bool> &placed,
                   uint64_t &bestCost, std::vector<uint64_t> &bestOrder,
                   uint64_t &count, std::vector<std::string> &outputs,
                   const Cost &cost) const;

    /**
     * @brief Finds the best order for a cost
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t>
    best(const Cost &cost) const;

  public:
    /**
     * @brief Constructor
     *
     * @param graph Hyper graph representing the scheduling problem
     * @param debug Enable debug output
     * @param objective Criterion to minimize
     * @param footprints Footprints of the nodes (Objective::TransferMemory)
     * @param memoryLimit Limit on the peak (Objective::TransferMemory)
     */
    BruteForceSolver(const hypergraph::HyperGraph &graph, bool debug = false,
                     Objective objective = Objective::Memory,
                     const PartitionFootprints *footprints = nullptr,
                     uint64_t memoryLimit = UINT64_MAX);

    virtual ~BruteForceSolver() = default;

    /**
     * @brief Solves the scheduling problem via brute-force enumeration
     *
     * @return Pair of (schedule, cost)
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> solve() const;
};

} // namespace bruteforce

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_BRUTEFORCESOLVER_H
