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
#include <utility>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace bruteforce {

/**
 * @brief Tries every order, to verify the other solvers on small graphs
 */
class BruteForceSolver {
  public:
    /**
     * @param objective Criterion to minimize
     * @param footprints Footprints of the nodes (Objective::TransferMemory)
     * @param memoryLimit Limit on the peak footprint
     * (Objective::TransferMemory); if no order meets it, the lowest peak any
     * order reaches
     */
    explicit BruteForceSolver(const hypergraph::HyperGraph &graph,
                              Objective objective = Objective::Memory,
                              const PartitionFootprints *footprints = nullptr,
                              uint64_t memoryLimit = UINT64_MAX);

    /**
     * @return The best order and its cost
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> solve() const;

  private:
    const hypergraph::HyperGraph &graph;
    Objective objective;
    const PartitionFootprints *footprints;
    uint64_t memoryLimit;
};

} // namespace bruteforce

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_BRUTEFORCESOLVER_H
