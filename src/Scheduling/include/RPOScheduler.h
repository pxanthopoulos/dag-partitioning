/**
 * @file RPOScheduler.h
 * @brief Reverse post-order heuristic scheduler
 */

#ifndef DAG_PARTITIONING_SCHEDULING_RPOSCHEDULER_H
#define DAG_PARTITIONING_SCHEDULING_RPOSCHEDULER_H

#include "HyperGraph.h"

#include <cstdint>
#include <utility>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace heuristic {

/**
 * @brief RPO (Reverse Post-Order) heuristic scheduler
 *
 * Computes a schedule using reverse post-order traversal of the DAG.
 * RPO tends to complete subtrees before moving to siblings, which
 * reduces the number of simultaneously live tensors.
 */
class RPOScheduler {
  protected:
    const hypergraph::HyperGraph &graph;
    bool debug = false;

  public:
    /**
     * @brief Constructor
     * @param graph Hyper graph representation
     * @param debug Enable debug output
     */
    RPOScheduler(const hypergraph::HyperGraph &graph, bool debug = false);

    virtual ~RPOScheduler() = default;

    /**
     * @brief Computes the RPO schedule
     * @return Pair of (schedule, peak memory)
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> solve() const;
};

} // namespace heuristic

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_RPOSCHEDULER_H
