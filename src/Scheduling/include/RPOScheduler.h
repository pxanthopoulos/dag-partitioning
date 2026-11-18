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
 * @brief Orders the nodes in reverse post-order of a depth-first search,
 * which tends to finish a subgraph before it starts the next, so that few
 * tensors are live at once
 */
class RPOScheduler {
  public:
    explicit RPOScheduler(const hypergraph::HyperGraph &graph,
                          bool debug = false);

    /**
     * @return The order and its peak memory
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> solve() const;

  private:
    const hypergraph::HyperGraph &graph;
    bool debug;
};

} // namespace heuristic

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_RPOSCHEDULER_H
