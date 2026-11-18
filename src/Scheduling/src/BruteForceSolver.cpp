/**
 * @file BruteForceSolver.cpp
 * @brief Exhaustive search over all orders, to verify the other solvers
 */

#include "BruteForceSolver.h"

#include <algorithm>
#include <functional>

namespace dag_partitioning {

namespace scheduling {

namespace bruteforce {

namespace {

/**
 * @brief The topological order of lowest cost, the first found on ties, with
 * the nodes tried in id order; a cost of UINT64_MAX means not allowed
 */
std::pair<std::vector<uint64_t>, uint64_t>
best(const hypergraph::HyperGraph &graph,
     const std::function<uint64_t(const std::vector<uint64_t> &)> &cost) {
    uint64_t n = graph.getSize();
    std::vector<uint64_t> order, bestOrder;
    std::vector<bool> placed(n, false);
    uint64_t bestCost = UINT64_MAX;
    auto ready = [&](uint64_t node) {
        for (uint64_t t : graph.getInputTensors()[node]) {
            if (!placed[graph.getTensorProducers()[t]])
                return false;
        }
        return true;
    };
    std::function<void()> extend = [&] {
        if (order.size() == n) {
            uint64_t orderCost = cost(order);
            if (orderCost < bestCost) {
                bestCost = orderCost;
                bestOrder = order;
            }
            return;
        }
        for (uint64_t node = 0; node < n; node++) {
            if (placed[node] || !ready(node))
                continue;
            placed[node] = true;
            order.push_back(node);
            extend();
            order.pop_back();
            placed[node] = false;
        }
    };
    extend();
    return {bestOrder, bestCost};
}

} // namespace

BruteForceSolver::BruteForceSolver(const hypergraph::HyperGraph &graph,
                                   Objective objective,
                                   const PartitionFootprints *footprints,
                                   uint64_t memoryLimit)
    : graph(graph), objective(objective), footprints(footprints),
      memoryLimit(memoryLimit) {}

std::pair<std::vector<uint64_t>, uint64_t> BruteForceSolver::solve() const {
    if (objective == Objective::Memory) {
        return best(graph, [&](const std::vector<uint64_t> &order) {
            return graph.computePeakMemory(order);
        });
    }

    // The lowest peak any order reaches can be above the limit
    uint64_t limit = memoryLimit;
    if (limit != UINT64_MAX) {
        uint64_t lowest = best(graph, [&](const std::vector<uint64_t> &order) {
                              return footprints->peak(order);
                          }).second;
        limit = std::max(limit, lowest);
    }
    return best(graph, [&](const std::vector<uint64_t> &order) {
        return limit == UINT64_MAX || footprints->peak(order) <= limit
                   ? graph.computeTransfers(order)
                   : UINT64_MAX;
    });
}

} // namespace bruteforce

} // namespace scheduling

} // namespace dag_partitioning
