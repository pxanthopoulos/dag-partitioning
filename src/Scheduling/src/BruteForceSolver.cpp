/**
 * @file BruteForceSolver.cpp
 * @brief Exhaustive search over all orders, to verify the other solvers
 */

#include "BruteForceSolver.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace dag_partitioning {

namespace scheduling {

namespace bruteforce {

BruteForceSolver::BruteForceSolver(const hypergraph::HyperGraph &graph,
                                   bool debug, Objective objective,
                                   const PartitionFootprints *footprints,
                                   uint64_t memoryLimit)
    : graph(graph), debug(debug), objective(objective), footprints(footprints),
      memoryLimit(memoryLimit) {}

bool BruteForceSolver::canPlace(uint64_t node,
                                const std::vector<bool> &placed) const {
    const auto &inputTensors = graph.getInputTensors();
    for (uint64_t t : inputTensors[node]) {
        if (!placed[graph.getTensorProducers()[t]])
            return false;
    }
    return true;
}

void BruteForceSolver::enumerate(std::vector<uint64_t> &current,
                                 std::vector<bool> &placed, uint64_t &bestCost,
                                 std::vector<uint64_t> &bestOrder,
                                 uint64_t &count,
                                 std::vector<std::string> &outputs,
                                 const Cost &cost) const {
    uint64_t n = graph.getSize();
    if (current.size() == n) {
        count++;
        uint64_t orderCost = cost(current);

        if (debug) {
            std::ostringstream oss;
            oss << "Order " << count - 1 << " (cost " << orderCost << "): ";
            for (const uint64_t &node : current) {
                oss << node << " ";
            }
            outputs.emplace_back(oss.str());
        }

        if (orderCost < bestCost) {
            bestCost = orderCost;
            bestOrder = current;
        }
        return;
    }

    for (uint64_t node = 0; node < n; node++) {
        if (!placed[node] && canPlace(node, placed)) {
            placed[node] = true;
            current.push_back(node);

            enumerate(current, placed, bestCost, bestOrder, count, outputs,
                      cost);

            current.pop_back();
            placed[node] = false;
        }
    }
}

std::pair<std::vector<uint64_t>, uint64_t> BruteForceSolver::solve() const {
    if (objective == Objective::Memory) {
        return best([&](const std::vector<uint64_t> &order) {
            return graph.computePeakMemory(order);
        });
    }
    if (objective != Objective::TransferMemory)
        throw std::invalid_argument("Unknown scheduling objective");

    // The lowest peak any order reaches can be above the limit
    uint64_t limit = memoryLimit;
    if (limit != UINT64_MAX) {
        uint64_t lowest = best([&](const std::vector<uint64_t> &order) {
                              return footprints->peak(order);
                          }).second;
        limit = std::max(limit, lowest);
    }
    return best([&](const std::vector<uint64_t> &order) {
        return limit == UINT64_MAX || footprints->peak(order) <= limit
                   ? graph.computeTransfers(order)
                   : UINT64_MAX;
    });
}

std::pair<std::vector<uint64_t>, uint64_t>
BruteForceSolver::best(const Cost &cost) const {
    uint64_t n = graph.getSize();
    std::vector<std::string> outputs;

    std::vector<uint64_t> current, bestOrder;
    std::vector<bool> placed(n, false);
    uint64_t bestCost = UINT64_MAX;
    uint64_t count = 0;

    enumerate(current, placed, bestCost, bestOrder, count, outputs, cost);

    if (debug) {
        std::cerr << "\nBrute-force results:\n";
        for (const auto &out : outputs) {
            std::cerr << out << std::endl;
        }
        std::cerr << "Best result: cost " << bestCost << " in order: ";
        for (const uint64_t &node : bestOrder) {
            std::cerr << node << " ";
        }
        std::cerr << std::endl;
    }

    return {bestOrder, bestCost};
}

} // namespace bruteforce

} // namespace scheduling

} // namespace dag_partitioning
