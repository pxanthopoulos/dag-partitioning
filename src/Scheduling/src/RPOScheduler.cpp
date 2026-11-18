/**
 * @file RPOScheduler.cpp
 * @brief Reverse post-order heuristic scheduler
 */

#include "RPOScheduler.h"

#include <iostream>
#include <stack>
#include <tuple>

namespace dag_partitioning {

namespace scheduling {

namespace heuristic {

RPOScheduler::RPOScheduler(const hypergraph::HyperGraph &graph, bool debug)
    : graph(graph), debug(debug) {}

std::pair<std::vector<uint64_t>, uint64_t> RPOScheduler::solve() const {
    uint64_t n = graph.getSize();
    std::vector<uint64_t> schedule;
    schedule.reserve(n);

    // Build adjacency list for DFS (node -> children/successors)
    std::vector<std::vector<uint64_t>> children(n);
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t input : graph.getInputTensors()[i]) {
            children[graph.getTensorProducers()[input]].push_back(i);
        }
    }

    // Find root nodes (nodes with no inputs)
    std::vector<uint64_t> roots;
    for (uint64_t i = 0; i < n; i++) {
        if (graph.getInputTensors()[i].empty()) {
            roots.push_back(i);
        }
    }

    // DFS to compute reverse post-order
    std::vector<bool> visited(n, false);
    std::vector<uint64_t> postOrder;
    postOrder.reserve(n);

    // Iterative DFS using explicit stack to avoid stack overflow
    std::stack<std::tuple<uint64_t, size_t, bool>> dfsStack;

    for (uint64_t root : roots) {
        if (visited[root])
            continue;

        dfsStack.push({root, 0, false});

        while (!dfsStack.empty()) {
            auto &[node, childIdx, postVisit] = dfsStack.top();

            if (postVisit) {
                postOrder.push_back(node);
                dfsStack.pop();
                continue;
            }

            if (!visited[node]) {
                visited[node] = true;
            }

            // Find next unvisited child
            while (childIdx < children[node].size() &&
                   visited[children[node][childIdx]]) {
                childIdx++;
            }

            if (childIdx < children[node].size()) {
                uint64_t child = children[node][childIdx];
                childIdx++;
                dfsStack.push({child, 0, false});
            } else {
                postVisit = true;
            }
        }
    }

    // Reverse post-order
    for (auto it = postOrder.rbegin(); it != postOrder.rend(); ++it) {
        schedule.push_back(*it);
    }

    if (debug) {
        std::cerr << "\nRPO schedule: ";
        for (uint64_t op : schedule) {
            std::cerr << op << " ";
        }
        std::cerr << std::endl;
    }

    // Compute peak memory
    uint64_t peakMemory = graph.computePeakMemory(schedule);

    if (debug) {
        std::cerr << "RPO peak memory: " << peakMemory << std::endl;
    }

    return {schedule, peakMemory};
}

} // namespace heuristic

} // namespace scheduling

} // namespace dag_partitioning
