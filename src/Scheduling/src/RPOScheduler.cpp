/**
 * @file RPOScheduler.cpp
 * @brief Reverse post-order heuristic scheduler
 */

#include "RPOScheduler.h"

#include <iostream>

namespace dag_partitioning {

namespace scheduling {

namespace heuristic {

RPOScheduler::RPOScheduler(const hypergraph::HyperGraph &graph, bool debug)
    : graph(graph), debug(debug) {}

std::pair<std::vector<uint64_t>, uint64_t> RPOScheduler::solve() const {
    uint64_t n = graph.getSize();
    std::vector<std::vector<uint64_t>> children(n);
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t input : graph.getInputTensors()[i])
            children[graph.getTensorProducers()[input]].push_back(i);
    }

    // Iterative depth-first search from every node without inputs; a stack
    // entry is a node and the index of its next child
    std::vector<bool> visited(n, false);
    std::vector<uint64_t> postOrder;
    std::vector<std::pair<uint64_t, size_t>> stack;
    for (uint64_t root = 0; root < n; root++) {
        if (!graph.getInputTensors()[root].empty() || visited[root])
            continue;
        visited[root] = true;
        stack.push_back({root, 0});
        while (!stack.empty()) {
            auto &[node, next] = stack.back();
            const auto &nodeChildren = children[node];
            while (next < nodeChildren.size() && visited[nodeChildren[next]])
                next++;
            if (next == nodeChildren.size()) {
                postOrder.push_back(node);
                stack.pop_back();
            } else {
                uint64_t child = nodeChildren[next++];
                visited[child] = true;
                stack.push_back({child, 0});
            }
        }
    }

    std::vector<uint64_t> schedule(postOrder.rbegin(), postOrder.rend());
    uint64_t peakMemory = graph.computePeakMemory(schedule);
    if (debug)
        std::cerr << "RPO peak memory: " << peakMemory << std::endl;
    return {schedule, peakMemory};
}

} // namespace heuristic

} // namespace scheduling

} // namespace dag_partitioning
