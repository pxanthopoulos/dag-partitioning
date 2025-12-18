/**
 * @file HyperGraph.cpp
 * @brief Implementation of hyper graph representation for memory-aware
 * scheduling
 */

#include "HyperGraph.h"

#include <queue>

namespace dag_partitioning {

namespace scheduling {

namespace hypergraph {

HyperGraph::HyperGraph(const core::Graph &graph) {
    size = graph.size;
    tensorSizes.resize(size);
    extraSizes.resize(size);
    inputTensors.resize(size);

    for (uint64_t node = 0; node < size; ++node) {
        tensorSizes[node] = 0;
        extraSizes[node] = graph.nodeWeights[node];

        for (const auto &[neighbor, edgeWeight] : graph.adj[node]) {
            inputTensors[neighbor].insert(node);
            tensorSizes[node] = std::max(tensorSizes[node], edgeWeight);
        }
    }

    computeTopology();
}

void HyperGraph::computeTopology() {
    computeAncestors();
    computeDescendants();
}

void HyperGraph::computeAncestors() {
    ancestors.clear();
    ancestors.resize(size);
    tensorUsers.clear();
    tensorUsers.resize(size);

    // Compute in-degree for topological sort
    std::vector<uint64_t> inDegree(size, 0);
    for (uint64_t i = 0; i < size; i++) {
        inDegree[i] = inputTensors[i].size();
    }

    // Topological sort using Kahn's algorithm
    std::queue<uint64_t> queue;
    for (uint64_t i = 0; i < size; i++) {
        if (inDegree[i] == 0) {
            queue.push(i);
        }
    }

    while (!queue.empty()) {
        uint64_t node = queue.front();
        queue.pop();

        // For this node, its ancestors are the union of:
        // 1. Its direct inputs
        // 2. All ancestors of its inputs
        for (uint64_t input : inputTensors[node]) {
            ancestors[node].insert(input);
            ancestors[node].insert(ancestors[input].begin(),
                                   ancestors[input].end());
            tensorUsers[input].insert(node);
        }

        // Update in-degrees and enqueue nodes
        for (uint64_t i = 0; i < size; i++) {
            if (inputTensors[i].count(node)) {
                inDegree[i]--;
                if (inDegree[i] == 0) {
                    queue.push(i);
                }
            }
        }
    }
}

void HyperGraph::computeDescendants() {
    descendants.clear();
    descendants.resize(size);

    // Build descendants from ancestors: if j is ancestor of i, then i is
    // descendant of j
    for (uint64_t i = 0; i < size; i++) {
        for (uint64_t ancestor : ancestors[i]) {
            descendants[ancestor].insert(i);
        }
    }
}

void HyperGraph::print(std::ostream &os) const {
    os << "Nodes: " << size << std::endl;

    for (uint64_t i = 0; i < size; i++) {
        os << "\nNode: " << i << std::endl;
        os << "    tensorSizes=" << tensorSizes[i]
           << ", extraSizes=" << extraSizes[i] << std::endl;

        os << "    inputs={";
        for (auto it = inputTensors[i].begin(); it != inputTensors[i].end();
             ++it) {
            if (it != inputTensors[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}, users={";
        for (auto it = tensorUsers[i].begin(); it != tensorUsers[i].end();
             ++it) {
            if (it != tensorUsers[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}" << std::endl;

        os << "    anc={";
        for (auto it = ancestors[i].begin(); it != ancestors[i].end(); ++it) {
            if (it != ancestors[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}, des={";
        for (auto it = descendants[i].begin(); it != descendants[i].end();
             ++it) {
            if (it != descendants[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}" << std::endl;
    }
}

uint64_t HyperGraph::getSize() const { return size; }

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getAncestors() const {
    return ancestors;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getDescendants() const {
    return descendants;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getTensorUsers() const {
    return tensorUsers;
}

const std::vector<uint64_t> &HyperGraph::getTensorSizes() const {
    return tensorSizes;
}

const std::vector<uint64_t> &HyperGraph::getExtraSizes() const {
    return extraSizes;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getInputTensors() const {
    return inputTensors;
}

std::ostream &operator<<(std::ostream &os, const HyperGraph &graph) {
    graph.print(os);
    return os;
}

} // namespace hypergraph

} // namespace scheduling

} // namespace dag_partitioning