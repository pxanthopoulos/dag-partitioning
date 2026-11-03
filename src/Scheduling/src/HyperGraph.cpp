/**
 * @file HyperGraph.cpp
 * @brief Implementation of hyper graph representation for memory-aware
 * scheduling
 */

#include "HyperGraph.h"

#include <algorithm>
#include <queue>

namespace dag_partitioning {

namespace scheduling {

namespace hypergraph {

HyperGraph::HyperGraph(std::vector<uint64_t> extraSizes,
                       const std::vector<Tensor> &tensors)
    : size(extraSizes.size()), extraSizes(std::move(extraSizes)) {
    outputTensors.resize(size);
    inputTensors.resize(size);
    tensorUsers.resize(tensors.size());

    for (uint64_t t = 0; t < tensors.size(); ++t) {
        tensorSizes.push_back(tensors[t].size);
        tensorProducers.push_back(tensors[t].producer);
        liveOut.push_back(tensors[t].liveOut);
        outputTensors[tensors[t].producer].push_back(t);
        for (uint64_t user : tensors[t].users) {
            tensorUsers[t].insert(user);
            inputTensors[user].insert(t);
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

    // Predecessors (producers of input tensors) and successors of every node
    std::vector<robin_hood::unordered_set<uint64_t>> predecessors(size);
    std::vector<std::vector<uint64_t>> successors(size);
    for (uint64_t i = 0; i < size; i++) {
        for (uint64_t t : inputTensors[i]) {
            uint64_t producer = tensorProducers[t];
            if (predecessors[i].insert(producer).second)
                successors[producer].push_back(i);
        }
    }

    // Compute in-degree for topological sort
    std::vector<uint64_t> inDegree(size, 0);
    for (uint64_t i = 0; i < size; i++) {
        inDegree[i] = predecessors[i].size();
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
        // 1. Its direct predecessors
        // 2. All ancestors of its predecessors
        for (uint64_t pred : predecessors[node]) {
            ancestors[node].insert(pred);
            ancestors[node].insert(ancestors[pred].begin(),
                                   ancestors[pred].end());
        }

        // Update in-degrees and enqueue nodes
        for (uint64_t succ : successors[node]) {
            inDegree[succ]--;
            if (inDegree[succ] == 0) {
                queue.push(succ);
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

uint64_t
HyperGraph::computePeakMemory(const std::vector<uint64_t> &order) const {
    uint64_t n = order.size();
    std::vector<uint64_t> pos(size);
    for (uint64_t step = 0; step < n; step++)
        pos[order[step]] = step;

    // Tensor t adds S_t to the steps [birth(t), death(t)]
    std::vector<int64_t> delta(n + 1, 0);
    for (uint64_t t = 0; t < tensorSizes.size(); t++) {
        uint64_t birth = pos[tensorProducers[t]];
        uint64_t death = liveOut[t] ? n - 1 : birth;
        for (uint64_t user : tensorUsers[t])
            death = std::max(death, pos[user]);
        delta[birth] += static_cast<int64_t>(tensorSizes[t]);
        delta[death + 1] -= static_cast<int64_t>(tensorSizes[t]);
    }

    uint64_t peak = 0;
    int64_t live = 0;
    for (uint64_t step = 0; step < n; step++) {
        live += delta[step];
        peak = std::max(peak,
                        static_cast<uint64_t>(live) + extraSizes[order[step]]);
    }
    return peak;
}

uint64_t
HyperGraph::computeTransfers(const std::vector<uint64_t> &order) const {
    uint64_t n = order.size();
    std::vector<uint64_t> pos(size);
    for (uint64_t step = 0; step < n; step++)
        pos[order[step]] = step;

    // The producer and the consumers of a tensor form runs of consecutive
    // steps (a live-out tensor is consumed at step n, after the last step).
    // Every run after the first is a gap: one fetch each, all served by one
    // offload
    uint64_t transfers = 0;
    std::vector<uint64_t> steps;
    for (uint64_t t = 0; t < tensorSizes.size(); t++) {
        steps.assign(1, pos[tensorProducers[t]]);
        for (uint64_t user : tensorUsers[t])
            steps.push_back(pos[user]);
        if (liveOut[t])
            steps.push_back(n);
        std::sort(steps.begin(), steps.end());
        uint64_t gaps = 0;
        for (uint64_t k = 1; k < steps.size(); k++)
            gaps += steps[k] != steps[k - 1] + 1;
        if (gaps != 0)
            transfers += tensorSizes[t] * (1 + gaps);
    }
    return transfers;
}

void HyperGraph::print(std::ostream &os) const {
    os << "Nodes: " << size << ", tensors: " << tensorSizes.size() << std::endl;

    for (uint64_t i = 0; i < size; i++) {
        os << "\nNode: " << i << std::endl;
        os << "    extraSizes=" << extraSizes[i] << std::endl;

        os << "    inputs={";
        for (auto it = inputTensors[i].begin(); it != inputTensors[i].end();
             ++it) {
            if (it != inputTensors[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}, outputs={";
        for (uint64_t k = 0; k < outputTensors[i].size(); ++k) {
            uint64_t t = outputTensors[i][k];
            if (k != 0)
                os << ", ";
            os << t << ":" << tensorSizes[t] << (liveOut[t] ? "(out)" : "")
               << "->[";
            for (auto it = tensorUsers[t].begin(); it != tensorUsers[t].end();
                 ++it) {
                if (it != tensorUsers[t].begin())
                    os << " ";
                os << *it;
            }
            os << "]";
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

uint64_t HyperGraph::getNumTensors() const { return tensorSizes.size(); }

const std::vector<uint64_t> &HyperGraph::getTensorProducers() const {
    return tensorProducers;
}

const std::vector<bool> &HyperGraph::getLiveOut() const { return liveOut; }

const std::vector<std::vector<uint64_t>> &HyperGraph::getOutputTensors() const {
    return outputTensors;
}

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