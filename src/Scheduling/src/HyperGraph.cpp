/**
 * @file HyperGraph.cpp
 * @brief Hyper graph of nodes and the tensors between them, for scheduling
 */

#include "HyperGraph.h"

#include <algorithm>
#include <queue>

namespace dag_partitioning {

namespace scheduling {

namespace hypergraph {

HyperGraph::HyperGraph(std::vector<uint64_t> extraSizes,
                       const std::vector<Tensor> &tensors,
                       std::vector<LiveIn> liveIns)
    : size(extraSizes.size()), liveIns(std::move(liveIns)),
      extraSizes(std::move(extraSizes)), ancestors(size), descendants(size),
      inputTensors(size), tensorUsers(tensors.size()) {
    for (uint64_t t = 0; t < tensors.size(); ++t) {
        tensorSizes.push_back(tensors[t].size);
        tensorProducers.push_back(tensors[t].producer);
        liveOut.push_back(tensors[t].liveOut);
        for (uint64_t user : tensors[t].users) {
            tensorUsers[t].insert(user);
            inputTensors[user].insert(t);
        }
    }

    // Ancestors in topological order (Kahn's algorithm): the predecessors
    // of a node and their ancestors
    std::vector<robin_hood::unordered_set<uint64_t>> predecessors(size);
    std::vector<std::vector<uint64_t>> successors(size);
    for (uint64_t i = 0; i < size; i++) {
        for (uint64_t t : inputTensors[i]) {
            if (predecessors[i].insert(tensorProducers[t]).second)
                successors[tensorProducers[t]].push_back(i);
        }
    }
    std::vector<uint64_t> inDegree(size);
    std::queue<uint64_t> ready;
    for (uint64_t i = 0; i < size; i++) {
        inDegree[i] = predecessors[i].size();
        if (inDegree[i] == 0)
            ready.push(i);
    }
    while (!ready.empty()) {
        uint64_t node = ready.front();
        ready.pop();
        for (uint64_t pred : predecessors[node]) {
            ancestors[node].insert(pred);
            ancestors[node].insert(ancestors[pred].begin(),
                                   ancestors[pred].end());
        }
        for (uint64_t succ : successors[node]) {
            if (--inDegree[succ] == 0)
                ready.push(succ);
        }
    }
    for (uint64_t i = 0; i < size; i++) {
        for (uint64_t ancestor : ancestors[i])
            descendants[ancestor].insert(i);
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
    // A live-in tensor adds its size to the steps [0, its last user]
    for (const LiveIn &liveIn : liveIns) {
        if (liveIn.users.empty())
            continue;
        uint64_t death = 0;
        for (uint64_t user : liveIn.users)
            death = std::max(death, pos[user]);
        delta[0] += static_cast<int64_t>(liveIn.size);
        delta[death + 1] -= static_cast<int64_t>(liveIn.size);
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

std::ostream &operator<<(std::ostream &os, const HyperGraph &graph) {
    os << graph.size << " nodes, extra sizes:";
    for (uint64_t extra : graph.extraSizes)
        os << " " << extra;
    os << "\n";
    for (uint64_t t = 0; t < graph.tensorSizes.size(); t++) {
        std::vector<uint64_t> users(graph.tensorUsers[t].begin(),
                                    graph.tensorUsers[t].end());
        std::sort(users.begin(), users.end());
        os << "tensor " << t << ": " << graph.tensorSizes[t] << " bytes, "
           << graph.tensorProducers[t] << " ->";
        for (uint64_t user : users)
            os << " " << user;
        os << (graph.liveOut[t] ? " (live-out)\n" : "\n");
    }
    for (const LiveIn &liveIn : graph.liveIns) {
        os << "live-in: " << liveIn.size << " bytes ->";
        for (uint64_t user : liveIn.users)
            os << " " << user;
        os << "\n";
    }
    return os;
}

} // namespace hypergraph

} // namespace scheduling

} // namespace dag_partitioning
