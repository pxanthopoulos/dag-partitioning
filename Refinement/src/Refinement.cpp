//
// Created by panagiotis on 22/12/2024.
//

#include "Refinement.h"

Refinement::Refinement(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t maxNumberOfPasses)
        : workingGraph(graph), initialBisectionInfo(initialBisectionInfo), maxNumberOfPasses(maxNumberOfPasses) {}

bool Refinement::checkValidBisection() const {
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, _]: neighbors) {
            if (initialBisectionInfo[i] && !initialBisectionInfo[neighborId]) return false;
        }
    }
    return true;
}

void Refinement::run() {
    uint64_t countPasses = 0;
    while (countPasses < maxNumberOfPasses) {
        assert(checkValidBisection() && "Bisection is invalid, it has edge from V1 to V0");
        if (!onePassRefinement()) break;
    }
}
