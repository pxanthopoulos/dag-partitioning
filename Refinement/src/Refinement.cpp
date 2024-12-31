//
// Created by panagiotis on 22/12/2024.
//

#include "Refinement.h"

Refinement::Refinement(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t &initialEdgeCut,
                       uint64_t maxNumberOfPasses,
                       double upperBoundPartWeight, double lowerBoundPartWeight)
        : workingGraph(graph), initialBisectionInfo(initialBisectionInfo), initialEdgeCut(initialEdgeCut),
          maxNumberOfPasses(maxNumberOfPasses),
          upperBoundPartWeight(upperBoundPartWeight), lowerBoundPartWeight(lowerBoundPartWeight) {}

bool Refinement::checkValidBisection() const {
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, _]: neighbors) {
            if (initialBisectionInfo[i] && !initialBisectionInfo[neighborId]) return false;
        }
    }
    return true;
}

std::pair<uint64_t, uint64_t> Refinement::calculatePartSizes() const {
    uint64_t sizeV0 = 0, sizeV1 = 0;
    for (uint64_t i = 0; i < initialBisectionInfo.size(); ++i) {
        if (!initialBisectionInfo[i]) sizeV0 += workingGraph.nodeWeights[i];
        else sizeV1 += workingGraph.nodeWeights[i];
    }
    return {sizeV0, sizeV1};
}

bool Refinement::checkBalance(uint64_t maxNodeWeight) const {
    auto [sizeV0, sizeV1] = calculatePartSizes();
    if ((double) sizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
        (double) sizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
        (double) sizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
        (double) sizeV1 > upperBoundPartWeight + (double) maxNodeWeight)
        return false;
    return true;
}

void Refinement::run() {
    uint64_t countPasses = 0;
    while (countPasses < maxNumberOfPasses) {
        assert(checkValidBisection() && "Bisection is invalid, it has edge from V1 to V0");
        if (!onePassRefinement()) break;
        countPasses++;
    }
    assert(checkBalance(workingGraph.maxNodeWeight()) && "Resulting partition is unbalanced");
}
