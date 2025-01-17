/**
 * @file Refinement.cpp
 * @brief Implementation of base refinement class methods
 */

#include <iostream>
#include "Refinement.h"

Refinement::Refinement(const Graph &graph, std::vector<uint8_t> &initialBisectionInfo, uint64_t &initialEdgeCut,
                       uint64_t maxNumberOfPasses, double upperBoundPartWeight, double lowerBoundPartWeight)
        : workingGraph(graph),
          initialBisectionInfo(initialBisectionInfo),
          initialEdgeCut(initialEdgeCut),
          maxNumberOfPasses(maxNumberOfPasses),
          upperBoundPartWeight(upperBoundPartWeight),
          lowerBoundPartWeight(lowerBoundPartWeight) {}

bool Refinement::checkValidBisection() const {
    // Check each edge to ensure no V1->V0 connections exist
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, _]: neighbors) {
            // If edge from V1 (true) to V0 (false) found, bisection is invalid
            if (initialBisectionInfo[i] == 1 && initialBisectionInfo[neighborId] == 0)
                return false;
        }
    }
    return true;
}

std::pair<uint64_t, uint64_t> Refinement::calculatePartSizes() const {
    uint64_t sizeV0 = 0, sizeV1 = 0;

    // Sum weights for each partition
    for (uint64_t i = 0; i < initialBisectionInfo.size(); ++i) {
        if (initialBisectionInfo[i] == 0)
            sizeV0 += workingGraph.nodeWeights[i];
        else
            sizeV1 += workingGraph.nodeWeights[i];
    }
    return {sizeV0, sizeV1};
}

bool Refinement::checkBalance(uint64_t maxNodeWeight) const {
    // Get current partition weights
    auto [sizeV0, sizeV1] = calculatePartSizes();

    // Check balance constraints with allowance for heaviest node
    if ((double) sizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
        (double) sizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
        (double) sizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
        (double) sizeV1 > upperBoundPartWeight + (double) maxNodeWeight)
        return false;

    return true;
}

bool Refinement::checkValidEdgeCut() {
    uint64_t edgeCut = 0;
    // Check each edge to ensure no V1->V0 connections exist
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, edgeWeight]: neighbors) {
            // If edge from V1 (true) to V0 (false) found, bisection is invalid
            if (initialBisectionInfo[i] != initialBisectionInfo[neighborId]) edgeCut += edgeWeight;
        }
    }
    if (initialEdgeCut != edgeCut) return false;
    return true;
}

void Refinement::run() {
    uint64_t countPasses = 0;

    // Continue refinement until max passes reached or no improvement possible
    while (countPasses < maxNumberOfPasses) {
        // Verify acyclicity maintained
        assert(checkValidBisection() && "Bisection is invalid, it has edge from V1 to V0");

        // Attempt one pass of refinement
        if (!onePassRefinement()) break;  // Stop if no improvement made

        countPasses++;
    }
    // Verify acyclicity maintained
    assert(checkValidBisection() && "Bisection is invalid, it has edge from V1 to V0");
    
    // Verify that the edge cut is consistent with the bisection info
    assert(checkValidEdgeCut());
}