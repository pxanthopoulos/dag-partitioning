/**
 * @file Refinement.cpp
 * @brief Implementation of base refinement class methods
 */

#include "Refinement.h"
#include <cassert>
#include <iostream>

Refinement::Refinement(const Graph &graph,
                       std::vector<uint8_t> &initialBisectionInfo,
                       uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
                       double upperBoundPartWeight, double lowerBoundPartWeight)
    : workingGraph(graph), initialBisectionInfo(initialBisectionInfo),
      initialEdgeCut(initialEdgeCut), maxNumberOfPasses(maxNumberOfPasses),
      upperBoundPartWeight(upperBoundPartWeight),
      lowerBoundPartWeight(lowerBoundPartWeight) {}

bool Refinement::checkValidBisection(const std::vector<uint8_t> &bisectionInfo,
                                     const Graph &graph) {
    // Check each edge to ensure no V1->V0 connections exist
    for (uint64_t i = 0; i < graph.size; ++i) {
        const auto &neighbors = graph.adj[i];
        for (const auto &[neighborId, _] : neighbors) {
            // If edge from V1 (true) to V0 (false) found, bisection is invalid
            if (bisectionInfo[i] == 1 && bisectionInfo[neighborId] == 0)
                return false;
        }
    }
    return true;
}

std::pair<uint64_t, uint64_t>
Refinement::calculatePartSizes(const std::vector<uint8_t> &bisectionInfo,
                               const Graph &graph) {
    uint64_t sizeV0 = 0, sizeV1 = 0;

    // Sum weights for each partition
    for (uint64_t i = 0; i < bisectionInfo.size(); ++i) {
        if (bisectionInfo[i] == 0)
            sizeV0 += graph.nodeWeights[i];
        else
            sizeV1 += graph.nodeWeights[i];
    }
    return {sizeV0, sizeV1};
}

bool Refinement::checkBalance(const std::vector<uint8_t> &bisectionInfo,
                              const Graph &graph, uint64_t maxNodeWeight,
                              double upperBoundPartWeight,
                              double lowerBoundPartWeight) {
    // Get current partition weights
    auto [sizeV0, sizeV1] = calculatePartSizes(bisectionInfo, graph);

    // Check balance constraints with allowance for heaviest node
    if ((double)sizeV0 < lowerBoundPartWeight - (double)maxNodeWeight ||
        (double)sizeV0 > upperBoundPartWeight + (double)maxNodeWeight ||
        (double)sizeV1 < lowerBoundPartWeight - (double)maxNodeWeight ||
        (double)sizeV1 > upperBoundPartWeight + (double)maxNodeWeight)
        return false;

    return true;
}

bool Refinement::checkBalance(uint64_t sizeV0, uint64_t sizeV1,
                              uint64_t maxNodeWeight,
                              double upperBoundPartWeight,
                              double lowerBoundPartWeight) {
    // Check balance constraints with allowance for heaviest node
    if ((double)sizeV0 < lowerBoundPartWeight - (double)maxNodeWeight ||
        (double)sizeV0 > upperBoundPartWeight + (double)maxNodeWeight ||
        (double)sizeV1 < lowerBoundPartWeight - (double)maxNodeWeight ||
        (double)sizeV1 > upperBoundPartWeight + (double)maxNodeWeight)
        return false;

    return true;
}

bool Refinement::checkValidEdgeCut(const std::vector<uint8_t> &bisectionInfo,
                                   const Graph &graph,
                                   uint64_t currentEdgeCut) {
    uint64_t edgeCut = 0;
    // Check each edge to ensure no V1->V0 connections exist
    for (uint64_t i = 0; i < graph.size; ++i) {
        const auto &neighbors = graph.adj[i];
        for (const auto &[neighborId, edgeWeight] : neighbors) {
            // If edge from V1 (true) to V0 (false) found, bisection is invalid
            if (bisectionInfo[i] != bisectionInfo[neighborId])
                edgeCut += edgeWeight;
        }
    }
    if (currentEdgeCut != edgeCut)
        return false;
    return true;
}

void Refinement::run() {
    uint64_t countPasses = 0;

    // Continue refinement until max passes reached or no improvement possible
    while (countPasses < maxNumberOfPasses) {
        // Verify acyclicity maintained
        assert(checkValidBisection(initialBisectionInfo, workingGraph) &&
               "Bisection is invalid, it has edge from V1 to V0");
        // Verify that the edge cut is consistent with the bisection info
        assert(checkValidEdgeCut(initialBisectionInfo, workingGraph,
                                 initialEdgeCut) &&
               "Computed edge cut is invalid");

        // Attempt one pass of refinement
        uint64_t initialEdgeCutOld = initialEdgeCut;
        if (!onePassRefinement())
            break; // Stop if no improvement made
        if ((double)initialEdgeCut >= 0.99 * (double)initialEdgeCutOld)
            break;

        countPasses++;
    }
    // Verify acyclicity maintained
    assert(checkValidBisection(initialBisectionInfo, workingGraph) &&
           "Bisection is invalid, it has edge from V1 to V0");
    // Verify that the edge cut is consistent with the bisection info
    assert(
        checkValidEdgeCut(initialBisectionInfo, workingGraph, initialEdgeCut) &&
        "Computed edge cut is invalid");
}