/**
 * @file Bisection.cpp
 * @brief Implementation of base bisection class methods
 */

#include "Bisection.h"

Bisection::Bisection(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight,
                     RefinementMethod refinementMethod, uint64_t refinementPasses)
        : workingGraph(graph), upperBoundPartWeight(upperBoundPartWeight), lowerBoundPartWeight(lowerBoundPartWeight),
          refinementMethod(refinementMethod), refinementPasses(refinementPasses) {}

bool Bisection::checkValidBisection(const std::vector<uint8_t> &bisection) const {
    // Check each vertex and its edges for acyclicity violation
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, _]: neighbors) {
            // If there's an edge from V1 (1) to V0 (0), bisection is invalid
            if (bisection[i] == 1 && bisection[neighborId] == 0) return false;
        }
    }
    return true;
}

uint64_t Bisection::computeEdgeCut(const std::vector<uint8_t> &bisection) const {
    uint64_t edgeCut = 0;

    // Check each edge to see if it crosses between partitions
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, edgeWeight]: neighbors) {
            // Add edge weight to cut if vertices are in different partitions
            if (bisection[i] != bisection[neighborId]) edgeCut += edgeWeight;
        }
    }
    return edgeCut;
}