/**
 * @file Bisection.cpp
 * @brief Implementation of base bisection class methods
 */

#include "Bisection.h"

Bisection::Bisection(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight,
                     RefinementMethod refinementMethod, uint64_t refinementPasses)
        : workingGraph(graph), upperBoundPartWeight(upperBoundPartWeight), lowerBoundPartWeight(lowerBoundPartWeight),
          refinementMethod(refinementMethod), refinementPasses(refinementPasses) {}

uint64_t Bisection::computeEdgeCut(const std::vector<uint8_t> &bisection, const Graph &graph) {
    uint64_t edgeCut = 0;

    // Check each edge to see if it crosses between partitions
    for (uint64_t i = 0; i < graph.size; ++i) {
        const auto &neighbors = graph.adj[i];
        for (const auto &[neighborId, edgeWeight]: neighbors) {
            // Add edge weight to cut if vertices are in different partitions
            if (bisection[i] != bisection[neighborId]) edgeCut += edgeWeight;
        }
    }
    return edgeCut;
}