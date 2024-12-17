//
// Created by panagiotis on 16/12/2024.
//

#include "Bisection.h"

Bisection::Bisection(Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight)
        : workingGraph(graph), upperBoundPartWeight(upperBoundPartWeight),
          lowerBoundPartWeight(lowerBoundPartWeight) {}

bool Bisection::checkValidBisection(const std::vector<bool> &bisection) {
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, _]: neighbors) {
            if (bisection[i] && !bisection[neighborId]) return false;
        }
    }
    return true;
}
