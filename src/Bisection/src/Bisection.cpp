/**
 * @file Bisection.cpp
 * @brief Implementation of base bisection class methods
 */

#include "Bisection.h"

Bisection::Bisection(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight,
                     RefinementMethod refinementMethod, uint64_t refinementPasses)
        : workingGraph(graph), upperBoundPartWeight(upperBoundPartWeight), lowerBoundPartWeight(lowerBoundPartWeight),
          refinementMethod(refinementMethod), refinementPasses(refinementPasses) {}

uint64_t Bisection::selectBestResult(const std::vector<std::tuple<uint64_t, uint8_t, double>> &results) {
    using Result = std::tuple<uint64_t, uint8_t, double>;  // (edge_cut, is_balanced, imbalance)
    std::vector<std::pair<uint64_t, Result>> balancedResults;
    std::vector<std::pair<uint64_t, Result>> unbalancedResults;

    // Separate balanced and unbalanced results, keeping track of original indices
    for (uint64_t i = 0; i < results.size(); ++i) {
        if (std::get<1>(results[i])) {  // is_balanced
            balancedResults.emplace_back(i, results[i]);
        } else {
            unbalancedResults.emplace_back(i, results[i]);
        }
    }

    // If we have balanced results
    if (!balancedResults.empty()) {
        // Check if all balanced results have zero edge cut
        bool allZero = std::all_of(balancedResults.begin(), balancedResults.end(),
                                   [](const auto &pair) { return std::get<0>(pair.second) == 0; });

        if (allZero) {
            // If all are zero, return index of first balanced result
            return balancedResults[0].first;
        } else {
            // Find the smallest non-zero edge cut among balanced results
            size_t bestIndex = balancedResults[0].first;
            uint64_t minNonZeroCut = std::numeric_limits<uint64_t>::max();

            for (const auto &[idx, result]: balancedResults) {
                uint64_t edgeCut = std::get<0>(result);
                if (edgeCut != 0 && edgeCut < minNonZeroCut) {
                    minNonZeroCut = edgeCut;
                    bestIndex = idx;
                }
            }
            return bestIndex;
        }
    }

    // If all results are unbalanced, find the least unbalanced one
    // with tie breaks favoring lower edge cuts
    size_t bestIndex = unbalancedResults[0].first;

    for (const auto &[idx, result]: unbalancedResults) {
        const auto &currentImbalance = std::get<2>(result);
        const auto &bestImbalance = std::get<2>(results[bestIndex]);

        if (currentImbalance < bestImbalance ||
            (currentImbalance == bestImbalance &&
             std::get<0>(result) < std::get<0>(results[bestIndex]))) {
            bestIndex = idx;
        }
    }

    return bestIndex;
}

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