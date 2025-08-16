/**
 * @file GreedyDirectedGraphGrowing.cpp
 * @brief Implementation of greedy directed graph growing bisection algorithm
 */

#include "GreedyDirectedGraphGrowing.h"
#include "Graph.h"
#include "Refinement.h"
#include "RefinementWrapper.h"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <queue>
#include <random>

namespace dag_partitioning {

namespace bisection {

GreedyDirectedGraphGrowing::GreedyDirectedGraphGrowing(
    const core::Graph &graph, double upperBoundPartWeight,
    double lowerBoundPartWeight, refinement::RefinementMethod refinementMethod,
    uint64_t refinementPasses)
    : Bisection(graph, upperBoundPartWeight, lowerBoundPartWeight,
                refinementMethod, refinementPasses) {}

/**
 * @brief Comparator for priority queue ordering
 *
 * Orders vertices by:
 * 1. In-edge weight (bigger weight -> bigger priority)
 * 2. Distance from first moved vertex (smaller distance -> bigger priority)
 * 3. Node ID (bigger ID -> bigger priority)
 */
struct HeapComparator {
    bool operator()(const std::tuple<uint64_t, uint64_t, uint64_t> &a,
                    const std::tuple<uint64_t, uint64_t, uint64_t> &b) {
        if (std::get<0>(a) != std::get<0>(b))
            return std::get<0>(a) < std::get<0>(b); // Edge weight priority
        if (std::get<1>(a) != std::get<1>(b))
            return std::get<1>(a) > std::get<1>(b); // Distance priority
        return std::get<2>(a) < std::get<2>(b);     // Node ID
    }
};

std::pair<std::vector<uint8_t>, uint64_t>
GreedyDirectedGraphGrowing::runOnNormalGraph() const {
    uint64_t bestEdgeCut = UINT64_MAX;
    auto bestImbalance = DBL_MAX;
    uint64_t currentEdgeCut = 0;
    uint64_t partWeightV0 = 0;
    uint64_t partWeightV1 = workingGraph.totalWeight;
    uint64_t bestMovePrefix = 0;
    std::vector<uint64_t> moveSequence;
    bool isValid = false;

    // Data structures for two-phase algorithm
    std::priority_queue<std::tuple<uint64_t, uint64_t, uint64_t>,
                        std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>,
                        HeapComparator>
        maxHeapPhase1;
    std::vector<uint64_t> localInDegree = workingGraph.inDegree;
    std::vector<uint64_t> sources;
    std::vector<uint64_t> weightedInDegree(workingGraph.size, 0);
    std::vector<uint64_t> weightedOutDegree(workingGraph.size, 0);
    uint64_t maxNodeWeight = workingGraph.maxNodeWeight;

    // Find sources and compute weighted degrees
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        if (localInDegree[i] == 0) {
            sources.emplace_back(i);
        }
        // Compute weighted in-degree and out-degree
        for (const auto &[neighborId, edgeWeight] : workingGraph.adj[i]) {
            weightedInDegree[neighborId] += edgeWeight;
            weightedOutDegree[i] += edgeWeight;
        }
    }

    // Randomly order source vertices
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(sources.begin(), sources.end(), gen);

    // Get distances from first source for tie-breaking
    std::vector<uint64_t> distancesFromFirstMovedNode =
        workingGraph.distancesFromNode(sources[0], false);
    for (const uint64_t &source : sources) {
        maxHeapPhase1.emplace(weightedInDegree[source],
                              distancesFromFirstMovedNode[source], source);
    }

    // Phase 1: Move vertices using weighted in-degree priority
    // Continue until V0 reaches 90% of max allowed weight
    while (!maxHeapPhase1.empty() &&
           ((double)partWeightV0 <=
            0.9 * (upperBoundPartWeight + (double)maxNodeWeight)) &&
           ((double)partWeightV1 >=
            1.1 * (lowerBoundPartWeight - (double)maxNodeWeight))) {
        const auto [_, distance, nodeId] = maxHeapPhase1.top();
        maxHeapPhase1.pop();

        // Move node to V0
        currentEdgeCut += weightedOutDegree[nodeId];
        currentEdgeCut -= weightedInDegree[nodeId];
        partWeightV0 += workingGraph.nodeWeights[nodeId];
        partWeightV1 -= workingGraph.nodeWeights[nodeId];
        moveSequence.emplace_back(nodeId);

        // Check neighbors for new movable vertices
        for (const auto &[neighborId, edgeWeight] : workingGraph.adj[nodeId]) {
            localInDegree[neighborId]--;
            assert(localInDegree[neighborId] >= 0 &&
                   "in-degree of node became negative");
            if (localInDegree[neighborId] == 0)
                maxHeapPhase1.emplace(weightedInDegree[neighborId],
                                      distancesFromFirstMovedNode[neighborId],
                                      neighborId);
        }

        double maxImbalance =
            std::max((double)partWeightV0 / upperBoundPartWeight,
                     (double)partWeightV1 / upperBoundPartWeight);

        // If no moves yet and balance is improved, update best solution
        if (bestMovePrefix == 0 && maxImbalance < bestImbalance) {
            bestEdgeCut = currentEdgeCut;
            bestImbalance = maxImbalance;
            bestMovePrefix = moveSequence.size();
        }

        // If balance is obtained
        if (((double)partWeightV0 <
             upperBoundPartWeight + (double)maxNodeWeight) &&
            ((double)partWeightV1 <
             upperBoundPartWeight + (double)maxNodeWeight)) {
            // If balance is obtained for the first time
            // Or edge cut is improved
            // Or edge cut is the same but balance is improved
            if (!isValid || (currentEdgeCut < bestEdgeCut ||
                             (currentEdgeCut == bestEdgeCut &&
                              maxImbalance < bestImbalance))) {
                bestEdgeCut = currentEdgeCut;
                bestImbalance = maxImbalance;
                isValid = true;
                bestMovePrefix = moveSequence.size();
            }
        }
    }

    // Prepare for Phase 2
    std::vector<uint64_t> readyNodes;
    while (!maxHeapPhase1.empty()) {
        const auto [_, distance, nodeId] = maxHeapPhase1.top();
        readyNodes.push_back(nodeId);
        maxHeapPhase1.pop();
    }

    // Phase 2: Use gain (in-degree - out-degree) as priority
    std::priority_queue<std::tuple<uint64_t, uint64_t, uint64_t>> maxHeapPhase2;
    for (uint64_t nodeId : readyNodes) {
        maxHeapPhase2.emplace(weightedInDegree[nodeId] -
                                  weightedOutDegree[nodeId],
                              workingGraph.nodeWeights[nodeId], nodeId);
    }

    // Continue until reaching balance constraints
    while (!maxHeapPhase2.empty() &&
           ((double)partWeightV0 <=
            upperBoundPartWeight + (double)maxNodeWeight) &&
           ((double)partWeightV1 >=
            lowerBoundPartWeight - (double)maxNodeWeight)) {
        const auto [_, distance, nodeId] = maxHeapPhase2.top();
        maxHeapPhase2.pop();

        // Move node to V0
        currentEdgeCut += weightedOutDegree[nodeId];
        currentEdgeCut -= weightedInDegree[nodeId];
        partWeightV0 += workingGraph.nodeWeights[nodeId];
        partWeightV1 -= workingGraph.nodeWeights[nodeId];
        moveSequence.emplace_back(nodeId);

        // Check neighbors for new movable vertices
        for (const auto &[neighborId, edgeWeight] : workingGraph.adj[nodeId]) {
            localInDegree[neighborId]--;
            assert(localInDegree[neighborId] >= 0 &&
                   "in-degree of node became negative");
            if (localInDegree[neighborId] == 0)
                maxHeapPhase2.emplace(weightedInDegree[neighborId] -
                                          weightedOutDegree[neighborId],
                                      workingGraph.nodeWeights[neighborId],
                                      neighborId);
        }

        // If balance is obtained
        if (((double)partWeightV0 <
             upperBoundPartWeight + (double)maxNodeWeight) &&
            ((double)partWeightV1 <
             upperBoundPartWeight + (double)maxNodeWeight)) {
            double maxImbalance =
                std::max((double)partWeightV0 / upperBoundPartWeight,
                         (double)partWeightV1 / upperBoundPartWeight);

            // If balance is obtained for the first time
            // Or edge cut is improved
            // Or edge cut is the same but balance is improved
            if (!isValid || (currentEdgeCut < bestEdgeCut ||
                             (currentEdgeCut == bestEdgeCut &&
                              maxImbalance < bestImbalance))) {
                bestEdgeCut = currentEdgeCut;
                bestImbalance = maxImbalance;
                isValid = true;
                bestMovePrefix = moveSequence.size();
            }
        }
    }

    assert(bestMovePrefix > 0 && "No possible moves found");

    // Initialize all nodes in V1
    std::vector<uint8_t> bisection(workingGraph.size, 1);
    for (uint64_t i = 0; i < bestMovePrefix; ++i)
        bisection[moveSequence[i]] = 0;

    return {bisection, bestEdgeCut};
}

std::pair<std::vector<uint8_t>, uint64_t>
GreedyDirectedGraphGrowing::runOnReverseGraph() const {
    // Similar to runOnNormalGraph but:
    // 1. Starts with all vertices in V0
    // 2. Uses out-degree instead of in-degree
    // 3. Moves vertices to V1
    // Structure and logic mirrors runOnNormalGraph with reversed directions
    uint64_t bestEdgeCut = UINT64_MAX;
    auto bestImbalance = DBL_MAX;
    uint64_t currentEdgeCut = 0;
    uint64_t partWeightV0 = workingGraph.totalWeight;
    uint64_t partWeightV1 = 0;
    uint64_t bestMovePrefix = 0;
    std::vector<uint64_t> moveSequence;
    bool isValid = false;

    std::priority_queue<std::tuple<uint64_t, uint64_t, uint64_t>,
                        std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>,
                        HeapComparator>
        maxHeapPhase1;
    std::vector<uint64_t> localOutDegree(workingGraph.size);
    std::vector<uint64_t> sinks;
    std::vector<uint64_t> weightedInDegree(workingGraph.size, 0);
    std::vector<uint64_t> weightedOutDegree(workingGraph.size, 0);
    uint64_t maxNodeWeight = workingGraph.maxNodeWeight;

    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        localOutDegree[i] = workingGraph.adj[i].size();
        if (localOutDegree[i] == 0) {
            sinks.emplace_back(i);
        }

        for (const auto &[neighborId, edgeWeight] : workingGraph.adj[i]) {
            weightedInDegree[neighborId] += edgeWeight;
            weightedOutDegree[i] += edgeWeight;
        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(sinks.begin(), sinks.end(), gen);

    std::vector<uint64_t> distancesFromFirstMovedNode =
        workingGraph.distancesFromNode(sinks[0], true);
    for (const uint64_t &sink : sinks) {
        maxHeapPhase1.emplace(weightedOutDegree[sink],
                              distancesFromFirstMovedNode[sink], sink);
    }

    while (!maxHeapPhase1.empty() &&
           ((double)partWeightV1 <=
            0.9 * (upperBoundPartWeight + (double)maxNodeWeight)) &&
           ((double)partWeightV0 >=
            1.1 * (lowerBoundPartWeight - (double)maxNodeWeight))) {
        const auto [_, distance, nodeId] = maxHeapPhase1.top();
        maxHeapPhase1.pop();

        currentEdgeCut -= weightedOutDegree[nodeId];
        currentEdgeCut += weightedInDegree[nodeId];
        partWeightV0 -= workingGraph.nodeWeights[nodeId];
        partWeightV1 += workingGraph.nodeWeights[nodeId];
        moveSequence.emplace_back(nodeId);

        for (const auto &[neighborId, edgeWeight] :
             workingGraph.revAdj[nodeId]) {
            localOutDegree[neighborId]--;
            assert(localOutDegree[neighborId] >= 0 &&
                   "out-degree of node became negative");
            if (localOutDegree[neighborId] == 0)
                maxHeapPhase1.emplace(weightedOutDegree[neighborId],
                                      distancesFromFirstMovedNode[neighborId],
                                      neighborId);
        }

        double maxImbalance =
            std::max((double)partWeightV0 / upperBoundPartWeight,
                     (double)partWeightV1 / upperBoundPartWeight);

        if (bestMovePrefix == 0 && maxImbalance < bestImbalance) {
            bestEdgeCut = currentEdgeCut;
            bestImbalance = maxImbalance;
            bestMovePrefix = moveSequence.size();
        }

        if (((double)partWeightV0 <
             upperBoundPartWeight + (double)maxNodeWeight) &&
            ((double)partWeightV1 <
             upperBoundPartWeight + (double)maxNodeWeight)) {
            if (!isValid || (currentEdgeCut < bestEdgeCut ||
                             (currentEdgeCut == bestEdgeCut &&
                              maxImbalance < bestImbalance))) {
                bestEdgeCut = currentEdgeCut;
                bestImbalance = maxImbalance;
                isValid = true;
                bestMovePrefix = moveSequence.size();
            }
        }
    }

    std::vector<uint64_t> readyNodes;
    while (!maxHeapPhase1.empty()) {
        const auto [_, distance, nodeId] = maxHeapPhase1.top();
        readyNodes.push_back(nodeId);
        maxHeapPhase1.pop();
    }
    std::priority_queue<std::tuple<uint64_t, uint64_t, uint64_t>> maxHeapPhase2;
    for (uint64_t nodeId : readyNodes) {
        maxHeapPhase2.emplace(weightedOutDegree[nodeId] -
                                  weightedInDegree[nodeId],
                              workingGraph.nodeWeights[nodeId], nodeId);
    }

    while (!maxHeapPhase2.empty() &&
           ((double)partWeightV1 <=
            upperBoundPartWeight + (double)maxNodeWeight) &&
           ((double)partWeightV0 >=
            lowerBoundPartWeight - (double)maxNodeWeight)) {
        const auto [_, distance, nodeId] = maxHeapPhase2.top();
        maxHeapPhase2.pop();

        currentEdgeCut -= weightedOutDegree[nodeId];
        currentEdgeCut += weightedInDegree[nodeId];
        partWeightV0 -= workingGraph.nodeWeights[nodeId];
        partWeightV1 += workingGraph.nodeWeights[nodeId];
        moveSequence.emplace_back(nodeId);

        for (const auto &[neighborId, edgeWeight] :
             workingGraph.revAdj[nodeId]) {
            localOutDegree[neighborId]--;
            assert(localOutDegree[neighborId] >= 0 &&
                   "out-degree of node became negative");
            if (localOutDegree[neighborId] == 0)
                maxHeapPhase2.emplace(weightedOutDegree[neighborId] -
                                          weightedInDegree[neighborId],
                                      workingGraph.nodeWeights[neighborId],
                                      neighborId);
        }

        if (((double)partWeightV0 <
             upperBoundPartWeight + (double)maxNodeWeight) &&
            ((double)partWeightV1 <
             upperBoundPartWeight + (double)maxNodeWeight)) {
            double maxImbalance =
                std::max((double)partWeightV0 / upperBoundPartWeight,
                         (double)partWeightV1 / upperBoundPartWeight);

            if (!isValid || (currentEdgeCut < bestEdgeCut ||
                             (currentEdgeCut == bestEdgeCut &&
                              maxImbalance < bestImbalance))) {
                bestEdgeCut = currentEdgeCut;
                bestImbalance = maxImbalance;
                isValid = true;
                bestMovePrefix = moveSequence.size();
            }
        }
    }

    assert(bestMovePrefix > 0 && "No possible moves found");

    // Initialize all nodes in V1
    std::vector<uint8_t> bisection(workingGraph.size, 0);
    for (uint64_t i = 0; i < bestMovePrefix; ++i)
        bisection[moveSequence[i]] = 1;

    return {bisection, bestEdgeCut};
}

std::pair<std::vector<uint8_t>, uint64_t>
GreedyDirectedGraphGrowing::run() const {
    std::vector<std::pair<std::vector<uint8_t>, uint64_t>> bisections;
    std::vector<std::tuple<uint64_t, uint8_t, double>> results;

    // Run both normal and reverse algorithms
    std::pair<std::vector<uint8_t>, uint64_t> bisectionNormal =
        runOnNormalGraph();
    refinement::refinementWrapper(workingGraph, bisectionNormal.first,
                                  bisectionNormal.second, refinementMethod,
                                  refinementPasses, upperBoundPartWeight,
                                  lowerBoundPartWeight);
    uint8_t isZero = (bisectionNormal.second == 0) ? 0 : 1;
    auto [sizeV0, sizeV1] = refinement::Refinement::calculatePartSizes(
        bisectionNormal.first, workingGraph);
    double imbalance =
        (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
         upperBoundPartWeight) *
        100;
    results.emplace_back(bisectionNormal.second, isZero, imbalance);
    bisections.emplace_back(bisectionNormal);

    std::pair<std::vector<uint8_t>, uint64_t> bisectionReverse =
        runOnReverseGraph();
    refinement::refinementWrapper(workingGraph, bisectionReverse.first,
                                  bisectionReverse.second, refinementMethod,
                                  refinementPasses, upperBoundPartWeight,
                                  lowerBoundPartWeight);
    isZero = (bisectionReverse.second == 0) ? 0 : 1;
    std::tie(sizeV0, sizeV1) = refinement::Refinement::calculatePartSizes(
        bisectionReverse.first, workingGraph);
    imbalance = (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
                 upperBoundPartWeight) *
                100;
    results.emplace_back(bisectionReverse.second, isZero, imbalance);
    bisections.emplace_back(bisectionReverse);

    // Verify both solutions maintain acyclicity and have valid edge cuts
    assert(refinement::Refinement::checkValidBisection(bisectionNormal.first,
                                                       workingGraph) &&
           "Bisection on normal graph is invalid");
    assert(refinement::Refinement::checkValidEdgeCut(
               bisectionNormal.first, workingGraph, bisectionNormal.second) &&
           "Edge cut on normal graph is invalid");
    assert(refinement::Refinement::checkValidBisection(bisectionReverse.first,
                                                       workingGraph) &&
           "Bisection on reverse graph is invalid");
    assert(refinement::Refinement::checkValidEdgeCut(
               bisectionReverse.first, workingGraph, bisectionReverse.second) &&
           "Edge cut on reverse graph is invalid");

    return bisections[selectBestResult(results)];
}

} // namespace bisection

} // namespace dag_partitioning