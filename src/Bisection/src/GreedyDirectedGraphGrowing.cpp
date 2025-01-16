/**
 * @file GreedyDirectedGraphGrowing.cpp
 * @brief Implementation of greedy directed graph growing bisection algorithm
 */

#include "GreedyDirectedGraphGrowing.h"
#include <queue>
#include <random>
#include <algorithm>
#include <cfloat>

GreedyDirectedGraphGrowing::GreedyDirectedGraphGrowing(const Graph &graph, double upperBoundPartWeight,
                                                       double lowerBoundPartWeight)
        : Bisection(graph, upperBoundPartWeight, lowerBoundPartWeight) {}

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
            return std::get<0>(a) < std::get<0>(b);  // Edge weight priority
        if (std::get<1>(a) != std::get<1>(b))
            return std::get<1>(a) > std::get<1>(b);  // Distance priority
        return std::get<2>(a) < std::get<2>(b);      // Node ID
    }
};

std::pair<std::vector<uint8_t>, uint64_t> GreedyDirectedGraphGrowing::runOnNormalGraph() const {
    // Initialize all vertices in V1
    std::vector<uint8_t> bisection(workingGraph.size, 1);
    std::vector<uint8_t> bestBisection;
    uint64_t bestEdgeCut = UINT64_MAX;
    auto bestImbalanceDiff = DBL_MAX;
    uint64_t currentEdgeCut = 0;
    uint64_t partWeightV0 = 0;
    uint64_t partWeightV1 = workingGraph.totalWeight;

    // Data structures for two-phase algorithm
    std::priority_queue<std::tuple<uint64_t, uint64_t, uint64_t>,
            std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>,
            HeapComparator> maxHeapPhase1;
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
        for (const auto &[reverseNeighborId, edgeWeight]: workingGraph.revAdj[i])
            weightedInDegree[i] += edgeWeight;
        for (const auto &[neighborId, edgeWeight]: workingGraph.adj[i])
            weightedOutDegree[i] += edgeWeight;
    }

    // Randomly order source vertices
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(sources.begin(), sources.end(), gen);

    // Get distances from first source for tie-breaking
    std::vector<uint64_t> distancesFromFirstMovedNode = workingGraph.distancesFromNode(sources[0], false);
    for (const uint64_t &source: sources) {
        maxHeapPhase1.emplace(weightedInDegree[source], distancesFromFirstMovedNode[source], source);
    }

    // Phase 1: Move vertices using weighted in-degree priority
    // Continue until V0 reaches 90% of max allowed weight
    while (!maxHeapPhase1.empty() &&
           ((double) partWeightV0 <= 0.9 * (upperBoundPartWeight + (double) maxNodeWeight)) &&
           ((double) partWeightV1 >= 1.1 * (lowerBoundPartWeight - (double) maxNodeWeight))) {
        const auto [_, distance, nodeId] = maxHeapPhase1.top();
        maxHeapPhase1.pop();

        // Move node to V0
        bisection[nodeId] = 0;
        currentEdgeCut += weightedOutDegree[nodeId];
        currentEdgeCut -= weightedInDegree[nodeId];
        partWeightV0 += workingGraph.nodeWeights[nodeId];
        partWeightV1 -= workingGraph.nodeWeights[nodeId];

        // Check neighbors for new movable vertices
        for (const auto &[neighborId, edgeWeight]: workingGraph.adj[nodeId]) {
            localInDegree[neighborId]--;
            assert(localInDegree[neighborId] >= 0 && "in-degree of node became negative");
            if (localInDegree[neighborId] == 0)
                maxHeapPhase1.emplace(weightedInDegree[neighborId],
                                      distancesFromFirstMovedNode[neighborId],
                                      neighborId);
        }

        // Update best solution if current is better
        if (((double) partWeightV0 < upperBoundPartWeight + (double) maxNodeWeight) &&
            ((double) partWeightV1 < upperBoundPartWeight + (double) maxNodeWeight)) {
            double maxImbalance = std::max((double) partWeightV0 / upperBoundPartWeight,
                                           (double) partWeightV1 / upperBoundPartWeight);

            if (currentEdgeCut < bestEdgeCut ||
                (currentEdgeCut == bestEdgeCut && std::abs(maxImbalance - 1.0) < bestImbalanceDiff)) {
                bestEdgeCut = currentEdgeCut;
                bestImbalanceDiff = maxImbalance;
                bestBisection = bisection;
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
    for (uint64_t nodeId: readyNodes) {
        maxHeapPhase2.emplace(weightedInDegree[nodeId] - weightedOutDegree[nodeId],
                              workingGraph.nodeWeights[nodeId],
                              nodeId);
    }

    // Continue until reaching balance constraints
    while (!maxHeapPhase2.empty() &&
           ((double) partWeightV0 <= upperBoundPartWeight + (double) maxNodeWeight) &&
           ((double) partWeightV1 >= lowerBoundPartWeight - (double) maxNodeWeight)) {
        const auto [_, distance, nodeId] = maxHeapPhase2.top();
        maxHeapPhase2.pop();

        // Move node to V0
        bisection[nodeId] = 0;
        currentEdgeCut += weightedOutDegree[nodeId];
        currentEdgeCut -= weightedInDegree[nodeId];
        partWeightV0 += workingGraph.nodeWeights[nodeId];
        partWeightV1 -= workingGraph.nodeWeights[nodeId];

        // Check neighbors for new movable vertices
        for (const auto &[neighborId, edgeWeight]: workingGraph.adj[nodeId]) {
            localInDegree[neighborId]--;
            assert(localInDegree[neighborId] >= 0 && "in-degree of node became negative");
            if (localInDegree[neighborId] == 0)
                maxHeapPhase2.emplace(weightedInDegree[neighborId] - weightedOutDegree[neighborId],
                                      workingGraph.nodeWeights[neighborId],
                                      neighborId);
        }

        // Update best solution if current is better
        if (((double) partWeightV0 < upperBoundPartWeight + (double) maxNodeWeight) &&
            ((double) partWeightV1 < upperBoundPartWeight + (double) maxNodeWeight)) {
            double maxImbalance = std::max((double) partWeightV0 / upperBoundPartWeight,
                                           (double) partWeightV1 / upperBoundPartWeight);

            if (currentEdgeCut < bestEdgeCut ||
                (currentEdgeCut == bestEdgeCut && std::abs(maxImbalance - 1.0) < bestImbalanceDiff)) {
                bestEdgeCut = currentEdgeCut;
                bestImbalanceDiff = maxImbalance;
                bestBisection = bisection;
            }
        }
    }

    return {bestBisection, bestEdgeCut};
}

std::pair<std::vector<uint8_t>, uint64_t> GreedyDirectedGraphGrowing::runOnReverseGraph() const {
    // Similar to runOnNormalGraph but:
    // 1. Starts with all vertices in V0
    // 2. Uses out-degree instead of in-degree
    // 3. Moves vertices to V1
    // Structure and logic mirrors runOnNormalGraph with reversed directions
    std::vector<uint8_t> bisection(workingGraph.size, 0);
    std::vector<uint8_t> bestBisection;
    uint64_t bestEdgeCut = UINT64_MAX;
    auto bestImbalanceDiff = DBL_MAX;
    uint64_t currentEdgeCut = 0;
    uint64_t partWeightV0 = workingGraph.totalWeight;
    uint64_t partWeightV1 = 0;
    std::priority_queue<std::tuple<uint64_t, uint64_t, uint64_t>,
            std::vector<std::tuple<uint64_t, uint64_t, uint64_t>>,
            HeapComparator> maxHeapPhase1;
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
        const auto &reverseNeighbors = workingGraph.revAdj[i];
        for (const auto &[reverseNeighborId, edgeWeight]: reverseNeighbors) weightedInDegree[i] += edgeWeight;
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, edgeWeight]: neighbors) weightedOutDegree[i] += edgeWeight;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(sinks.begin(), sinks.end(), gen);

    std::vector<uint64_t> distancesFromFirstMovedNode = workingGraph.distancesFromNode(sinks[0], true);
    for (const uint64_t &sink: sinks) {
        maxHeapPhase1.emplace(weightedOutDegree[sink], distancesFromFirstMovedNode[sink], sink);
    }

    while (!maxHeapPhase1.empty() && ((double) partWeightV1 <= 0.9 * (upperBoundPartWeight + (double) maxNodeWeight)) &&
           ((double) partWeightV0 >= 1.1 * (lowerBoundPartWeight - (double) maxNodeWeight))) {
        const auto [_, distance, nodeId] = maxHeapPhase1.top();
        maxHeapPhase1.pop();

        bisection[nodeId] = 1;
        currentEdgeCut -= weightedOutDegree[nodeId];
        currentEdgeCut += weightedInDegree[nodeId];
        partWeightV0 -= workingGraph.nodeWeights[nodeId];
        partWeightV1 += workingGraph.nodeWeights[nodeId];

        for (const auto &[neighborId, edgeWeight]: workingGraph.revAdj[nodeId]) {
            localOutDegree[neighborId]--;
            assert(localOutDegree[neighborId] >= 0 && "out-degree of node became negative");
            if (localOutDegree[neighborId] == 0)
                maxHeapPhase1.emplace(weightedOutDegree[neighborId], distancesFromFirstMovedNode[neighborId],
                                      neighborId);
        }

        if (((double) partWeightV0 < upperBoundPartWeight + (double) maxNodeWeight) &&
            ((double) partWeightV1 < upperBoundPartWeight + (double) maxNodeWeight)) {
            double maxImbalance = std::max((double) partWeightV0 / upperBoundPartWeight,
                                           (double) partWeightV1 / upperBoundPartWeight);

            if (currentEdgeCut < bestEdgeCut ||
                (currentEdgeCut == bestEdgeCut && std::abs(maxImbalance - 1.0) < bestImbalanceDiff)) {
                bestEdgeCut = currentEdgeCut;
                bestImbalanceDiff = maxImbalance;
                bestBisection = bisection;
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
    for (uint64_t nodeId: readyNodes) {
        maxHeapPhase2.emplace(weightedOutDegree[nodeId] - weightedInDegree[nodeId], workingGraph.nodeWeights[nodeId],
                              nodeId);
    }

    while (!maxHeapPhase2.empty() && ((double) partWeightV1 <= upperBoundPartWeight + (double) maxNodeWeight) &&
           ((double) partWeightV0 >= lowerBoundPartWeight - (double) maxNodeWeight)) {
        const auto [_, distance, nodeId] = maxHeapPhase2.top();
        maxHeapPhase2.pop();

        bisection[nodeId] = 1;
        currentEdgeCut -= weightedOutDegree[nodeId];
        currentEdgeCut += weightedInDegree[nodeId];
        partWeightV0 -= workingGraph.nodeWeights[nodeId];
        partWeightV1 += workingGraph.nodeWeights[nodeId];

        for (const auto &[neighborId, edgeWeight]: workingGraph.revAdj[nodeId]) {
            localOutDegree[neighborId]--;
            assert(localOutDegree[neighborId] >= 0 && "out-degree of node became negative");
            if (localOutDegree[neighborId] == 0)
                maxHeapPhase2.emplace(weightedOutDegree[neighborId] - weightedInDegree[neighborId],
                                      workingGraph.nodeWeights[neighborId],
                                      neighborId);
        }

        if (((double) partWeightV0 < upperBoundPartWeight + (double) maxNodeWeight) &&
            ((double) partWeightV1 < upperBoundPartWeight + (double) maxNodeWeight)) {
            double maxImbalance = std::max((double) partWeightV0 / upperBoundPartWeight,
                                           (double) partWeightV1 / upperBoundPartWeight);

            if (currentEdgeCut < bestEdgeCut ||
                (currentEdgeCut == bestEdgeCut && std::abs(maxImbalance - 1.0) < bestImbalanceDiff)) {
                bestEdgeCut = currentEdgeCut;
                bestImbalanceDiff = maxImbalance;
                bestBisection = bisection;
            }
        }
    }

    return {bestBisection, bestEdgeCut};
}

std::pair<std::vector<uint8_t>, uint64_t> GreedyDirectedGraphGrowing::run() const {
    // Run both normal and reverse algorithms
    std::pair<std::vector<uint8_t>, uint64_t> bisectionNormal = runOnNormalGraph();
    std::pair<std::vector<uint8_t>, uint64_t> bisectionReverse = runOnReverseGraph();

    // Verify both solutions maintain acyclicity
    assert(checkValidBisection(bisectionNormal.first) == true && "bisection on normal graph is invalid");
    assert(checkValidBisection(bisectionReverse.first) == true && "bisection on reverse graph is invalid");

    // Return the solution with smaller edge cut
    return (bisectionNormal.second < bisectionReverse.second ?
            std::move(bisectionNormal) : std::move(bisectionReverse));
}