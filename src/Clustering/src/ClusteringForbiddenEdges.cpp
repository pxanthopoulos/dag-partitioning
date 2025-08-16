/**
 * @file ClusteringForbiddenEdges.cpp
 * @brief Implementation of the forbidden edges clustering algorithm
 */

#include "ClusteringForbiddenEdges.h"
#include <algorithm>
#include <cmath>
#include <cstddef>

ClusteringForbiddenEdges::ClusteringForbiddenEdges(const Graph &graph,
                                                   uint64_t maxRounds,
                                                   uint64_t minVertices,
                                                   double vertexRatio)
    : Clustering(graph, maxRounds, minVertices, vertexRatio) {}

std::vector<std::pair<uint64_t, uint64_t>>
ClusteringForbiddenEdges::findValidNeighbors(
    uint64_t node,
    const std::vector<std::tuple<uint64_t, uint64_t, bool>> &neighbors,
    const std::vector<uint64_t> &leaders,
    const std::vector<uint64_t> &clusterWeights,
    const std::vector<uint64_t> &topLevels,
    const std::vector<uint64_t> &numberOfBadNeighbors,
    const std::vector<uint64_t> &leaderOfBadNeighbors,
    const std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>>
        &leadersToMinMaxTopValues) const {
    std::vector<std::pair<uint64_t, uint64_t>> validNeighbors;

    // If node has two bad neighbors, it must remain singleton (Theorem 4.2)
    if (numberOfBadNeighbors[node] == 2)
        return {};

    // Special case: node has exactly one bad neighbor cluster
    if (numberOfBadNeighbors[node] == 1) {
        uint64_t leaderOfBadNeighbor = leaderOfBadNeighbors[node];

        // Check cluster size constraint (10% of total graph weight)
        if ((clusterWeights[leaderOfBadNeighbor] +
             workingGraph.nodeWeights[node]) >
            (uint64_t)ceil((double)workingGraph.totalWeight * 0.1))
            return {};

        // Check top-level compatibility with the cluster
        auto &value = leadersToMinMaxTopValues.at(leaderOfBadNeighbor);
        if (topLevels[node] != value.first && topLevels[node] != value.second)
            return {};

        // Can only join the bad neighbor cluster
        for (const auto &[neighborId, edgeWeight, isSuccessor] : neighbors) {
            if (leaders[neighborId] == leaderOfBadNeighbor)
                return {{neighborId, edgeWeight}};
        }
    }

    // Check all neighbors for valid clustering candidates
    for (const auto &[neighborId, edgeWeight, isSuccessor] : neighbors) {
        uint64_t leaderOfNeighbor = leaders[neighborId];
        uint64_t weightOfNeighborsCluster = clusterWeights[leaderOfNeighbor];

        // Check cluster size constraint (10% of total graph weight)
        if ((weightOfNeighborsCluster + workingGraph.nodeWeights[node]) >
            (uint64_t)ceil((double)workingGraph.totalWeight * 0.1))
            continue;

        // Check top-level difference constraint (Theorem 4.2)
        uint64_t diff = (topLevels[node] > topLevels[neighborId])
                            ? (topLevels[node] - topLevels[neighborId])
                            : (topLevels[neighborId] - topLevels[node]);
        if (diff > 1)
            continue;

        // Skip neighbors that have bad neighbors themselves
        if (numberOfBadNeighbors[neighborId] != 0)
            continue;

        validNeighbors.emplace_back(neighborId, edgeWeight);
    }

    return validNeighbors;
}

uint64_t ClusteringForbiddenEdges::findBestNeighbor(
    const std::vector<std::pair<uint64_t, uint64_t>> &validNeighbors) {
    uint64_t bestNeighborId;
    uint64_t maxEdgeWeight = UINT64_MAX;

    // Select neighbor with minimum edge weight
    for (const auto &[validNeighborId, validNeighborEdgeWeight] :
         validNeighbors) {
        if (validNeighborEdgeWeight < maxEdgeWeight) {
            bestNeighborId = validNeighborId;
            maxEdgeWeight = validNeighborEdgeWeight;
        }
    }
    return bestNeighborId;
}

std::pair<std::vector<uint64_t>, uint64_t>
ClusteringForbiddenEdges::oneRoundClustering() const {
    uint64_t newSize = workingGraph.size;
    std::vector<uint8_t> isNotSingleton(workingGraph.size, 0);
    std::vector<uint64_t> leaders(workingGraph.size);
    std::vector<uint64_t> clusterWeights(workingGraph.size);

    // Initialize each node as its own cluster
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        leaders[i] = i;
        clusterWeights[i] = workingGraph.nodeWeights[i];
    }

    // Data structures for tracking bad neighbors (Theorem 4.2)
    std::vector<uint64_t> numberOfBadNeighbors(workingGraph.size, 0);
    std::vector<uint64_t> leaderOfBadNeighbors(workingGraph.size, UINT64_MAX);
    std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>>
        leadersToMinMaxTopValues;

    // Get topological order and top-level values
    const auto [topologicalOrder, topLevels] =
        workingGraph.topologicalSortAndTopLevels();

    // Process nodes in topological order
    for (uint64_t node : topologicalOrder) {
        if (isNotSingleton[node] == 1)
            continue;

        // Find valid clustering candidates
        std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors =
            workingGraph.getNeighbors(node);
        std::vector<std::pair<uint64_t, uint64_t>> validNeighbors =
            findValidNeighbors(node, neighbors, leaders, clusterWeights,
                               topLevels, numberOfBadNeighbors,
                               leaderOfBadNeighbors, leadersToMinMaxTopValues);

        if (validNeighbors.empty())
            continue;

        // Merge node into best neighbor's cluster
        newSize--;
        uint64_t bestNeighbor = findBestNeighbor(validNeighbors);
        uint64_t leaderOfBestNeighbor = leaders[bestNeighbor];
        leaders[node] = leaderOfBestNeighbor;
        clusterWeights[leaderOfBestNeighbor] += workingGraph.nodeWeights[node];

        // Update min/max top-level values for the cluster
        auto [it, inserted] = leadersToMinMaxTopValues.try_emplace(
            leaderOfBestNeighbor,
            std::min(topLevels[node], topLevels[bestNeighbor]),
            std::max(topLevels[node], topLevels[bestNeighbor]));

        if (!inserted) {
            it->second = {std::min(topLevels[node], it->second.first),
                          std::max(topLevels[node], it->second.second)};
        }

        // Update bad neighbor information for node's neighbors (only those with
        // top level diff <= 1)
        for (const auto &[neighborId, edgeWeight, isSuccessor] : neighbors) {
            uint64_t diff = (topLevels[node] > topLevels[neighborId])
                                ? (topLevels[node] - topLevels[neighborId])
                                : (topLevels[neighborId] - topLevels[node]);
            if (diff > 1)
                continue;

            if (numberOfBadNeighbors[neighborId] == 0) {
                numberOfBadNeighbors[neighborId] = 1;
                leaderOfBadNeighbors[neighborId] = leaderOfBestNeighbor;
            } else if ((numberOfBadNeighbors[neighborId] == 1) &&
                       (leaderOfBadNeighbors[neighborId] !=
                        leaderOfBestNeighbor)) {
                numberOfBadNeighbors[neighborId] =
                    2; // Node now has 2 bad neighbors
            }
        }

        // If best neighbor was singleton, update its neighbors too (only those
        // with top level diff <= 1)
        if (isNotSingleton[bestNeighbor] == 0) {
            std::vector<std::tuple<uint64_t, uint64_t, bool>>
                neighborsOfBestNeighbor =
                    workingGraph.getNeighbors(bestNeighbor);
            for (const auto &[neighborId, edgeWeight, isSuccessor] :
                 neighborsOfBestNeighbor) {
                uint64_t diff =
                    (topLevels[bestNeighbor] > topLevels[neighborId])
                        ? (topLevels[bestNeighbor] - topLevels[neighborId])
                        : (topLevels[neighborId] - topLevels[bestNeighbor]);
                if (diff > 1)
                    continue;

                if (numberOfBadNeighbors[neighborId] == 0) {
                    numberOfBadNeighbors[neighborId] = 1;
                    leaderOfBadNeighbors[neighborId] = leaderOfBestNeighbor;
                } else if ((numberOfBadNeighbors[neighborId] == 1) &&
                           (leaderOfBadNeighbors[neighborId] !=
                            leaderOfBestNeighbor)) {
                    numberOfBadNeighbors[neighborId] = 2;
                }
            }

            isNotSingleton[bestNeighbor] = 1;
        }
        isNotSingleton[node] = 1;

        // Stop if reached minimum vertices target
        if (newSize == minVertices)
            break;
    }

    return {leaders, newSize};
}