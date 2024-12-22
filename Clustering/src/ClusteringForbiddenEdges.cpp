//
// Created by panagiotis on 6/12/2024.
//

#include "ClusteringForbiddenEdges.h"
#include <cstddef>
#include <cmath>
#include <algorithm>

ClusteringForbiddenEdges::ClusteringForbiddenEdges(const Graph &graph, uint64_t maxRounds, uint64_t minVertices)
        : Clustering(graph, maxRounds, minVertices) {}

std::vector<std::pair<uint64_t, uint64_t>>
ClusteringForbiddenEdges::findValidNeighbors(uint64_t node,
                                             const std::vector<std::tuple<uint64_t, uint64_t, bool>> &neighbors,
                                             const std::vector<uint64_t> &leaders,
                                             const std::vector<uint64_t> &clusterWeights,
                                             const std::vector<uint64_t> &topLevels,
                                             const std::vector<uint64_t> &numberOfBadNeighbors,
                                             const std::vector<uint64_t> &leaderOfBadNeighbors,
                                             const std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>> &leadersToMinMaxTopValues) const {
    std::vector<std::pair<uint64_t, uint64_t>> validNeighbors;

    if (numberOfBadNeighbors[node] == 2) return {};

    if (numberOfBadNeighbors[node] == 1) {
        uint64_t leaderOfBadNeighbor = leaderOfBadNeighbors[node];
        if ((clusterWeights[leaderOfBadNeighbor] + workingGraph.nodeWeights[node]) >
            (uint64_t) ceil((double) workingGraph.totalWeight * 0.1))
            return {};
        auto &value = leadersToMinMaxTopValues.at(leaderOfBadNeighbor);
        if (topLevels[node] != value.first && topLevels[node] != value.second) return {};
        for (const auto &[neighborId, edgeWeight, isSuccessor]: neighbors) {
            if (leaders[neighborId] == leaderOfBadNeighbor) return {{neighborId, edgeWeight}};
        }
    }

    for (const auto &[neighborId, edgeWeight, isSuccessor]: neighbors) {
        uint64_t leaderOfNeighbor = leaders[neighborId];
        uint64_t weightOfNeighborsCluster = clusterWeights[leaderOfNeighbor];
        if ((weightOfNeighborsCluster + workingGraph.nodeWeights[node]) >
            (uint64_t) ceil((double) workingGraph.totalWeight * 0.1))
            continue;

        uint64_t diff = (topLevels[node] > topLevels[neighborId]) ? (topLevels[node] - topLevels[neighborId]) : (
                topLevels[neighborId] - topLevels[node]);
        if (diff > 1) continue;

        validNeighbors.emplace_back(neighborId, edgeWeight);
    }

    return validNeighbors;
}

uint64_t
ClusteringForbiddenEdges::findBestNeighbor(const std::vector<std::pair<uint64_t, uint64_t>> &validNeighbors) {
    uint64_t bestNeighborId;
    uint64_t maxEdgeWeight = UINT64_MAX;
    for (const auto &[validNeighborId, validNeighborEdgeWeight]: validNeighbors) {
        if (validNeighborEdgeWeight < maxEdgeWeight) {
            bestNeighborId = validNeighborId;
            maxEdgeWeight = validNeighborEdgeWeight;
        }
    }
    return bestNeighborId;
}

std::pair<std::vector<uint64_t>, uint64_t> ClusteringForbiddenEdges::oneRoundClustering() const {
    uint64_t newSize = workingGraph.size;
    std::vector<bool> isNotSingleton(workingGraph.size, false);
    std::vector<uint64_t> leaders(workingGraph.size);
    std::vector<uint64_t> clusterWeights(workingGraph.size);
    for (size_t i = 0; i < workingGraph.size; ++i) {
        leaders[i] = i;
        clusterWeights[i] = workingGraph.nodeWeights[i];
    }
    std::vector<uint64_t> numberOfBadNeighbors(workingGraph.size, 0);
    std::vector<uint64_t> leaderOfBadNeighbors(workingGraph.size, UINT64_MAX);
    std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>> leadersToMinMaxTopValues;

    const auto [topologicalOrder, topLevels] = workingGraph.topologicalSortAndTopLevels();

    for (uint64_t node: topologicalOrder) {
        if (isNotSingleton[node]) continue;

        std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors = workingGraph.getNeighbors(node);
        std::vector<std::pair<uint64_t, uint64_t>> validNeighbors = findValidNeighbors(node, neighbors,
                                                                                       leaders,
                                                                                       clusterWeights,
                                                                                       topLevels,
                                                                                       numberOfBadNeighbors,
                                                                                       leaderOfBadNeighbors,
                                                                                       leadersToMinMaxTopValues);

        if (validNeighbors.empty()) continue;

        newSize--;
        uint64_t bestNeighbor = findBestNeighbor(validNeighbors);
        uint64_t leaderOfBestNeighbor = leaders[bestNeighbor];
        leaders[node] = leaderOfBestNeighbor;
        clusterWeights[leaderOfBestNeighbor] += workingGraph.nodeWeights[node];

        auto [it, inserted] = leadersToMinMaxTopValues.try_emplace(
                leaderOfBestNeighbor,
                std::min(topLevels[node], topLevels[bestNeighbor]),
                std::max(topLevels[node], topLevels[bestNeighbor])
        );

        if (!inserted) {
            it->second = {
                    std::min(topLevels[node], it->second.first),
                    std::max(topLevels[node], it->second.second)
            };
        }

        for (const auto &[neighborId, edgeWeight, isSuccessor]: neighbors) {
            uint64_t diff = (topLevels[node] > topLevels[neighborId]) ? (topLevels[node] - topLevels[neighborId])
                                                                      : (
                                    topLevels[neighborId] - topLevels[node]);
            if (diff > 1) continue;
            if (numberOfBadNeighbors[neighborId] == 0) {
                numberOfBadNeighbors[neighborId] = 1;
                leaderOfBadNeighbors[neighborId] = leaderOfBestNeighbor;
            } else if ((numberOfBadNeighbors[neighborId] == 1) &&
                       (leaderOfBadNeighbors[neighborId] != leaderOfBestNeighbor)) {
                numberOfBadNeighbors[neighborId] = 2;
            }
        }

        if (!isNotSingleton[bestNeighbor]) {
            std::vector<std::tuple<uint64_t, uint64_t, bool>> neighborsOfBestNeighbor = workingGraph.getNeighbors(
                    bestNeighbor);
            for (const auto &[neighborId, edgeWeight, isSuccessor]: neighborsOfBestNeighbor) {
                uint64_t diff = (topLevels[bestNeighbor] > topLevels[neighborId]) ? (topLevels[bestNeighbor] -
                                                                                     topLevels[neighborId])
                                                                                  : (
                                        topLevels[neighborId] - topLevels[bestNeighbor]);
                if (diff > 1) continue;
                if (numberOfBadNeighbors[neighborId] == 0) {
                    numberOfBadNeighbors[neighborId] = 1;
                    leaderOfBadNeighbors[neighborId] = leaderOfBestNeighbor;
                } else if ((numberOfBadNeighbors[neighborId] == 1) &&
                           (leaderOfBadNeighbors[neighborId] != leaderOfBestNeighbor)) {
                    numberOfBadNeighbors[neighborId] = 2;
                }
            }

            isNotSingleton[bestNeighbor] = true;
        }
        isNotSingleton[node] = true;
        if (newSize == minVertices) break;
    }

    return {leaders, newSize};
}
