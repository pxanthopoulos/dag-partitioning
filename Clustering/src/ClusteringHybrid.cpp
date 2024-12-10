//
// Created by panagiotis on 10/12/2024.
//

#include "ClusteringHybrid.h"
#include <cmath>

bool ClusteringHybrid::checkLargeDegrees(uint64_t from, uint64_t to) const {
    auto maxDegree = (uint64_t) (sqrt((double) workingGraph.size) / 10.0);
    if (workingGraph.adj[from].size() > maxDegree) return true;
    if (workingGraph.inDegree[to] > maxDegree) return true;
    return false;
}

void ClusteringHybrid::bookKeepingForForbiddenEdges(const std::vector<uint64_t> &topLevels, uint64_t node,
                                                    const std::vector<std::tuple<uint64_t, uint64_t, bool>> &sortedNeighbors,
                                                    uint64_t neighborId, uint64_t leaderOfNeighbor,
                                                    const std::vector<bool> &markup,
                                                    const std::vector<bool> &markdown,
                                                    std::vector<uint64_t> &numberOfBadNeighbors,
                                                    std::vector<uint64_t> &leaderOfBadNeighbors) const {
    for (const auto &[neighborIdInternal, edgeWeightInternal, isSuccessorInternal]: sortedNeighbors) {
        uint64_t diffInternal = (topLevels[node] > topLevels[neighborIdInternal]) ? (topLevels[node] -
                                                                                     topLevels[neighborIdInternal])
                                                                                  : (
                                        topLevels[neighborIdInternal] - topLevels[node]);
        if (diffInternal > 1) continue;
        if (numberOfBadNeighbors[neighborIdInternal] == 0) {
            numberOfBadNeighbors[neighborIdInternal] = 1;
            leaderOfBadNeighbors[neighborIdInternal] = leaderOfNeighbor;
        } else if ((numberOfBadNeighbors[neighborIdInternal] == 1) &&
                   (leaderOfBadNeighbors[neighborIdInternal] != leaderOfNeighbor)) {
            numberOfBadNeighbors[neighborIdInternal] = 2;
        }
    }

    if (!markup[neighborId] && !markdown[neighborId]) {
        std::vector<std::tuple<uint64_t, uint64_t, bool>> neighborsOfMergedNeighbor = workingGraph.getNeighbors(
                neighborId);
        for (const auto &[neighborIdInternal, edgeWeightInternal, isSuccessorInternal]: neighborsOfMergedNeighbor) {
            uint64_t diffInternal = (topLevels[neighborId] > topLevels[neighborIdInternal]) ? (
                    topLevels[neighborId] - topLevels[neighborIdInternal]) : (
                                            topLevels[neighborIdInternal] - topLevels[neighborId]);
            if (diffInternal > 1) continue;
            if (numberOfBadNeighbors[neighborIdInternal] == 0) {
                numberOfBadNeighbors[neighborIdInternal] = 1;
                leaderOfBadNeighbors[neighborIdInternal] = leaderOfNeighbor;
            } else if ((numberOfBadNeighbors[neighborIdInternal] == 1) &&
                       (leaderOfBadNeighbors[neighborIdInternal] != leaderOfNeighbor)) {
                numberOfBadNeighbors[neighborIdInternal] = 2;
            }
        }
    }
}

std::pair<std::vector<uint64_t>, uint64_t> ClusteringHybrid::oneRoundClustering() const {
    uint64_t newSize = workingGraph.size;
    std::vector<uint64_t> leaders(workingGraph.size);
    std::vector<uint64_t> clusterWeights(workingGraph.size);
    for (size_t i = 0; i < workingGraph.size; ++i) {
        leaders[i] = i;
        clusterWeights[i] = workingGraph.nodeWeights[i];
    }
    std::vector<bool> markup(workingGraph.size, false);
    std::vector<bool> markdown(workingGraph.size, false);
    std::vector<uint64_t> numberOfBadNeighbors(workingGraph.size, 0);
    std::vector<uint64_t> leaderOfBadNeighbors(workingGraph.size, UINT64_MAX);

    const auto [topologicalOrder, topLevels] = workingGraph.topologicalSortAndTopLevels();

    for (uint64_t node: topologicalOrder) {
        if (markup[node] || markdown[node]) continue;

        std::vector<std::tuple<uint64_t, uint64_t, bool>> sortedNeighbors = workingGraph.getNeighborsSortedByEdgeWeightAsc(
                node);
        for (const auto &[neighborId, edgeWeight, isSuccessor]: sortedNeighbors) {
            uint64_t diff = (topLevels[node] > topLevels[neighborId]) ? (topLevels[node] - topLevels[neighborId]) : (
                    topLevels[neighborId] - topLevels[node]);
            if (diff > 1) continue;

            uint64_t leaderOfNeighbor = leaders[neighborId];
            uint64_t weightOfNeighborsCluster = clusterWeights[leaderOfNeighbor];
            if ((weightOfNeighborsCluster + workingGraph.nodeWeights[node]) >
                (uint64_t) ceil((double) workingGraph.totalWeight * 0.1))
                continue;
            if (isSuccessor) {
                if (markup[neighborId]) continue;

                bool largeDegrees = checkLargeDegrees(node, neighborId);

                if (!largeDegrees && detectCycle(node, neighborId, topLevels, leaders))
                    continue;

                if (numberOfBadNeighbors[node] == 2) continue;
                if (numberOfBadNeighbors[node] == 1) {
                    uint64_t leaderOfBadNeighbor = leaderOfBadNeighbors[node];
                    if (leaderOfNeighbor != leaderOfBadNeighbor) continue;
                }
                leaders[node] = leaderOfNeighbor;
                newSize--;
                clusterWeights[leaderOfNeighbor] += workingGraph.nodeWeights[node];
                bookKeepingForForbiddenEdges(topLevels, node, sortedNeighbors, neighborId, leaderOfNeighbor, markup,
                                             markdown, numberOfBadNeighbors, leaderOfBadNeighbors);
                markup[node] = markdown[neighborId] = true;
            } else {
                if (markdown[neighborId]) continue;

                bool largeDegrees = checkLargeDegrees(neighborId, node);

                if (!largeDegrees && detectCycle(neighborId, node, topLevels, leaders))
                    continue;

                if (numberOfBadNeighbors[node] == 2) continue;
                if (numberOfBadNeighbors[node] == 1) {
                    uint64_t leaderOfBadNeighbor = leaderOfBadNeighbors[node];
                    if (leaderOfNeighbor != leaderOfBadNeighbor) continue;
                }
                leaders[node] = leaderOfNeighbor;
                newSize--;
                clusterWeights[leaderOfNeighbor] += workingGraph.nodeWeights[node];

                bookKeepingForForbiddenEdges(topLevels, node, sortedNeighbors, neighborId, leaderOfNeighbor, markup,
                                             markdown, numberOfBadNeighbors, leaderOfBadNeighbors);

                markdown[node] = markup[neighborId] = true;
            }
            break;
        }
    }

    return {leaders, newSize};
}
