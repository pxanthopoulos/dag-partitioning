//
// Created by panagiotis on 8/12/2024.
//

#include "ClusteringCycleDetection.h"
#include <queue>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include "llvm/ADT/STLExtras.h"

ClusteringCycleDetection::ClusteringCycleDetection(const Graph &graph, uint64_t maxRounds, uint64_t minVertices)
        : Clustering(graph, maxRounds, minVertices) {}

void ClusteringCycleDetection::hardCheckCycle(const std::vector<uint64_t> &leaders, uint64_t newSize) const {
    uint64_t maxNewNodeId = -1;
    std::unordered_map<uint64_t, uint64_t> leadersToNewNodeIds;
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> newNodes(newSize);
    std::unordered_set<uint64_t> seenLeaders;
    for (uint64_t nodeId = 0; nodeId < workingGraph.size; ++nodeId) {
        uint64_t leader = leaders[nodeId];
        if (seenLeaders.find(leader) == seenLeaders.end()) {
            seenLeaders.insert(leader);
            maxNewNodeId++;
            leadersToNewNodeIds[leader] = maxNewNodeId;
        }

        newNodes[leadersToNewNodeIds[leader]].first.emplace_back(nodeId);
        newNodes[leadersToNewNodeIds[leader]].second += workingGraph.nodeWeights[nodeId];
    }

    Graph newGraph(newSize);
    std::vector<std::unordered_map<uint64_t, uint64_t>> newAdj(newSize);
    for (uint64_t i = 0; i < newSize; ++i) {
        newGraph.addNode(i, newNodes[i].second);
        for (uint64_t containedNode: newNodes[i].first) {
            for (const auto &edge: workingGraph.adj[containedNode]) {
                if (leadersToNewNodeIds[leaders[edge.first]] == i) continue;
                auto &subAdj = newAdj[i];
                auto it = subAdj.find(leadersToNewNodeIds[leaders[edge.first]]);
                if (it != subAdj.end()) {
                    it->second = it->second + edge.second;
                } else {
                    subAdj[leadersToNewNodeIds[leaders[edge.first]]] = edge.second;
                }
            }
        }
    }


    for (uint64_t i = 0; i < newSize; ++i) {
        auto &subAdj = newAdj[i];
        for (const auto &[to, weight]: subAdj) {
            newGraph.addEdge(i, to, weight);
        }
    }
    assert(!newGraph.hasCycle() && "Hard check for cycle failed");
}

uint64_t
ClusteringCycleDetection::findMinimumTopLevelInCluster(uint64_t node, const std::vector<uint64_t> &topLevels,
                                                       const std::vector<uint64_t> &leaders) {
    uint64_t nodeLeader = leaders[node];
    uint64_t minimumTopLevelInCluster = UINT64_MAX;
    for (auto [i, leader]: llvm::enumerate(leaders)) {
        if (leader == nodeLeader) {
            if (topLevels[i] < minimumTopLevelInCluster) minimumTopLevelInCluster = topLevels[i];
        }
    }
    return minimumTopLevelInCluster;
}

bool
ClusteringCycleDetection::detectCycle(uint64_t from, uint64_t to, const std::vector<uint64_t> &topLevels,
                                      const std::vector<uint64_t> &leaders) const {
    uint64_t minimumTopLevelInCluster = findMinimumTopLevelInCluster(to, topLevels, leaders);
    std::vector<bool> visited(workingGraph.size, false);
    std::queue<uint64_t> q;

    visited[from] = true;
    q.push(from);

    while (!q.empty()) {
        uint64_t u = q.front();
        visited[u] = true;
        if (leaders[u] == leaders[to]) return true;

        q.pop();

        for (const auto &[successorId, edgeWeight]: workingGraph.adj[u]) {
            if (successorId == to && u == from) continue;
            uint64_t diff = (topLevels[successorId] > minimumTopLevelInCluster) ? (topLevels[successorId] -
                                                                                   minimumTopLevelInCluster) : (
                                    minimumTopLevelInCluster - topLevels[successorId]);
            if (!visited[successorId] && diff <= 1) {
                q.push(successorId);
            }
        }
        for (const auto &[predecessorId, edgeWeight]: workingGraph.revAdj[u]) {
            uint64_t diff = (topLevels[predecessorId] > minimumTopLevelInCluster) ? (topLevels[predecessorId] -
                                                                                     minimumTopLevelInCluster) : (
                                    minimumTopLevelInCluster - topLevels[predecessorId]);
            if (!visited[predecessorId] && diff <= 1 && leaders[predecessorId] == leaders[u]) {
                q.push(predecessorId);
            }
        }

    }
    return false;
}

std::pair<std::vector<uint64_t>, uint64_t> ClusteringCycleDetection::oneRoundClustering() const {
    uint64_t newSize = workingGraph.size;
    std::vector<bool> markup(workingGraph.size, false);
    std::vector<bool> markdown(workingGraph.size, false);
    std::vector<uint64_t> leaders(workingGraph.size);
    std::vector<uint64_t> clusterWeights(workingGraph.size);
    for (size_t i = 0; i < workingGraph.size; ++i) {
        leaders[i] = i;
        clusterWeights[i] = workingGraph.nodeWeights[i];
    }

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
                if (markup[neighborId] ||
                    detectCycle(node, neighborId, topLevels, leaders))
                    continue;
                leaders[node] = leaderOfNeighbor;
                markup[node] = markdown[neighborId] = true;
                newSize--;
                clusterWeights[leaderOfNeighbor] += workingGraph.nodeWeights[node];
            } else {
                if (markdown[neighborId] ||
                    detectCycle(neighborId, node, topLevels, leaders))
                    continue;
                leaders[node] = leaderOfNeighbor;
                markdown[node] = markup[neighborId] = true;
                newSize--;
                clusterWeights[leaderOfNeighbor] += workingGraph.nodeWeights[node];
            }
            break;
        }
        if (newSize == minVertices) break;
    }

    return {leaders, newSize};
}
