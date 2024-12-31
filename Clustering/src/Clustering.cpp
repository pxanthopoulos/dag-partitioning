//
// Created by panagiotis on 8/12/2024.
//

#include "Clustering.h"
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <numeric>

Clustering::Clustering(Graph graph, uint64_t maxRounds, uint64_t minVertices) :
        workingGraph(std::move(graph)), maxRounds(maxRounds), minVertices(minVertices) {}

bool Clustering::updateGraphAndClusters(const std::vector<uint64_t> &leaders, uint64_t newSize) {
    assert(leaders.size() == workingGraph.size && "Leader value must be specified for all nodes");
    if (newSize == workingGraph.size) return false;

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

    workingGraph = newGraph;
    std::vector<uint64_t> clustering(leaders.size());
    for (uint64_t i = 0; i < clustering.size(); ++i) clustering[i] = leadersToNewNodeIds[leaders[i]];
    intermediateGraphsAndClusters.emplace(workingGraph, clustering);

    if (newSize == minVertices) return false;
    return true;
}

std::stack<std::pair<Graph, std::vector<uint64_t>>> Clustering::run() {
    if (workingGraph.size <= minVertices) {
        std::vector<uint64_t> clustering(workingGraph.size);
        iota(clustering.begin(), clustering.end(), 0);
        intermediateGraphsAndClusters.emplace(workingGraph, clustering);
        return intermediateGraphsAndClusters;
    }
    uint64_t countRounds = 0;
    while (true) {
        const auto &pair = oneRoundClustering();
        if (!updateGraphAndClusters(pair.first, pair.second)) break;
        if (++countRounds == maxRounds) break;
    }
    return intermediateGraphsAndClusters;
}
