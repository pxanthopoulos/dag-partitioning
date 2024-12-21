//
// Created by panagiotis on 8/12/2024.
//

#include "Clustering.h"
#include <map>
#include <unordered_set>

Clustering::Clustering(const Graph &graph) :
        workingGraph(graph), clusters(graph.size) {
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        clusters[i] = i;
    }
}

bool Clustering::updateGraphAndClusters(const std::vector<uint64_t> &leaders, uint64_t newSize) {
    assert(leaders.size() == workingGraph.size && "Leader value must be specified for all nodes");
    if (newSize == workingGraph.size) return false;

    uint64_t maxNewNodeId = -1;
    std::map<uint64_t, uint64_t> leadersToNewNodeIds;
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
    for (unsigned long &cluster: clusters) cluster = leadersToNewNodeIds[leaders[cluster]];

    return true;
}

void Clustering::run() {
    while (true) {
        const auto &pair = oneRoundClustering();
        if (!updateGraphAndClusters(pair.first, pair.second)) break;
    }
}
