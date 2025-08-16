/**
 * @file Clustering.cpp
 * @brief Implementation of the Clustering base class
 */

#include "Clustering.h"

#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <utility>

// Constructor initializes the working graph and clustering parameters
Clustering::Clustering(Graph graph, uint64_t maxRounds, uint64_t minVertices,
                       double vertexRatio)
    : workingGraph(std::move(graph)), maxRounds(maxRounds),
      minVertices(minVertices), vertexRatio(vertexRatio) {}

bool Clustering::updateGraphAndClusters(const std::vector<uint64_t> &leaders,
                                        uint64_t newSize) {
    assert(leaders.size() == workingGraph.size &&
           "Leader value must be specified for all nodes");

    // If no clustering occurred, return false to stop
    if (newSize == workingGraph.size)
        return false;

    // Create mapping from leader nodes to new node IDs
    uint64_t maxNewNodeId = -1;
    std::unordered_map<uint64_t, uint64_t> leadersToNewNodeIds;
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> newNodes(newSize);
    std::unordered_set<uint64_t> seenLeaders;

    // First pass: Assign new IDs to leaders and group nodes
    for (uint64_t nodeId = 0; nodeId < workingGraph.size; ++nodeId) {
        uint64_t leader = leaders[nodeId];
        if (seenLeaders.find(leader) == seenLeaders.end()) {
            seenLeaders.insert(leader);
            maxNewNodeId++;
            leadersToNewNodeIds[leader] = maxNewNodeId;
        }

        // Add node to its cluster and accumulate weights
        newNodes[leadersToNewNodeIds[leader]].first.emplace_back(nodeId);
        newNodes[leadersToNewNodeIds[leader]].second +=
            workingGraph.nodeWeights[nodeId];
    }

    // Create new graph with merged nodes
    Graph newGraph(newSize);
    std::vector<std::unordered_map<uint64_t, uint64_t>> newAdj(newSize);

    // Second pass: Build new graph structure
    for (uint64_t i = 0; i < newSize; ++i) {
        // Add node with combined weight of its cluster
        newGraph.addNode(i, newNodes[i].second);

        // Process edges for all nodes in this cluster
        for (uint64_t containedNode : newNodes[i].first) {
            for (const auto &edge : workingGraph.adj[containedNode]) {
                // Skip self-loops within same cluster
                if (leadersToNewNodeIds[leaders[edge.first]] == i)
                    continue;

                // Combine parallel edges between clusters
                auto &subAdj = newAdj[i];
                auto it = subAdj.find(leadersToNewNodeIds[leaders[edge.first]]);
                if (it != subAdj.end()) {
                    it->second = it->second + edge.second;
                } else {
                    subAdj[leadersToNewNodeIds[leaders[edge.first]]] =
                        edge.second;
                }
            }
        }
    }

    // Add combined edges to new graph
    for (uint64_t i = 0; i < newSize; ++i) {
        auto &subAdj = newAdj[i];
        for (const auto &[to, weight] : subAdj) {
            newGraph.addEdge(i, to, weight);
        }
    }

    // Update working graph and save intermediate result
    workingGraph = newGraph;
    std::vector<uint64_t> clustering(leaders.size());
    for (uint64_t i = 0; i < clustering.size(); ++i)
        clustering[i] = leadersToNewNodeIds[leaders[i]];
    intermediateGraphsAndClusters.emplace(workingGraph, clustering);

    // Stop if minimum size reached or passed
    if (newSize <= minVertices)
        return false;
    if ((double)newSize / (double)leaders.size() > vertexRatio)
        return false;
    return true;
}

std::stack<std::pair<Graph, std::vector<uint64_t>>> Clustering::run() {
    // Handle case where graph is already small enough (return empty stack)
    if (workingGraph.size <= minVertices) {
        return intermediateGraphsAndClusters;
    }

    // Main clustering loop
    uint64_t countRounds = 0;
    while (true) {
        // Perform one round of clustering
        const auto &pair = oneRoundClustering();

        // Update graph and check if the min vertices stopping criteria is met
        if (!updateGraphAndClusters(pair.first, pair.second))
            break;

        // Check if maximum rounds reached
        if (++countRounds == maxRounds)
            break;
    }

    return intermediateGraphsAndClusters;
}