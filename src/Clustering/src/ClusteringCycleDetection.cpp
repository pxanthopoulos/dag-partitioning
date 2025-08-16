/**
 * @file ClusteringCycleDetection.cpp
 * @brief Implementation of clustering with cycle detection algorithm
 */

#include "ClusteringCycleDetection.h"
#include <cassert>
#include <cmath>
#include <queue>
#include <unordered_map>
#include <unordered_set>

ClusteringCycleDetection::ClusteringCycleDetection(const Graph &graph,
                                                   uint64_t maxRounds,
                                                   uint64_t minVertices,
                                                   double vertexRatio)
    : Clustering(graph, maxRounds, minVertices, vertexRatio) {}

void ClusteringCycleDetection::hardCheckCycle(
    const std::vector<uint64_t> &leaders, uint64_t newSize) const {
    // Map original leader IDs to consecutive new IDs for the coarsened graph
    uint64_t maxNewNodeId = -1;
    std::unordered_map<uint64_t, uint64_t> leadersToNewNodeIds;
    std::vector<std::pair<std::vector<uint64_t>, uint64_t>> newNodes(newSize);
    std::unordered_set<uint64_t> seenLeaders;

    // First pass: Create mapping and collect nodes for each cluster
    for (uint64_t nodeId = 0; nodeId < workingGraph.size; ++nodeId) {
        uint64_t leader = leaders[nodeId];
        if (seenLeaders.find(leader) == seenLeaders.end()) {
            seenLeaders.insert(leader);
            maxNewNodeId++;
            leadersToNewNodeIds[leader] = maxNewNodeId;
        }

        newNodes[leadersToNewNodeIds[leader]].first.emplace_back(nodeId);
        newNodes[leadersToNewNodeIds[leader]].second +=
            workingGraph.nodeWeights[nodeId];
    }

    // Construct the coarsened graph
    Graph newGraph(newSize);
    std::vector<std::unordered_map<uint64_t, uint64_t>> newAdj(newSize);

    // Add nodes and combine edges
    for (uint64_t i = 0; i < newSize; ++i) {
        // Add node with combined cluster weight
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

    // Add combined edges to coarsened graph
    for (uint64_t i = 0; i < newSize; ++i) {
        auto &subAdj = newAdj[i];
        for (const auto &[to, weight] : subAdj) {
            newGraph.addEdge(i, to, weight);
        }
    }

    // Verify acyclicity
    assert(!newGraph.hasCycle() && "Hard check for cycle failed");
}

uint64_t ClusteringCycleDetection::findMinimumTopLevelInCluster(
    uint64_t node, const std::vector<uint64_t> &topLevels,
    const std::vector<uint8_t> &markup, const std::vector<uint8_t> &markdown) {
    if (markup[node] == 1)
        return topLevels[node];
    if (markdown[node] == 1)
        return topLevels[node] - 1;
    return topLevels[node];
}

bool ClusteringCycleDetection::detectCycle(
    uint64_t from, uint64_t to, const std::vector<uint64_t> &topLevels,
    const std::vector<uint64_t> &leaders, const std::vector<uint8_t> &markup,
    const std::vector<uint8_t> &markdown) const {
    // Get minimum top-level value in target cluster
    uint64_t minimumTopLevelInCluster =
        findMinimumTopLevelInCluster(to, topLevels, markup, markdown);
    std::vector<uint8_t> visited(workingGraph.size, 0);
    std::queue<uint64_t> q;

    // Start BFS from the node being added
    visited[from] = 1;
    q.push(from);

    // BFS traversal checking for paths back to target cluster
    while (!q.empty()) {
        uint64_t u = q.front();
        visited[u] = 1;

        // If we reach any node in the target cluster, we found a cycle
        if (leaders[u] == leaders[to])
            return true;

        q.pop();

        // Check successors - only follow edges where level difference <= 1
        for (const auto &[successorId, edgeWeight] : workingGraph.adj[u]) {
            // Skip direct edge being tested
            if (successorId == to && u == from)
                continue;

            // Only traverse to nodes with top-level close to cluster minimum
            uint64_t diff =
                (topLevels[successorId] > minimumTopLevelInCluster)
                    ? (topLevels[successorId] - minimumTopLevelInCluster)
                    : (minimumTopLevelInCluster - topLevels[successorId]);
            if (visited[successorId] == 0 && diff <= 1) {
                q.push(successorId);
            }
        }

        // Check predecessors - only follow edges within same cluster and level
        // difference <= 1
        for (const auto &[predecessorId, edgeWeight] : workingGraph.revAdj[u]) {
            uint64_t diff =
                (topLevels[predecessorId] > minimumTopLevelInCluster)
                    ? (topLevels[predecessorId] - minimumTopLevelInCluster)
                    : (minimumTopLevelInCluster - topLevels[predecessorId]);
            if (visited[predecessorId] == 0 && diff <= 1 &&
                leaders[predecessorId] == leaders[u]) {
                q.push(predecessorId);
            }
        }
    }
    return false;
}

std::pair<std::vector<uint64_t>, uint64_t>
ClusteringCycleDetection::oneRoundClustering() const {
    uint64_t newSize = workingGraph.size;
    std::vector<uint8_t> markup(workingGraph.size,
                                0); // Nodes at level t in their cluster
    std::vector<uint8_t> markdown(workingGraph.size,
                                  0); // Nodes at level t+1 in their cluster
    std::vector<uint64_t> leaders(workingGraph.size);
    std::vector<uint64_t> clusterWeights(workingGraph.size);

    // Initialize each node as its own cluster
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        leaders[i] = i;
        clusterWeights[i] = workingGraph.nodeWeights[i];
    }

    // Get topological order and top-level values
    const auto [topologicalOrder, topLevels] =
        workingGraph.topologicalSortAndTopLevels();

    // Process nodes in topological order
    for (uint64_t node : topologicalOrder) {
        // Skip nodes already in clusters
        if (markup[node] == 1 || markdown[node] == 1)
            continue;

        // Get neighbors sorted by edge weight
        std::vector<std::tuple<uint64_t, uint64_t, bool>> sortedNeighbors =
            workingGraph.getNeighborsSortedByEdgeWeightAsc(node);

        // Try to merge with each neighbor in order
        for (const auto &[neighborId, edgeWeight, isSuccessor] :
             sortedNeighbors) {
            // Check top-level difference constraint
            uint64_t diff = (topLevels[node] > topLevels[neighborId])
                                ? (topLevels[node] - topLevels[neighborId])
                                : (topLevels[neighborId] - topLevels[node]);
            if (diff > 1)
                continue;

            // Check cluster size constraint (10% of total weight)
            uint64_t leaderOfNeighbor = leaders[neighborId];
            uint64_t weightOfNeighborsCluster =
                clusterWeights[leaderOfNeighbor];
            if ((weightOfNeighborsCluster + workingGraph.nodeWeights[node]) >
                (uint64_t)ceil((double)workingGraph.totalWeight * 0.1))
                continue;

            if (isSuccessor) {
                // Edge from node to neighbor
                if (markup[neighborId] == 1 ||
                    detectCycle(node, neighborId, topLevels, leaders, markup,
                                markdown))
                    continue;
                leaders[node] = leaderOfNeighbor;
                markup[node] = markdown[neighborId] = 1;
                newSize--;
                clusterWeights[leaderOfNeighbor] +=
                    workingGraph.nodeWeights[node];
            } else {
                // Edge from neighbor to node
                if (markdown[neighborId] == 1 ||
                    detectCycle(neighborId, node, topLevels, leaders, markup,
                                markdown))
                    continue;
                leaders[node] = leaderOfNeighbor;
                markdown[node] = markup[neighborId] = 1;
                newSize--;
                clusterWeights[leaderOfNeighbor] +=
                    workingGraph.nodeWeights[node];
            }
            break;
        }

        // Stop if reached minimum vertices target
        if (newSize == minVertices)
            break;
    }

    return {leaders, newSize};
}