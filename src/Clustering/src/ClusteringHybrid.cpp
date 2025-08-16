/**
 * @file ClusteringHybrid.cpp
 * @brief Implementation of hybrid clustering combining cycle detection and
 * forbidden edges
 */

#include "ClusteringHybrid.h"
#include <cmath>

ClusteringHybrid::ClusteringHybrid(const Graph &graph, uint64_t maxRounds,
                                   uint64_t minVertices, double vertexRatio)
    : ClusteringCycleDetection(graph, maxRounds, minVertices, vertexRatio) {}

bool ClusteringHybrid::checkLargeDegrees(uint64_t from, uint64_t to) const {
    // Calculate degree threshold as sqrt(|V|)/10 as specified in the paper
    auto maxDegree = (uint64_t)(sqrt((double)workingGraph.size) / 10.0);

    // Check if either endpoint exceeds threshold
    if (workingGraph.adj[from].size() > maxDegree)
        return true; // Large out-degree
    if (workingGraph.inDegree[to] > maxDegree)
        return true; // Large in-degree
    return false;
}

void ClusteringHybrid::bookKeepingForForbiddenEdges(
    const std::vector<uint64_t> &topLevels, uint64_t node,
    const std::vector<std::tuple<uint64_t, uint64_t, bool>> &sortedNeighbors,
    uint64_t neighborId, uint64_t leaderOfNeighbor,
    const std::vector<uint8_t> &markup, const std::vector<uint8_t> &markdown,
    std::vector<uint64_t> &numberOfBadNeighbors,
    std::vector<uint64_t> &leaderOfBadNeighbors) const {
    // Update bad neighbor information for all neighbors of the node being
    // merged (only those with top level diff <= 1)
    for (const auto &[neighborIdInternal, edgeWeightInternal,
                      isSuccessorInternal] : sortedNeighbors) {
        // Only consider neighbors with top-level difference <= 1
        uint64_t diffInternal =
            (topLevels[node] > topLevels[neighborIdInternal])
                ? (topLevels[node] - topLevels[neighborIdInternal])
                : (topLevels[neighborIdInternal] - topLevels[node]);
        if (diffInternal > 1)
            continue;

        // Update bad neighbor tracking
        if (numberOfBadNeighbors[neighborIdInternal] == 0) {
            // First bad neighbor cluster
            numberOfBadNeighbors[neighborIdInternal] = 1;
            leaderOfBadNeighbors[neighborIdInternal] = leaderOfNeighbor;
        } else if ((numberOfBadNeighbors[neighborIdInternal] == 1) &&
                   (leaderOfBadNeighbors[neighborIdInternal] !=
                    leaderOfNeighbor)) {
            // Second bad neighbor cluster - node becomes unmatchable
            numberOfBadNeighbors[neighborIdInternal] = 2;
        }
    }

    // If neighbor was a singleton, update its neighbors too (only those with
    // top level diff <= 1)
    if (markup[neighborId] == 0 && markdown[neighborId] == 0) {
        std::vector<std::tuple<uint64_t, uint64_t, bool>>
            neighborsOfMergedNeighbor = workingGraph.getNeighbors(neighborId);

        for (const auto &[neighborIdInternal, edgeWeightInternal,
                          isSuccessorInternal] : neighborsOfMergedNeighbor) {
            uint64_t diffInternal =
                (topLevels[neighborId] > topLevels[neighborIdInternal])
                    ? (topLevels[neighborId] - topLevels[neighborIdInternal])
                    : (topLevels[neighborIdInternal] - topLevels[neighborId]);
            if (diffInternal > 1)
                continue;

            // Update bad neighbor tracking for neighbor's neighbors
            if (numberOfBadNeighbors[neighborIdInternal] == 0) {
                numberOfBadNeighbors[neighborIdInternal] = 1;
                leaderOfBadNeighbors[neighborIdInternal] = leaderOfNeighbor;
            } else if ((numberOfBadNeighbors[neighborIdInternal] == 1) &&
                       (leaderOfBadNeighbors[neighborIdInternal] !=
                        leaderOfNeighbor)) {
                numberOfBadNeighbors[neighborIdInternal] = 2;
            }
        }
    }
}

std::pair<std::vector<uint64_t>, uint64_t>
ClusteringHybrid::oneRoundClustering() const {
    uint64_t newSize = workingGraph.size;
    std::vector<uint64_t> leaders(workingGraph.size);
    std::vector<uint64_t> clusterWeights(workingGraph.size);

    // Initialize each node as its own cluster
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        leaders[i] = i;
        clusterWeights[i] = workingGraph.nodeWeights[i];
    }

    // Initialize tracking data structures for both strategies
    std::vector<uint8_t> markup(workingGraph.size, 0);   // For cycle detection
    std::vector<uint8_t> markdown(workingGraph.size, 0); // For cycle detection
    std::vector<uint64_t> numberOfBadNeighbors(workingGraph.size,
                                               0); // For forbidden edges
    std::vector<uint64_t> leaderOfBadNeighbors(
        workingGraph.size, UINT64_MAX); // For forbidden edges

    // Get topological order and top-level values
    const auto [topologicalOrder, topLevels] =
        workingGraph.topologicalSortAndTopLevels();

    // Process nodes in topological order
    for (uint64_t node : topologicalOrder) {
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
                if (markup[neighborId] == 1)
                    continue;

                // Check if should use forbidden edges strategy
                bool largeDegrees = checkLargeDegrees(node, neighborId);

                if (!largeDegrees) {
                    // Use cycle detection
                    if (detectCycle(node, neighborId, topLevels, leaders,
                                    markup, markdown))
                        continue;
                } else {
                    // Use forbidden edges rules
                    if (numberOfBadNeighbors[node] == 2)
                        continue;
                    if (numberOfBadNeighbors[node] == 1) {
                        uint64_t leaderOfBadNeighbor =
                            leaderOfBadNeighbors[node];
                        if (leaderOfNeighbor != leaderOfBadNeighbor)
                            continue;
                    }
                    if (numberOfBadNeighbors[neighborId] != 0)
                        continue;
                }

                // Merge nodes
                leaders[node] = leaderOfNeighbor;
                newSize--;
                clusterWeights[leaderOfNeighbor] +=
                    workingGraph.nodeWeights[node];

                // Update forbidden edges tracking
                bookKeepingForForbiddenEdges(
                    topLevels, node, sortedNeighbors, neighborId,
                    leaderOfNeighbor, markup, markdown, numberOfBadNeighbors,
                    leaderOfBadNeighbors);

                // Update cycle detection tracking
                markup[node] = markdown[neighborId] = 1;
            } else {
                // Edge from neighbor to node (similar logic with reversed
                // roles)
                if (markdown[neighborId] == 1)
                    continue;

                bool largeDegrees = checkLargeDegrees(neighborId, node);

                if (!largeDegrees) {
                    if (detectCycle(neighborId, node, topLevels, leaders,
                                    markup, markdown))
                        continue;
                } else {
                    if (numberOfBadNeighbors[node] == 2)
                        continue;
                    if (numberOfBadNeighbors[node] == 1) {
                        uint64_t leaderOfBadNeighbor =
                            leaderOfBadNeighbors[node];
                        if (leaderOfNeighbor != leaderOfBadNeighbor)
                            continue;
                    }
                    if (numberOfBadNeighbors[neighborId] != 0)
                        continue;
                }

                leaders[node] = leaderOfNeighbor;
                newSize--;
                clusterWeights[leaderOfNeighbor] +=
                    workingGraph.nodeWeights[node];

                bookKeepingForForbiddenEdges(
                    topLevels, node, sortedNeighbors, neighborId,
                    leaderOfNeighbor, markup, markdown, numberOfBadNeighbors,
                    leaderOfBadNeighbors);

                markdown[node] = markup[neighborId] = 1;
            }
            break;
        }

        // Stop if reached minimum vertices target
        if (newSize == minVertices)
            break;
    }

    return {leaders, newSize};
}