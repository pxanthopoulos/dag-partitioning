/**
 * @file RecursivePartitioner.cpp
 * @brief Implementation of recursive k-way graph partitioning
 *
 * Uses recursive bisection to create k partitions while maintaining:
 * - Acyclicity in each partition
 * - Balance constraints across all partitions
 * - Minimized edge cut between partitions
 */

#include "RecursivePartitioner.h"
#include "Graph.h"
#include "MultilevelBisectioner.h"

#include <algorithm>
#include <cassert>

namespace dag_partitioning {

namespace driver {

RecursivePartitioner::RecursivePartitioner(
    const core::Graph &graph, uint64_t partitions,
    clustering::ClusteringMethod clusteringMethod, uint64_t maxClusteringRounds,
    uint64_t minClusteringVertices, double clusteringVertexRatio,
    bisection::BisectionMethod bisectionMethod, double imbalanceRatio,
    refinement::RefinementMethod refinementMethod, uint64_t refinementPasses)
    : workingGraph(graph), partitions(partitions),
      clusteringMethod(clusteringMethod),
      maxClusteringRounds(maxClusteringRounds),
      minClusteringVertices(minClusteringVertices),
      clusteringVertexRatio(clusteringVertexRatio),
      bisectionMethod(bisectionMethod), imbalanceRatio(imbalanceRatio),
      refinementMethod(refinementMethod), refinementPasses(refinementPasses) {
    assert(partitions <= workingGraph.size &&
           "Cannot create more partitions than the number of nodes");
}

std::tuple<core::Graph, std::unordered_map<uint64_t, uint64_t>, core::Graph,
           std::unordered_map<uint64_t, uint64_t>>
RecursivePartitioner::createSubgraphs(
    const std::vector<uint8_t> &bisection) const {
    // Track sizes for resulting subgraphs
    uint64_t subGraphSize0 = 0;
    for (const auto &val : bisection) {
        if (val == 0)
            subGraphSize0++;
    }

    uint64_t subGraphSize1 = bisection.size() - subGraphSize0;

    // Create subgraphs with appropriate sizes
    core::Graph subGraph0(subGraphSize0);
    core::Graph subGraph1(subGraphSize1);

    uint64_t newNodeId0 = 0, newNodeId1 = 0;

    // Maps that map new node IDs to old node IDs
    std::unordered_map<uint64_t, uint64_t> map0, map1;

    // Map that maps old node IDs to new node IDs
    std::unordered_map<uint64_t, uint64_t> map;

    for (uint64_t i = 0; i < bisection.size(); ++i) {
        if (bisection[i] == 0) {
            subGraph0.addNode(newNodeId0, workingGraph.nodeWeights[i]);
            map0[newNodeId0] = i;
            map[i] = newNodeId0;
            newNodeId0++;
        } else {
            subGraph1.addNode(newNodeId1, workingGraph.nodeWeights[i]);
            map1[newNodeId1] = i;
            map[i] = newNodeId1;
            newNodeId1++;
        }
    }

    for (uint64_t i = 0; i < bisection.size(); ++i) {
        const auto &neighbors = workingGraph.adj[i];
        if (bisection[i] == 0) {
            for (const auto &[neighborId, edgeWeight] : neighbors) {
                if (bisection[neighborId] == 0)
                    subGraph0.addEdge(map[i], map[neighborId], edgeWeight);
            }
        } else {
            for (const auto &[neighborId, edgeWeight] : neighbors) {
                if (bisection[neighborId] == 1)
                    subGraph1.addEdge(map[i], map[neighborId], edgeWeight);
            }
        }
    }

    return {subGraph0, map0, subGraph1, map1};
}

std::pair<std::vector<uint64_t>, uint64_t> RecursivePartitioner::run() const {
    std::vector<uint64_t> partitionMapping(workingGraph.size, 0);
    uint64_t totalEdgeCut = 0;

    // Base cases
    if (workingGraph.totalWeight == 0) {
        return {partitionMapping,
                totalEdgeCut}; // All vertices stay in partition 0
    }
    if (partitions == 1) {
        return {partitionMapping,
                totalEdgeCut}; // All vertices stay in partition 0
    }
    if (partitions == 2) {
        // Direct bisection
        MultilevelBisectioner bisectioner(
            workingGraph, clusteringMethod, maxClusteringRounds,
            minClusteringVertices, clusteringVertexRatio, bisectionMethod,
            imbalanceRatio, refinementMethod, refinementPasses);
        auto [bisectionInfo, edgeCut] = bisectioner.run();

        // Convert bool vector to partition numbers (0 and 1)
        for (uint64_t i = 0; i < bisectionInfo.size(); i++) {
            partitionMapping[i] = bisectionInfo[i];
        }
        return {partitionMapping, edgeCut};
    }

    // Recursive case
    // 1. First bisect the graph
    MultilevelBisectioner bisectioner(
        workingGraph, clusteringMethod, maxClusteringRounds,
        minClusteringVertices, clusteringVertexRatio, bisectionMethod,
        imbalanceRatio, refinementMethod, refinementPasses);
    auto [bisectionInfo, edgeCut] = bisectioner.run();
    totalEdgeCut += edgeCut;

    // 2. Create subgraphs based on bisection
    // If bisection places all nodes in the same part, return like when
    // requesting 1 partition
    if (std::all_of(
            bisectionInfo.begin(), bisectionInfo.end(),
            [first = bisectionInfo[0]](bool val) { return val == first; }))
        return {partitionMapping, totalEdgeCut};

    const auto &[subGraph0, map0, subGraph1, map1] =
        createSubgraphs(bisectionInfo);

    // 3. Calculate number of parts for each subgraph
    uint64_t partsLeft = partitions / 2;
    uint64_t partsRight = partitions - partsLeft;

    // If either subgraph has fewer vertices than allocated parts,
    // give it exactly one part per vertex and give the remaining parts to the
    // other subgraph
    if (subGraph0.size < partsLeft) {
        partsLeft = subGraph0.size;
        partsRight = partitions - partsLeft;
    } else if (subGraph1.size < partsRight) {
        partsRight = subGraph1.size;
        partsLeft = partitions - partsRight;
    }

    // 4. Recursively partition each subgraph
    RecursivePartitioner left_partitioner(
        subGraph0, partsLeft, clusteringMethod, maxClusteringRounds,
        minClusteringVertices, clusteringVertexRatio, bisectionMethod,
        imbalanceRatio, refinementMethod, refinementPasses);
    auto [leftMapping, leftEdgeCut] = left_partitioner.run();

    RecursivePartitioner right_partitioner(
        subGraph1, partsRight, clusteringMethod, maxClusteringRounds,
        minClusteringVertices, clusteringVertexRatio, bisectionMethod,
        imbalanceRatio, refinementMethod, refinementPasses);
    auto [rightMapping, rightEdgeCut] = right_partitioner.run();

    // 5. Combine results into final mapping
    for (uint64_t i = 0; i < leftMapping.size(); ++i) {
        partitionMapping[map0.at(i)] = leftMapping[i];
    }
    for (uint64_t i = 0; i < rightMapping.size(); ++i) {
        partitionMapping[map1.at(i)] =
            rightMapping[i] + partsLeft; // Offset right partitions
    }

    return {partitionMapping, totalEdgeCut + leftEdgeCut + rightEdgeCut};
}

} // namespace driver

} // namespace dag_partitioning
