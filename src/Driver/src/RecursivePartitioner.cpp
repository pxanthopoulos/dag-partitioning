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
#include "Bisection.h"
#include "Clustering.h"
#include "Graph.h"
#include "MultilevelBisectioner.h"
#include "Refinement.h"

#include <algorithm>
#include <cassert>
#include <set>
#include <stdexcept>
#include <thread>

#ifdef DAG_PARTITIONING_OPENMP_ENABLED
#include <omp.h>
#endif

namespace dag_partitioning {

namespace driver {

RecursivePartitioner::RecursivePartitioner(
    const core::Graph &graph, uint64_t partitions,
    clustering::ClusteringMethod clusteringMethod, uint64_t maxClusteringRounds,
    uint64_t minClusteringVertices, double clusteringVertexRatio,
    bisection::BisectionMethod bisectionMethod, double imbalanceRatio,
    refinement::RefinementMethod refinementMethod, uint64_t refinementPasses,
    bool enableParallel, uint64_t minSizeForParallel, uint64_t maxParallelDepth)
    : workingGraph(graph), partitions(partitions),
      clusteringMethod(clusteringMethod),
      maxClusteringRounds(maxClusteringRounds),
      minClusteringVertices(minClusteringVertices),
      clusteringVertexRatio(clusteringVertexRatio),
      bisectionMethod(bisectionMethod), imbalanceRatio(imbalanceRatio),
      refinementMethod(refinementMethod), refinementPasses(refinementPasses),
      enableParallel(enableParallel), minSizeForParallel(minSizeForParallel),
      maxParallelDepth(maxParallelDepth) {
    if (partitions > workingGraph.size) {
        throw std::invalid_argument("Cannot create more partitions (" +
                                    std::to_string(partitions) +
                                    ") than the number of nodes (" +
                                    std::to_string(workingGraph.size) + ")");
    }
}

RecursivePartitioner::RecursivePartitioner(
    const core::Graph &graph, uint64_t partitions, std::string clusteringMethod,
    uint64_t maxClusteringRounds, uint64_t minClusteringVertices,
    double clusteringVertexRatio, std::string bisectionMethod,
    double imbalanceRatio, std::string refinementMethod,
    uint64_t refinementPasses, bool enableParallel, uint64_t minSizeForParallel,
    uint64_t maxParallelDepth)
    : workingGraph(graph), partitions(partitions),
      maxClusteringRounds(maxClusteringRounds),
      minClusteringVertices(minClusteringVertices),
      clusteringVertexRatio(clusteringVertexRatio),
      imbalanceRatio(imbalanceRatio), refinementPasses(refinementPasses),
      enableParallel(enableParallel), minSizeForParallel(minSizeForParallel),
      maxParallelDepth(maxParallelDepth) {
    if (partitions > workingGraph.size) {
        throw std::invalid_argument("Cannot create more partitions (" +
                                    std::to_string(partitions) +
                                    ") than the number of nodes (" +
                                    std::to_string(workingGraph.size) + ")");
    }

    if (clusteringMethod == "FORB") {
        this->clusteringMethod = clustering::ClusteringMethod::FORB;
    } else if (clusteringMethod == "CYC") {
        this->clusteringMethod = clustering::ClusteringMethod::CYC;
    } else if (clusteringMethod == "HYB") {
        this->clusteringMethod = clustering::ClusteringMethod::HYB;
    } else {
        this->clusteringMethod = clustering::ClusteringMethod::HYB;
    }

    if (bisectionMethod == "GGG") {
        this->bisectionMethod = bisection::BisectionMethod::GGG;
    } else if (bisectionMethod == "UNDIRMETIS") {
        this->bisectionMethod = bisection::BisectionMethod::UNDIRMETIS;
    } else if (bisectionMethod == "UNDIRSCOTCH") {
        this->bisectionMethod = bisection::BisectionMethod::UNDIRSCOTCH;
    } else if (bisectionMethod == "UNDIRBOTH") {
        this->bisectionMethod = bisection::BisectionMethod::UNDIRBOTH;
    } else {
        this->bisectionMethod = bisection::BisectionMethod::UNDIRBOTH;
    }

    if (refinementMethod == "BOUNDARYFM") {
        this->refinementMethod = refinement::RefinementMethod::BOUNDARYFM;
    } else if (refinementMethod == "BOUNDARYKL") {
        this->refinementMethod = refinement::RefinementMethod::BOUNDARYKL;
    } else if (refinementMethod == "MIXED") {
        this->refinementMethod = refinement::RefinementMethod::MIXED;
    } else {
        this->refinementMethod = refinement::RefinementMethod::MIXED;
    }
}

std::tuple<core::Graph, robin_hood::unordered_map<uint64_t, uint64_t>,
           core::Graph, robin_hood::unordered_map<uint64_t, uint64_t>>
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
    robin_hood::unordered_map<uint64_t, uint64_t> map0, map1;
    map0.reserve(subGraphSize0);
    map1.reserve(subGraphSize1);

    // Map that maps old node IDs to new node IDs
    robin_hood::unordered_map<uint64_t, uint64_t> map;
    map.reserve(bisection.size());

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

std::pair<std::vector<uint64_t>, uint64_t>
RecursivePartitioner::run(uint64_t currentDepth) const {
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

    auto subgraphs = createSubgraphs(bisectionInfo);
    const auto &subGraph0 = std::get<0>(subgraphs);
    const auto &map0 = std::get<1>(subgraphs);
    const auto &subGraph1 = std::get<2>(subgraphs);
    const auto &map1 = std::get<3>(subgraphs);

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

    // 4. Recursively partition each subgraph (parallel if conditions are met)
    std::vector<uint64_t> leftMapping;
    leftMapping.reserve(subGraph0.size);
    uint64_t leftEdgeCut = 0;
    std::vector<uint64_t> rightMapping;
    rightMapping.reserve(subGraph1.size);
    uint64_t rightEdgeCut = 0;

    // Check if we should parallelize based on size, settings, and depth
    bool shouldParallelize = enableParallel &&
                             currentDepth < maxParallelDepth &&
                             (subGraph0.size >= minSizeForParallel ||
                              subGraph1.size >= minSizeForParallel);

    if (shouldParallelize) {
#ifdef DAG_PARTITIONING_OPENMP_ENABLED
// Use OpenMP tasks for better task scheduling
#pragma omp parallel
        {
#pragma omp single nowait
            {
                int leftPriority = subGraph0.size > subGraph1.size ? 1 : 0;
                int rightPriority = subGraph1.size > subGraph0.size ? 1 : 0;

// Create task for left subgraph
#pragma omp task priority(leftPriority) shared(leftMapping, leftEdgeCut)
                {
                    RecursivePartitioner left_partitioner(
                        subGraph0, partsLeft, clusteringMethod,
                        maxClusteringRounds, minClusteringVertices,
                        clusteringVertexRatio, bisectionMethod, imbalanceRatio,
                        refinementMethod, refinementPasses, enableParallel,
                        minSizeForParallel, maxParallelDepth);
                    auto result = left_partitioner.run(currentDepth + 1);
                    leftMapping = std::move(result.first);
                    leftEdgeCut = result.second;
                }

// Create task for right subgraph
#pragma omp task priority(rightPriority) shared(rightMapping, rightEdgeCut)
                {
                    RecursivePartitioner right_partitioner(
                        subGraph1, partsRight, clusteringMethod,
                        maxClusteringRounds, minClusteringVertices,
                        clusteringVertexRatio, bisectionMethod, imbalanceRatio,
                        refinementMethod, refinementPasses, enableParallel,
                        minSizeForParallel, maxParallelDepth);
                    auto result = right_partitioner.run(currentDepth + 1);
                    rightMapping = std::move(result.first);
                    rightEdgeCut = result.second;
                }

// Wait for both tasks to complete
#pragma omp taskwait
            }
        }
#else
        // Fallback to std::thread when OpenMP is not available
        std::thread leftThread([&subGraph0, partsLeft, this, &leftMapping,
                                &leftEdgeCut, currentDepth]() {
            RecursivePartitioner left_partitioner(
                subGraph0, partsLeft, clusteringMethod, maxClusteringRounds,
                minClusteringVertices, clusteringVertexRatio, bisectionMethod,
                imbalanceRatio, refinementMethod, refinementPasses,
                enableParallel, minSizeForParallel, maxParallelDepth);
            auto result = left_partitioner.run(currentDepth + 1);
            leftMapping = std::move(result.first);
            leftEdgeCut = result.second;
        });

        std::thread rightThread([&subGraph1, partsRight, this, &rightMapping,
                                 &rightEdgeCut, currentDepth]() {
            RecursivePartitioner right_partitioner(
                subGraph1, partsRight, clusteringMethod, maxClusteringRounds,
                minClusteringVertices, clusteringVertexRatio, bisectionMethod,
                imbalanceRatio, refinementMethod, refinementPasses,
                enableParallel, minSizeForParallel, maxParallelDepth);
            auto result = right_partitioner.run(currentDepth + 1);
            rightMapping = std::move(result.first);
            rightEdgeCut = result.second;
        });

        leftThread.join();
        rightThread.join();
#endif
    } else {
        // Sequential execution for small subgraphs or when parallelization is
        // disabled
        RecursivePartitioner left_partitioner(
            subGraph0, partsLeft, clusteringMethod, maxClusteringRounds,
            minClusteringVertices, clusteringVertexRatio, bisectionMethod,
            imbalanceRatio, refinementMethod, refinementPasses, enableParallel,
            minSizeForParallel, maxParallelDepth);
        auto leftResult = left_partitioner.run(currentDepth + 1);
        leftMapping = std::move(leftResult.first);
        leftEdgeCut = leftResult.second;

        RecursivePartitioner right_partitioner(
            subGraph1, partsRight, clusteringMethod, maxClusteringRounds,
            minClusteringVertices, clusteringVertexRatio, bisectionMethod,
            imbalanceRatio, refinementMethod, refinementPasses, enableParallel,
            minSizeForParallel, maxParallelDepth);
        auto rightResult = right_partitioner.run(currentDepth + 1);
        rightMapping = std::move(rightResult.first);
        rightEdgeCut = rightResult.second;
    }

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

std::pair<std::vector<uint64_t>, uint64_t> RecursivePartitioner::resetPartitionNumbers(
    std::pair<std::vector<uint64_t>, uint64_t> partitionResult) {
    auto &partitionMapping = partitionResult.first;

    // Find unique partition IDs
    std::set<uint64_t> uniquePartitions(partitionMapping.begin(),
                                        partitionMapping.end());

    // Create mapping from old partition IDs to new contiguous IDs
    std::unordered_map<uint64_t, uint64_t> idMapping;
    uint64_t newId = 0;
    for (const auto &oldId : uniquePartitions) {
        idMapping[oldId] = newId++;
    }

    // Update partition mapping to use new contiguous IDs
    for (auto &partId : partitionMapping) {
        partId = idMapping[partId];
    }

    return partitionResult;
}

} // namespace driver

} // namespace dag_partitioning
