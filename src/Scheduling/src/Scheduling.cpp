/**
 * @file Scheduling.cpp
 * @brief Implementation of base scheduling class methods
 */

#include "Scheduling.h"
#include "Graph.h"
#include "robin_hood.h"

#include <cassert>
#include <iostream>

namespace dag_partitioning {

namespace scheduling {

Scheduler::Scheduler(const core::Graph &originalGraph,
                     std::vector<uint64_t> &partitionMapping)
    : originalGraph(originalGraph), partitionMapping(partitionMapping) {}

uint64_t Scheduler::packBestFit(std::vector<Tensor> &tensors) {
    auto overlaps = [](const Tensor &a, const Tensor &b) {
        return !(a.death < b.birth || b.death < a.birth);
    };

    auto verify = [&overlaps](const std::vector<Tensor> &tensors) {
        for (size_t i = 0; i < tensors.size(); ++i) {
            for (size_t j = i + 1; j < tensors.size(); ++j) {
                if (overlaps(tensors[i], tensors[j])) {
                    uint64_t end1 = tensors[i].offset + tensors[i].size;
                    uint64_t end2 = tensors[j].offset + tensors[j].size;
                    if (!(end1 <= tensors[j].offset ||
                          end2 <= tensors[i].offset)) {
                        std::cerr << "CONFLICT: tensor " << i << " and " << j
                                  << "\n";
                        return false;
                    }
                }
            }
        }
        return true;
    };

    if (tensors.empty())
        return 0;

    for (size_t i = 0; i < tensors.size(); ++i) {
        Tensor &t = tensors[i];

        std::vector<std::pair<uint64_t, uint64_t>> occupied;
        for (size_t j = 0; j < i; ++j) { // all previously placed
            if (overlaps(t, tensors[j])) {
                occupied.push_back(
                    {tensors[j].offset, tensors[j].offset + tensors[j].size});
            }
        }

        // Sort and merge overlapping intervals
        std::sort(occupied.begin(), occupied.end());
        std::vector<std::pair<uint64_t, uint64_t>> merged;
        for (auto &iv : occupied) {
            if (merged.empty() || merged.back().second < iv.first) {
                merged.push_back(iv);
            } else {
                merged.back().second =
                    std::max(merged.back().second, iv.second);
            }
        }

        // Best-fit: find smallest gap >= t.size
        uint64_t bestOffset = 0;
        uint64_t bestWaste = UINT64_MAX;

        if (merged.empty()) {
            bestOffset = 0;
        } else {
            // Gap before first interval [0, merged[0].first)
            if (merged[0].first >= t.size &&
                merged[0].first - t.size < bestWaste) {
                bestWaste = merged[0].first - t.size;
                bestOffset = 0;
            }

            // Gaps between intervals
            for (size_t i = 0; i + 1 < merged.size(); ++i) {
                uint64_t gapSize = merged[i + 1].first - merged[i].second;
                if (gapSize >= t.size && gapSize - t.size < bestWaste) {
                    bestWaste = gapSize - t.size;
                    bestOffset = merged[i].second;
                }
            }

            // If no gap fits, place after last interval
            if (bestWaste == UINT64_MAX) {
                bestOffset = merged.back().second;
            }
        }

        t.offset = bestOffset;
    }

    uint64_t peak = 0;
    for (const Tensor &t : tensors) {
        peak = std::max(peak, t.offset + t.size);
    }

    assert(verify(tensors) && "Memory packing verification failed");

    return peak;
}

std::vector<uint64_t> Scheduler::calculatePartitionWeights() const {
    std::vector<uint64_t> topologicalOrder = originalGraph.topologicalSort();

    // Create global position mapping
    std::vector<uint64_t> globalNodeToPosition(originalGraph.size);
    for (uint64_t pos = 0; pos < topologicalOrder.size(); ++pos) {
        globalNodeToPosition[topologicalOrder[pos]] = pos;
    }

    // Find number of partitions
    uint64_t numPartitions = 0;
    for (uint64_t partId : partitionMapping) {
        numPartitions = std::max(numPartitions, partId + 1);
    }

    // Group nodes by partition and sort by global topological position
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> partitionNodes(
        numPartitions);
    for (uint64_t node = 0; node < originalGraph.size; ++node) {
        uint64_t partition = partitionMapping[node];
        partitionNodes[partition].emplace_back(globalNodeToPosition[node],
                                               node);
    }

    // Sort nodes within each partition by their global topological position
    for (auto &nodes : partitionNodes) {
        std::sort(nodes.begin(), nodes.end());
    }

    // Create local position mapping: nodeToLocalPosition[node] = local position
    // within its partition
    std::vector<uint64_t> nodeToLocalPosition(originalGraph.size);
    for (uint64_t partition = 0; partition < numPartitions; ++partition) {
        for (uint64_t localPos = 0; localPos < partitionNodes[partition].size();
             ++localPos) {
            uint64_t node = partitionNodes[partition][localPos].second;
            nodeToLocalPosition[node] = localPos;
        }
    }

    // Collect tensors per partition
    std::vector<std::vector<Tensor>> partitionTensors(numPartitions);

    for (uint64_t producer = 0; producer < originalGraph.size; ++producer) {
        uint64_t producerPartition = partitionMapping[producer];

        // Find all consumers of this producer within the same partition
        uint64_t lastConsumerPosition =
            nodeToLocalPosition[producer]; // Birth time
        bool hasInternalConsumers = false;

        for (const auto &[consumer, edgeWeight] : originalGraph.adj[producer]) {
            // Only consider consumers in the same partition
            if (partitionMapping[consumer] == producerPartition) {
                hasInternalConsumers = true;
                lastConsumerPosition = std::max(lastConsumerPosition,
                                                nodeToLocalPosition[consumer]);
            }
        }

        // If there are internal consumers, add tensor to partition
        if (hasInternalConsumers) {
            Tensor tensor;
            tensor.producer = producer;
            tensor.size = originalGraph.nodeWeights[producer];
            tensor.birth = nodeToLocalPosition[producer];
            tensor.death = lastConsumerPosition;
            tensor.offset = 0; // Will be assigned during packing
            partitionTensors[producerPartition].push_back(tensor);
        }
    }

    // Perform best-fit packing for each partition
    std::vector<uint64_t> partitionWeights(numPartitions, 0);

    for (uint64_t partition = 0; partition < numPartitions; ++partition) {
        auto &tensors = partitionTensors[partition];

        if (tensors.empty()) {
            continue;
        }

        // Sort tensors by birth time, then by death time (for deterministic
        // packing)
        std::sort(tensors.begin(), tensors.end(),
                  [](const Tensor &a, const Tensor &b) {
                      if (a.birth != b.birth)
                          return a.birth < b.birth;
                      return a.death < b.death;
                  });

        // Best-fit packing algorithm
        uint64_t peakMemory = packBestFit(tensors);

        partitionWeights[partition] = peakMemory;

        // Pretty print the packing
        std::cout << "\nPartition " << partition << " memory packing:\n";
        std::cout << "Peak memory: " << peakMemory << "\n";
        for (const auto &tensor : tensors) {
            std::cout << "  Tensor (producer=" << tensor.producer
                      << ", size=" << tensor.size << ", lifetime=["
                      << tensor.birth << "," << tensor.death << "]"
                      << ") -> offset=" << tensor.offset << "\n";
        }
    }

    return partitionWeights;
}

void Scheduler::buildCoarseGraph() {
    // Find the node IDs of the new graph
    robin_hood::unordered_set<uint64_t> partitionIds(partitionMapping.begin(),
                                                     partitionMapping.end());

    size_t coarseGraphSize = partitionIds.size();
    assert(coarseGraphSize != 0 && "Size of coarse graph must be > 0");

    coarseGraph = std::make_unique<core::Graph>(coarseGraphSize);

    std::vector<uint64_t> partitionWeights = calculatePartitionWeights();

    uint64_t nodeId = 0;
    partitionsToNewNodes.assign(coarseGraphSize, 0);
    newNodesToPartitions.assign(coarseGraphSize, 0);
    for (const uint64_t &partitionId : partitionIds) {
        newNodesToPartitions[nodeId] = partitionId;
        partitionsToNewNodes[partitionId] = nodeId;
        coarseGraph->addNode(nodeId++, partitionWeights[partitionId]);
    }

    assert(nodeId == coarseGraphSize && "IDs are not contiguous");

    // For each edge in the original graph, create an edge between partitions
    // if needed, with total size of edges between partitions
    std::vector<std::unordered_map<uint64_t, uint64_t>> newEdges(
        coarseGraphSize);
    for (uint64_t from = 0; from < originalGraph.size; ++from) {
        for (const auto &[to, edgeWeight] : originalGraph.adj[from]) {
            uint64_t fromPartition = partitionMapping[from];
            uint64_t fromPartitionNewNode = partitionsToNewNodes[fromPartition];
            uint64_t toPartition = partitionMapping[to];
            uint64_t toPartitionNewNode = partitionsToNewNodes[toPartition];

            if (fromPartition != toPartition) {
                if (newEdges[fromPartitionNewNode].find(toPartitionNewNode) !=
                    newEdges[fromPartitionNewNode].end()) {
                    newEdges[fromPartitionNewNode][toPartitionNewNode] +=
                        edgeWeight;
                } else {
                    newEdges[fromPartitionNewNode].emplace(toPartitionNewNode,
                                                           edgeWeight);
                }
            }
        }
    }

    for (uint64_t newNodeId = 0; newNodeId < newEdges.size(); ++newNodeId) {
        for (const auto &[neighborId, edgeWeight] : newEdges[newNodeId]) {
            coarseGraph->addEdge(newNodeId, neighborId, edgeWeight);
        }
    }
}

void Scheduler::run() { buildCoarseGraph(); }

const std::vector<uint64_t> &Scheduler::getNewNodesToPartitions() const {
    return newNodesToPartitions;
};

const std::vector<uint64_t> &Scheduler::getPartitionsToNewNodes() const {
    return partitionsToNewNodes;
};

const std::unique_ptr<core::Graph> &Scheduler::getCoarseGraph() const {
    return coarseGraph;
}

} // namespace scheduling

} // namespace dag_partitioning