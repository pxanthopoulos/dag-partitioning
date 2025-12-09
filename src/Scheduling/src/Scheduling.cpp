/**
 * @file Scheduling.cpp
 * @brief Implementation of base scheduling class methods
 */

#include "Scheduling.h"
#include "Graph.h"
#include "robin_hood.h"

#include <cassert>
#include <iostream>
#include <queue>

namespace dag_partitioning {

namespace scheduling {

namespace packing {

uint64_t packBestFit(std::vector<Tensor> &tensors) {
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

} // namespace packing

namespace ilp {

struct ILPGraph {
    uint64_t num_nodes;
    std::vector<uint64_t> tensor_size; // S_i: output tensor size
    std::vector<uint64_t> extra_size;  // ES_i: extra memory during execution
    std::vector<robin_hood::unordered_set<uint64_t>>
        input_tensors; // IN_i: input tensors

    // Computed topology
    std::vector<robin_hood::unordered_set<uint64_t>> ancestors;
    std::vector<robin_hood::unordered_set<uint64_t>> descendants;
    std::vector<robin_hood::unordered_set<uint64_t>> tensor_users;

    void computeTopology() {
        computeAncestors();
        computeDescendants();
    }

    void computeAncestors() {
        ancestors.clear();
        ancestors.resize(num_nodes);
        tensor_users.clear();
        tensor_users.resize(num_nodes);

        // Compute in-degree for topological sort
        std::vector<uint64_t> in_degree(num_nodes, 0);
        for (uint64_t i = 0; i < num_nodes; i++) {
            in_degree[i] = input_tensors[i].size();
        }

        // Topological sort using Kahn's algorithm
        std::queue<uint64_t> queue;
        for (uint64_t i = 0; i < num_nodes; i++) {
            if (in_degree[i] == 0) {
                queue.push(i);
            }
        }

        while (!queue.empty()) {
            uint64_t node = queue.front();
            queue.pop();

            // For this node, its ancestors are the union of:
            // 1. Its direct inputs
            // 2. All ancestors of its inputs
            for (uint64_t input : input_tensors[node]) {
                ancestors[node].insert(input);
                ancestors[node].insert(ancestors[input].begin(),
                                       ancestors[input].end());
                tensor_users[input].insert(node);
            }

            // Update in-degrees and enqueue nodes
            for (uint64_t i = 0; i < num_nodes; i++) {
                if (input_tensors[i].count(node)) {
                    in_degree[i]--;
                    if (in_degree[i] == 0) {
                        queue.push(i);
                    }
                }
            }
        }
    }

    void computeDescendants() {
        descendants.clear();
        descendants.resize(num_nodes);

        // Build descendants from ancestors: if j is ancestor of i, then i is
        // descendant of j
        for (uint64_t i = 0; i < num_nodes; i++) {
            for (uint64_t ancestor : ancestors[i]) {
                descendants[ancestor].insert(i);
            }
        }
    }

    void print(std::ostream &os) const {
        os << "Nodes: " << num_nodes << std::endl;

        for (uint64_t i = 0; i < num_nodes; i++) {
            os << "\nNode: " << i << std::endl;
            os << "    tensor_size=" << tensor_size[i]
               << ", extra_size=" << extra_size[i] << std::endl;

            os << "    inputs={";
            for (auto it = input_tensors[i].begin();
                 it != input_tensors[i].end(); ++it) {
                if (it != input_tensors[i].begin())
                    os << ", ";
                os << *it;
            }
            os << "}, users={";
            for (auto it = tensor_users[i].begin(); it != tensor_users[i].end();
                 ++it) {
                if (it != tensor_users[i].begin())
                    os << ", ";
                os << *it;
            }
            os << "}" << std::endl;

            os << "    anc={";
            for (auto it = ancestors[i].begin(); it != ancestors[i].end();
                 ++it) {
                if (it != ancestors[i].begin())
                    os << ", ";
                os << *it;
            }
            os << "}, des={";
            for (auto it = descendants[i].begin(); it != descendants[i].end();
                 ++it) {
                if (it != descendants[i].begin())
                    os << ", ";
                os << *it;
            }
            os << "}" << std::endl;
        }
    }

    friend std::ostream &operator<<(std::ostream &os, const ILPGraph &graph) {
        graph.print(os);
        return os;
    }
};

ILPGraph buildILPGraph(const core::Graph &graph) {
    ILPGraph ilpGraph;
    ilpGraph.num_nodes = graph.size;
    ilpGraph.tensor_size.resize(graph.size);
    ilpGraph.extra_size.resize(graph.size);
    ilpGraph.input_tensors.resize(graph.size);

    for (uint64_t node = 0; node < graph.size; ++node) {
        ilpGraph.tensor_size[node] = 0;
        ilpGraph.extra_size[node] = graph.nodeWeights[node];

        for (const auto &[neighbor, edgeWeight] : graph.adj[node]) {
            ilpGraph.input_tensors[neighbor].insert(node);
            ilpGraph.tensor_size[node] =
                std::max(ilpGraph.tensor_size[node], edgeWeight);
        }
    }

    ilpGraph.computeTopology();

    return ilpGraph;
}

} // namespace ilp

Scheduler::Scheduler(const core::Graph &originalGraph,
                     std::vector<uint64_t> &partitionMapping)
    : originalGraph(originalGraph), partitionMapping(partitionMapping) {}

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
    std::vector<std::vector<packing::Tensor>> partitionTensors(numPartitions);

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
            packing::Tensor tensor;
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
                  [](const packing::Tensor &a, const packing::Tensor &b) {
                      if (a.birth != b.birth)
                          return a.birth < b.birth;
                      return a.death < b.death;
                  });

        // Best-fit packing algorithm
        uint64_t peakMemory = packing::packBestFit(tensors);

        partitionWeights[partition] = peakMemory;
    }

    return partitionWeights;
}

void Scheduler::buildCoarseGraph() {
    // Find the number of partitions (assuming contiguous IDs from 0 to k-1)
    uint64_t numPartitions =
        *std::max_element(partitionMapping.begin(), partitionMapping.end()) + 1;
    assert(numPartitions != 0 && "Number of partitions must be > 0");

    coarseGraph = std::make_unique<core::Graph>(numPartitions);

    std::vector<uint64_t> partitionWeights = calculatePartitionWeights();

    // Add nodes to coarse graph - partition ID is the node ID
    for (uint64_t partitionId = 0; partitionId < numPartitions; ++partitionId) {
        coarseGraph->addNode(partitionId, partitionWeights[partitionId]);
    }

    // For each edge in the original graph, create an edge between partitions
    // if needed, with total size of edges between partitions
    std::vector<std::unordered_map<uint64_t, uint64_t>> newEdges(numPartitions);
    for (uint64_t from = 0; from < originalGraph.size; ++from) {
        for (const auto &[to, edgeWeight] : originalGraph.adj[from]) {
            uint64_t fromPartition = partitionMapping[from];
            uint64_t toPartition = partitionMapping[to];

            if (fromPartition != toPartition) {
                if (newEdges[fromPartition].find(toPartition) !=
                    newEdges[fromPartition].end()) {
                    newEdges[fromPartition][toPartition] += edgeWeight;
                } else {
                    newEdges[fromPartition].emplace(toPartition, edgeWeight);
                }
            }
        }
    }

    for (uint64_t partitionId = 0; partitionId < newEdges.size();
         ++partitionId) {
        for (const auto &[neighborId, edgeWeight] : newEdges[partitionId]) {
            coarseGraph->addEdge(partitionId, neighborId, edgeWeight);
        }
    }
}

void Scheduler::run() {
    buildCoarseGraph();
    std::cout << *coarseGraph;
    ilp::ILPGraph ilpGraph = ilp::buildILPGraph(*coarseGraph);
    std::cout << ilpGraph;
}

const std::unique_ptr<core::Graph> &Scheduler::getCoarseGraph() const {
    return coarseGraph;
}

} // namespace scheduling

} // namespace dag_partitioning