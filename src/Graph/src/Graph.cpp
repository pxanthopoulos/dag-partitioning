/**
 * @file Graph.cpp
 * @brief Implementation of the Graph class methods
 */

#include "Graph.h"

#include <algorithm>
#include <ostream>
#include <queue>
#include <stack>
#include <stdexcept>

namespace dag_partitioning {

namespace core {

// Constructor initializes all data structures with given size
Graph::Graph(uint64_t size)
    : size(size), adj(size), revAdj(size), nodeWeights(size, 0), totalWeight(0),
      inDegree(size, 0), maxNodeWeight(0) {}

// Add a node with given weight to the graph
void Graph::addNode(uint64_t id, uint64_t weight) {
    if (id >= size) {
        throw std::out_of_range("Node ID " + std::to_string(id) +
                                " exceeds graph size " + std::to_string(size));
    }
    nodeWeights[id] = weight;
    totalWeight += weight;

    if (weight > maxNodeWeight)
        maxNodeWeight = weight;
}

// Add a weighted edge between two nodes
void Graph::addEdge(uint64_t from, uint64_t to, uint64_t weight) {
    if (from >= size) {
        throw std::out_of_range("Source node ID " + std::to_string(from) +
                                " exceeds graph size " + std::to_string(size));
    }
    if (to >= size) {
        throw std::out_of_range("Target node ID " + std::to_string(to) +
                                " exceeds graph size " + std::to_string(size));
    }
    // Add forward edge
    adj[from].emplace_back(to, weight);
    // Add reverse edge for backward traversal
    revAdj[to].emplace_back(from, weight);
    // Update in-degree count
    inDegree[to]++;
}

// Get all neighbors (both incoming and outgoing) of a node
std::vector<std::tuple<uint64_t, uint64_t, bool>>
Graph::getNeighbors(uint64_t node) const {
    if (node >= size) {
        throw std::out_of_range("Node ID " + std::to_string(node) +
                                " exceeds graph size " + std::to_string(size));
    }
    std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors;
    neighbors.reserve(adj[node].size() + revAdj[node].size());

    // Add outgoing neighbors (marked with true)
    for (const auto &[outNeighbor, edgeWeight] : adj[node]) {
        neighbors.emplace_back(outNeighbor, edgeWeight, true);
    }

    // Add incoming neighbors (marked with false)
    for (const auto &[inNeighbor, edgeWeight] : revAdj[node]) {
        neighbors.emplace_back(inNeighbor, edgeWeight, false);
    }

    return neighbors;
}

// Get neighbors sorted by edge weight in ascending order
std::vector<std::tuple<uint64_t, uint64_t, bool>>
Graph::getNeighborsSortedByEdgeWeightAsc(uint64_t node) const {
    if (node >= size) {
        throw std::out_of_range("Node ID " + std::to_string(node) +
                                " exceeds graph size " + std::to_string(size));
    }
    std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors =
        getNeighbors(node);

    // Sort by edge weight (second element of tuple)
    std::sort(neighbors.begin(), neighbors.end(),
              [](const std::tuple<uint64_t, uint64_t, bool> &a,
                 const std::tuple<uint64_t, uint64_t, bool> &b) {
                  return std::get<1>(a) < std::get<1>(b);
              });

    return neighbors;
}

// Check for cycles using iterative DFS from a start node
bool Graph::iterativeDfsHasCycle(uint64_t start) const {
    if (start >= size) {
        throw std::out_of_range("Start node ID " + std::to_string(start) +
                                " exceeds graph size " + std::to_string(size));
    }
    std::vector<uint8_t> visited(size, 0);
    std::vector<uint8_t> inStack(size, 0); // Track nodes in current DFS path
    std::stack<std::pair<uint64_t, uint64_t>> stack;

    stack.emplace(start, 0);
    visited[start] = 1;
    inStack[start] = 1;

    while (!stack.empty()) {
        auto [node, neighborIndex] = stack.top();

        // If all neighbors of current node are processed
        if (neighborIndex >= adj[node].size()) {
            inStack[node] = 0; // Remove from current path
            stack.pop();
            continue;
        }

        // Move to next neighbor
        stack.top().second++;
        uint64_t next = adj[node][neighborIndex].first;

        if (visited[next] == 0) {
            // Process unvisited neighbor
            visited[next] = 1;
            inStack[next] = 1;
            stack.emplace(next, 0);
        } else if (inStack[next] == 1) {
            // Back edge found - cycle detected
            return true;
        }
    }
    return false;
}

// Check if the entire graph has any cycles
bool Graph::hasCycle() const {
    std::vector<uint8_t> visited(adj.size(), 0);

    // Try to find cycles starting from each unvisited node
    for (uint64_t i = 0; i < adj.size(); i++) {
        if (visited[i] == 0 && iterativeDfsHasCycle(i)) {
            return true;
        }
    }
    return false;
}

// Perform topological sort using Kahn's algorithm
std::vector<uint64_t> Graph::topologicalSort() const {
    std::vector<uint64_t> topologicalOrder;
    topologicalOrder.reserve(size);
    std::vector<uint64_t> localInDegree = inDegree;

    // Start with nodes having no incoming edges
    std::queue<uint64_t> q;
    for (uint64_t i = 0; i < size; i++) {
        if (localInDegree[i] == 0) {
            q.push(i);
        }
    }

    // Process nodes in topological order
    while (!q.empty()) {
        uint64_t curr = q.front();
        q.pop();
        topologicalOrder.push_back(curr);

        // Update in-degrees and add new nodes with zero in-degree
        for (const auto &[next, weight] : adj[curr]) {
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
        }
    }

    if (topologicalOrder.size() != size) {
        throw std::runtime_error(
            "Graph contains cycles - cannot perform topological sort");
    }
    return topologicalOrder;
}

// Compute top levels (longest path lengths from roots)
std::vector<uint64_t> Graph::computeTopLevels() const {
    std::vector<uint64_t> localInDegree = inDegree;
    std::vector<uint64_t> topLevels(size, 0);

    // Start with nodes having no incoming edges
    std::queue<uint64_t> q;
    for (uint64_t i = 0; i < size; i++) {
        if (localInDegree[i] == 0) {
            q.push(i);
        }
    }

    while (!q.empty()) {
        uint64_t curr = q.front();
        q.pop();

        // Update levels of neighbors
        for (const auto &[next, weight] : adj[curr]) {
            topLevels[next] = std::max(topLevels[next], topLevels[curr] + 1);
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
            if (localInDegree[next] < 0) {
                throw std::runtime_error("Graph contains cycles - invalid "
                                         "in-degree during computation");
            }
        }
    }

    return topLevels;
}

// Compute top levels using pre-computed topological order
std::vector<uint64_t>
Graph::computeTopLevels(const std::vector<uint64_t> &topologicalOrder) const {
    std::vector<uint64_t> topLevels(size, 0);

    // Process nodes in topological order
    for (uint64_t u : topologicalOrder) {
        if (u >= size) {
            throw std::out_of_range(
                "Node ID " + std::to_string(u) +
                " in topological order exceeds graph size " +
                std::to_string(size));
        }
        for (const auto &[v, weight] : adj[u]) {
            topLevels[v] = std::max(topLevels[v], topLevels[u] + 1);
        }
    }

    return topLevels;
}

// Compute both topological sort and top levels in a single pass
std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
Graph::topologicalSortAndTopLevels() const {
    std::vector<uint64_t> localInDegree = inDegree;
    std::vector<uint64_t> topologicalOrder;
    topologicalOrder.reserve(size);
    std::vector<uint64_t> topLevels(size, 0);

    // Initialize with nodes having no incoming edges
    std::queue<uint64_t> q;
    for (uint64_t i = 0; i < size; i++) {
        if (localInDegree[i] == 0) {
            q.push(i);
            topLevels[i] = 0;
        }
    }

    while (!q.empty()) {
        uint64_t curr = q.front();
        q.pop();
        topologicalOrder.push_back(curr);

        // Update both top levels and in-degrees
        for (const auto &[next, weight] : adj[curr]) {
            topLevels[next] = std::max(topLevels[next], topLevels[curr] + 1);
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
            if (localInDegree[next] < 0) {
                throw std::runtime_error("Graph contains cycles - invalid "
                                         "in-degree during computation");
            }
        }
    }

    if (topologicalOrder.size() != size) {
        throw std::runtime_error(
            "Graph contains cycles - cannot perform topological sort");
    }
    return {topologicalOrder, topLevels};
}

// Compute shortest paths using Dijkstra's algorithm
std::vector<uint64_t> Graph::distancesFromNode(uint64_t startNode,
                                               bool reverseGraph) const {
    if (startNode >= size) {
        throw std::out_of_range("Start node ID " + std::to_string(startNode) +
                                " exceeds graph size " + std::to_string(size));
    }
    std::vector<uint64_t> distances(size, UINT64_MAX);

    // Priority queue for Dijkstra's algorithm
    std::priority_queue<std::pair<uint64_t, uint64_t>,
                        std::vector<std::pair<uint64_t, uint64_t>>,
                        std::greater<>>
        pq;

    distances[startNode] = 0;
    pq.emplace(startNode, 0);

    while (!pq.empty()) {
        uint64_t currentNode = pq.top().first;
        uint64_t currentDist = pq.top().second;
        pq.pop();

        // Process edges in forward or reverse direction
        for (const auto &edge :
             (reverseGraph ? revAdj[currentNode] : adj[currentNode])) {
            uint64_t neighborId = edge.first;
            uint64_t edgeWeight = edge.second;

            uint64_t newDist = currentDist + edgeWeight;

            // Update distance if shorter path found
            if (newDist < distances[neighborId]) {
                distances[neighborId] = newDist;
                pq.emplace(neighborId, newDist);
            }
        }
    }

    return distances;
}

std::vector<uint64_t> Graph::partitionTopologicalOrder(
    const std::vector<uint64_t> &partitionInfo) const {
    if (partitionInfo.size() != size) {
        throw std::invalid_argument(
            "Partition info size (" + std::to_string(partitionInfo.size()) +
            ") must match graph size (" + std::to_string(size) + ")");
    }

    // Step 1: Create a graph of partitions
    // Find the number of partitions
    uint64_t numPartitions = 0;
    for (uint64_t p : partitionInfo) {
        numPartitions = std::max(numPartitions, p + 1);
    }

    // Create a condensed graph where nodes are partitions
    std::vector<std::vector<uint64_t>> coarseGraph(numPartitions);
    std::vector<uint64_t> coarseInDegree(numPartitions, 0);

    // For each edge in the original graph, add an edge between partitions if
    // needed
    for (uint64_t from = 0; from < size; ++from) {
        for (const auto &[to, _] : adj[from]) {
            uint64_t fromPartition = partitionInfo[from];
            uint64_t toPartition = partitionInfo[to];

            // If edge crosses partition boundaries, add it to partition graph
            if (fromPartition != toPartition) {
                // Check if this edge already exists in the partition graph
                bool edgeExists = false;
                for (uint64_t existingTo : coarseGraph[fromPartition]) {
                    if (existingTo == toPartition) {
                        edgeExists = true;
                        break;
                    }
                }

                if (!edgeExists) {
                    coarseGraph[fromPartition].push_back(toPartition);
                    coarseInDegree[toPartition]++;
                }
            }
        }
    }

    // Step 2: Topologically sort the partition graph
    std::vector<uint64_t> partitionTopologicalOrder;
    partitionTopologicalOrder.reserve(numPartitions);

    // Start with nodes having no incoming edges
    std::queue<uint64_t> q;
    for (uint64_t i = 0; i < numPartitions; i++) {
        if (coarseInDegree[i] == 0) {
            q.push(i);
        }
    }

    // Process nodes in topological order
    while (!q.empty()) {
        uint64_t curr = q.front();
        q.pop();
        partitionTopologicalOrder.push_back(curr);

        // Update in-degrees and add new nodes with zero in-degree
        for (uint64_t next : coarseGraph[curr]) {
            coarseInDegree[next]--;
            if (coarseInDegree[next] == 0) {
                q.push(next);
            }
        }
    }

    if (partitionTopologicalOrder.size() != numPartitions) {
        throw std::runtime_error(
            "Partition graph contains cycles - cannot determine ordering");
    }

    return partitionTopologicalOrder;
}

// Print graph information to output stream
void Graph::print(std::ostream &os) const {
    if (adj.size() != nodeWeights.size()) {
        throw std::runtime_error(
            "Graph internal inconsistency: adjacency list size (" +
            std::to_string(adj.size()) + ") != node weights size (" +
            std::to_string(nodeWeights.size()) + ")");
    }
    if (size != adj.size()) {
        throw std::runtime_error(
            "Graph internal inconsistency: declared size (" +
            std::to_string(size) + ") != adjacency list size (" +
            std::to_string(adj.size()) + ")");
    }

    // Print each node's information
    for (uint64_t i = 0; i < size; ++i) {
        const auto &adjList = adj[i];
        uint64_t weight = nodeWeights[i];
        os << "Node: " << i << "\nWeight: " << weight << "\n";
        os << "Neighbors: ";
        for (const auto &[neighborId, edgeWeight] : adjList) {
            os << neighborId << "(" << edgeWeight << "), ";
        }
        os << "\n\n";
    }
}

// Stream output operator
std::ostream &operator<<(std::ostream &os, const Graph &graph) {
    graph.print(os);
    return os;
}

} // namespace core

} // namespace dag_partitioning