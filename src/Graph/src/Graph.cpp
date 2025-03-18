/**
 * @file Graph.cpp
 * @brief Implementation of the Graph class methods
 */

#include "Graph.h"
#include <fstream>
#include <regex>
#include <unordered_map>
#include <cassert>
#include <queue>

// Constructor initializes all data structures with given size
Graph::Graph(uint64_t size)
        : size(size), adj(size), revAdj(size), nodeWeights(size, 0), totalWeight(0), inDegree(size, 0),
          maxNodeWeight(0) {
}

// Add a node with given weight to the graph
void Graph::addNode(uint64_t id, uint64_t weight) {
    assert(id < size && "Node ID must be smaller than graph size");
    nodeWeights[id] = weight;
    totalWeight += weight;

    if (weight > maxNodeWeight) maxNodeWeight = weight;
}

// Add a weighted edge between two nodes
void Graph::addEdge(uint64_t from, uint64_t to, uint64_t weight) {
    assert(from < size && "Node ID `from` must be smaller than graph size");
    assert(to < size && "Node ID `to` must be smaller than graph size");
    // Add forward edge
    adj[from].emplace_back(to, weight);
    // Add reverse edge for backward traversal
    revAdj[to].emplace_back(from, weight);
    // Update in-degree count
    inDegree[to]++;
}

// Get all neighbors (both incoming and outgoing) of a node
std::vector<std::tuple<uint64_t, uint64_t, bool>> Graph::getNeighbors(uint64_t node) const {
    assert(node < size && "node id must be smaller than graph size");
    std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors;

    // Add outgoing neighbors (marked with true)
    for (const auto &[outNeighbor, edgeWeight]: adj[node]) {
        neighbors.emplace_back(outNeighbor, edgeWeight, true);
    }

    // Add incoming neighbors (marked with false)
    for (const auto &[inNeighbor, edgeWeight]: revAdj[node]) {
        neighbors.emplace_back(inNeighbor, edgeWeight, false);
    }

    return neighbors;
}

// Get neighbors sorted by edge weight in ascending order
std::vector<std::tuple<uint64_t, uint64_t, bool>> Graph::getNeighborsSortedByEdgeWeightAsc(uint64_t node) const {
    assert(node < size && "node id must be smaller than graph size");
    std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors = getNeighbors(node);

    // Sort by edge weight (second element of tuple)
    std::sort(neighbors.begin(), neighbors.end(),
              [](const std::tuple<uint64_t, uint64_t, bool> &a, const std::tuple<uint64_t, uint64_t, bool> &b) {
                  return std::get<1>(a) < std::get<1>(b);
              });

    return neighbors;
}

// Check for cycles using iterative DFS from a start node
bool Graph::iterativeDfsHasCycle(uint64_t start) const {
    assert(start < size && "start node id must be smaller than graph size");
    std::vector<uint8_t> visited(size, 0);
    std::vector<uint8_t> inStack(size, 0);  // Track nodes in current DFS path
    std::stack<std::pair<uint64_t, uint64_t>> stack;

    stack.emplace(start, 0);
    visited[start] = 1;
    inStack[start] = 1;

    while (!stack.empty()) {
        auto [node, neighborIndex] = stack.top();

        // If all neighbors of current node are processed
        if (neighborIndex >= adj[node].size()) {
            inStack[node] = 0;  // Remove from current path
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
        for (const auto &[next, weight]: adj[curr]) {
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
        }
    }

    assert(topologicalOrder.size() == size && "Graph is cyclic");
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
        for (const auto &[next, weight]: adj[curr]) {
            topLevels[next] = std::max(topLevels[next], topLevels[curr] + 1);
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
            assert(localInDegree[next] >= 0 && "Graph is cyclic");
        }
    }

    return topLevels;
}

// Compute top levels using pre-computed topological order
std::vector<uint64_t> Graph::computeTopLevels(const std::vector<uint64_t> &topologicalOrder) const {
    std::vector<uint64_t> topLevels(size, 0);

    // Process nodes in topological order
    for (uint64_t u: topologicalOrder) {
        assert(u < size && "node id in topological order must be smaller than graph size");
        for (const auto &[v, weight]: adj[u]) {
            topLevels[v] = std::max(topLevels[v], topLevels[u] + 1);
        }
    }

    return topLevels;
}

// Compute both topological sort and top levels in a single pass
std::pair<std::vector<uint64_t>, std::vector<uint64_t>> Graph::topologicalSortAndTopLevels() const {
    std::vector<uint64_t> localInDegree = inDegree;
    std::vector<uint64_t> topologicalOrder;
    std::vector<uint64_t> topLevels(size, 0);
    topologicalOrder.reserve(size);

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
        for (const auto &[next, weight]: adj[curr]) {
            topLevels[next] = std::max(topLevels[next], topLevels[curr] + 1);
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
            assert(localInDegree[next] >= 0 && "Graph is cyclic");
        }
    }

    assert(topologicalOrder.size() == size && "Graph is cyclic");
    return {topologicalOrder, topLevels};
}

// Compute shortest paths using Dijkstra's algorithm
std::vector<uint64_t> Graph::distancesFromNode(uint64_t startNode, bool reverseGraph) const {
    assert(startNode < size && "start node id must be smaller than graph size");
    std::vector<uint64_t> distances(size, UINT64_MAX);

    // Priority queue for Dijkstra's algorithm
    std::priority_queue<
            std::pair<uint64_t, uint64_t>,
            std::vector<std::pair<uint64_t, uint64_t>>,
            std::greater<>
    > pq;

    distances[startNode] = 0;
    pq.emplace(startNode, 0);

    while (!pq.empty()) {
        uint64_t currentNode = pq.top().first;
        uint64_t currentDist = pq.top().second;
        pq.pop();

        // Process edges in forward or reverse direction
        for (const auto &edge: (reverseGraph ? revAdj[currentNode] : adj[currentNode])) {
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

std::vector<uint64_t> Graph::groupedTopSortPositions(const std::vector<uint64_t> &partitionInfo) const {
    assert(partitionInfo.size() == size && "Partitioning vector size must match graph size");

    // Step 1: Create a graph of partitions
    // Find the number of partitions
    uint64_t numPartitions = 0;
    for (uint64_t p: partitionInfo) {
        numPartitions = std::max(numPartitions, p + 1);
    }

    // Create a condensed graph where nodes are partitions
    std::vector<std::vector<uint64_t>> coarseGraph(numPartitions);
    std::vector<uint64_t> coarseInDegree(numPartitions, 0);

    // For each edge in the original graph, add an edge between partitions if needed
    for (uint64_t from = 0; from < size; ++from) {
        for (const auto &[to, _]: adj[from]) {
            uint64_t fromPartition = partitionInfo[from];
            uint64_t toPartition = partitionInfo[to];

            // If edge crosses partition boundaries, add it to partition graph
            if (fromPartition != toPartition) {
                // Check if this edge already exists in the partition graph
                bool edgeExists = false;
                for (uint64_t existingTo: coarseGraph[fromPartition]) {
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
    std::vector<uint64_t> coarseTopologicalOrder;
    coarseTopologicalOrder.reserve(numPartitions);


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
        coarseTopologicalOrder.push_back(curr);

        // Update in-degrees and add new nodes with zero in-degree
        for (uint64_t next: coarseGraph[curr]) {
            coarseInDegree[next]--;
            if (coarseInDegree[next] == 0) {
                q.push(next);
            }
        }
    }

    assert(coarseTopologicalOrder.size() == numPartitions && "Coarse graph is cyclic");

    // Step 3: Calculate the position of each partition in the topological order
    std::vector<uint64_t> coarseTopSortPositions(numPartitions);
    for (uint64_t i = 0; i < numPartitions; ++i) {
        coarseTopSortPositions[coarseTopologicalOrder[i]] = i;
    }

    // Step 4: Assign the same topological position to all nodes in the same partition
    std::vector<uint64_t> result(size);
    for (uint64_t i = 0; i < size; ++i) {
        uint64_t partition = partitionInfo[i];
        result[i] = coarseTopSortPositions[partition];
    }

    return result;
}

// Print graph information to output stream
void Graph::print(std::ostream &os) const {
    assert(adj.size() == nodeWeights.size() &&
           "Adjacency-list 2D vectors' size and nodeWeights vectors' size must be equal");
    assert(size == adj.size() && "Adjacency-list 2D vectors' size must be the same as the graph size");

    // Print each node's information
    for (uint64_t i = 0; i < size; ++i) {
        const auto &adjList = adj[i];
        uint64_t weight = nodeWeights[i];
        os << "Node: " << i << "\nWeight: " << weight << "\n";
        os << "Neighbors: ";
        for (const auto &[neighborId, edgeWeight]: adjList) {
            os << neighborId << "(" << edgeWeight << "), ";
        }
        os << "\n\n";
    }
}

// Export graph to DOT format file
void Graph::printToDot(const std::string &dotFilename) const {
    std::ofstream dotFile(dotFilename);

    // Write graph size as comment
    dotFile << "// size=" << size << "\n";
    dotFile << "digraph cfg {\n";

    // Write node weights
    for (uint64_t i = 0; i < size; ++i) {
        dotFile << i << "[weight=" << nodeWeights[i] << "];\n";
    }

    // Write edges and their weights
    for (uint64_t i = 0; i < size; ++i) {
        for (const auto &[neighborId, edgeWeight]: adj[i]) {
            dotFile << i << "->" << neighborId << "[weight=" << edgeWeight << "];\n";
        }
    }
    dotFile << "}\n";
}

// Read graph from DOT format file
Graph readDotFile(const std::string &dotFilename, const std::string &mappingFilename) {
    std::ifstream dotFile(dotFilename);
    assert(dotFile.is_open() && "DOT file could not be opened!");
    std::string line;
    std::regex sizeRegex(R"(//\s*size\s*=\s*(\d+))");

    // Extract graph size from comment
    uint64_t graphSize = UINT64_MAX;
    while (std::getline(dotFile, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, sizeRegex)) {
            graphSize = std::stoi(matches[1]);
            break;
        }
    }

    assert(graphSize != UINT64_MAX && "Failed to extract graph size");

    Graph g(graphSize);

    // Parse nodes and edges
    std::unordered_map<std::string, uint64_t> nodeMap;
    std::regex nodeRegex(R"(([a-zA-Z0-9_]+)\[weight=(\d+)\];)");
    std::regex edgeRegex(R"(([a-zA-Z0-9_]+)->([a-zA-Z0-9_]+)\[weight=(\d+)\];)");
    uint64_t nodeId = 0;

    while (std::getline(dotFile, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, edgeRegex)) {
            // Parse edge
            uint64_t from = nodeMap[matches[1].str()];
            uint64_t to = nodeMap[matches[2].str()];
            uint64_t nodeWeight = std::stoull(matches[3].str());
            g.addEdge(from, to, nodeWeight);
        } else if (std::regex_search(line, matches, nodeRegex)) {
            // Parse node
            uint64_t edgeWeight = std::stoull(matches[2].str());
            g.addNode(nodeId, edgeWeight);
            nodeMap[matches[1].str()] = nodeId;
            nodeId++;
        }
    }

    assert(g.adj.size() == graphSize && "given graph size must match the computed one");

    // Write node mapping
    std::ofstream mapFile(mappingFilename);
    for (const auto &[name, id]: nodeMap) {
        mapFile << name << " -> " << id << "\n";
    }

    return g;
}