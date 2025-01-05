/**
 * @file Graph.h
 * @brief Header file for the Graph class representing a weighted directed graph (DAG)
 */

#ifndef DAG_PARTITIONING_GRAPH_H
#define DAG_PARTITIONING_GRAPH_H

#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include "llvm/Support/raw_ostream.h"

/**
 * @class Graph
 * @brief Represents a weighted directed graph (DAG) with node weights and edge weights
 *
 * This class provides functionality for creating and manipulating directed acyclic graphs
 * with weighted nodes and edges. It supports operations like cycle detection, topological
 * sorting, and distance calculations.
 */
class Graph {
public:
    uint64_t size;                                                    // Number of nodes in the graph
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> adj;     // Forward adjacency list (node_id, edge_weight)
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> revAdj;  // Reverse adjacency list for backward traversal
    std::vector<uint64_t> nodeWeights;                               // Weights of nodes
    uint64_t totalWeight;                                            // Sum of all node weights
    std::vector<uint64_t> inDegree;                                 // Number of incoming edges for each node

    /**
     * @brief Constructs a graph with the specified number of nodes
     * @param size Number of nodes in the graph
     */
    explicit Graph(uint64_t size);

    /**
     * @brief Adds a node with specified weight to the graph
     * @param id Node identifier
     * @param weight Weight of the node
     */
    void addNode(uint64_t id, uint64_t weight);

    /**
     * @brief Adds a weighted edge between two nodes
     * @param from Source node ID
     * @param to Destination node ID
     * @param weight Edge weight
     */
    void addEdge(uint64_t from, uint64_t to, uint64_t weight);

    /**
     * @brief Gets both incoming and outgoing neighbors of a node
     * @param node Node ID to get neighbors for
     * @return Vector of tuples containing (neighbor_id, edge_weight, is_outgoing)
     */
    [[nodiscard]] std::vector<std::tuple<uint64_t, uint64_t, bool>> getNeighbors(uint64_t node) const;

    /**
     * @brief Gets neighbors sorted by edge weight in ascending order
     * @param node Node ID to get sorted neighbors for
     * @return Sorted vector of tuples (neighbor_id, edge_weight, is_outgoing)
     */
    [[nodiscard]] std::vector<std::tuple<uint64_t, uint64_t, bool>>
    getNeighborsSortedByEdgeWeightAsc(uint64_t node) const;

    /**
     * @brief Checks if the graph contains any cycles
     * @return true if graph has cycles, false otherwise
     */
    [[nodiscard]] bool hasCycle() const;

    /**
     * @brief Helper function for cycle detection using iterative DFS
     * @param start Starting node for DFS
     * @return true if cycle is found from start node, false otherwise
     */
    [[nodiscard]] bool iterativeDfsHasCycle(uint64_t start) const;

    /**
     * @brief Performs topological sorting of the graph
     * @return Vector of node IDs in topological order
     */
    [[nodiscard]] std::vector<uint64_t> topologicalSort() const;

    /**
     * @brief Computes the top level (longest path length from any root) for each node
     * @return Vector of top levels for each node
     */
    [[nodiscard]] std::vector<uint64_t> computeTopLevels() const;

    /**
     * @brief Computes top levels using a pre-computed topological order
     * @param topologicalOrder Vector containing nodes in topological order
     * @return Vector of top levels for each node
     */
    [[nodiscard]] std::vector<uint64_t> computeTopLevels(const std::vector<uint64_t> &topologicalOrder) const;

    /**
     * @brief Computes both topological sort and top levels in a single pass
     * @return Pair of vectors containing topological order and top levels
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, std::vector<uint64_t>> topologicalSortAndTopLevels() const;

    /**
     * @brief Computes shortest distances from a start node using Dijkstra's algorithm
     * @param startNode Starting node for distance calculation
     * @param reverseGraph If true, computes distances using reverse edges
     * @return Vector of distances from start node to all other nodes
     */
    [[nodiscard]] std::vector<uint64_t> distancesFromNode(uint64_t startNode, bool reverseGraph = false) const;

    /**
     * @brief Gets the maximum node weight in the graph
     * @return Maximum node weight
     */
    [[nodiscard]] uint64_t maxNodeWeight() const;

    /**
     * @brief Prints graph information to the specified output stream
     * @param os Output stream to print to
     */
    void print(llvm::raw_ostream &os) const;

    /**
     * @brief Exports the graph to a DOT file format
     * @param dotFilename Name of the output DOT file
     */
    void printToDot(const std::string &dotFilename) const;
};

/**
 * @brief Reads a graph from a DOT file and creates node ID mapping
 * @param dotFilename Input DOT file path
 * @param mappingFilename Output file for node name to ID mapping
 * @return Constructed Graph object
 */
Graph readDotFile(const std::string &dotFilename, const std::string &mappingFilename);

#endif //DAG_PARTITIONING_GRAPH_H