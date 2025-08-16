/**
 * @file Clustering.h
 * @brief Abstract base class for graph clustering algorithms
 *
 * This class provides the framework for implementing different graph clustering
 * algorithms. It maintains the working state of the graph during clustering
 * and tracks intermediate results.
 */

#ifndef DAG_PARTITIONING_CLUSTERING_H
#define DAG_PARTITIONING_CLUSTERING_H

#include "Graph.h"

#include <stack>
#include <vector>

class Clustering {
  protected:
    Graph workingGraph; // Current state of the graph being clustered
    std::stack<std::pair<Graph, std::vector<uint64_t>>>
        intermediateGraphsAndClusters; // History of clustering steps
    uint64_t maxRounds;   // Maximum number of clustering rounds to perform
    uint64_t minVertices; // Minimum number of vertices to stop clustering
    double vertexRatio;   // If, after a clustering round, this ratio is not
                          // surpassed, stop clustering

    /**
     * @brief Pure virtual method to implement one round of clustering
     * @return Pair containing:
     *         - Vector mapping each node to its cluster leader
     *         - Number of clusters (new graph size) after this round
     */
    [[nodiscard]] virtual std::pair<std::vector<uint64_t>, uint64_t>
    oneRoundClustering() const = 0;

    /**
     * @brief Updates the working graph based on clustering results
     * @param leaders Vector mapping each node to its cluster leader
     * @param newSize Number of nodes in the new clustered graph
     * @return true if clustering should continue, false if minimum vertices
     * stopping criteria met
     */
    [[nodiscard]] bool
    updateGraphAndClusters(const std::vector<uint64_t> &leaders,
                           uint64_t newSize);

    /**
     * @brief Protected constructor for derived classes
     * @param graph Original graph to be clustered
     * @param maxRounds Maximum number of clustering iterations
     * @param minVertices Minimum number of vertices to stop at
     * @param vertexRatio If, after a clustering round, this ratio is not
     * surpassed, stop clustering
     */
    Clustering(Graph graph, uint64_t maxRounds, uint64_t minVertices,
               double vertexRatio);

  public:
    /**
     * @brief Virtual destructor for proper cleanup in derived classes
     */
    virtual ~Clustering() = default;

    /**
     * @brief Executes the clustering algorithm, while checking for the stopping
     * criteria
     * @return Stack of intermediate clustering results, with each entry
     * containing:
     *         - The graph state at that step
     *         - Vector mapping original nodes to cluster IDs
     */
    [[nodiscard]] std::stack<std::pair<Graph, std::vector<uint64_t>>> run();
};

#endif // DAG_PARTITIONING_CLUSTERING_H