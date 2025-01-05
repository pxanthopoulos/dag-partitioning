/**
 * @file ClusteringHybrid.h
 * @brief Implementation of hybrid acyclic clustering as described in Section 4.1.3
 *
 * This class implements a hybrid approach that combines cycle detection and forbidden
 * edges strategies. It uses cycle detection by default but switches to forbidden edges
 * rules for vertices with large degrees to avoid the quadratic runtime behavior.
 */

#ifndef DAG_PARTITIONING_CLUSTERINGHYBRID_H
#define DAG_PARTITIONING_CLUSTERINGHYBRID_H

#include "ClusteringCycleDetection.h"

class ClusteringHybrid : public ClusteringCycleDetection {
private:
    /**
     * @brief Determines if edge endpoints have large degrees
     *
     * Implements the degree threshold check from section 4.1.3:
     * "We define a limit on the degree of a vertex (typically sqrt(|V|)/10)"
     *
     * @param from Source node of edge being considered
     * @param to Target node of edge being considered
     * @return true if either endpoint exceeds degree threshold
     */
    [[nodiscard]] bool checkLargeDegrees(uint64_t from, uint64_t to) const;

    /**
     * @brief Updates forbidden edge tracking data when using forbidden edges strategy
     *
     * This function maintains the bad neighbor tracking data structures described
     * in section 4.1.1, for when the algorithm switches to forbidden edges mode
     * for high-degree vertices.
     *
     * @param topLevels Top-level values for all nodes
     * @param node Current node being processed
     * @param sortedNeighbors Neighbors sorted by edge weight
     * @param neighborId ID of neighbor being merged with
     * @param leaderOfNeighbor Cluster leader of the neighbor
     * @param markup Tracks nodes at level t in their cluster
     * @param markdown Tracks nodes at level t+1 in their cluster
     * @param numberOfBadNeighbors Count of bad neighbor clusters per node
     * @param leaderOfBadNeighbors Leader of first bad neighbor cluster per node
     */
    void bookKeepingForForbiddenEdges(const std::vector<uint64_t> &topLevels,
                                      uint64_t node,
                                      const std::vector<std::tuple<uint64_t, uint64_t, bool>> &sortedNeighbors,
                                      uint64_t neighborId,
                                      uint64_t leaderOfNeighbor,
                                      const std::vector<bool> &markup,
                                      const std::vector<bool> &markdown,
                                      std::vector<uint64_t> &numberOfBadNeighbors,
                                      std::vector<uint64_t> &leaderOfBadNeighbors) const;

    /**
     * @brief Implements one round of hybrid clustering
     *
     * Processes nodes in topological order, using cycle detection for normal vertices
     * and forbidden edges rules for high-degree vertices. Combines data structures
     * from both approaches to maintain consistency.
     *
     * @return Pair containing:
     *         - Vector mapping nodes to cluster leaders
     *         - Number of clusters after this round
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const override;

public:
    /**
     * @brief Constructs the hybrid clustering algorithm instance
     * @param graph Original graph to be clustered
     * @param maxRounds Maximum number of clustering rounds
     * @param minVertices Minimum number of vertices to stop at
     */
    ClusteringHybrid(const Graph &graph, uint64_t maxRounds, uint64_t minVertices);
};

#endif //DAG_PARTITIONING_CLUSTERINGHYBRID_H