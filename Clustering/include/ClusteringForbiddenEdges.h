/**
 * @file ClusteringForbiddenEdges.h
 * @brief Implementation of acyclic clustering with forbidden edges as described in "Acyclic Partitioning of DAGs"
 *
 * This class implements the clustering algorithm described in Section 4.1.1 of the paper,
 * which uses static information about top-level differences and a concept of "bad neighbors"
 * to ensure the resulting coarsened graph remains acyclic.
 */

#ifndef DAG_PARTITIONING_CLUSTERINGFORBIDDENEDGES_H
#define DAG_PARTITIONING_CLUSTERINGFORBIDDENEDGES_H

#include "Clustering.h"
#include <vector>
#include <cstdint>
#include <utility>

class ClusteringForbiddenEdges : public Clustering {
private:
    /**
     * @brief Finds valid neighboring nodes that can form a cluster with the given node
     *
     * Implements the neighbor selection criteria from Theorem 4.2:
     * 1. Top-level difference must be ≤ 1
     * 2. A node cannot join a cluster if it has more than one "bad neighbor" cluster
     * 3. Total cluster weight must not exceed 10% of graph weight
     *
     * @param node Current node being processed
     * @param neighbors List of node's neighbors with their edge weights and direction
     * @param leaders Current cluster assignments
     * @param clusterWeights Weights of each cluster
     * @param topLevels Top-level values of each node
     * @param numberOfBadNeighbors Number of bad neighbor clusters for each node
     * @param leaderOfBadNeighbors Leader of the bad neighbor cluster for each node
     * @param leadersToMinMaxTopValues Min/max top-level values for each cluster
     * @return Vector of valid neighbor nodes and their edge weights
     */
    [[nodiscard]] std::vector<std::pair<uint64_t, uint64_t>>
    findValidNeighbors(uint64_t node,
                       const std::vector<std::tuple<uint64_t, uint64_t, bool>> &neighbors,
                       const std::vector<uint64_t> &leaders,
                       const std::vector<uint64_t> &clusterWeights,
                       const std::vector<uint64_t> &topLevels,
                       const std::vector<uint64_t> &numberOfBadNeighbors,
                       const std::vector<uint64_t> &leaderOfBadNeighbors,
                       const std::unordered_map<uint64_t, std::pair<uint64_t, uint64_t>> &leadersToMinMaxTopValues) const;

    /**
     * @brief Selects the best neighbor from valid candidates based on edge weight
     *
     * From the list of valid neighbors that satisfy Theorem 4.2's conditions,
     * selects the one with the minimum edge weight as described in Algorithm 1.
     *
     * @param validNeighbors Vector of valid neighbor nodes and their edge weights
     * @return ID of the selected best neighbor
     */
    [[nodiscard]] static uint64_t findBestNeighbor(const std::vector<std::pair<uint64_t, uint64_t>> &);

    /**
     * @brief Implements one round of the clustering algorithm
     *
     * Implementation of Algorithm 1 from the paper. Processes nodes in topological order,
     * maintaining invariants about bad neighbors and top-level differences to ensure
     * the resulting coarsened graph remains acyclic.
     *
     * @return Pair containing:
     *         - Vector mapping each node to its cluster leader
     *         - Number of clusters after this round
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const override;

public:
    /**
     * @brief Constructs the clustering algorithm instance
     * @param graph Original graph to be clustered
     * @param maxRounds Maximum number of clustering rounds
     * @param minVertices Minimum number of vertices to stop at
     */
    ClusteringForbiddenEdges(const Graph &graph, uint64_t maxRounds, uint64_t minVertices);
};

#endif //DAG_PARTITIONING_CLUSTERINGFORBIDDENEDGES_H