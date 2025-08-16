/**
 * @file ClusteringCycleDetection.h
 * @brief Implementation of acyclic clustering with cycle detection as described
 * in Section 4.1.2
 *
 * This class implements a less restrictive clustering algorithm than the
 * forbidden edges approach. It maintains the top-level difference constraint
 * but uses explicit cycle detection instead of forbidden edge rules to ensure
 * acyclicity.
 */

#ifndef DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H
#define DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H

#include "Clustering.h"

namespace dag_partitioning {

namespace clustering {

class ClusteringCycleDetection : public Clustering {
  private:
    /**
     * @brief Performs a thorough cycle check on the coarsened graph
     *
     * This is a verification function that constructs the coarsened graph and
     * checks for cycles using full graph traversal. Used as a safety check to
     * ensure the faster cycle detection didn't miss anything.
     *
     * @param leaders Mapping of nodes to their cluster leaders
     * @param newSize Number of clusters in the coarsened graph
     */
    void hardCheckCycle(const std::vector<uint64_t> &leaders,
                        uint64_t newSize) const;

    /**
     * @brief Finds the minimum top-level value among all nodes in a cluster
     *
     * Used by the cycle detection algorithm as described in the paper:
     * "Let t denote the minimum top level in C"
     *
     * @param node A node in the cluster
     * @param topLevels Top-level values for all nodes
     * @param leaders Cluster assignments
     * @param markup Vector that signifies if a clustered node has the small top
     * value in its cluster
     * @param markdown Vector that signifies if a clustered node has the high
     * top value in its cluster
     * @return Minimum top-level value in the cluster
     */
    [[nodiscard]] static uint64_t
    findMinimumTopLevelInCluster(uint64_t node,
                                 const std::vector<uint64_t> &topLevels,
                                 const std::vector<uint8_t> &markup,
                                 const std::vector<uint8_t> &markdown);

    /**
     * @brief Implementation of one round of clustering using cycle detection
     *
     * Implements Algorithm 2 from the paper, processing nodes in topological
     * order and using detectCycle() to check if merging nodes would create
     * cycles.
     *
     * @return Pair containing:
     *         - Vector mapping nodes to cluster leaders
     *         - Number of clusters after this round
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t>
    oneRoundClustering() const override;

  protected:
    /**
     * @brief Checks if merging nodes `from` and `to` (directed graph) would
     * create a cycle
     *
     * Implements the cycle detection algorithm described in the paper:
     * Only needs to check nodes with top-level values t or t+1 where t is the
     * minimum top-level in the target cluster. Uses BFS traversal.
     *
     * @param from Node being considered for addition to cluster
     * @param to Node already in the target cluster
     * @param topLevels Top-level values for all nodes
     * @param leaders Current cluster assignments
     * @param markup Vector that signifies if a clustered node has the small top
     * value in its cluster
     * @param markdown Vector that signifies if a clustered node has the high
     * top value in its cluster
     * @return true if adding the node would create a cycle
     */
    [[nodiscard]] bool detectCycle(uint64_t from, uint64_t to,
                                   const std::vector<uint64_t> &topLevels,
                                   const std::vector<uint64_t> &leaders,
                                   const std::vector<uint8_t> &markup,
                                   const std::vector<uint8_t> &markdown) const;

  public:
    /**
     * @brief Constructs the clustering algorithm instance
     * @param graph Original graph to be clustered
     * @param maxRounds Maximum number of clustering rounds
     * @param minVertices Minimum number of vertices to stop at
     * @param vertexRatio If, after a clustering round, this ratio is not
     * surpassed, stop clustering
     */
    ClusteringCycleDetection(const core::Graph &graph, uint64_t maxRounds,
                             uint64_t minVertices, double vertexRatio);
};

} // namespace clustering

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H