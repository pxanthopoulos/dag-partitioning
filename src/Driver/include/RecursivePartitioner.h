/**
 * @file RecursivePartitioner.h
 * @brief Recursive k-way partitioning implementation using recursive bisection
 *
 * Implements k-way partitioning by recursively applying bisection to create
 * hierarchical partitions. Each bisection uses multilevel methods while
 * maintaining acyclicity constraints throughout the process.
 */

#ifndef DAG_PARTITIONING_RECURSIVEPARTITIONER_H
#define DAG_PARTITIONING_RECURSIVEPARTITIONER_H

#include "robin_hood.h"
#include <cstdint>
#include <string>
#include <vector>

namespace dag_partitioning {

namespace core {
class Graph;
}
namespace refinement {
enum class RefinementMethod;
}

namespace clustering {
enum class ClusteringMethod;
}

namespace bisection {
enum class BisectionMethod;
}

namespace driver {

class RecursivePartitioner {
  private:
    const core::Graph &workingGraph; // Graph to be partitioned
    uint64_t partitions;             // Number of partitions to create
    clustering::ClusteringMethod
        clusteringMethod;           // Selected clustering strategy
    uint64_t maxClusteringRounds;   // Maximum coarsening levels
    uint64_t minClusteringVertices; // Stop coarsening at this size
    double clusteringVertexRatio; // If, after a clustering round, this ratio is
                                  // not surpassed, stop clustering
    bisection::BisectionMethod bisectionMethod; // Selected bisection strategy
    double imbalanceRatio; // Maximum allowed partition imbalance
    refinement::RefinementMethod
        refinementMethod;      // Selected refinement strategy
    uint64_t refinementPasses; // Maximum refinement passes per level

    // Parallelization options
    bool enableParallel;         // Enable parallel recursive calls
    uint64_t minSizeForParallel; // Minimum subgraph size to parallelize
    uint64_t maxParallelDepth;   // Maximum recursion depth for parallelization

    /**
     * @brief Creates two subgraphs by splitting the input graph based on a
     * bisection
     *
     * Takes a partition vector and splits the graph into two subgraphs while
     * maintaining internal edges and removing cross-partition edges. Also,
     * returns one map for each subgraph that maps new node IDs to old node IDs
     *
     * @param bisection Vector where bisection[i] indicates which partition
     * (0=V0, 1=V1) vertex i belongs to
     * @return The two resulting subgraphs (V0, V1) and the two maps
     */
    [[nodiscard]] std::tuple<
        core::Graph, robin_hood::unordered_map<uint64_t, uint64_t>, core::Graph,
        robin_hood::unordered_map<uint64_t, uint64_t>>
    createSubgraphs(const std::vector<uint8_t> &bisection) const;

    /**
     * @brief Constructs recursive partitioner with selected strategies
     * @param graph Graph to be partitioned
     * @param partitions Number of partitions to create
     * @param clusteringMethod Coarsening strategy
     * @param maxClusteringRounds Maximum coarsening levels
     * @param minClusteringVertices Minimum coarse graph size
     * @param clusteringVertexRatio If, after a clustering round, this ratio is
     * not surpassed, stop clustering
     * @param bisectionMethod Initial bisection strategy
     * @param imbalanceRatio Maximum allowed imbalance
     * @param refinementMethod Refinement strategy
     * @param refinementPasses Maximum refinement passes per level
     * @param enableParallel Enable parallel recursive calls
     * @param minSizeForParallel Minimum subgraph size to parallelize
     * @param maxParallelDepth Maximum recursion depth for parallelization
     */
    RecursivePartitioner(const core::Graph &graph, uint64_t partitions,
                         clustering::ClusteringMethod clusteringMethod,
                         uint64_t maxClusteringRounds,
                         uint64_t minClusteringVertices,
                         double clusteringVertexRatio,
                         bisection::BisectionMethod bisectionMethod,
                         double imbalanceRatio,
                         refinement::RefinementMethod refinementMethod,
                         uint64_t refinementPasses, bool enableParallel = true,
                         uint64_t minSizeForParallel = 100,
                         uint64_t maxParallelDepth = 10);

  public:
    /**
     * @brief Public constructor overload with strings
     * @param graph Graph to be partitioned
     * @param partitions Number of partitions to create
     * @param clusteringMethod Coarsening strategy
     * @param maxClusteringRounds Maximum coarsening levels
     * @param minClusteringVertices Minimum coarse graph size
     * @param clusteringVertexRatio If, after a clustering round, this ratio is
     * not surpassed, stop clustering
     * @param bisectionMethod Initial bisection strategy
     * @param imbalanceRatio Maximum allowed imbalance
     * @param refinementMethod Refinement strategy
     * @param refinementPasses Maximum refinement passes per level
     * @param enableParallel Enable parallel recursive calls
     * @param minSizeForParallel Minimum subgraph size to parallelize
     * @param maxParallelDepth Maximum recursion depth for parallelization
     */
    RecursivePartitioner(const core::Graph &graph, uint64_t partitions,
                         std::string clusteringMethod,
                         uint64_t maxClusteringRounds,
                         uint64_t minClusteringVertices,
                         double clusteringVertexRatio,
                         std::string bisectionMethod, double imbalanceRatio,
                         std::string refinementMethod,
                         uint64_t refinementPasses, bool enableParallel = true,
                         uint64_t minSizeForParallel = 100,
                         uint64_t maxParallelDepth = 10);

    /**
     * @brief Executes recursive partitioning process
     *
     * Recursively divides the graph into k parts by:
     * 1. Base cases:
     *    - k=1: Return all vertices in partition 0
     *    - k=2: Use multilevel bisection directly
     * 2. Recursive case:
     *    - Bisect graph into V0, V1
     *    - Split k parts between subgraphs
     *    - Recursively partition each subgraph
     *    - Combine results with appropriate offsets
     *
     * @param currentDepth Current recursion depth (default 0 for public calls)
     * @return Pair containing:
     *         - Vector mapping each vertex to its partition (0 to k-1)
     *         - Total edge cut weight across all partitions
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t>
    run(uint64_t currentDepth = 0) const;
};

} // namespace driver

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_RECURSIVEPARTITIONER_H
