/**
 * @file MultilevelBisectioner.h
 * @brief Implementation of multilevel directed graph partitioning
 *
 * Orchestrates the three main phases of multilevel partitioning:
 * 1. Coarsening: Cluster graph while maintaining acyclicity
 * 2. Initial Partitioning: Bisect coarsest graph
 * 3. Uncoarsening + Refinement: Project and refine through all levels
 */

#ifndef DAG_PARTITIONING_MULTILEVELBISECTIONER_H
#define DAG_PARTITIONING_MULTILEVELBISECTIONER_H

#include <cstdint>
#include <stack>
#include <vector>

class Graph;
enum class RefinementMethod;

/**
 * @brief Available clustering methods for coarsening phase
 */
enum class ClusteringMethod {
    FORB, // Forbidden edges (Section 4.1.1)
    CYC,  // Cycle detection (Section 4.1.2)
    HYB   // Hybrid approach (Section 4.1.3)
};

/**
 * @brief Available methods for initial bisection
 */
enum class BisectionMethod {
    GGG,         // Greedy directed graph growing (Section 4.2.1)
    UNDIRMETIS,  // Undirected METIS + fix acyclicity (Section 4.2.2)
    UNDIRSCOTCH, // Undirected Scotch + fix acyclicity (Section 4.2.2)
    UNDIRBOTH    // Try both METIS and Scotch
};

class MultilevelBisectioner {
  private:
    const Graph &workingGraph;         // Graph to be partitioned
    ClusteringMethod clusteringMethod; // Selected clustering strategy
    uint64_t maxClusteringRounds;      // Maximum coarsening levels
    uint64_t minClusteringVertices;    // Stop coarsening at this size
    double clusteringVertexRatio; // If, after a clustering round, this ratio is
                                  // not surpassed, stop clustering
    BisectionMethod bisectionMethod;   // Selected bisection strategy
    double imbalanceRatio;             // Maximum allowed partition imbalance
    RefinementMethod refinementMethod; // Selected refinement strategy
    uint64_t refinementPasses;         // Maximum refinement passes per level

  public:
    /**
     * @brief Constructs multilevel bisectioner with selected strategies
     * @param graph Graph to be partitioned
     * @param clusteringMethod Coarsening strategy
     * @param maxClusteringRounds Maximum coarsening levels
     * @param minClusteringVertices Minimum coarse graph size
     * @param clusteringVertexRatio If, after a clustering round, this ratio is
     * not surpassed, stop clustering
     * @param bisectionMethod Initial bisection strategy
     * @param imbalanceRatio Maximum allowed imbalance
     * @param refinementMethod Refinement strategy
     * @param refinementPasses Maximum refinement passes per level
     */
    MultilevelBisectioner(const Graph &graph, ClusteringMethod clusteringMethod,
                          uint64_t maxClusteringRounds,
                          uint64_t minClusteringVertices,
                          double clusteringVertexRatio,
                          BisectionMethod bisectionMethod,
                          double imbalanceRatio,
                          RefinementMethod refinementMethod,
                          uint64_t refinementPasses);

    /**
     * @brief Executes coarsening phase
     * @return Stack of intermediate graphs and mappings
     */
    [[nodiscard]] std::stack<std::pair<Graph, std::vector<uint64_t>>>
    runClustering() const;

    /**
     * @brief Computes initial bisection of coarsest graph
     * @param graph Graph to bisect
     * @return Pair of bisection vector and edge cut weight
     */
    [[nodiscard]] std::pair<std::vector<uint8_t>, uint64_t>
    runBisection(const Graph &graph) const;

    /**
     * @brief Projects bisection to finer level using node mapping
     * @param bisectionInfoPair Current bisection to be projected
     * @param mapping Node mapping from finer to coarser level (clustering)
     */
    static void projectBisection(
        std::pair<std::vector<uint8_t>, uint64_t> &bisectionInfoPair,
        const std::vector<uint64_t> &mapping);

    /**
     * @brief Refines bisection at current level
     * @param graph Current level graph
     * @param bisectionInfoPair Bisection to be refined
     */
    void runRefinement(
        const Graph &graph,
        std::pair<std::vector<uint8_t>, uint64_t> &bisectionInfoPair) const;

    /**
     * @brief Executes complete multilevel bisection process
     * @return Final bisection and edge cut weight
     */
    [[nodiscard]] std::pair<std::vector<uint8_t>, uint64_t> run() const;
};

#endif // DAG_PARTITIONING_MULTILEVELBISECTIONER_H