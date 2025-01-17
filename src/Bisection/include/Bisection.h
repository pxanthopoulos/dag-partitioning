/**
 * @file Bisection.h
 * @brief Base class for graph bisection algorithms with acyclicity constraints
 *
 * This class provides the framework for implementing different graph bisection
 * strategies for the initial partitioning phase. Bisections must maintain the
 * acyclicity property: no edges from V1 to V0 are allowed.
 */

#ifndef DAG_PARTITIONING_BISECTION_H
#define DAG_PARTITIONING_BISECTION_H

#include "Graph.h"
#include "Refinement.h"

class Bisection {
protected:
    const Graph &workingGraph;                  // Graph to be bisected
    double upperBoundPartWeight;         // Maximum allowed weight for each partition
    double lowerBoundPartWeight;         // Minimum allowed weight for each partition
    RefinementMethod refinementMethod;   // Refinement method
    uint64_t refinementPasses;           // Number of refinement passes

    /**
     * @brief Protected constructor for derived classes
     * @param graph Graph to be bisected
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum required partition weight
     * @param refinementMethod Refinement method
     * @param refinementPasses Number of refinement passes
     */
    Bisection(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight,
              RefinementMethod refinementMethod, uint64_t refinementPasses);

public:
    /**
     * @brief Virtual destructor for proper cleanup in derived classes
     */
    virtual ~Bisection() = default;

    /**
     * @brief Computes total weight of edges crossing the bisection
     *
     * Sums up weights of all edges that connect vertices in different
     * partitions (edge cut).
     *
     * @param bisection Vector where 1 indicates V1, 0 indicates V0
     * @param graph The corresponding graph
     * @return Total weight of edges crossing between partitions
     */
    [[nodiscard]] static uint64_t computeEdgeCut(const std::vector<uint8_t> &bisection, const Graph &graph);

    /**
     * @brief Pure virtual method to execute the bisection algorithm
     * @return Pair containing:
     *         - Vector indicating partition assignment (0=V0, 1=V1)
     *         - Total edge cut weight
     */
    [[nodiscard]] virtual std::pair<std::vector<uint8_t>, uint64_t> run() const = 0;
};

#endif //DAG_PARTITIONING_BISECTION_H