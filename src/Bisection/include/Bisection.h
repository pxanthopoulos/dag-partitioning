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

class Bisection {
protected:
    const Graph &workingGraph;                  // Graph to be bisected
    double upperBoundPartWeight;         // Maximum allowed weight for each partition
    double lowerBoundPartWeight;         // Minimum allowed weight for each partition

    /**
     * @brief Verifies that a bisection maintains acyclicity
     *
     * Checks that there are no edges from V1 (1) to V0 (0).
     * A valid bisection must not have any backward edges from the second
     * partition to the first.
     *
     * @param bisection Vector where true indicates V1, false indicates V0
     * @return true if bisection is acyclic
     */
    [[nodiscard]] virtual bool checkValidBisection(const std::vector<uint8_t> &bisection) const;

    /**
     * @brief Computes total weight of edges crossing the bisection
     *
     * Sums up weights of all edges that connect vertices in different
     * partitions (edge cut).
     *
     * @param bisection Vector where true indicates V1, false indicates V0
     * @return Total weight of edges crossing between partitions
     */
    [[nodiscard]] virtual uint64_t computeEdgeCut(const std::vector<uint8_t> &bisection) const;

    /**
     * @brief Protected constructor for derived classes
     * @param graph Graph to be bisected
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum required partition weight
     */
    Bisection(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight);

public:
    /**
     * @brief Virtual destructor for proper cleanup in derived classes
     */
    virtual ~Bisection() = default;

    /**
     * @brief Pure virtual method to execute the bisection algorithm
     * @return Pair containing:
     *         - Vector indicating partition assignment (0=V0, 1=V1)
     *         - Total edge cut weight
     */
    [[nodiscard]] virtual std::pair<std::vector<uint8_t>, uint64_t> run() const = 0;
};

#endif //DAG_PARTITIONING_BISECTION_H