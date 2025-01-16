/**
 * @file Refinement.h
 * @brief Base class for partition refinement algorithms
 *
 * Provides framework for iteratively improving an initial partition while
 * maintaining acyclicity and balance constraints. Derived classes implement
 * specific refinement strategies.
 */

#ifndef DAG_PARTITIONING_REFINEMENT_H
#define DAG_PARTITIONING_REFINEMENT_H

#include "Graph.h"

class Refinement {
protected:
    const Graph &workingGraph;               // Graph being partitioned
    std::vector<uint8_t> &initialBisectionInfo; // Current partition assignment (modified by refinement)
    uint64_t &initialEdgeCut;               // Current edge cut weight (modified by refinement)
    uint64_t maxNumberOfPasses;             // Maximum refinement iterations
    double upperBoundPartWeight;            // Maximum allowed partition weight
    double lowerBoundPartWeight;            // Minimum allowed partition weight

    /**
     * @brief Checks if current bisection maintains acyclicity
     *
     * Verifies that there are no edges from V1 (1) to V0 (0).
     * Must be maintained throughout refinement.
     *
     * @return true if bisection is acyclic
     */
    [[nodiscard]] bool checkValidBisection() const;

    /**
     * @brief Computes current weight of both partitions
     * @return Pair containing (V0 weight, V1 weight)
     */
    [[nodiscard]] std::pair<uint64_t, uint64_t> calculatePartSizes() const;

    /**
     * @brief Checks if partition weights are within balance constraints
     *
     * Allows for some imbalance up to the weight of the heaviest node,
     * as moving even one node could be necessary to maintain acyclicity.
     *
     * @param maxNodeWeight Weight of heaviest node in graph
     * @return true if partition weights are valid
     */
    [[nodiscard]] bool checkBalance(uint64_t maxNodeWeight) const;

    /**
     * @brief Checks if the computed edge cut is consistent with the bisection info
     *
     * @return true if edge cut is correct
     */
    [[nodiscard]] bool checkValidEdgeCut();

    /**
     * @brief Pure virtual method for one refinement pass
     *
     * Derived classes implement specific refinement strategy here.
     * Should maintain acyclicity and balance constraints.
     *
     * @return true if refinement pass made any changes at all
     */
    virtual bool onePassRefinement() = 0;

    /**
     * @brief Protected constructor for derived classes
     * @param graph Graph to be refined
     * @param initialBisectionInfo Current partition assignments
     * @param initialEdgeCut Current edge cut weight
     * @param maxNumberOfPasses Maximum refinement iterations
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum allowed partition weight
     */
    Refinement(const Graph &graph, std::vector<uint8_t> &initialBisectionInfo, uint64_t &initialEdgeCut,
               uint64_t maxNumberOfPasses, double upperBoundPartWeight, double lowerBoundPartWeight);

public:
    /**
     * @brief Virtual destructor for proper cleanup
     */
    virtual ~Refinement() = default;

    /**
     * @brief Executes refinement process
     *
     * Iteratively calls onePassRefinement() until either:
     * 1. Maximum passes reached
     * 2. No further improvement possible
     * Maintains invariants throughout process.
     */
    void run();
};

#endif //DAG_PARTITIONING_REFINEMENT_H