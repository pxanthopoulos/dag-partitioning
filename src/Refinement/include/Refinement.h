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

#include <cstdint>
#include <vector>

namespace dag_partitioning {

namespace core {
class Graph;
}

namespace refinement {

class Refinement {
  protected:
    const core::Graph &workingGraph;            // Graph being partitioned
    std::vector<uint8_t> &initialBisectionInfo; // Current partition assignment
                                                // (modified by refinement)
    uint64_t
        &initialEdgeCut; // Current edge cut weight (modified by refinement)
    uint64_t maxNumberOfPasses;  // Maximum refinement iterations
    double upperBoundPartWeight; // Maximum allowed partition weight
    double lowerBoundPartWeight; // Minimum allowed partition weight

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
    Refinement(const core::Graph &graph,
               std::vector<uint8_t> &initialBisectionInfo,
               uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
               double upperBoundPartWeight, double lowerBoundPartWeight);

  public:
    /**
     * @brief Virtual destructor for proper cleanup
     */
    virtual ~Refinement() = default;

    /**
     * @brief Checks if current bisection maintains acyclicity
     *
     * Verifies that there are no edges from V1 (1) to V0 (0).
     * Must be maintained throughout refinement.
     *
     * @param bisectionInfo Bisection information
     * @param graph Corresponding graph
     * @return true if bisection is acyclic
     */
    [[nodiscard]] static bool
    checkValidBisection(const std::vector<uint8_t> &bisectionInfo,
                        const core::Graph &graph);

    /**
     * @brief Computes current weight of both partitions
     * @param bisectionInfo Bisection information
     * @param graph Corresponding graph
     * @return Pair containing (V0 weight, V1 weight)
     */
    [[nodiscard]] static std::pair<uint64_t, uint64_t>
    calculatePartSizes(const std::vector<uint8_t> &bisectionInfo,
                       const core::Graph &graph);

    /**
     * @brief Checks if partition weights are within balance constraints
     *
     * Allows for some imbalance up to the weight of the heaviest node,
     * as moving even one node could be necessary to maintain acyclicity.
     *
     * @param bisectionInfo Bisection information
     * @param graph Corresponding graph
     * @param maxNodeWeight Weight of heaviest node in graph
     * @param upperBoundPartWeight Upper bound for partition total weight
     * @param lowerBoundPartWeight Lower bound for partition total weight
     * @return true if partition weights are valid
     */
    [[nodiscard]] static bool
    checkBalance(const std::vector<uint8_t> &bisectionInfo,
                 const core::Graph &graph, uint64_t maxNodeWeight,
                 double upperBoundPartWeight, double lowerBoundPartWeight);

    /**
     * @brief Checks if partition weights are within balance constraints
     *
     * Allows for some imbalance up to the weight of the heaviest node,
     * as moving even one node could be necessary to maintain acyclicity.
     *
     * @param sizeV0 Size of V0 part
     * @param sizeV1 Size of V1 part
     * @param maxNodeWeight Weight of heaviest node in graph
     * @param upperBoundPartWeight Upper bound for partition total weight
     * @param lowerBoundPartWeight Lower bound for partition total weight
     * @return true if partition weights are valid
     */
    [[nodiscard]] static bool checkBalance(uint64_t sizeV0, uint64_t sizeV1,
                                           uint64_t maxNodeWeight,
                                           double upperBoundPartWeight,
                                           double lowerBoundPartWeight);

    /**
     * @brief Checks if the computed edge cut is consistent with the bisection
     * info
     *
     * @param bisectionInfo Bisection information
     * @param graph Corresponding graph
     * @param currentEdgeCut Current edge cut
     * @return true if edge cut is correct
     */
    [[nodiscard]] static bool
    checkValidEdgeCut(const std::vector<uint8_t> &bisectionInfo,
                      const core::Graph &graph, uint64_t currentEdgeCut);

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

/**
 * @brief Available refinement methods for refinement phase
 */
enum class RefinementMethod {
    BOUNDARYFM, // Boundary FM adaptation (Section 4.3)
    BOUNDARYKL, // Boundary KL adaptation
    MIXED       // Mixed (1 pass of KL followed by one pass of FM)
};

} // namespace refinement

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_REFINEMENT_H