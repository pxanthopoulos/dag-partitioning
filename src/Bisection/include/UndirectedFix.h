/**
 * @file UndirectedFix.h
 * @brief Implementation of undirected bisection with acyclicity fixing as described in Section 4.2.2
 *
 * This class implements a two-step approach:
 * 1. Use standard undirected graph partitioners (METIS/Scotch) to get initial bisection
 * 2. Fix acyclicity by propagating partition assignments up or down the graph
 */

#ifndef DAG_PARTITIONING_UNDIRECTEDFIX_H
#define DAG_PARTITIONING_UNDIRECTEDFIX_H

#include "Bisection.h"
#include <stdio.h>
#include "scotch.h"
#include "metis.h"

class UndirectedFix : public Bisection {
private:
    bool useMetis;   // Whether to try METIS partitioning
    bool useScotch;  // Whether to try Scotch partitioning

    /**
     * @brief Computes total number of edges in the graph
     * @return Edge count including both forward and reverse edges
     */
    [[nodiscard]] int64_t computeNumberOfEdges() const;

    /**
     * @brief Converts graph to Compressed Sparse Row format for external partitioners
     *
     * Creates CSR representation treating the graph as undirected by combining
     * forward and reverse edges.
     *
     * @param edgeNumber Total number of edges
     * @param nodeNeighborsOffset CSR row pointers
     * @param nodeNeighbors CSR column indices
     * @param edgeWeights CSR edge weights
     * @param nodeWeights Node weights
     */
    void graphToCSRFormat(int64_t edgeNumber,
                          std::vector<int64_t> &nodeNeighborsOffset,
                          std::vector<int64_t> &nodeNeighbors,
                          std::vector<int64_t> &edgeWeights,
                          std::vector<int64_t> &nodeWeights) const;

    /**
     * @brief Gets initial bisection using Scotch partitioner
     * @return Boolean vector indicating partition assignment
     */
    [[nodiscard]] std::vector<uint8_t> getUndirectedBisectionScotch() const;

    /**
     * @brief Gets initial bisection using METIS partitioner
     * @return Boolean vector indicating partition assignment
     */
    [[nodiscard]] std::vector<uint8_t> getUndirectedBisectionMetis() const;

    /**
     * @brief Fixes acyclicity by propagating partition 0 assignment upward
     *
     * Implements Algorithm 3 from the paper: moves ancestors of V0 vertices
     * to V0 in reverse topological order.
     *
     * @param undirectedBisection Bisection to be fixed
     */
    void fixAcyclicityUp(std::vector<uint8_t> &undirectedBisection) const;

    /**
     * @brief Fixes acyclicity by propagating partition 1 assignment downward
     *
     * Implements Algorithm 4 from the paper: moves descendants of V1 vertices
     * to V1 in topological order.
     *
     * @param undirectedBisection Bisection to be fixed
     */
    void fixAcyclicityDown(std::vector<uint8_t> &undirectedBisection) const;

    /**
     * @brief Returns the best result based on some criteria
     *
     * From the balanced results, from the ones with non-zero edge cut, take the one with the smallest edge cut.
     * From the balanced results, if all have zero edge cut, take whichever.
     * If all are unbalanced, take the least unbalanced one with tie breaks favouring the lowest edge cut.
     *
     * @param results Vector of results
     * @return The index of the best result
     */
    static uint64_t selectBestResult(const std::vector<std::tuple<uint64_t, uint8_t, double>> &results);

    /**
     * @brief Runs the complete undirected bisection and fixing process
     *
     * For the METIS partitioner:
     * 1. Gets initial undirected bisection
     * 2. Tries both partition assignments (P0=V0,P1=V1 and P1=V0,P0=V1)
     * 3. Fixes acyclicity both upward and downward
     * 4. Returns best result among all attempts
     *
     * @return Pair containing:
     *         - Acyclic bisection vector
     *         - Edge cut weight
     */
    [[nodiscard]] std::pair<std::vector<uint8_t>, uint64_t> runMetis() const;

    /**
     * @brief Runs the complete undirected bisection and fixing process
     *
     * For the Scotch partitioner:
     * 1. Gets initial undirected bisection
     * 2. Tries both partition assignments (P0=V0,P1=V1 and P1=V0,P0=V1)
     * 3. Fixes acyclicity both upward and downward
     * 4. Returns best result among all attempts
     *
     * @return Pair containing:
     *         - Acyclic bisection vector
     *         - Edge cut weight
     */
    [[nodiscard]] std::pair<std::vector<uint8_t>, uint64_t> runScotch() const;

    /**
     * @brief Runs the 2 partitioners and obtains best result
     *
     * Based on the enabled flags, runs the 2 partitioners (METIS & Scotch)
     * and returns the best result, giving priority to non-zero edge cuts.
     *
     * @return Pair containing:
     *         - Acyclic bisection vector
     *         - Edge cut weight
     */
    [[nodiscard]] std::pair<std::vector<uint8_t>, uint64_t> run() const override;

public:
    /**
     * @brief Constructs the undirected fix bisection algorithm
     * @param graph Graph to be bisected
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum required partition weight
     * @param refinementMethod Refinement method
     * @param refinementPasses Number of refinement passes
     * @param useMetis Whether to try METIS partitioning
     * @param useScotch Whether to try Scotch partitioning
     */
    UndirectedFix(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight,
                  RefinementMethod refinementMethod, uint64_t refinementPasses, bool useMetis, bool useScotch);
};

#endif //DAG_PARTITIONING_UNDIRECTEDFIX_H