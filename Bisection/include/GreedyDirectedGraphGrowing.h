/**
 * @file GreedyDirectedGraphGrowing.h
 * @brief Implementation of greedy directed graph growing bisection as described in Section 4.2.1
 *
 * This class implements a two-phase greedy algorithm that moves vertices between partitions
 * while maintaining acyclicity. It tries both directions (V1->V0 and V0->V1) and selects
 * the better result.
 */

#ifndef DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H
#define DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H

#include "Bisection.h"

class GreedyDirectedGraphGrowing : public Bisection {
private:
    /**
     * @brief Runs the greedy algorithm in normal direction (V1->V0)
     *
     * Starts with all vertices in V1 and moves vertices to V0. Uses two phases:
     * 1. First phase uses incoming edge weights as priority
     * 2. Second phase uses gain (incoming - outgoing weights) as priority
     *
     * @return Pair containing:
     *         - Bisection vector (false=V0, true=V1)
     *         - Total edge cut weight
     */
    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> runOnNormalGraph() const;

    /**
     * @brief Runs the greedy algorithm in reverse direction (V0->V1)
     *
     * Similar to normal direction but starts with all vertices in V0 and moves to V1,
     * considering out-neighbors instead of in-neighbors.
     *
     * @return Pair containing:
     *         - Bisection vector (false=V0, true=V1)
     *         - Total edge cut weight
     */
    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> runOnReverseGraph() const;

    /**
     * @brief Runs both normal and reverse algorithms and selects better result
     *
     * As described in the paper: "returns the better of these two alternatives"
     * based on edge cut weight.
     *
     * @return Best bisection result from either direction
     */
    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> run() const override;

public:
    /**
     * @brief Constructs the greedy bisection algorithm instance
     * @param graph Graph to be bisected
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum required partition weight
     */
    GreedyDirectedGraphGrowing(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight);
};

#endif //DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H