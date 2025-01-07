/**
 * @file BoundaryKL.h
 * @brief Implementation of boundary Kernighan-Lin refinement for directed graphs
 *
 * Essentially the same as FM but chooses pairs of nodes (one in V0 and one in V1)
 * to swap. Picks the pair with the best total gain (for balanced partitions) or
 * the pair with the best balance improvement (for unbalanced partitions) at each step.
 */

#ifndef DAG_PARTITIONING_BOUNDARYKL_H
#define DAG_PARTITIONING_BOUNDARYKL_H

#include <list>
#include "Refinement.h"

class BoundaryKL : public Refinement {
private:
    /**
     * @brief Implements one pass of boundary KL refinement
     *
     * Repeatedly finds pairs of vertices (one from V0, one from V1) to swap.
     * For balanced partitions, selects pair with highest total gain.
     * For unbalanced partitions, selects pair that most improves balance.
     * Each vertex can only move once per pass. At end of pass,
     * realizes most profitable prefix of pairs.
     *
     * @return true if improvements were made in this pass
     */
    [[nodiscard]] bool onePassRefinement() override;

    /**
     * @brief Identifies initially movable vertices and adds them to appropriate lists
     *
     * A vertex is movable if:
     * - From V0: all out-neighbors in V1 or no out-neighbors
     * - From V1: all in-neighbors in V0 or no in-neighbors
     * Gain is calculated as sum of outgoing weights minus sum of incoming weights
     *
     * @param listV0 List for V0 vertices
     * @param listV1 List for V1 vertices
     * @param inList Tracks which vertices are in lists
     */
    void insertMovableNodesIntoLists(std::list<std::pair<int64_t, uint64_t>> &listV0,
                                     std::list<std::pair<int64_t, uint64_t>> &listV1,
                                     std::vector<bool> &inList) const;

    /**
     * @brief Updates movable vertices after a move
     *
     * When a vertex moves, some neighbors may become movable (added to lists). The lists are kept sorted with the additions.
     *
     * @param currentBisectionInfo Current partition assignments
     * @param listV0 List for V0 vertices
     * @param listV1 List for V1 vertices
     * @param inList Tracks which vertices are in lists
     * @param movedNodeId Node that was just moved
     * @param checkOutNeighbors Whether to check outgoing (true) or incoming (false) neighbors
     */
    void insertMovableNeighborsIntoLists(const std::vector<bool> &currentBisectionInfo,
                                         std::list<std::pair<int64_t, uint64_t>> &listV0,
                                         std::list<std::pair<int64_t, uint64_t>> &listV1,
                                         std::vector<bool> &inList, uint64_t movedNodeId,
                                         bool checkOutNeighbors) const;

    /**
     * @brief Finds the best swappable pair from the movable elements, if the partition is already balanced
     *
     * Checks all pairs and identifies the best (maximum gain) and actually swappable (no edge between the nodes) pair
     * from the movable elements, while ensuring the swap does not break the balance. If none is found, the function returns a flag.
     *
     * @param listV0 List for V0 vertices
     * @param listV1 List for V1 vertices
     * @param moved Tracks which have already been moved (or marked unmovable)
     * @param sizeV0 Size of V0
     * @param sizeV1 Size of V1
     * @return Returns a tuple. The first element informs if a swappable pair was found. The other 2 elements contain the swappable pair. The last element is the gain.
     */
    [[nodiscard]] std::tuple<bool, uint64_t, uint64_t, int64_t>
    findBestMovablePairBalanced(std::list<std::pair<int64_t, uint64_t>> &listV0,
                                std::list<std::pair<int64_t, uint64_t>> &listV1,
                                std::vector<bool> &moved, uint64_t &sizeV0,
                                uint64_t &sizeV1) const;

    /**
     * @brief Finds the best swappable pair from the movable elements, if the partition is unbalanced
     *
     * Checks all pairs and identifies the best (maximum balance improvement) and actually swappable (no edge between the nodes) pair
     * from the movable elements, while ensuring the swap does not cause unbalance the other way. If none is found, the function returns a flag.
     *
     * @param listV0 List for V0 vertices
     * @param listV1 List for V1 vertices
     * @param moved Tracks which have already been moved (or marked unmovable)
     * @param sizeV0 Size of V0
     * @param sizeV1 Size of V1
     * @return Returns a tuple. The first element informs if a swappable pair was found. The other 2 elements contain the swappable pair. The last element is the gain.
     */
    [[nodiscard]] std::tuple<bool, uint64_t, uint64_t, int64_t>
    findBestMovablePairUnbalanced(std::list<std::pair<int64_t, uint64_t>> &listV0,
                                  std::list<std::pair<int64_t, uint64_t>> &listV1,
                                  std::vector<bool> &moved, uint64_t &sizeV0,
                                  uint64_t &sizeV1) const;

public:
    /**
     * @brief Constructs the boundary KL refinement algorithm
     * @param graph Graph being refined
     * @param initialBisectionInfo Current partition assignments
     * @param initialEdgeCut Current edge cut weight
     * @param maxNumberOfPasses Maximum refinement iterations
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum allowed partition weight
     */
    BoundaryKL(const Graph &graph, std::vector<bool> &initialBisectionInfo,
               uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
               double upperBoundPartWeight, double lowerBoundPartWeight);
};

#endif //DAG_PARTITIONING_BOUNDARYKL_H