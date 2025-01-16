/**
 * @file BoundaryFM.h
 * @brief Implementation of boundary Fiduccia-Mattheyses refinement for directed graphs
 *
 * Adapts the FM algorithm for directed graphs by defining vertex movability based on
 * neighbor conditions and using static gain calculations. A vertex can move:
 * - V0->V1: if all out-neighbors in V1 or no out-neighbors
 * - V1->V0: if all in-neighbors in V0 or no in-neighbors
 */

#ifndef DAG_PARTITIONING_BOUNDARYFM_H
#define DAG_PARTITIONING_BOUNDARYFM_H

#include "Refinement.h"
#include <queue>

class BoundaryFM : virtual public Refinement {
protected:
    /**
     * @brief Implements one pass of boundary FM refinement
     *
     * Repeatedly selects highest gain movable vertex and tentatively moves it.
     * At each round, only vertices from the heavier part are considered for moving.
     * At end of pass, realizes most profitable prefix of moves.
     * Each vertex can move at most once per pass.
     *
     * @return true if improvements were made in this pass
     */
    [[nodiscard]] bool onePassRefinement() override;

private:
    /**
     * @brief Identifies initially movable vertices and adds them to appropriate heaps
     *
     * A vertex is movable if:
     * - From V0: all out-neighbors in V1 or no out-neighbors
     * - From V1: all in-neighbors in V0 or no in-neighbors
     * Gain is calculated as sum of outgoing weights minus sum of incoming weights
     *
     * @param heapV0 Priority queue for V0 vertices
     * @param heapV1 Priority queue for V1 vertices
     * @param inHeap Tracks which vertices are in heaps
     */
    void insertMovableNodesIntoHeaps(std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                     std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                     std::vector<uint8_t> &inHeap) const;

    /**
     * @brief Updates movable vertices after a move
     *
     * When a vertex moves, some neighbors may become movable (added to heap)
     *
     * @param currentBisectionInfo Current partition assignments
     * @param heapV0 Priority queue for V0 vertices
     * @param heapV1 Priority queue for V1 vertices
     * @param inHeap Tracks which vertices are in heaps
     * @param movedNodeId Node that was just moved
     * @param checkOutNeighbors Whether to check outgoing (true) or incoming (false) neighbors
     */
    void insertMovableNeighborsIntoHeaps(const std::vector<uint8_t> &currentBisectionInfo,
                                         std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                         std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                         std::vector<uint8_t> &inHeap, uint64_t movedNodeId,
                                         bool checkOutNeighbors) const;

    /**
     * @brief Checks if move maintains or improves partition balance
     *
     * If partition is already balanced, move must maintain balance.
     * If unbalanced, move must improve balance without going too far the other way.
     *
     * @param moveFromV0 Whether move is from V0->V1 (true) or V1->V0 (false)
     * @param movedNodeId Node being moved
     * @param maxNodeWeight Weight of heaviest node
     * @param sizeV0 Current V0 partition weight
     * @param sizeV1 Current V1 partition weight
     * @param isBalanced Whether partition is currently balanced
     * @return true if move is acceptable for balance
     */
    [[nodiscard]] bool isBalanceImprovedOrMaintained(bool moveFromV0, uint64_t movedNodeId,
                                                     uint64_t maxNodeWeight, uint64_t sizeV0,
                                                     uint64_t sizeV1, bool isBalanced) const;

public:
    /**
     * @brief Constructs the boundary FM refinement algorithm
     * @param graph Graph being refined
     * @param initialBisectionInfo Current partition assignments
     * @param initialEdgeCut Current edge cut weight
     * @param maxNumberOfPasses Maximum refinement iterations
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum allowed partition weight
     */
    BoundaryFM(const Graph &graph, std::vector<uint8_t> &initialBisectionInfo,
               uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
               double upperBoundPartWeight, double lowerBoundPartWeight);
};

#endif //DAG_PARTITIONING_BOUNDARYFM_H