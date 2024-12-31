//
// Created by panagiotis on 23/12/2024.
//

#ifndef DAG_PARTITIONING_BOUNDARYFM_H
#define DAG_PARTITIONING_BOUNDARYFM_H


#include "Refinement.h"
#include <queue>

class BoundaryFM : public Refinement {
private:
    bool onePassRefinement() override;

    void insertMovableNodesIntoHeaps(std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                     std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                     std::vector<bool> &inHeap) const;

    void insertMovableNeighborsIntoHeaps(const std::vector<bool> &currentBisectionInfo,
                                         std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                         std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                         std::vector<bool> &inHeap, uint64_t movedNodeId, bool checkOutNeighbors) const;

    bool isBalanceImprovedOrMaintained(bool moveFromV0, uint64_t movedNodeId, uint64_t maxNodeWeight, uint64_t sizeV0,
                                       uint64_t sizeV1, bool isBalanced);

public:
    BoundaryFM(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t &initialEdgeCut,
               uint64_t maxNumberOfPasses,
               double upperBoundPartWeight, double lowerBoundPartWeight);
};


#endif //DAG_PARTITIONING_BOUNDARYFM_H
