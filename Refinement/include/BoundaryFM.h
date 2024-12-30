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

    void insertMovableNodesIntoHeaps(const std::vector<bool> &initialBisectionInfo,
                                     std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                     std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                     std::vector<bool> &inHeap);

    void insertMovableNeighborsIntoHeaps(const std::vector<bool> &currentBisectionInfo,
                                         std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                         std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                         std::vector<bool> &inHeap, uint64_t movedNodeId, bool checkOutNeighbors);

public:
    BoundaryFM(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t maxNumberOfPasses);
};


#endif //DAG_PARTITIONING_BOUNDARYFM_H
