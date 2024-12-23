//
// Created by panagiotis on 23/12/2024.
//

#ifndef DAG_PARTITIONING_BOUNDARYFM_H
#define DAG_PARTITIONING_BOUNDARYFM_H


#include "Refinement.h"

class BoundaryFM : public Refinement {
private:
    bool onePassRefinement() override;

public:
    BoundaryFM(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t maxNumberOfPasses);
};


#endif //DAG_PARTITIONING_BOUNDARYFM_H
