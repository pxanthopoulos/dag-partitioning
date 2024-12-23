//
// Created by panagiotis on 22/12/2024.
//

#ifndef DAG_PARTITIONING_REFINEMENT_H
#define DAG_PARTITIONING_REFINEMENT_H


#include "Graph.h"

class Refinement {
protected:
    const Graph &workingGraph;
    std::vector<bool> &initialBisectionInfo;
    uint64_t maxNumberOfPasses;

    [[nodiscard]] bool checkValidBisection() const;

    virtual bool onePassRefinement() = 0;

    Refinement(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t maxNumberOfPasses);

public:
    virtual ~Refinement() = default;

    void run();
};


#endif //DAG_PARTITIONING_REFINEMENT_H
