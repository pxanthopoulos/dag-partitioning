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
    double upperBoundPartWeight, lowerBoundPartWeight;

    [[nodiscard]] bool checkValidBisection() const;

    [[nodiscard]] std::pair<uint64_t, uint64_t> calculatePartSizes() const;

    [[nodiscard]] bool checkBalance(uint64_t maxNodeWeight) const;

    virtual bool onePassRefinement() = 0;

    Refinement(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t maxNumberOfPasses,
               double upperBoundPartWeight, double lowerBoundPartWeight);

public:
    virtual ~Refinement() = default;

    void run();
};


#endif //DAG_PARTITIONING_REFINEMENT_H
