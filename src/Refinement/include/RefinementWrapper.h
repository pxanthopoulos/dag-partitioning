//
// Created by panagiotis on 17/1/2025.
//

#ifndef DAG_PARTITIONING_REFINEMENTWRAPPER_H
#define DAG_PARTITIONING_REFINEMENTWRAPPER_H

#include "Graph.h"
#include "Refinement.h"

void refinementWrapper(const Graph &graph, std::vector<uint8_t> &bisectionInfo, uint64_t &edgeCut,
                       RefinementMethod refinementMethod, uint64_t refinementPasses, double upperBoundPartWeight,
                       double lowerBoundPartWeight);


#endif //DAG_PARTITIONING_REFINEMENTWRAPPER_H
