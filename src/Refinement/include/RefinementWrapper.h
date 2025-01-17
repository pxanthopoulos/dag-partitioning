//
// Created by panagiotis on 17/1/2025.
//

#ifndef DAG_PARTITIONING_REFINEMENTWRAPPER_H
#define DAG_PARTITIONING_REFINEMENTWRAPPER_H

#include "Graph.h"
#include "Refinement.h"

void refinementWrapper(const Graph &graph, std::pair<std::vector<uint8_t>, uint64_t> &bisectionInfoPair,
                       RefinementMethod refinementMethod, uint64_t refinementPasses, double upperBoundPartWeight,
                       double lowerBoundPartWeight);


#endif //DAG_PARTITIONING_REFINEMENTWRAPPER_H
