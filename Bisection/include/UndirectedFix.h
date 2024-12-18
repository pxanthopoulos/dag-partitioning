//
// Created by panagiotis on 17/12/2024.
//

#ifndef DAG_PARTITIONING_UNDIRECTEDFIX_H
#define DAG_PARTITIONING_UNDIRECTEDFIX_H


#include "Bisection.h"
#include <stdio.h>
#include "scotch.h"
#include "metis.h"

class UndirectedFix : Bisection {
public:
    bool useMetis;
    bool useScotch;

    explicit UndirectedFix(Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight, bool useMetis, bool useScotch);

    [[nodiscard]] int64_t computeNumberOfEdges() const;

    void graphToCSRFormat(int64_t edgeNumber, std::vector<int64_t> &nodeNeighborsOffset, std::vector<int64_t> &nodeNeighbors, std::vector<int64_t> &edgeWeights, std::vector<int64_t> &nodeWeights) const;

    [[nodiscard]] std::vector<bool> getUndirectedPartitionScotch();

    [[nodiscard]] std::vector<bool> getUndirectedPartitionMetis();

    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> run() const override;
};


#endif //DAG_PARTITIONING_UNDIRECTEDFIX_H
