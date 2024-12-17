//
// Created by panagiotis on 17/12/2024.
//

#ifndef DAG_PARTITIONING_UNDIRECTEDFIX_H
#define DAG_PARTITIONING_UNDIRECTEDFIX_H


#include "Bisection.h"

class UndirectedFix : Bisection {
public:
    explicit UndirectedFix(Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight);

    void getUndirectedPartition();

    std::pair<std::vector<bool>, uint64_t> run() override;
};


#endif //DAG_PARTITIONING_UNDIRECTEDFIX_H
