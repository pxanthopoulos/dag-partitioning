//
// Created by panagiotis on 10/12/2024.
//

#ifndef DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H
#define DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H


#include "Bisection.h"

class GreedyDirectedGraphGrowing : Bisection {
public:
    explicit GreedyDirectedGraphGrowing(Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight);

    std::pair<std::vector<bool>, uint64_t> runOnNormalGraph();

    std::pair<std::vector<bool>, uint64_t> runOnReverseGraph();

    std::pair<std::vector<bool>, uint64_t> run() override;
};


#endif //DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H