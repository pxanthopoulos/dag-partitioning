//
// Created by panagiotis on 16/12/2024.
//

#ifndef DAG_PARTITIONING_BISECTION_H
#define DAG_PARTITIONING_BISECTION_H


#include "Graph.h"

class Bisection {
public:
    Graph workingGraph;
    double upperBoundPartWeight, lowerBoundPartWeight;

    explicit Bisection(Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight);

    virtual bool checkValidBisection(const std::vector<bool> &bisection);

    virtual std::pair<std::vector<bool>, uint64_t> run() = 0;
};


#endif //DAG_PARTITIONING_BISECTION_H
