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

    [[nodiscard]] virtual bool checkValidBisection(const std::vector<bool> &bisection) const;

    [[nodiscard]] virtual uint64_t computeEdgeCut(const std::vector<bool> &bisection) const;

    [[nodiscard]] virtual std::pair<std::vector<bool>, uint64_t> run() const = 0;
};


#endif //DAG_PARTITIONING_BISECTION_H
