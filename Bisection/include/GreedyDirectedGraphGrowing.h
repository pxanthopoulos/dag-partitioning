//
// Created by panagiotis on 10/12/2024.
//

#ifndef DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H
#define DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H


#include "Bisection.h"

class GreedyDirectedGraphGrowing : public Bisection {
private:
    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> runOnNormalGraph() const;

    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> runOnReverseGraph() const;

public:
    GreedyDirectedGraphGrowing(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight);

    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> run() const override;
};


#endif //DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H