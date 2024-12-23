//
// Created by panagiotis on 17/12/2024.
//

#ifndef DAG_PARTITIONING_UNDIRECTEDFIX_H
#define DAG_PARTITIONING_UNDIRECTEDFIX_H


#include "Bisection.h"
#include <stdio.h>
#include "scotch.h"
#include "metis.h"

class UndirectedFix : public Bisection {
private:
    bool useMetis;
    bool useScotch;

    [[nodiscard]] int64_t computeNumberOfEdges() const;

    void
    graphToCSRFormat(int64_t edgeNumber, std::vector<int64_t> &nodeNeighborsOffset, std::vector<int64_t> &nodeNeighbors,
                     std::vector<int64_t> &edgeWeights, std::vector<int64_t> &nodeWeights) const;

    [[nodiscard]] std::vector<bool> getUndirectedBisectionScotch() const;

    [[nodiscard]] std::vector<bool> getUndirectedBisectionMetis() const;

    void fixAcyclicityUp(std::vector<bool> &undirectedBisection) const;

    void fixAcyclicityDown(std::vector<bool> &undirectedBisection) const;

    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> run() const override;

public:
    UndirectedFix(const Graph &graph, double upperBoundPartWeight, double lowerBoundPartWeight, bool useMetis,
                  bool useScotch);
};


#endif //DAG_PARTITIONING_UNDIRECTEDFIX_H
