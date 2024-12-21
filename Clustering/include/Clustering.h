//
// Created by panagiotis on 8/12/2024.
//

#ifndef DAG_PARTITIONING_CLUSTERING_H
#define DAG_PARTITIONING_CLUSTERING_H


#include "Graph.h"
#include <stack>
#include <vector>

class Clustering {
protected:
    Graph workingGraph;
    std::stack<std::pair<Graph, std::vector<uint64_t>>> intermediateGraphsAndClusters;
    uint64_t maxRounds;
    uint64_t minVertices;

    [[nodiscard]] virtual std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const = 0;

    [[nodiscard]] bool updateGraphAndClusters(const std::vector<uint64_t> &leaders, uint64_t newSize);

public:
    explicit Clustering(Graph graph, uint64_t maxRounds, uint64_t minVertices);

    virtual ~Clustering() = default;

    [[nodiscard]] std::stack<std::pair<Graph, std::vector<uint64_t>>> run();
};


#endif //DAG_PARTITIONING_CLUSTERING_H
