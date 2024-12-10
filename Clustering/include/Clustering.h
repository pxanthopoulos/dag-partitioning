//
// Created by panagiotis on 8/12/2024.
//

#ifndef DAG_PARTITIONING_CLUSTERING_H
#define DAG_PARTITIONING_CLUSTERING_H


#include "Graph.h"
#include <vector>

class Clustering {
public:
    Graph workingGraph;
    std::vector<uint64_t> clusters;

    explicit Clustering(const Graph &graph);

    virtual ~Clustering() = default;

    [[nodiscard]] virtual std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const = 0;

    [[nodiscard]] bool updateGraphAndClusters(const std::vector<uint64_t> &leaders, uint64_t newSize);

    void runClustering();
};


#endif //DAG_PARTITIONING_CLUSTERING_H
