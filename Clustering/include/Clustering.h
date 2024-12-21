//
// Created by panagiotis on 8/12/2024.
//

#ifndef DAG_PARTITIONING_CLUSTERING_H
#define DAG_PARTITIONING_CLUSTERING_H


#include "Graph.h"
#include <vector>

class Clustering {
protected:
    Graph workingGraph;
    std::vector<uint64_t> clusters;

    [[nodiscard]] virtual std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const = 0;

    [[nodiscard]] bool updateGraphAndClusters(const std::vector<uint64_t> &leaders, uint64_t newSize);

public:
    explicit Clustering(const Graph &graph);

    virtual ~Clustering() = default;

    void run();
};


#endif //DAG_PARTITIONING_CLUSTERING_H
