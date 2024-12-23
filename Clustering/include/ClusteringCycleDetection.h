//
// Created by panagiotis on 8/12/2024.
//

#ifndef DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H
#define DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H


#include "Clustering.h"
#include <chrono>

class ClusteringCycleDetection : virtual public Clustering {
public:
    ClusteringCycleDetection(const Graph &graph, uint64_t maxRounds, uint64_t minVertices);

private:
    void hardCheckCycle(const std::vector<uint64_t> &leaders, uint64_t newSize) const;

    [[nodiscard]] static uint64_t findMinimumTopLevelInCluster(uint64_t node, const std::vector<uint64_t> &topLevels,
                                                               const std::vector<uint64_t> &leaders);

protected:
    [[nodiscard]] bool detectCycle(uint64_t from, uint64_t to, const std::vector<uint64_t> &topLevels,
                                   const std::vector<uint64_t> &leaders) const;

    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const override;
};


#endif //DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H
