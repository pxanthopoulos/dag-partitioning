//
// Created by panagiotis on 8/12/2024.
//

#ifndef DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H
#define DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H


#include "Clustering.h"
#include <chrono>

class ClusteringCycleDetection : public Clustering {
private:
    void hardCheckCycle(const std::vector<uint64_t> &leaders, uint64_t newSize) const;

    [[nodiscard]] static uint64_t findMinimumTopLevelInCluster(uint64_t node, const std::vector<uint64_t> &topLevels,
                                                               const std::vector<uint64_t> &leaders);

    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const override;

protected:
    [[nodiscard]] bool detectCycle(uint64_t from, uint64_t to, const std::vector<uint64_t> &topLevels,
                                   const std::vector<uint64_t> &leaders) const;

public:
    ClusteringCycleDetection(const Graph &graph, uint64_t maxRounds, uint64_t minVertices);
};


#endif //DAG_PARTITIONING_CLUSTERINGCYCLEDETECTION_H
