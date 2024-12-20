//
// Created by panagiotis on 10/12/2024.
//

#ifndef DAG_PARTITIONING_CLUSTERINGHYBRID_H
#define DAG_PARTITIONING_CLUSTERINGHYBRID_H


#include "ClusteringCycleDetection.h"

class ClusteringHybrid : public ClusteringCycleDetection {
public:
    explicit ClusteringHybrid(const Graph &graph);

    [[nodiscard]] bool checkLargeDegrees(uint64_t from, uint64_t to) const;

    void bookKeepingForForbiddenEdges(const std::vector<uint64_t> &topLevels, uint64_t node,
                                      const std::vector<std::tuple<uint64_t, uint64_t, bool>> &sortedNeighbors,
                                      uint64_t neighborId, uint64_t leaderOfNeighbor,
                                      const std::vector<bool> &markup,
                                      const std::vector<bool> &markdown,
                                      std::vector<uint64_t> &numberOfBadNeighbors,
                                      std::vector<uint64_t> &leaderOfBadNeighbors) const;

protected:
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const override;
};


#endif //DAG_PARTITIONING_CLUSTERINGHYBRID_H
