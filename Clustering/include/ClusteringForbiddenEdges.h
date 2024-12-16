//
// Created by panagiotis on 6/12/2024.
//

#ifndef DAG_PARTITIONING_CLUSTERINGFORBIDDENEDGES_H
#define DAG_PARTITIONING_CLUSTERINGFORBIDDENEDGES_H


#include "Clustering.h"
#include <vector>
#include <cstdint>
#include <utility>

class ClusteringForbiddenEdges : virtual public Clustering {
public:
    explicit ClusteringForbiddenEdges(const Graph &graph);

    [[nodiscard]] std::vector<std::pair<uint64_t, uint64_t>>
    findValidNeighbors(uint64_t node,
                       const std::vector<std::tuple<uint64_t, uint64_t, bool>> &neighbors,
                       const std::vector<uint64_t> &leaders,
                       const std::vector<uint64_t> &clusterWeights,
                       const std::vector<uint64_t> &topLevels,
                       const std::vector<uint64_t> &numberOfBadNeighbors,
                       const std::vector<uint64_t> &leaderOfBadNeighbors) const;

    [[nodiscard]] static uint64_t findBestNeighbor(const std::vector<std::pair<uint64_t, uint64_t>> &);

    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> oneRoundClustering() const override;
};


#endif //DAG_PARTITIONING_CLUSTERINGFORBIDDENEDGES_H
