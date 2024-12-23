//
// Created by panagiotis on 21/12/2024.
//

#ifndef DAG_PARTITIONING_MULTILEVELBISECTIONER_H
#define DAG_PARTITIONING_MULTILEVELBISECTIONER_H


#include <stack>
#include "Graph.h"

enum class ClusteringMethod {
    FORB,
    CYC,
    HYB
};

enum class BisectionMethod {
    GGG,
    UNDIRMETIS,
    UNDIRSCOTCH,
    UNDIRBOTH
};

class MultilevelBisectioner {
private:
    Graph workingGraph;
    ClusteringMethod clusteringMethod;
    uint64_t maxClusteringRounds;
    uint64_t minClusteringVertices;
    BisectionMethod bisectionMethod;
    double imbalanceRatio;
public:
    MultilevelBisectioner(Graph graph, ClusteringMethod clusteringMethod,
                          uint64_t maxClusteringRounds,
                          uint64_t minClusteringVertices,
                          BisectionMethod bisectionMethod,
                          double imbalanceRatio);

    [[nodiscard]] std::stack<std::pair<Graph, std::vector<uint64_t>>> runClustering() const;

    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> runBisection(const Graph &graph) const;

    static void projectBisection(std::pair<std::vector<bool>, uint64_t> &bisectionInfo,
                                 const std::vector<uint64_t> &mapping);

    static void dummyRefinement(const Graph &graph, std::pair<std::vector<bool>, uint64_t> &bisectionInfo);

    [[nodiscard]] std::pair<std::vector<bool>, uint64_t> run() const;
};


#endif //DAG_PARTITIONING_MULTILEVELBISECTIONER_H
