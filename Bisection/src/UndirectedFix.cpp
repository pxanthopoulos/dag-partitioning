//
// Created by panagiotis on 17/12/2024.
//

#include <iostream>
#include "UndirectedFix.h"

UndirectedFix::UndirectedFix(Graph &graph, double upperBoundPartWeight,
                             double lowerBoundPartWeight, bool useMetis, bool useScotch) : Bisection(graph,
                                                                      upperBoundPartWeight,
                                                                      lowerBoundPartWeight), useMetis(useMetis), useScotch(useScotch) {
    assert((useMetis || useScotch) && "Specify at least 1 undirected partitioning algorithm");
}

int64_t UndirectedFix::computeNumberOfEdges() const {
    int64_t edgeNumber = 0;
    for (const auto & neighbors : workingGraph.adj) edgeNumber += static_cast<int64_t>(neighbors.size());
    return edgeNumber;
}

void UndirectedFix::graphToScotchGraph(int64_t edgeNumber, std::vector<int64_t> &nodeNeighborsOffset, std::vector<int64_t> &nodeNeighbors, std::vector<int64_t> &edgeWeights, std::vector<int64_t> &nodeWeights) const {
    int64_t pos = 0;
    int64_t cumulativeNeighbors = 0;
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        nodeNeighborsOffset[i] = cumulativeNeighbors;
        nodeWeights[i] = static_cast<int64_t>(workingGraph.nodeWeights[i]);
        const auto & neighbors = workingGraph.adj[i];
        for (const auto & [neighborId, edgeWeight] : neighbors) {
            assert((pos < 2 * edgeNumber) && "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = static_cast<int64_t>(neighborId);
            edgeWeights[pos] = static_cast<int64_t>(edgeWeight);
            pos++;
        }
        const auto & reverseNeighbors = workingGraph.revAdj[i];
        for (const auto & [reverseNeighborId, reverseEdgeWeight] : reverseNeighbors) {
            assert((pos < 2 * edgeNumber) && "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = static_cast<int64_t>(reverseNeighborId);
            edgeWeights[pos] = static_cast<int64_t>(reverseEdgeWeight);
            pos++;
        }
        cumulativeNeighbors += static_cast<int64_t>(neighbors.size());
        cumulativeNeighbors += static_cast<int64_t>(reverseNeighbors.size());
    }
    nodeNeighborsOffset[workingGraph.size] = cumulativeNeighbors;
    assert(pos == 2 * edgeNumber && "pos must be equal to 2 * #Edges");
}

std::vector<int64_t> UndirectedFix::getUndirectedPartitionScotch() {
    SCOTCH_Graph scotchGraph;
    SCOTCH_Strat scotchStrategy;
    SCOTCH_graphInit(&scotchGraph);
    SCOTCH_stratInit(&scotchStrategy);

    int64_t edgeNumber = computeNumberOfEdges();
    std::vector<int64_t> nodeNeighborsOffset(workingGraph.size + 1);
    std::vector<int64_t> nodeNeighbors(2 * edgeNumber);
    std::vector<int64_t> edgeWeights(2 * edgeNumber);
    std::vector<int64_t> nodeWeights(workingGraph.size);
    graphToScotchGraph(edgeNumber, nodeNeighborsOffset, nodeNeighbors, edgeWeights, nodeWeights);

    int ret = SCOTCH_graphBuild(&scotchGraph, 0, static_cast<int64_t>(workingGraph.size), nodeNeighborsOffset.data(), nodeNeighborsOffset.data() + 1, nodeWeights.data(), nullptr, edgeNumber * 2, nodeNeighbors.data(), edgeWeights.data());
    assert(ret == 0 && "scotch graph build failed");

    ret = SCOTCH_graphCheck(&scotchGraph);
    assert(ret == 0 && "scotch graph check failed");

    double imbalanceRatio = upperBoundPartWeight * 2 / (double) workingGraph.totalWeight;
    ret = SCOTCH_stratGraphMapBuild(&scotchStrategy, SCOTCH_STRATQUALITY | SCOTCH_STRATBALANCE, 2, imbalanceRatio - 1);
    assert(ret == 0 && "scotch strategy build failed");

    std::vector<int64_t> partitions(workingGraph.size);
    ret = SCOTCH_graphPart (&scotchGraph, 2, &scotchStrategy, partitions.data());
    assert(ret == 0 && "scotch partition failed");

    SCOTCH_graphExit(&scotchGraph);
    SCOTCH_stratExit(&scotchStrategy);

    return partitions;
}

std::pair<std::vector<bool>, uint64_t> UndirectedFix::run() const {
    return {};
}
