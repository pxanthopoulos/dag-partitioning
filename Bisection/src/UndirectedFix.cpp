//
// Created by panagiotis on 17/12/2024.
//

#include <iostream>
#include <cmath>
#include "UndirectedFix.h"

UndirectedFix::UndirectedFix(Graph &graph, double upperBoundPartWeight,
                             double lowerBoundPartWeight, bool useMetis, bool useScotch) : Bisection(graph,
                                                                      upperBoundPartWeight,
                                                                      lowerBoundPartWeight), useMetis(useMetis), useScotch(useScotch) {
    assert((useMetis || useScotch) && "Specify at least 1 undirected partitioning algorithm");
}

int64_t UndirectedFix::computeNumberOfEdges() const {
    int64_t edgeNumber = 0;
    for (const auto & neighbors : workingGraph.adj) edgeNumber += (int64_t)neighbors.size();
    return edgeNumber;
}

void UndirectedFix::graphToCSRFormat(int64_t edgeNumber, std::vector<int64_t> &nodeNeighborsOffset, std::vector<int64_t> &nodeNeighbors, std::vector<int64_t> &edgeWeights, std::vector<int64_t> &nodeWeights) const {
    int64_t pos = 0;
    int64_t cumulativeNeighbors = 0;
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        nodeNeighborsOffset[i] = cumulativeNeighbors;
        nodeWeights[i] = (int64_t)workingGraph.nodeWeights[i];
        const auto & neighbors = workingGraph.adj[i];
        for (const auto & [neighborId, edgeWeight] : neighbors) {
            assert((pos < 2 * edgeNumber) && "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = (int64_t)neighborId;
            edgeWeights[pos] = (int64_t)edgeWeight;
            pos++;
        }
        const auto & reverseNeighbors = workingGraph.revAdj[i];
        for (const auto & [reverseNeighborId, reverseEdgeWeight] : reverseNeighbors) {
            assert((pos < 2 * edgeNumber) && "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = (int64_t)reverseNeighborId;
            edgeWeights[pos] = (int64_t)reverseEdgeWeight;
            pos++;
        }
        cumulativeNeighbors += (int64_t)neighbors.size();
        cumulativeNeighbors += (int64_t)reverseNeighbors.size();
    }
    nodeNeighborsOffset[workingGraph.size] = cumulativeNeighbors;
    assert(pos == 2 * edgeNumber && "pos must be equal to 2 * #Edges");
}

std::vector<bool> UndirectedFix::getUndirectedPartitionScotch() {
    SCOTCH_Graph scotchGraph;
    SCOTCH_Strat scotchStrategy;
    SCOTCH_graphInit(&scotchGraph);
    SCOTCH_stratInit(&scotchStrategy);

    int64_t edgeNumber = computeNumberOfEdges();
    std::vector<int64_t> nodeNeighborsOffset(workingGraph.size + 1);
    std::vector<int64_t> nodeNeighbors(2 * edgeNumber);
    std::vector<int64_t> edgeWeights(2 * edgeNumber);
    std::vector<int64_t> nodeWeights(workingGraph.size);
    graphToCSRFormat(edgeNumber, nodeNeighborsOffset, nodeNeighbors, edgeWeights, nodeWeights);

    int ret = SCOTCH_graphBuild(&scotchGraph, 0, (int64_t)workingGraph.size, nodeNeighborsOffset.data(), nodeNeighborsOffset.data() + 1, nodeWeights.data(), nullptr, edgeNumber * 2, nodeNeighbors.data(), edgeWeights.data());
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

    std::vector<bool> result(partitions.size());
    std::transform(partitions.begin(), partitions.end(), result.begin(),
                   [](int64_t x) { return x != 0; });
    return result;
}

std::vector<bool> UndirectedFix::getUndirectedPartitionMetis() {
    int64_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_CONTIG] = 1;
    options[METIS_OPTION_NUMBERING] = 0;
    double imbalanceRatio = upperBoundPartWeight * 2 / (double) workingGraph.totalWeight;
    options[METIS_OPTION_UFACTOR] = (int64_t) ceil((imbalanceRatio - 1) * 1000);

    int64_t edgeNumber = computeNumberOfEdges();
    std::vector<int64_t> nodeNeighborsOffset(workingGraph.size + 1);
    std::vector<int64_t> nodeNeighbors(2 * edgeNumber);
    std::vector<int64_t> edgeWeights(2 * edgeNumber);
    std::vector<int64_t> nodeWeights(workingGraph.size);
    graphToCSRFormat(edgeNumber, nodeNeighborsOffset, nodeNeighbors, edgeWeights, nodeWeights);

    int64_t one = 1;
    int64_t two = 2;
    int64_t edgeCut;
    auto size = (int64_t)workingGraph.size;
    double targetWeights[2] = {0.5, 0.5};
    int64_t partitions[workingGraph.size];

    int ret = METIS_PartGraphKway(&size, &one, nodeNeighborsOffset.data(), nodeNeighbors.data(),
                              nodeWeights.data(), nullptr, edgeWeights.data(), &two, targetWeights,
                              nullptr, options, &edgeCut, partitions);
    assert(ret == METIS_OK && "metis partition failed");

    std::vector<bool> result(size);
    std::transform(partitions, partitions + size, result.begin(), [](int64_t x) { return x != 0; });
    return result;
}

std::pair<std::vector<bool>, uint64_t> UndirectedFix::run() const {
    return {};
}
