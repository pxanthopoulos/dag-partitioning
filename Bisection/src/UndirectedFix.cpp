/**
 * @file UndirectedFix.cpp
 * @brief Implementation of undirected bisection with acyclicity fixing
 */

#include <cmath>
#include <llvm/ADT/STLExtras.h>
#include "UndirectedFix.h"

UndirectedFix::UndirectedFix(const Graph &graph, double upperBoundPartWeight,
                             double lowerBoundPartWeight, bool useMetis, bool useScotch)
        : Bisection(graph, upperBoundPartWeight, lowerBoundPartWeight),
          useMetis(useMetis),
          useScotch(useScotch) {
    // Require at least one partitioning method
    assert((useMetis || useScotch) && "Specify at least 1 undirected partitioning algorithm");
}

// Count edges in both directions for CSR format
int64_t UndirectedFix::computeNumberOfEdges() const {
    int64_t edgeNumber = 0;
    for (const auto &neighbors: workingGraph.adj)
        edgeNumber += (int64_t) neighbors.size();
    return edgeNumber;
}

// Convert DAG to undirected CSR format for external partitioners
void UndirectedFix::graphToCSRFormat(int64_t edgeNumber,
                                     std::vector<int64_t> &nodeNeighborsOffset,
                                     std::vector<int64_t> &nodeNeighbors,
                                     std::vector<int64_t> &edgeWeights,
                                     std::vector<int64_t> &nodeWeights) const {
    int64_t pos = 0;
    int64_t cumulativeNeighbors = 0;

    // Process each vertex
    for (uint64_t i = 0; i < workingGraph.size; ++i) {
        nodeNeighborsOffset[i] = cumulativeNeighbors;
        nodeWeights[i] = (int64_t) workingGraph.nodeWeights[i];

        // Add forward edges
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, edgeWeight]: neighbors) {
            assert((pos < 2 * edgeNumber) && "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = (int64_t) neighborId;
            edgeWeights[pos] = (int64_t) edgeWeight;
            pos++;
        }

        // Add reverse edges to make graph undirected
        const auto &reverseNeighbors = workingGraph.revAdj[i];
        for (const auto &[reverseNeighborId, reverseEdgeWeight]: reverseNeighbors) {
            assert((pos < 2 * edgeNumber) && "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = (int64_t) reverseNeighborId;
            edgeWeights[pos] = (int64_t) reverseEdgeWeight;
            pos++;
        }

        cumulativeNeighbors += (int64_t) (neighbors.size() + reverseNeighbors.size());
    }

    nodeNeighborsOffset[workingGraph.size] = cumulativeNeighbors;
    assert(pos == 2 * edgeNumber && "pos must be equal to 2 * #Edges");
}

// Get initial bisection using Scotch partitioner
std::vector<bool> UndirectedFix::getUndirectedBisectionScotch() const {
    // Initialize Scotch structures
    SCOTCH_Graph scotchGraph;
    SCOTCH_Strat scotchStrategy;
    SCOTCH_graphInit(&scotchGraph);
    SCOTCH_stratInit(&scotchStrategy);

    // Create CSR format data
    int64_t edgeNumber = computeNumberOfEdges();
    std::vector<int64_t> nodeNeighborsOffset(workingGraph.size + 1);
    std::vector<int64_t> nodeNeighbors(2 * edgeNumber);
    std::vector<int64_t> edgeWeights(2 * edgeNumber);
    std::vector<int64_t> nodeWeights(workingGraph.size);
    graphToCSRFormat(edgeNumber, nodeNeighborsOffset, nodeNeighbors, edgeWeights, nodeWeights);

    // Build and verify Scotch graph
    int ret = SCOTCH_graphBuild(&scotchGraph, 0, (int64_t) workingGraph.size,
                                nodeNeighborsOffset.data(), nodeNeighborsOffset.data() + 1,
                                nodeWeights.data(), nullptr, edgeNumber * 2,
                                nodeNeighbors.data(), edgeWeights.data());
    assert(ret == 0 && "scotch graph build failed");

    ret = SCOTCH_graphCheck(&scotchGraph);
    assert(ret == 0 && "scotch graph check failed");

    // Set partitioning strategy with balance constraint
    double imbalanceRatio = upperBoundPartWeight * 2 / (double) workingGraph.totalWeight;
    ret = SCOTCH_stratGraphMapBuild(&scotchStrategy,
                                    SCOTCH_STRATQUALITY | SCOTCH_STRATBALANCE,
                                    2, imbalanceRatio - 1);
    assert(ret == 0 && "scotch strategy build failed");

    // Compute partitioning
    std::vector<int64_t> bisectionData(workingGraph.size);
    ret = SCOTCH_graphPart(&scotchGraph, 2, &scotchStrategy, bisectionData.data());
    assert(ret == 0 && "scotch partition failed");

    // Clean up
    SCOTCH_graphExit(&scotchGraph);
    SCOTCH_stratExit(&scotchStrategy);

    // Convert to boolean representation
    std::vector<bool> result(bisectionData.size());
    std::transform(bisectionData.begin(), bisectionData.end(), result.begin(),
                   [](int64_t x) { return x != 0; });
    return result;
}

// Get initial bisection using METIS partitioner
std::vector<bool> UndirectedFix::getUndirectedBisectionMetis() const {
    // Configure METIS options
    int64_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;

    // Set balance constraint
    double imbalanceRatio = upperBoundPartWeight * 2 / (double) workingGraph.totalWeight;
    options[METIS_OPTION_UFACTOR] = (int64_t) ceil((imbalanceRatio - 1) * 1000);

    // Create CSR format data
    int64_t edgeNumber = computeNumberOfEdges();
    std::vector<int64_t> nodeNeighborsOffset(workingGraph.size + 1);
    std::vector<int64_t> nodeNeighbors(2 * edgeNumber);
    std::vector<int64_t> edgeWeights(2 * edgeNumber);
    std::vector<int64_t> nodeWeights(workingGraph.size);
    graphToCSRFormat(edgeNumber, nodeNeighborsOffset, nodeNeighbors, edgeWeights, nodeWeights);

    // Set up METIS parameters
    int64_t one = 1;
    int64_t two = 2;
    int64_t edgeCut;
    auto size = (int64_t) workingGraph.size;
    double targetWeights[2] = {0.5, 0.5};
    int64_t bisectionData[workingGraph.size];

    // Compute partitioning
    int ret = METIS_PartGraphKway(&size, &one, nodeNeighborsOffset.data(), nodeNeighbors.data(),
                                  nodeWeights.data(), nullptr, edgeWeights.data(), &two, targetWeights,
                                  nullptr, options, &edgeCut, bisectionData);
    assert(ret == METIS_OK && "metis partition failed");

    // Convert to boolean representation
    std::vector<bool> result(size);
    std::transform(bisectionData, bisectionData + size, result.begin(),
                   [](int64_t x) { return x != 0; });
    return result;
}

// Fix cycles by moving ancestors to V0 (Algorithm 3)
void UndirectedFix::fixAcyclicityUp(std::vector<bool> &undirectedBisection) const {
    std::vector<uint64_t> topologicalOrder = workingGraph.topologicalSort();

    // Process vertices in reverse topological order
    for (const uint64_t &nodeId: llvm::reverse(topologicalOrder)) {
        if (!undirectedBisection[nodeId]) {  // If node is in V0
            // Move all predecessors to V0
            for (const auto &[predecessorId, edgeWeight]: workingGraph.revAdj[nodeId])
                undirectedBisection[predecessorId] = false;
        }
    }
}

// Fix cycles by moving descendants to V1 (Algorithm 4)
void UndirectedFix::fixAcyclicityDown(std::vector<bool> &undirectedBisection) const {
    std::vector<uint64_t> topologicalOrder = workingGraph.topologicalSort();

    // Process vertices in topological order
    for (const uint64_t &nodeId: topologicalOrder) {
        if (undirectedBisection[nodeId]) {  // If node is in V1
            // Move all successors to V1
            for (const auto &[successorId, edgeWeight]: workingGraph.adj[nodeId])
                undirectedBisection[successorId] = true;
        }
    }
}

// Main algorithm that tries all fixing combinations
std::pair<std::vector<bool>, uint64_t> UndirectedFix::run() const {
    uint64_t bestEdgeCut = UINT64_MAX, edgeCut;
    std::vector<bool> bestBisection;

    // Try METIS partitioning if enabled
    if (useMetis) {
        // Get initial undirected bisection
        std::vector<bool> undirectedBisection = getUndirectedBisectionMetis();
        std::vector<bool> reverseUndirectedBisection(undirectedBisection.size());

        // Create reversed assignment
        std::transform(undirectedBisection.begin(), undirectedBisection.end(),
                       reverseUndirectedBisection.begin(), [](bool x) { return !x; });
        std::vector<bool> copy = undirectedBisection;

        // Try Approach 1: Fix up from original assignment
        fixAcyclicityUp(undirectedBisection);
        edgeCut = computeEdgeCut(undirectedBisection);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = undirectedBisection;
        }

        // Try Approach 2: Fix down from original assignment
        fixAcyclicityDown(copy);
        edgeCut = computeEdgeCut(copy);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = copy;
        }

        // Try Approach 3: Fix up from reversed assignment
        copy = reverseUndirectedBisection;
        fixAcyclicityUp(reverseUndirectedBisection);
        edgeCut = computeEdgeCut(reverseUndirectedBisection);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = reverseUndirectedBisection;
        }

        // Try Approach 4: Fix down from reversed assignment
        fixAcyclicityDown(copy);
        edgeCut = computeEdgeCut(copy);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = copy;
        }
    }

    // Repeat process with Scotch if enabled
    if (useScotch) {
        std::vector<bool> undirectedBisection = getUndirectedBisectionScotch();
        std::vector<bool> reverseUndirectedBisection(undirectedBisection.size());
        std::transform(undirectedBisection.begin(), undirectedBisection.end(), reverseUndirectedBisection.begin(),
                       [](bool x) { return !x; });
        std::vector<bool> copy = undirectedBisection;

        fixAcyclicityUp(undirectedBisection);
        edgeCut = computeEdgeCut(undirectedBisection);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = undirectedBisection;
        }

        fixAcyclicityDown(copy);
        edgeCut = computeEdgeCut(copy);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = copy;
        }

        copy = reverseUndirectedBisection;
        fixAcyclicityUp(reverseUndirectedBisection);
        edgeCut = computeEdgeCut(reverseUndirectedBisection);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = reverseUndirectedBisection;
        }

        fixAcyclicityDown(copy);
        edgeCut = computeEdgeCut(copy);
        if (edgeCut < bestEdgeCut && edgeCut != 0) {
            bestEdgeCut = edgeCut;
            bestBisection = copy;
        }
    }
    assert(bestEdgeCut != UINT64_MAX && "No valid edgeCut computed");
    assert(checkValidBisection(bestBisection) && "Best bisection is cyclic");
    return {bestBisection, bestEdgeCut};
}