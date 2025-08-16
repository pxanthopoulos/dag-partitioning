/**
 * @file UndirectedFix.cpp
 * @brief Implementation of undirected bisection with acyclicity fixing
 */

#include "UndirectedFix.h"
#include "BoundaryFM.h"
#include "Graph.h"
#include "Refinement.h"
#include "RefinementWrapper.h"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace dag_partitioning {

namespace bisection {

UndirectedFix::UndirectedFix(const core::Graph &graph,
                             double upperBoundPartWeight,
                             double lowerBoundPartWeight,
                             refinement::RefinementMethod refinementMethod,
                             uint64_t refinementPasses, bool useMetis,
                             bool useScotch)
    : Bisection(graph, upperBoundPartWeight, lowerBoundPartWeight,
                refinementMethod, refinementPasses),
      useMetis(useMetis), useScotch(useScotch) {
    // Require at least one partitioning method
    assert((useMetis || useScotch) &&
           "Specify at least 1 undirected partitioning algorithm");
}

// Count edges in both directions for CSR format
int64_t UndirectedFix::computeNumberOfEdges() const {
    int64_t edgeNumber = 0;
    for (const auto &neighbors : workingGraph.adj)
        edgeNumber += (int64_t)neighbors.size();
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
        nodeWeights[i] = (int64_t)workingGraph.nodeWeights[i];

        // Add forward edges
        const auto &neighbors = workingGraph.adj[i];
        for (const auto &[neighborId, edgeWeight] : neighbors) {
            assert((pos < 2 * edgeNumber) &&
                   "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = (int64_t)neighborId;
            edgeWeights[pos] = (int64_t)edgeWeight;
            pos++;
        }

        // Add reverse edges to make graph undirected
        const auto &reverseNeighbors = workingGraph.revAdj[i];
        for (const auto &[reverseNeighborId, reverseEdgeWeight] :
             reverseNeighbors) {
            assert((pos < 2 * edgeNumber) &&
                   "pos must be inside the [0, 2 * #Edges) interval");
            nodeNeighbors[pos] = (int64_t)reverseNeighborId;
            edgeWeights[pos] = (int64_t)reverseEdgeWeight;
            pos++;
        }

        cumulativeNeighbors +=
            (int64_t)(neighbors.size() + reverseNeighbors.size());
    }

    nodeNeighborsOffset[workingGraph.size] = cumulativeNeighbors;
    assert(pos == 2 * edgeNumber && "pos must be equal to 2 * #Edges");
}

// Get initial bisection using Scotch partitioner
std::vector<uint8_t> UndirectedFix::getUndirectedBisectionScotch() const {
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
    graphToCSRFormat(edgeNumber, nodeNeighborsOffset, nodeNeighbors,
                     edgeWeights, nodeWeights);

    // Build and verify Scotch graph
    int ret = SCOTCH_graphBuild(
        &scotchGraph, 0, (int64_t)workingGraph.size, nodeNeighborsOffset.data(),
        nodeNeighborsOffset.data() + 1, nodeWeights.data(), nullptr,
        edgeNumber * 2, nodeNeighbors.data(), edgeWeights.data());
    assert(ret == 0 && "scotch graph build failed");

    ret = SCOTCH_graphCheck(&scotchGraph);
    assert(ret == 0 && "scotch graph check failed");

    // Set partitioning strategy with balance constraint
    double imbalanceRatio =
        upperBoundPartWeight * 2 / (double)workingGraph.totalWeight;
    ret = SCOTCH_stratGraphMapBuild(&scotchStrategy, SCOTCH_STRATQUALITY, 2,
                                    imbalanceRatio - 1);
    assert(ret == 0 && "scotch strategy build failed");

    // Compute partitioning
    std::vector<int64_t> bisectionData(workingGraph.size);
    ret = SCOTCH_graphPart(&scotchGraph, 2, &scotchStrategy,
                           bisectionData.data());
    assert(ret == 0 && "scotch partition failed");

    // Clean up
    SCOTCH_graphExit(&scotchGraph);
    SCOTCH_stratExit(&scotchStrategy);

    // Convert to boolean representation
    std::vector<uint8_t> result(bisectionData.size());
    std::transform(bisectionData.begin(), bisectionData.end(), result.begin(),
                   [](int64_t x) { return (uint8_t)x; });
    return result;
}

// Get initial bisection using METIS partitioner
std::vector<uint8_t> UndirectedFix::getUndirectedBisectionMetis() const {
    // Configure METIS options
    int64_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;

    // Set balance constraint
    double imbalanceRatio =
        upperBoundPartWeight * 2 / (double)workingGraph.totalWeight;
    options[METIS_OPTION_UFACTOR] = (int64_t)ceil((imbalanceRatio - 1) * 1000);

    // Create CSR format data
    int64_t edgeNumber = computeNumberOfEdges();
    std::vector<int64_t> nodeNeighborsOffset(workingGraph.size + 1);
    std::vector<int64_t> nodeNeighbors(2 * edgeNumber);
    std::vector<int64_t> edgeWeights(2 * edgeNumber);
    std::vector<int64_t> nodeWeights(workingGraph.size);
    graphToCSRFormat(edgeNumber, nodeNeighborsOffset, nodeNeighbors,
                     edgeWeights, nodeWeights);

    // Set up METIS parameters
    int64_t one = 1;
    int64_t two = 2;
    int64_t edgeCut;
    auto size = (int64_t)workingGraph.size;
    double targetWeights[2] = {0.5, 0.5};
    int64_t bisectionData[workingGraph.size];

    // Compute partitioning
    int ret = METIS_PartGraphKway(
        &size, &one, nodeNeighborsOffset.data(), nodeNeighbors.data(),
        nodeWeights.data(), nullptr, edgeWeights.data(), &two, targetWeights,
        nullptr, options, &edgeCut, bisectionData);
    assert(ret == METIS_OK && "metis partition failed");

    // Convert to boolean representation
    std::vector<uint8_t> result(size);
    std::transform(bisectionData, bisectionData + size, result.begin(),
                   [](int64_t x) { return (uint8_t)x; });
    return result;
}

// Fix cycles by moving ancestors to V0 (Algorithm 3)
void UndirectedFix::fixAcyclicityUp(
    std::vector<uint8_t> &undirectedBisection) const {
    std::vector<uint64_t> topologicalOrder = workingGraph.topologicalSort();

    // Process vertices in reverse topological order
    for (auto it = topologicalOrder.rbegin(); it != topologicalOrder.rend();
         ++it) {
        const uint64_t &nodeId = *it;
        if (undirectedBisection[nodeId] == 0) { // If node is in V0
            // Move all predecessors to V0
            for (const auto &[predecessorId, edgeWeight] :
                 workingGraph.revAdj[nodeId])
                undirectedBisection[predecessorId] = 0;
        }
    }
}

// Fix cycles by moving descendants to V1 (Algorithm 4)
void UndirectedFix::fixAcyclicityDown(
    std::vector<uint8_t> &undirectedBisection) const {
    std::vector<uint64_t> topologicalOrder = workingGraph.topologicalSort();

    // Process vertices in topological order
    for (const uint64_t &nodeId : topologicalOrder) {
        if (undirectedBisection[nodeId] == 1) { // If node is in V1
            // Move all successors to V1
            for (const auto &[successorId, edgeWeight] :
                 workingGraph.adj[nodeId])
                undirectedBisection[successorId] = 1;
        }
    }
}

std::pair<std::vector<uint8_t>, uint64_t> UndirectedFix::runMetis() const {
    std::vector<std::pair<std::vector<uint8_t>, uint64_t>> bisections;
    std::vector<std::tuple<uint64_t, uint8_t, double>> results;

    // Get initial undirected bisection
    std::vector<uint8_t> undirectedBisection = getUndirectedBisectionMetis();
    std::vector<uint8_t> reverseUndirectedBisection(undirectedBisection.size());

    // Create reversed assignment
    std::transform(undirectedBisection.begin(), undirectedBisection.end(),
                   reverseUndirectedBisection.begin(),
                   [](uint8_t x) { return 1 - x; });
    std::vector<uint8_t> copy = undirectedBisection;

    // Try Approach 1: Fix down from original assignment (always store result as
    // best)
    fixAcyclicityDown(undirectedBisection);
    uint64_t edgeCut = computeEdgeCut(undirectedBisection, workingGraph);
    // One round of FM
    refinement::BoundaryFM boundaryFM1 =
        refinement::BoundaryFM(workingGraph, undirectedBisection, edgeCut, 1,
                               upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM1.run();
    // One refinement call
    refinement::refinementWrapper(workingGraph, undirectedBisection, edgeCut,
                                  refinementMethod, refinementPasses,
                                  upperBoundPartWeight, lowerBoundPartWeight);
    edgeCut = computeEdgeCut(undirectedBisection, workingGraph);
    uint8_t isZero = (edgeCut == 0) ? 0 : 1;
    auto [sizeV0, sizeV1] = refinement::Refinement::calculatePartSizes(
        undirectedBisection, workingGraph);
    double imbalance =
        (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
         upperBoundPartWeight) *
        100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(undirectedBisection, edgeCut);

    // Try Approach 2: Fix up from original assignment
    fixAcyclicityUp(copy);
    edgeCut = computeEdgeCut(copy, workingGraph);
    refinement::BoundaryFM boundaryFM2 =
        refinement::BoundaryFM(workingGraph, copy, edgeCut, 1,
                               upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM2.run();
    refinement::refinementWrapper(workingGraph, copy, edgeCut, refinementMethod,
                                  refinementPasses, upperBoundPartWeight,
                                  lowerBoundPartWeight);
    edgeCut = computeEdgeCut(copy, workingGraph);
    isZero = (edgeCut == 0) ? 0 : 1;
    std::tie(sizeV0, sizeV1) =
        refinement::Refinement::calculatePartSizes(copy, workingGraph);
    imbalance = (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
                 upperBoundPartWeight) *
                100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(copy, edgeCut);

    // Try Approach 3: Fix down from reversed assignment
    copy = reverseUndirectedBisection;
    fixAcyclicityDown(reverseUndirectedBisection);
    edgeCut = computeEdgeCut(reverseUndirectedBisection, workingGraph);
    refinement::BoundaryFM boundaryFM3 = refinement::BoundaryFM(
        workingGraph, reverseUndirectedBisection, edgeCut, 1,
        upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM3.run();
    refinement::refinementWrapper(workingGraph, reverseUndirectedBisection,
                                  edgeCut, refinementMethod, refinementPasses,
                                  upperBoundPartWeight, lowerBoundPartWeight);
    edgeCut = computeEdgeCut(reverseUndirectedBisection, workingGraph);
    isZero = (edgeCut == 0) ? 0 : 1;
    std::tie(sizeV0, sizeV1) = refinement::Refinement::calculatePartSizes(
        reverseUndirectedBisection, workingGraph);
    imbalance = (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
                 upperBoundPartWeight) *
                100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(reverseUndirectedBisection, edgeCut);

    // Try Approach 4: Fix up from reversed assignment
    fixAcyclicityUp(copy);
    edgeCut = computeEdgeCut(copy, workingGraph);
    refinement::BoundaryFM boundaryFM4 =
        refinement::BoundaryFM(workingGraph, copy, edgeCut, 1,
                               upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM4.run();
    refinement::refinementWrapper(workingGraph, copy, edgeCut, refinementMethod,
                                  refinementPasses, upperBoundPartWeight,
                                  lowerBoundPartWeight);

    edgeCut = computeEdgeCut(copy, workingGraph);
    isZero = (edgeCut == 0) ? 0 : 1;
    std::tie(sizeV0, sizeV1) =
        refinement::Refinement::calculatePartSizes(copy, workingGraph);
    imbalance = (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
                 upperBoundPartWeight) *
                100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(copy, edgeCut);

    auto [bestBisection, bestEdgeCut] = bisections[selectBestResult(results)];

    assert(refinement::Refinement::checkValidBisection(bestBisection,
                                                       workingGraph) &&
           "Best bisection is cyclic");
    assert(refinement::Refinement::checkValidEdgeCut(
               bestBisection, workingGraph, bestEdgeCut) &&
           "Computed edge cut is invalid");
    return {bestBisection, bestEdgeCut};
}

std::pair<std::vector<uint8_t>, uint64_t> UndirectedFix::runScotch() const {
    std::vector<std::pair<std::vector<uint8_t>, uint64_t>> bisections;
    std::vector<std::tuple<uint64_t, uint8_t, double>> results;

    // Get initial undirected bisection
    std::vector<uint8_t> undirectedBisection = getUndirectedBisectionScotch();
    std::vector<uint8_t> reverseUndirectedBisection(undirectedBisection.size());

    // Create reversed assignment
    std::transform(undirectedBisection.begin(), undirectedBisection.end(),
                   reverseUndirectedBisection.begin(),
                   [](uint8_t x) { return 1 - x; });
    std::vector<uint8_t> copy = undirectedBisection;

    // Try Approach 1: Fix down from original assignment (always store result as
    // best)
    fixAcyclicityDown(undirectedBisection);
    uint64_t edgeCut = computeEdgeCut(undirectedBisection, workingGraph);
    // One round of FM
    refinement::BoundaryFM boundaryFM1 =
        refinement::BoundaryFM(workingGraph, undirectedBisection, edgeCut, 1,
                               upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM1.run();
    // One refinement call
    refinement::refinementWrapper(workingGraph, undirectedBisection, edgeCut,
                                  refinementMethod, refinementPasses,
                                  upperBoundPartWeight, lowerBoundPartWeight);
    edgeCut = computeEdgeCut(undirectedBisection, workingGraph);
    uint8_t isZero = (edgeCut == 0) ? 0 : 1;
    auto [sizeV0, sizeV1] = refinement::Refinement::calculatePartSizes(
        undirectedBisection, workingGraph);
    double imbalance =
        (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
         upperBoundPartWeight) *
        100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(undirectedBisection, edgeCut);

    // Try Approach 2: Fix up from original assignment
    fixAcyclicityUp(copy);
    edgeCut = computeEdgeCut(copy, workingGraph);
    refinement::BoundaryFM boundaryFM2 =
        refinement::BoundaryFM(workingGraph, copy, edgeCut, 1,
                               upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM2.run();
    refinement::refinementWrapper(workingGraph, copy, edgeCut, refinementMethod,
                                  refinementPasses, upperBoundPartWeight,
                                  lowerBoundPartWeight);
    edgeCut = computeEdgeCut(copy, workingGraph);
    isZero = (edgeCut == 0) ? 0 : 1;
    std::tie(sizeV0, sizeV1) =
        refinement::Refinement::calculatePartSizes(copy, workingGraph);
    imbalance = (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
                 upperBoundPartWeight) *
                100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(copy, edgeCut);

    // Try Approach 3: Fix down from reversed assignment
    copy = reverseUndirectedBisection;
    fixAcyclicityDown(reverseUndirectedBisection);
    edgeCut = computeEdgeCut(reverseUndirectedBisection, workingGraph);
    refinement::BoundaryFM boundaryFM3 = refinement::BoundaryFM(
        workingGraph, reverseUndirectedBisection, edgeCut, 1,
        upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM3.run();
    refinement::refinementWrapper(workingGraph, reverseUndirectedBisection,
                                  edgeCut, refinementMethod, refinementPasses,
                                  upperBoundPartWeight, lowerBoundPartWeight);
    edgeCut = computeEdgeCut(reverseUndirectedBisection, workingGraph);
    isZero = (edgeCut == 0) ? 0 : 1;
    std::tie(sizeV0, sizeV1) = refinement::Refinement::calculatePartSizes(
        reverseUndirectedBisection, workingGraph);
    imbalance = (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
                 upperBoundPartWeight) *
                100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(reverseUndirectedBisection, edgeCut);

    // Try Approach 4: Fix up from reversed assignment
    fixAcyclicityUp(copy);
    edgeCut = computeEdgeCut(copy, workingGraph);
    refinement::BoundaryFM boundaryFM4 =
        refinement::BoundaryFM(workingGraph, copy, edgeCut, 1,
                               upperBoundPartWeight, lowerBoundPartWeight);
    boundaryFM4.run();
    refinement::refinementWrapper(workingGraph, copy, edgeCut, refinementMethod,
                                  refinementPasses, upperBoundPartWeight,
                                  lowerBoundPartWeight);

    edgeCut = computeEdgeCut(copy, workingGraph);
    isZero = (edgeCut == 0) ? 0 : 1;
    std::tie(sizeV0, sizeV1) =
        refinement::Refinement::calculatePartSizes(copy, workingGraph);
    imbalance = (((double)std::max(sizeV0, sizeV1) - upperBoundPartWeight) /
                 upperBoundPartWeight) *
                100;
    results.emplace_back(edgeCut, isZero, imbalance);
    bisections.emplace_back(copy, edgeCut);

    auto [bestBisection, bestEdgeCut] = bisections[selectBestResult(results)];

    assert(refinement::Refinement::checkValidBisection(bestBisection,
                                                       workingGraph) &&
           "Best bisection is cyclic");
    assert(refinement::Refinement::checkValidEdgeCut(
               bestBisection, workingGraph, bestEdgeCut) &&
           "Computed edge cut is invalid");
    return {bestBisection, bestEdgeCut};
}

// Main algorithm that tries all fixing combinations
std::pair<std::vector<uint8_t>, uint64_t> UndirectedFix::run() const {
    if (!useMetis)
        return runScotch();
    if (!useScotch)
        return runMetis();

    // Both enabled - compare results
    const auto resultMetis = runMetis();
    const auto resultScotch = runScotch();

    return (resultMetis.second == 0)    ? std::move(resultScotch)
           : (resultScotch.second == 0) ? std::move(resultMetis)
           : (resultMetis.second < resultScotch.second)
               ? std::move(resultMetis)
               : std::move(resultScotch);
}

} // namespace bisection

} // namespace dag_partitioning