/**
 * @file MultilevelBisectioner.cpp
 * @brief Implementation of multilevel graph partitioning orchestration
 */

#include "MultilevelBisectioner.h"
#include "ClusteringCycleDetection.h"
#include "ClusteringForbiddenEdges.h"
#include "ClusteringHybrid.h"
#include "GreedyDirectedGraphGrowing.h"
#include "RefinementWrapper.h"
#include "UndirectedFix.h"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <utility>

namespace dag_partitioning {

namespace driver {

MultilevelBisectioner::MultilevelBisectioner(
    const core::Graph &graph, clustering::ClusteringMethod clusteringMethod,
    uint64_t maxClusteringRounds, uint64_t minClusteringVertices,
    double clusteringVertexRatio, bisection::BisectionMethod bisectionMethod,
    double imbalanceRatio, refinement::RefinementMethod refinementMethod,
    uint64_t refinementPasses)
    : workingGraph(graph), clusteringMethod(clusteringMethod),
      maxClusteringRounds(maxClusteringRounds),
      minClusteringVertices(minClusteringVertices),
      clusteringVertexRatio(clusteringVertexRatio),
      bisectionMethod(bisectionMethod), imbalanceRatio(imbalanceRatio),
      refinementMethod(refinementMethod), refinementPasses(refinementPasses) {}

std::stack<std::pair<core::Graph, std::vector<uint64_t>>>
MultilevelBisectioner::runClustering() const {
    // Create appropriate clustering algorithm based on selected method
    std::unique_ptr<clustering::Clustering> clustering;
    switch (clusteringMethod) {
    case clustering::ClusteringMethod::FORB:
        clustering = std::make_unique<clustering::ClusteringForbiddenEdges>(
            workingGraph, maxClusteringRounds, minClusteringVertices,
            clusteringVertexRatio);
        break;
    case clustering::ClusteringMethod::CYC:
        clustering = std::make_unique<clustering::ClusteringCycleDetection>(
            workingGraph, maxClusteringRounds, minClusteringVertices,
            clusteringVertexRatio);
        break;
    case clustering::ClusteringMethod::HYB:
        clustering = std::make_unique<clustering::ClusteringHybrid>(
            workingGraph, maxClusteringRounds, minClusteringVertices,
            clusteringVertexRatio);
        break;
    default:
        throw std::invalid_argument("Unknown clustering type");
    }

    // Execute coarsening and verify constraints
    std::stack<std::pair<core::Graph, std::vector<uint64_t>>>
        intermediateClusters = clustering->run();
    assert(intermediateClusters.size() <= maxClusteringRounds &&
           "maximum number of clustering rounds exceeded");
    return intermediateClusters;
}

std::pair<std::vector<uint8_t>, uint64_t>
MultilevelBisectioner::runBisection(const core::Graph &graph) const {
    // Calculate partition weight bounds based on imbalance ratio
    double lowerBoundPartWeight = 1.0;
    double upperBoundPartWeight =
        imbalanceRatio * ((double)graph.totalWeight / 2.0);

    // Create appropriate bisection algorithm
    std::unique_ptr<bisection::Bisection> bisection;
    switch (bisectionMethod) {
    case bisection::BisectionMethod::GGG:
        bisection = std::make_unique<bisection::GreedyDirectedGraphGrowing>(
            graph, upperBoundPartWeight, lowerBoundPartWeight, refinementMethod,
            refinementPasses);
        break;
    case bisection::BisectionMethod::UNDIRMETIS:
        bisection = std::make_unique<bisection::UndirectedFix>(
            graph, upperBoundPartWeight, lowerBoundPartWeight, refinementMethod,
            refinementPasses, true, false);
        break;
    case bisection::BisectionMethod::UNDIRSCOTCH:
        bisection = std::make_unique<bisection::UndirectedFix>(
            graph, upperBoundPartWeight, lowerBoundPartWeight, refinementMethod,
            refinementPasses, false, true);
        break;
    case bisection::BisectionMethod::UNDIRBOTH:
        bisection = std::make_unique<bisection::UndirectedFix>(
            graph, upperBoundPartWeight, lowerBoundPartWeight, refinementMethod,
            refinementPasses, true, true);
        break;
    default:
        throw std::invalid_argument("Unknown bisection type");
    }

    return bisection->run();
}

void MultilevelBisectioner::projectBisection(
    std::pair<std::vector<uint8_t>, uint64_t> &bisectionInfoPair,
    const std::vector<uint64_t> &mapping) {
    // Project partition assignments to finer level using mapping
    std::vector<uint8_t> newBisection(mapping.size());
    for (uint64_t i = 0; i < mapping.size(); ++i) {
        newBisection[i] = bisectionInfoPair.first[mapping[i]];
    }
    bisectionInfoPair.first = newBisection;
}

void MultilevelBisectioner::runRefinement(
    const core::Graph &graph,
    std::pair<std::vector<uint8_t>, uint64_t> &bisectionInfoPair) const {
    // Calculate partition weight bounds
    double lowerBoundPartWeight = 1.0;
    double upperBoundPartWeight =
        imbalanceRatio * ((double)graph.totalWeight / 2.0);

    refinement::refinementWrapper(graph, bisectionInfoPair.first,
                                  bisectionInfoPair.second, refinementMethod,
                                  refinementPasses, upperBoundPartWeight,
                                  lowerBoundPartWeight);
}

std::pair<std::vector<uint8_t>, uint64_t> MultilevelBisectioner::run() const {
    // Phase 1: Coarsening
    std::stack<std::pair<core::Graph, std::vector<uint64_t>>>
        intermediateClusters = runClustering();

    // If no coarsening occurred, partition original graph and refine it
    if (intermediateClusters.empty()) {
        std::pair<std::vector<uint8_t>, uint64_t> bisectionInfoPair =
            runBisection(workingGraph);
        runRefinement(workingGraph, bisectionInfoPair);
        return bisectionInfoPair;
    }

    // Phase 2: Initial Partitioning
    auto [coarsestGraph, coarsestMapping] = intermediateClusters.top();
    assert(coarsestGraph.size >= minClusteringVertices &&
           "minimum number of vertices for intermediate graphs violated");

    std::pair<std::vector<uint8_t>, uint64_t> bisectionInfoPair =
        runBisection(coarsestGraph);
    runRefinement(coarsestGraph, bisectionInfoPair);
    intermediateClusters.pop();

    // Phase 3: Uncoarsening and Refinement
    auto &currentMapping = coarsestMapping;
    while (!intermediateClusters.empty()) {
        // Project to next level (1 level coarser than the current)
        projectBisection(bisectionInfoPair, currentMapping);
        auto [intermediateGraph, intermediateMapping] =
            intermediateClusters.top();
        assert(intermediateGraph.size >= minClusteringVertices &&
               "minimum number of vertices for intermediate graphs violated");
        currentMapping = intermediateMapping;

        // Refine at current level
        runRefinement(intermediateGraph, bisectionInfoPair);
        intermediateClusters.pop();
    }

    // Final projection and refinement
    projectBisection(bisectionInfoPair, currentMapping);
    runRefinement(workingGraph, bisectionInfoPair);
    return bisectionInfoPair;
}

} // namespace driver

} // namespace dag_partitioning