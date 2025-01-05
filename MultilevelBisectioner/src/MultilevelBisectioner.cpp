/**
 * @file MultilevelBisectioner.cpp
 * @brief Implementation of multilevel graph partitioning orchestration
 */

#include "MultilevelBisectioner.h"
#include "ClusteringForbiddenEdges.h"
#include "ClusteringCycleDetection.h"
#include "ClusteringHybrid.h"
#include "GreedyDirectedGraphGrowing.h"
#include "UndirectedFix.h"
#include "BoundaryFM.h"
#include <utility>

MultilevelBisectioner::MultilevelBisectioner(Graph graph, ClusteringMethod clusteringMethod,
                                             uint64_t maxClusteringRounds,
                                             uint64_t minClusteringVertices,
                                             BisectionMethod bisectionMethod,
                                             double imbalanceRatio,
                                             RefinementMethod refinementMethod,
                                             uint64_t refinementPasses)
        : workingGraph(std::move(graph)),
          clusteringMethod(clusteringMethod),
          maxClusteringRounds(maxClusteringRounds),
          minClusteringVertices(minClusteringVertices),
          bisectionMethod(bisectionMethod),
          imbalanceRatio(imbalanceRatio),
          refinementMethod(refinementMethod),
          refinementPasses(refinementPasses) {}

std::stack<std::pair<Graph, std::vector<uint64_t>>> MultilevelBisectioner::runClustering() const {
    // Create appropriate clustering algorithm based on selected method
    std::unique_ptr<Clustering> clustering;
    switch (clusteringMethod) {
        case ClusteringMethod::FORB:
            clustering = std::make_unique<ClusteringForbiddenEdges>(
                    workingGraph, maxClusteringRounds, minClusteringVertices);
            break;
        case ClusteringMethod::CYC:
            clustering = std::make_unique<ClusteringCycleDetection>(
                    workingGraph, maxClusteringRounds, minClusteringVertices);
            break;
        case ClusteringMethod::HYB:
            clustering = std::make_unique<ClusteringHybrid>(
                    workingGraph, maxClusteringRounds, minClusteringVertices);
            break;
        default:
            throw std::invalid_argument("Unknown clustering type");
    }

    // Execute coarsening and verify constraints
    std::stack<std::pair<Graph, std::vector<uint64_t>>> intermediateClusters = clustering->run();
    assert(intermediateClusters.size() <= maxClusteringRounds &&
           "maximum number of clustering rounds exceeded");
    return intermediateClusters;
}

std::pair<std::vector<bool>, uint64_t> MultilevelBisectioner::runBisection(const Graph &graph) const {
    // Calculate partition weight bounds based on imbalance ratio
    double lowerBoundPartWeight = 1.0;
    double upperBoundPartWeight = imbalanceRatio * ((double) workingGraph.totalWeight / 2.0);

    // Create appropriate bisection algorithm
    std::unique_ptr<Bisection> bisection;
    switch (bisectionMethod) {
        case BisectionMethod::GGG:
            bisection = std::make_unique<GreedyDirectedGraphGrowing>(
                    graph, upperBoundPartWeight, lowerBoundPartWeight);
            break;
        case BisectionMethod::UNDIRMETIS:
            bisection = std::make_unique<UndirectedFix>(
                    graph, upperBoundPartWeight, lowerBoundPartWeight, true, false);
            break;
        case BisectionMethod::UNDIRSCOTCH:
            bisection = std::make_unique<UndirectedFix>(
                    graph, upperBoundPartWeight, lowerBoundPartWeight, false, true);
            break;
        case BisectionMethod::UNDIRBOTH:
            bisection = std::make_unique<UndirectedFix>(
                    graph, upperBoundPartWeight, lowerBoundPartWeight, true, true);
            break;
        default:
            throw std::invalid_argument("Unknown bisection type");
    }

    return bisection->run();
}

void MultilevelBisectioner::projectBisection(
        std::pair<std::vector<bool>, uint64_t> &bisectionInfo,
        const std::vector<uint64_t> &mapping) {
    // Project partition assignments to finer level using mapping
    std::vector<bool> newBisection(mapping.size());
    for (uint64_t i = 0; i < mapping.size(); ++i) {
        newBisection[i] = bisectionInfo.first[mapping[i]];
    }
    bisectionInfo.first = newBisection;
}

void MultilevelBisectioner::runRefinement(
        const Graph &graph,
        std::pair<std::vector<bool>, uint64_t> &bisectionInfo) const {
    // Calculate partition weight bounds
    double lowerBoundPartWeight = 1.0;
    double upperBoundPartWeight = imbalanceRatio * ((double) workingGraph.totalWeight / 2.0);

    // Create appropriate refinement algorithm
    std::unique_ptr<Refinement> refinement;
    switch (refinementMethod) {
        case RefinementMethod::BOUNDARYFM:
            refinement = std::make_unique<BoundaryFM>(
                    graph, bisectionInfo.first, bisectionInfo.second,
                    refinementPasses, upperBoundPartWeight, lowerBoundPartWeight);
            break;
        case RefinementMethod::BOUNDARYKL:
            throw std::invalid_argument("boundary kl not implemented yet");
            break;
        default:
            throw std::invalid_argument("Unknown refinement type");
    }
    refinement->run();
}

std::pair<std::vector<bool>, uint64_t> MultilevelBisectioner::run() const {
    // Phase 1: Coarsening
    std::stack<std::pair<Graph, std::vector<uint64_t>>> intermediateClusters = runClustering();

    // If no coarsening occurred, partition original graph and refine it
    if (intermediateClusters.empty()) {
        std::pair<std::vector<bool>, uint64_t> bisectionInfo = runBisection(workingGraph);
        runRefinement(workingGraph, bisectionInfo);
        return bisectionInfo;
    }

    // Phase 2: Initial Partitioning
    auto [coarsestGraph, coarsestMapping] = intermediateClusters.top();
    assert(coarsestGraph.size >= minClusteringVertices &&
           "minimum number of vertices for intermediate graphs violated");

    std::pair<std::vector<bool>, uint64_t> bisectionInfo = runBisection(coarsestGraph);
    runRefinement(coarsestGraph, bisectionInfo);
    intermediateClusters.pop();

    // Phase 3: Uncoarsening and Refinement
    auto &currentMapping = coarsestMapping;
    while (!intermediateClusters.empty()) {
        // Project to next level (1 level coarser than the current)
        projectBisection(bisectionInfo, currentMapping);
        auto [intermediateGraph, intermediateMapping] = intermediateClusters.top();
        assert(intermediateGraph.size >= minClusteringVertices &&
               "minimum number of vertices for intermediate graphs violated");
        currentMapping = intermediateMapping;

        // Refine at current level
        runRefinement(intermediateGraph, bisectionInfo);
        intermediateClusters.pop();
    }

    // Final projection and refinement
    projectBisection(bisectionInfo, currentMapping);
    runRefinement(workingGraph, bisectionInfo);
    return bisectionInfo;
}