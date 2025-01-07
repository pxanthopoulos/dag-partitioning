#include "Graph.h"
#include "MultilevelBisectioner.h"
#include <iostream>
#include <chrono>

#include <vector>
#include <set>
#include <cmath>
#include <iostream>

std::tuple<bool, uint64_t, uint64_t, uint64_t>
checkValidBisectionAndCompute(const std::vector<bool> &bisection, const Graph &graph) {
    uint64_t edgeCut = 0, sizeV0 = 0, sizeV1 = 0;
    for (uint64_t i = 0; i < graph.size; ++i) {
        const auto &neighbors = graph.adj[i];
        for (const auto &[neighborId, edgeWeight]: neighbors) {
            if (bisection[i] && !bisection[neighborId]) return {false, edgeCut, sizeV0, sizeV1};
            if (bisection[i] != bisection[neighborId]) edgeCut += edgeWeight;
        }
        if (!bisection[i]) sizeV0 += graph.nodeWeights[i];
        else sizeV1 += graph.nodeWeights[i];
    }
    return {true, edgeCut, sizeV0, sizeV1};
}

std::vector<std::pair<uint64_t, double>> generate_bisections(const Graph &graph) {
    std::vector<std::pair<uint64_t, double>> pairs;

    for (uint64_t i = 0; i < (1 << graph.size); i++) {
        std::vector<bool> bisection;
        uint64_t mask = i;
        for (int j = 0; j < graph.size; j++) {
            bisection.push_back(mask & 1);
            mask >>= 1;
        }
        const auto &res = checkValidBisectionAndCompute(bisection, graph);
        if (std::get<0>(res)) {
            uint64_t sizeV0 = std::get<2>(res);
            uint64_t sizeV1 = std::get<3>(res);
            uint64_t maxNodeWeight = graph.maxNodeWeight();
            double lowerBoundPartWeight = 1.0;
            double upperBoundPartWeight = 1.2 * ((double) graph.totalWeight / 2.0);
            if (!((double) sizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
                  (double) sizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
                  (double) sizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
                  (double) sizeV1 > upperBoundPartWeight + (double) maxNodeWeight)) {
                uint64_t edgeCut = std::get<1>(res);
                double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
                imbalance = round(imbalance * 10) / 10;
                pairs.emplace_back(edgeCut, imbalance);
            }
        }
    }

    sort(pairs.begin(), pairs.end(), [](auto &a, auto &b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    });

    return pairs;
}

std::pair<uint64_t, uint64_t> calculatePartSizes(const std::vector<bool> &bisection, const Graph &graph) {
    uint64_t sizeV0 = 0, sizeV1 = 0;
    for (uint64_t i = 0; i < bisection.size(); ++i) {
        if (!bisection[i]) sizeV0 += graph.nodeWeights[i];
        else sizeV1 += graph.nodeWeights[i];
    }
    return {sizeV0, sizeV1};
}

void printResults(const std::vector<bool> &partition, const Graph &graph,
                  uint64_t cutSize, const std::string &methodStr) {
    auto [sizeV0, sizeV1] = calculatePartSizes(partition, graph);
    double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
    imbalance = round(imbalance * 10) / 10;
    std::cout << methodStr << "," << cutSize << "," << imbalance << "\n";
}

struct MethodCombination {
    ClusteringMethod clustering;
    BisectionMethod bisection;
    RefinementMethod refinement;
    std::string clusteringStr;
    std::string bisectionStr;
    std::string refinementStr;
};

int main() {
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/2mm_10_20_30_40.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/test/dag.dot",
                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");

    /*
    auto pairs = generate_bisections(graph);
    auto it = lower_bound(pairs.begin(), pairs.end(),
                          std::make_pair(vec.second, imbalance),
                          [](auto &a, auto &b) {
                              if (a.first != b.first) return a.first < b.first;
                              return a.second < b.second;
                          });
    long location = it - pairs.begin();
    */

    std::vector<MethodCombination> methods = {
            {ClusteringMethod::FORB, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYFM, "FORB", "UNDIR", "FM"},
            {ClusteringMethod::CYC, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYFM, "CYC", "UNDIR", "FM"},
            {ClusteringMethod::HYB, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYFM, "HYB", "UNDIR", "FM"},
            {ClusteringMethod::FORB, BisectionMethod::GGG, RefinementMethod::BOUNDARYFM, "FORB", "GGG", "FM"},
            {ClusteringMethod::CYC, BisectionMethod::GGG, RefinementMethod::BOUNDARYFM, "CYC", "GGG", "FM"},
            {ClusteringMethod::HYB, BisectionMethod::GGG, RefinementMethod::BOUNDARYFM, "HYB", "GGG", "FM"},
            {ClusteringMethod::FORB, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYKL, "FORB", "UNDIR", "KL"},
            {ClusteringMethod::CYC, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYKL, "CYC", "UNDIR", "KL"},
            {ClusteringMethod::HYB, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYKL, "HYB", "UNDIR", "KL"},
            {ClusteringMethod::FORB, BisectionMethod::GGG, RefinementMethod::BOUNDARYKL, "FORB", "GGG", "KL"},
            {ClusteringMethod::CYC, BisectionMethod::GGG, RefinementMethod::BOUNDARYKL, "CYC", "GGG", "KL"},
            {ClusteringMethod::HYB, BisectionMethod::GGG, RefinementMethod::BOUNDARYKL, "HYB", "GGG", "KL"}
    };

    // Common parameters for all combinations
    const int maxLevel = 20;
    const int minSize = 1;
    const double maxImbalance = 1.2;
    const int maxPasses = 10;


    // Main testing loop
    for (const auto &method: methods) {
        MultilevelBisectioner bisectioner(
                graph,
                method.clustering,
                maxLevel,
                minSize,
                method.bisection,
                maxImbalance,
                method.refinement,
                maxPasses
        );

        auto [partition, cutSize] = bisectioner.run();
        std::string methodStr = method.clusteringStr + "," +
                                method.bisectionStr + "," +
                                method.refinementStr;
        printResults(partition, graph, cutSize, methodStr);
    }

    return 0;
}
