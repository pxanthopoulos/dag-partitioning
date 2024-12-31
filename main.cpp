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

int main() {
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/2mm_10_20_30_40.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/test/dag.dot",
                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");

//    auto pairs = generate_bisections(graph);

    MultilevelBisectioner mb_forb_undir(graph, ClusteringMethod::FORB, 20, 1, BisectionMethod::UNDIRBOTH, 1.2,
                                        RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_cyc_undir(graph, ClusteringMethod::CYC, 20, 1, BisectionMethod::UNDIRBOTH, 1.2,
                                       RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_hyb_undir(graph, ClusteringMethod::HYB, 20, 1, BisectionMethod::UNDIRBOTH, 1.2,
                                       RefinementMethod::BOUNDARYFM, 10);

    MultilevelBisectioner mb_forb_ggg(graph, ClusteringMethod::FORB, 20, 1, BisectionMethod::GGG, 1.2,
                                      RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_cyc_ggg(graph, ClusteringMethod::CYC, 20, 1, BisectionMethod::GGG, 1.2,
                                     RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_hyb_ggg(graph, ClusteringMethod::HYB, 20, 1, BisectionMethod::GGG, 1.2,
                                     RefinementMethod::BOUNDARYFM, 10);

    auto vec = mb_forb_undir.run();
    {
        auto [sizeV0, sizeV1] = calculatePartSizes(vec.first, graph);
        double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
        imbalance = round(imbalance * 10) / 10;
//        auto it = lower_bound(pairs.begin(), pairs.end(),
//                              std::make_pair(vec.second, imbalance),
//                              [](auto &a, auto &b) {
//                                  if (a.first != b.first) return a.first < b.first;
//                                  return a.second < b.second;
//                              });
//        long location = it - pairs.begin();
        std::cout << "FORB,UNDIR," << vec.second << "," << imbalance << "\n";
    }

    vec = mb_cyc_undir.run();
    {
        auto [sizeV0, sizeV1] = calculatePartSizes(vec.first, graph);
        double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
        imbalance = round(imbalance * 10) / 10;
//        auto it = lower_bound(pairs.begin(), pairs.end(),
//                              std::make_pair(vec.second, imbalance),
//                              [](auto &a, auto &b) {
//                                  if (a.first != b.first) return a.first < b.first;
//                                  return a.second < b.second;
//                              });
//        long location = it - pairs.begin();
        std::cout << "CYC,UNDIR," << vec.second << "," << imbalance << "\n";
    }

    vec = mb_hyb_undir.run();
    {
        auto [sizeV0, sizeV1] = calculatePartSizes(vec.first, graph);
        double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
        imbalance = round(imbalance * 10) / 10;
//        auto it = lower_bound(pairs.begin(), pairs.end(),
//                              std::make_pair(vec.second, imbalance),
//                              [](auto &a, auto &b) {
//                                  if (a.first != b.first) return a.first < b.first;
//                                  return a.second < b.second;
//                              });
//        long location = it - pairs.begin();
        std::cout << "HYB,UNDIR," << vec.second << "," << imbalance << "\n";
    }

    vec = mb_forb_ggg.run();
    {
        auto [sizeV0, sizeV1] = calculatePartSizes(vec.first, graph);
        double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
        imbalance = round(imbalance * 10) / 10;
//        auto it = lower_bound(pairs.begin(), pairs.end(),
//                              std::make_pair(vec.second, imbalance),
//                              [](auto &a, auto &b) {
//                                  if (a.first != b.first) return a.first < b.first;
//                                  return a.second < b.second;
//                              });
//        long location = it - pairs.begin();
        std::cout << "FORB,GGG," << vec.second << "," << imbalance << "\n";
    }

    vec = mb_cyc_ggg.run();
    {
        auto [sizeV0, sizeV1] = calculatePartSizes(vec.first, graph);
        double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
        imbalance = round(imbalance * 10) / 10;
//        auto it = lower_bound(pairs.begin(), pairs.end(),
//                              std::make_pair(vec.second, imbalance),
//                              [](auto &a, auto &b) {
//                                  if (a.first != b.first) return a.first < b.first;
//                                  return a.second < b.second;
//                              });
//        long location = it - pairs.begin();
        std::cout << "CYC,GGG," << vec.second << "," << imbalance << "\n";
    }

    vec = mb_hyb_ggg.run();
    {
        auto [sizeV0, sizeV1] = calculatePartSizes(vec.first, graph);
        double imbalance = ((double) std::max(sizeV0, sizeV1) / (double) graph.totalWeight) * 100;
        imbalance = round(imbalance * 10) / 10;
//        auto it = lower_bound(pairs.begin(), pairs.end(),
//                              std::make_pair(vec.second, imbalance),
//                              [](auto &a, auto &b) {
//                                  if (a.first != b.first) return a.first < b.first;
//                                  return a.second < b.second;
//                              });
//        long location = it - pairs.begin();
        std::cout << "HYB,GGG," << vec.second << "," << imbalance << "\n";
    }

    return 0;
}
