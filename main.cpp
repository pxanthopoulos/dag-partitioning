#include "Graph.h"
#include "RecursivePartitioner.h"

#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <thread>

uint64_t calculateMaxPartSize(const std::vector<uint64_t> &partition, const Graph &graph, uint64_t partitions) {
    std::vector<uint64_t> partSizes(partitions, 0);
    for (uint64_t i = 0; i < partition.size(); ++i) {
        partSizes[partition[i]] += graph.nodeWeights[i];
    }
    return *std::max_element(partSizes.begin(), partSizes.end());
}

void
printResults(const std::vector<uint64_t> &partition, const Graph &graph, uint64_t cutSize, const std::string &methodStr,
             uint64_t duration, uint64_t partitions) {
    uint64_t maxPartSize = calculateMaxPartSize(partition, graph, partitions);
    double imbalance = ((double) maxPartSize / (double) graph.totalWeight) * 100;
    imbalance = round(imbalance * 10) / 10;
    std::cout << methodStr << "," << cutSize << "," << imbalance << "," << duration << "\n";
}

struct MethodCombination {
    ClusteringMethod clustering;
    BisectionMethod bisection;
    RefinementMethod refinement;
    std::string clusteringStr;
    std::string bisectionStr;
    std::string refinementStr;
};

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <# of partitions> <dot file path>" << std::endl;
        return 1;
    }

    uint64_t partitions;
    try {
        partitions = std::stoull(argv[1]);
    } catch (const std::exception &e) {
        std::cerr << "Error: Please provide a valid positive integer" << std::endl;
        return 1;
    }

    std::string dotPath = argv[2];
    Graph graph = readDotFile(dotPath, dotPath + ".node-mappings.txt");

    std::vector<MethodCombination> methods = {
            {ClusteringMethod::FORB, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYFM, "FORB", "UNDIR", "FM"},
            {ClusteringMethod::CYC,  BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYFM, "CYC",  "UNDIR", "FM"},
            {ClusteringMethod::HYB,  BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYFM, "HYB",  "UNDIR", "FM"},
            {ClusteringMethod::FORB, BisectionMethod::GGG,       RefinementMethod::BOUNDARYFM, "FORB", "GGG",   "FM"},
            {ClusteringMethod::CYC,  BisectionMethod::GGG,       RefinementMethod::BOUNDARYFM, "CYC",  "GGG",   "FM"},
            {ClusteringMethod::HYB,  BisectionMethod::GGG,       RefinementMethod::BOUNDARYFM, "HYB",  "GGG",   "FM"},
            {ClusteringMethod::FORB, BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYKL, "FORB", "UNDIR", "KL"},
            {ClusteringMethod::CYC,  BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYKL, "CYC",  "UNDIR", "KL"},
            {ClusteringMethod::HYB,  BisectionMethod::UNDIRBOTH, RefinementMethod::BOUNDARYKL, "HYB",  "UNDIR", "KL"},
            {ClusteringMethod::FORB, BisectionMethod::GGG,       RefinementMethod::BOUNDARYKL, "FORB", "GGG",   "KL"},
            {ClusteringMethod::CYC,  BisectionMethod::GGG,       RefinementMethod::BOUNDARYKL, "CYC",  "GGG",   "KL"},
            {ClusteringMethod::HYB,  BisectionMethod::GGG,       RefinementMethod::BOUNDARYKL, "HYB",  "GGG",   "KL"},
            {ClusteringMethod::FORB, BisectionMethod::UNDIRBOTH, RefinementMethod::MIXED,      "FORB", "UNDIR", "MIX"},
            {ClusteringMethod::CYC,  BisectionMethod::UNDIRBOTH, RefinementMethod::MIXED,      "CYC",  "UNDIR", "MIX"},
            {ClusteringMethod::HYB,  BisectionMethod::UNDIRBOTH, RefinementMethod::MIXED,      "HYB",  "UNDIR", "MIX"},
            {ClusteringMethod::FORB, BisectionMethod::GGG,       RefinementMethod::MIXED,      "FORB", "GGG",   "MIX"},
            {ClusteringMethod::CYC,  BisectionMethod::GGG,       RefinementMethod::MIXED,      "CYC",  "GGG",   "MIX"},
            {ClusteringMethod::HYB,  BisectionMethod::GGG,       RefinementMethod::MIXED,      "HYB",  "GGG",   "MIX"}
    };

    const uint64_t maxLevel = 20;
    const uint64_t minSize = 1;
    const double maxImbalance = 1.01;
    const uint64_t maxPasses = 10;

    struct ThreadResult {
        std::vector<uint64_t> partition;
        uint64_t cutSize{};
        uint64_t duration{};
    };

    std::vector<std::thread> threads;
    std::vector<ThreadResult> results(methods.size());

    for (uint64_t i = 0; i < methods.size(); i++) {
        const auto &method = methods[i];
        threads.emplace_back([&graph, &method, &results, i, maxImbalance, &partitions]() {
            RecursivePartitioner partitioner(
                    graph,
                    partitions,
                    method.clustering,
                    maxLevel,
                    minSize,
                    method.bisection,
                    maxImbalance,
                    method.refinement,
                    maxPasses
            );

            auto start = std::chrono::high_resolution_clock::now();
            auto [partition, cutSize] = partitioner.run();
            auto end = std::chrono::high_resolution_clock::now();
            uint64_t duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

            results[i] = ThreadResult{std::move(partition), cutSize, duration};
        });
    }

    for (auto &thread: threads) {
        thread.join();
    }

    for (uint64_t i = 0; i < methods.size(); i++) {
        const auto &method = methods[i];
        const auto &result = results[i];

        std::string methodStr = method.clusteringStr + "," +
                                method.bisectionStr + "," +
                                method.refinementStr;
        printResults(result.partition, graph, result.cutSize, methodStr, result.duration, partitions);
    }

    return 0;
}
