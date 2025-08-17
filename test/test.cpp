#include "Graph.h"
#include "RecursivePartitioner.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <thread>
#include <vector>

uint64_t calculateMaxPartSize(const std::vector<uint64_t> &partition,
                              const dag_partitioning::core::Graph &graph,
                              uint64_t partitions) {
    std::vector<uint64_t> partSizes(partitions, 0);
    for (uint64_t i = 0; i < partition.size(); ++i) {
        partSizes[partition[i]] += graph.nodeWeights[i];
    }
    return *std::max_element(partSizes.begin(), partSizes.end());
}

void printResults(const std::vector<uint64_t> &partition,
                  const dag_partitioning::core::Graph &graph, uint64_t cutSize,
                  const std::string &methodStr, uint64_t duration,
                  uint64_t partitions) {
    uint64_t maxPartSize = calculateMaxPartSize(partition, graph, partitions);
    double imbalance =
        std::abs(((double)maxPartSize / (double)graph.totalWeight) * 100 -
                 ((double)100 / (double)partitions));
    imbalance = round(imbalance * 10) / 10;
    std::cout << methodStr << "," << cutSize << "," << imbalance << ","
              << duration << "\n";
}

struct MethodCombination {
    std::string clusteringStr;
    std::string bisectionStr;
    std::string refinementStr;
};

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <# of partitions> <dot file path>" << std::endl;
        return 1;
    }

    uint64_t partitions;
    try {
        partitions = std::stoull(argv[1]);
    } catch (const std::exception &e) {
        std::cerr << "Error: Please provide a valid positive integer"
                  << std::endl;
        return 1;
    }

    std::string dotPath = argv[2];
    dag_partitioning::core::Graph graph = dag_partitioning::core::readDotFile(
        dotPath, dotPath + ".node-mappings.txt");

    std::vector<MethodCombination> methods = {
        {"FORB", "UNDIRBOTH", "BOUNDARYFM"},
        {"CYC", "UNDIRBOTH", "BOUNDARYFM"},
        {"HYB", "UNDIRBOTH", "BOUNDARYFM"},
        {"FORB", "GGG", "BOUNDARYFM"},
        {"CYC", "GGG", "BOUNDARYFM"},
        {"HYB", "GGG", "BOUNDARYFM"},
        {"FORB", "UNDIRBOTH", "BOUNDARYKL"},
        {"CYC", "UNDIRBOTH", "BOUNDARYKL"},
        {"HYB", "UNDIRBOTH", "BOUNDARYKL"},
        {"FORB", "GGG", "BOUNDARYKL"},
        {"CYC", "GGG", "BOUNDARYKL"},
        {"HYB", "GGG", "BOUNDARYKL"},
        {"FORB", "UNDIRBOTH", "MIXED"},
        {"CYC", "UNDIRBOTH", "MIXED"},
        {"HYB", "UNDIRBOTH", "MIXED"},
        {"FORB", "GGG", "MIXED"},
        {"CYC", "GGG", "MIXED"},
        {"HYB", "GGG", "MIXED"}};

    const uint64_t maxLevel = 20;
    const uint64_t minSize = 50 * partitions;
    const double vertRatio = 0.9;
    const double maxImbalance = 1.1;
    const uint64_t maxPasses = 10;

    for (uint64_t i = 0; i < methods.size(); i++) {
        const auto &method = methods[i];
        dag_partitioning::driver::RecursivePartitioner partitioner(
            graph, partitions, method.clusteringStr, maxLevel, minSize,
            vertRatio, method.bisectionStr, maxImbalance, method.refinementStr,
            maxPasses);

        auto start = std::chrono::high_resolution_clock::now();
        auto [partition, cutSize] = partitioner.run();
        auto end = std::chrono::high_resolution_clock::now();
        uint64_t duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                .count();

        std::string methodStr = method.clusteringStr + "," +
                                method.bisectionStr + "," +
                                method.refinementStr;
        printResults(partition, graph, cutSize, methodStr, duration,
                     partitions);
    }

    return 0;
}
