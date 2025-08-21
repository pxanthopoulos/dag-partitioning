#include "Graph.h"
#include "RecursivePartitioner.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
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
    if (argc < 4 || argc > 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <# of partitions> <dot file path> <enable_parallel> "
                     "[min_size_for_parallel] [max_parallel_depth]"
                  << std::endl;
        std::cerr << "  enable_parallel: 1 for parallel, 0 for sequential"
                  << std::endl;
        std::cerr << "  min_size_for_parallel: minimum subgraph size to "
                     "parallelize (default: 100)"
                  << std::endl;
        std::cerr << "  max_parallel_depth: maximum recursion depth for "
                     "parallelization (default: 10)"
                  << std::endl;
        return 1;
    }

    uint64_t partitions;
    try {
        partitions = std::stoull(argv[1]);
    } catch (const std::exception &e) {
        std::cerr
            << "Error: Please provide a valid positive integer for partitions"
            << std::endl;
        return 1;
    }

    std::string dotPath = argv[2];

    // Parse enable parallel flag (required)
    bool enableParallel;
    try {
        int parallelFlag = std::stoi(argv[3]);
        if (parallelFlag != 0 && parallelFlag != 1) {
            std::cerr << "Error: enable_parallel must be 0 or 1" << std::endl;
            return 1;
        }
        enableParallel = (parallelFlag == 1);
    } catch (const std::exception &e) {
        std::cerr << "Error: Please provide a valid integer for "
                     "enable_parallel (0 or 1)"
                  << std::endl;
        return 1;
    }

    uint64_t minSizeForParallel = 100;
    if (argc >= 5) {
        try {
            minSizeForParallel = std::stoull(argv[4]);
        } catch (const std::exception &e) {
            std::cerr << "Error: Please provide a valid positive integer for "
                         "min_size_for_parallel"
                      << std::endl;
            return 1;
        }
    }

    uint64_t maxParallelDepth = 10;
    if (argc == 6) {
        try {
            maxParallelDepth = std::stoull(argv[5]);
        } catch (const std::exception &e) {
            std::cerr << "Error: Please provide a valid positive integer for "
                         "max_parallel_depth"
                      << std::endl;
            return 1;
        }
    }
    dag_partitioning::core::Graph graph = dag_partitioning::core::readDotFile(
        dotPath, dotPath + ".node-mappings.txt");

    std::vector<MethodCombination> methods = {
        {"FORB", "UNDIRMETIS", "BOUNDARYFM"},
        {"CYC", "UNDIRMETIS", "BOUNDARYFM"},
        {"HYB", "UNDIRMETIS", "BOUNDARYFM"},
        {"FORB", "GGG", "BOUNDARYFM"},
        {"CYC", "GGG", "BOUNDARYFM"},
        {"HYB", "GGG", "BOUNDARYFM"},
        {"FORB", "UNDIRMETIS", "BOUNDARYKL"},
        {"CYC", "UNDIRMETIS", "BOUNDARYKL"},
        {"HYB", "UNDIRMETIS", "BOUNDARYKL"},
        {"FORB", "GGG", "BOUNDARYKL"},
        {"CYC", "GGG", "BOUNDARYKL"},
        {"HYB", "GGG", "BOUNDARYKL"},
        {"FORB", "UNDIRMETIS", "MIXED"},
        {"CYC", "UNDIRMETIS", "MIXED"},
        {"HYB", "UNDIRMETIS", "MIXED"},
        {"FORB", "GGG", "MIXED"},
        {"CYC", "GGG", "MIXED"},
        {"HYB", "GGG", "MIXED"}};

    const uint64_t maxLevel = 20;
    const uint64_t minSize = 50 * partitions;
    const double vertRatio = 0.9;
    const double maxImbalance = 1.1;
    const uint64_t maxPasses = 10;

    std::cout << "CLUSTERING,BISECTION,REFINEMENT,CUTSIZE,IMBALANCE,DURATION"
              << std::endl;

    for (uint64_t i = 0; i < methods.size(); i++) {
        const auto &method = methods[i];
        dag_partitioning::driver::RecursivePartitioner partitioner(
            graph, partitions, method.clusteringStr, maxLevel, minSize,
            vertRatio, method.bisectionStr, maxImbalance, method.refinementStr,
            maxPasses, enableParallel, minSizeForParallel, maxParallelDepth);

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
