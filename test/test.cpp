#include "Graph.h"
#include "RecursivePartitioner.h"
#include "Scheduling.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <string>
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
                  const std::string &methodStr, uint64_t partitioningDuration,
                  uint64_t partitions, uint64_t schedulingDuration) {
    uint64_t maxPartSize = calculateMaxPartSize(partition, graph, partitions);
    double imbalance =
        std::abs(((double)maxPartSize / (double)graph.totalWeight) * 100 -
                 ((double)100 / (double)partitions));
    imbalance = round(imbalance * 10) / 10;
    std::cout << methodStr << "," << cutSize << "," << imbalance << ","
              << partitioningDuration << "," << schedulingDuration << "\n";
}

int main(int argc, char **argv) {
    if (argc < 7 || argc > 10) {
        std::cerr << "Usage: " << argv[0]
                  << " <# of partitions> <dot file path> <clustering method> "
                     "<bisection method> "
                     "<refinement method> <enable parallel> "
                     "[min size for parallel] [max parallel depth] [enable "
                     "scheduling]"
                  << std::endl;
        std::cerr << "  clustering method: FORB, CYC, or HYB" << std::endl;
        std::cerr
            << "  bisection method: GGG, UNDIRSCOTCH, UNDIRMETIS or UNDIRBOTH"
            << std::endl;
        std::cerr << "  refinement method: BOUNDARYFM, BOUNDARYKL, or MIXED"
                  << std::endl;
        std::cerr << "  enable parallel: 1 for parallel, 0 for sequential"
                  << std::endl;
        std::cerr << "  min size for parallel: minimum subgraph size to "
                     "parallelize (default: 100)"
                  << std::endl;
        std::cerr << "  max parallel depth: maximum recursion depth for "
                     "parallelization (default: 10)"
                  << std::endl;
        std::cerr << "  enable scheduling: 1 to enable scheduling, 0 to skip "
                     "(default: 0)"
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
    std::string clusteringMethod = argv[3];
    std::string bisectionMethod = argv[4];
    std::string refinementMethod = argv[5];

    // Validate clustering method
    if (clusteringMethod != "FORB" && clusteringMethod != "CYC" &&
        clusteringMethod != "HYB") {
        std::cerr << "Error: clustering method must be FORB, CYC, or HYB"
                  << std::endl;
        return 1;
    }

    // Validate bisection method
    if (bisectionMethod != "GGG" && bisectionMethod != "UNDIRMETIS" &&
        bisectionMethod != "UNDIRSCOTCH" && bisectionMethod != "UNDIRBOTH") {
        std::cerr << "Error: bisection method must be GGG, UNDIRMETIS, "
                     "UNDIRSCOTCH, or UNDIRBOTH"
                  << std::endl;
        return 1;
    }

    // Validate refinement method
    if (refinementMethod != "BOUNDARYFM" && refinementMethod != "BOUNDARYKL" &&
        refinementMethod != "MIXED") {
        std::cerr << "Error: refinement method must be BOUNDARYFM, BOUNDARYKL, "
                     "or MIXED"
                  << std::endl;
        return 1;
    }

    // Parse enable parallel flag (required)
    bool enableParallel;
    try {
        int parallelFlag = std::stoi(argv[6]);
        if (parallelFlag != 0 && parallelFlag != 1) {
            std::cerr << "Error: enable parallel must be 0 or 1" << std::endl;
            return 1;
        }
        enableParallel = (parallelFlag == 1);
    } catch (const std::exception &e) {
        std::cerr << "Error: Please provide a valid integer for "
                     "enable parallel (0 or 1)"
                  << std::endl;
        return 1;
    }

    uint64_t minSizeForParallel = 100;
    if (argc >= 8) {
        try {
            minSizeForParallel = std::stoull(argv[7]);
        } catch (const std::exception &e) {
            std::cerr << "Error: Please provide a valid positive integer for "
                         "min size for parallel"
                      << std::endl;
            return 1;
        }
    }

    uint64_t maxParallelDepth = 10;
    if (argc >= 9) {
        try {
            maxParallelDepth = std::stoull(argv[8]);
        } catch (const std::exception &e) {
            std::cerr << "Error: Please provide a valid positive integer for "
                         "max parallel depth"
                      << std::endl;
            return 1;
        }
    }

    // Parse enable scheduling flag (optional, defaults to false)
    bool enableScheduling = false;
    if (argc == 10) {
        try {
            int schedulingFlag = std::stoi(argv[9]);
            if (schedulingFlag != 0 && schedulingFlag != 1) {
                std::cerr << "Error: enable scheduling must be 0 or 1"
                          << std::endl;
                return 1;
            }
            enableScheduling = (schedulingFlag == 1);
        } catch (const std::exception &e) {
            std::cerr << "Error: Please provide a valid integer for "
                         "enable scheduling (0 or 1)"
                      << std::endl;
            return 1;
        }
    }
    dag_partitioning::core::Graph graph = dag_partitioning::core::readDotFile(
        dotPath, dotPath + ".node-mappings.txt");

    const uint64_t maxLevel = 20;
    const uint64_t minSize = 50 * partitions;
    const double vertRatio = 0.9;
    const double maxImbalance = 1.1;
    const uint64_t maxPasses = 10;

    dag_partitioning::driver::RecursivePartitioner partitioner(
        graph, partitions, clusteringMethod, maxLevel, minSize, vertRatio,
        bisectionMethod, maxImbalance, refinementMethod, maxPasses,
        enableParallel, minSizeForParallel, maxParallelDepth);

    auto start = std::chrono::high_resolution_clock::now();
    auto [partitionMapping, cutSize] =
        dag_partitioning::driver::RecursivePartitioner::resetPartitionNumbers(
            partitioner.run());
    auto end = std::chrono::high_resolution_clock::now();
    uint64_t duration =
        std::chrono::duration_cast<std::chrono::microseconds>(end - start)
            .count();

    std::string methodStr =
        clusteringMethod + "," + bisectionMethod + "," + refinementMethod;

    uint64_t schedulingDuration = 0;

    if (enableScheduling) {
        dag_partitioning::scheduling::Scheduler scheduler(
            graph, partitionMapping, false, false);
        start = std::chrono::high_resolution_clock::now();
        auto [schedule, peakMemory] = scheduler.run();
        end = std::chrono::high_resolution_clock::now();
        schedulingDuration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                .count();
    }

    printResults(partitionMapping, graph, cutSize, methodStr, duration,
                 partitions, schedulingDuration);

    return 0;
}
