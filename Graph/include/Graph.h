//
// Created by panagiotis on 5/12/2024.
//

#ifndef DAG_PARTITIONING_GRAPH_H
#define DAG_PARTITIONING_GRAPH_H


#include <vector>
#include <string>
#include <utility>
#include <tuple>
#include "llvm/Support/raw_ostream.h"

class Graph {
public:
    uint64_t size;
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> adj;
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> revAdj;
    std::vector<uint64_t> nodeWeights;
    uint64_t totalWeight;
    std::vector<uint64_t> inDegree;

    explicit Graph(uint64_t size);

    void addNode(uint64_t id, uint64_t weight);

    void addEdge(uint64_t from, uint64_t to, uint64_t weight);

    [[nodiscard]] std::vector<std::tuple<uint64_t, uint64_t, bool>> getNeighbors(uint64_t node) const;

    [[nodiscard]] std::vector<std::tuple<uint64_t, uint64_t, bool>> getNeighborsSortedByEdgeWeight(uint64_t node) const;

    [[nodiscard]] bool hasCycle() const;

    [[nodiscard]] bool iterativeDfsHasCycle(uint64_t start) const;

    [[nodiscard]] std::vector<uint64_t> topologicalSort() const;

    [[nodiscard]] std::vector<uint64_t> computeTopLevels() const;

    [[nodiscard]] std::vector<uint64_t> computeTopLevels(const std::vector<uint64_t> &topologicalOrder) const;

    [[nodiscard]] std::pair<std::vector<uint64_t>, std::vector<uint64_t>> topologicalSortAndTopLevels() const;

    void print(llvm::raw_ostream &os) const;

    void printToDot(const std::string &dotFilename) const;
};

Graph readDotFile(const std::string &dotFilename, const std::string &mappingFilename);


#endif //DAG_PARTITIONING_GRAPH_H
