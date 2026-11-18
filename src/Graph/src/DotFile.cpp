/**
 * @file DotFile.cpp
 * @brief Reading and writing graphs in the DOT format
 */

#include "Graph.h"

#include "robin_hood.h"
#include <fstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace dag_partitioning {

namespace core {

// Export graph to DOT format file
void Graph::printToDot(const std::string &dotFilename) const {
    std::ofstream dotFile(dotFilename);

    // Write graph size as comment
    dotFile << "// size=" << size << "\n";
    dotFile << "digraph cfg {\n";

    // Write node weights
    for (uint64_t i = 0; i < size; ++i) {
        dotFile << i << "[weight=" << nodeWeights[i] << "];\n";
    }

    // Write edges and their weights
    for (uint64_t i = 0; i < size; ++i) {
        for (const auto &[neighborId, edgeWeight] : adj[i]) {
            dotFile << i << "->" << neighborId << "[weight=" << edgeWeight
                    << "];\n";
        }
    }
    dotFile << "}\n";
}

namespace {

// The text of line in [begin, end) without surrounding blanks
std::string trimmed(const std::string &line, size_t begin, size_t end) {
    size_t first = line.find_first_not_of(" \t", begin);
    if (first == std::string::npos || first >= end)
        return "";
    size_t last = line.find_last_not_of(" \t", end - 1);
    return line.substr(first, last + 1 - first);
}

// The weight in "[weight=123]" after position bracket, if any
bool readWeight(const std::string &line, size_t bracket, uint64_t &weight) {
    size_t weightPos = line.find("weight=", bracket);
    if (weightPos == std::string::npos)
        return false;
    size_t weightStart = weightPos + 7; // length of "weight="
    size_t weightEnd = line.find_first_of("];", weightStart);
    if (weightEnd == std::string::npos)
        return false;
    weight = std::stoull(line.substr(weightStart, weightEnd - weightStart));
    return true;
}

} // namespace

// Read graph from DOT format file: nodes are "name[weight=123];" and edges
// "from->to[weight=123];", one per line; other lines are skipped
Graph readDotFile(const std::string &dotFilename,
                  const std::string &mappingFilename) {
    std::ifstream dotFile(dotFilename);
    if (!dotFile.is_open()) {
        throw std::runtime_error("Could not open DOT file: " + dotFilename);
    }

    struct Edge {
        std::string from, to;
        uint64_t weight;
    };
    std::vector<std::pair<std::string, uint64_t>> nodes;
    std::vector<Edge> edges;
    std::string line;
    while (std::getline(dotFile, line)) {
        // Skip blank lines and comments
        size_t start = line.find_first_not_of(" \t");
        if (start == std::string::npos || line[start] == '/')
            continue;

        size_t arrowPos = line.find("->");
        size_t bracketPos = line.find(
            "[", arrowPos == std::string::npos ? start : arrowPos + 2);
        uint64_t weight;
        if (bracketPos == std::string::npos ||
            !readWeight(line, bracketPos, weight))
            continue;
        if (arrowPos != std::string::npos) {
            edges.push_back({trimmed(line, start, arrowPos),
                             trimmed(line, arrowPos + 2, bracketPos), weight});
        } else {
            std::string name = trimmed(line, start, bracketPos);
            if (!name.empty())
                nodes.emplace_back(std::move(name), weight);
        }
    }

    Graph g(nodes.size());
    robin_hood::unordered_map<std::string, uint64_t> nodeMap;
    nodeMap.reserve(nodes.size());
    for (uint64_t nodeId = 0; nodeId < nodes.size(); nodeId++) {
        g.addNode(nodeId, nodes[nodeId].second);
        nodeMap[nodes[nodeId].first] = nodeId;
    }
    for (const Edge &edge : edges) {
        auto from = nodeMap.find(edge.from);
        auto to = nodeMap.find(edge.to);
        if (from == nodeMap.end() || to == nodeMap.end()) {
            throw std::runtime_error("DOT file parsing error: edge " +
                                     edge.from + "->" + edge.to +
                                     " has an undeclared node");
        }
        g.addEdge(from->second, to->second, edge.weight);
    }

    // Write node mapping
    std::ofstream mapFile(mappingFilename);
    if (mapFile.is_open()) {
        for (const auto &[name, id] : nodeMap) {
            mapFile << name << " -> " << id << '\n';
        }
    }

    return g;
}

} // namespace core

} // namespace dag_partitioning
