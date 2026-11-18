/**
 * @file DotFile.cpp
 * @brief Reading and writing graphs in the DOT format
 */

#include "Graph.h"

#include "robin_hood.h"
#include <fstream>
#include <stdexcept>
#include <string>

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

// Read graph from DOT format file
Graph readDotFile(const std::string &dotFilename,
                  const std::string &mappingFilename) {
    std::ifstream dotFile(dotFilename);
    if (!dotFile.is_open()) {
        throw std::runtime_error("Could not open DOT file: " + dotFilename);
    }

    std::string line;
    line.reserve(256); // Pre-allocate line buffer for typical line lengths

    uint64_t graphSize = 0;
    while (std::getline(dotFile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line.find_first_not_of(" \t") == '/' ||
            line.find_first_not_of(" \t") == std::string::npos)
            continue;

        // Skip edge lines
        size_t arrowPos = line.find("->");
        if (arrowPos != std::string::npos)
            continue;

        size_t semicolon = line.find(";");
        if (semicolon == std::string::npos)
            continue;

        graphSize++;
    }

    Graph g(graphSize);

    robin_hood::unordered_map<std::string, uint64_t> nodeMap;
    nodeMap.reserve(graphSize);

    uint64_t nodeId = 0;

    dotFile.clear();
    dotFile.seekg(0, std::ios::beg);

    while (std::getline(dotFile, line)) {
        // Skip empty lines and comments
        if (line.empty() || line[line.find_first_not_of(" \t")] == '/' ||
            line.find_first_not_of(" \t") == std::string::npos)
            continue;

        // Check for edge: "node1->node2[weight=123];"
        size_t arrowPos = line.find("->");
        if (arrowPos != std::string::npos) {
            // Parse edge
            size_t fromStart = line.find_first_not_of(" \t");
            size_t fromEnd = line.find_last_not_of(" \t", arrowPos - 1) + 1;

            size_t toStart = line.find_first_not_of(" \t", arrowPos + 2);
            size_t bracketPos = line.find("[", toStart);
            if (bracketPos == std::string::npos)
                continue;

            size_t toEnd = line.find_last_not_of(" \t", bracketPos - 1) + 1;

            // Extract weight
            size_t weightPos = line.find("weight=", bracketPos);
            if (weightPos == std::string::npos)
                continue;

            size_t weightStart = weightPos + 7; // length of "weight="
            size_t weightEnd = line.find_first_of("];", weightStart);
            if (weightEnd == std::string::npos)
                continue;

            std::string fromNode = line.substr(fromStart, fromEnd - fromStart);
            std::string toNode = line.substr(toStart, toEnd - toStart);
            std::string weightStr =
                line.substr(weightStart, weightEnd - weightStart);

            uint64_t from = nodeMap[fromNode];
            uint64_t to = nodeMap[toNode];
            uint64_t edgeWeight = std::stoull(weightStr);
            g.addEdge(from, to, edgeWeight);
        }
        // Check for node: "node_name[weight=123];"
        else {
            size_t bracketPos = line.find("[");
            if (bracketPos == std::string::npos)
                continue;

            size_t nodeStart = line.find_first_not_of(" \t");
            size_t nodeEnd = line.find_last_not_of(" \t", bracketPos - 1) + 1;

            if (nodeStart >= nodeEnd)
                continue;

            // Extract weight
            size_t weightPos = line.find("weight=", bracketPos);
            if (weightPos == std::string::npos)
                continue;

            size_t weightStart = weightPos + 7; // length of "weight="
            size_t weightEnd = line.find_first_of("];", weightStart);
            if (weightEnd == std::string::npos)
                continue;

            std::string nodeName = line.substr(nodeStart, nodeEnd - nodeStart);

            std::string weightStr =
                line.substr(weightStart, weightEnd - weightStart);

            uint64_t nodeWeight = std::stoull(weightStr);
            g.addNode(nodeId, nodeWeight);
            nodeMap[nodeName] = nodeId;
            nodeId++;
        }
    }

    if (g.adj.size() != graphSize) {
        throw std::runtime_error(
            "DOT file parsing error: computed graph size (" +
            std::to_string(g.adj.size()) + ") != declared size (" +
            std::to_string(graphSize) + ")");
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
