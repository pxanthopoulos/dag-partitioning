//
// Created by panagiotis on 5/12/2024.
//

#include "Graph.h"
#include <fstream>
#include <regex>
#include <map>
#include <cassert>
#include <queue>
#include <iostream>

Graph::Graph(uint64_t size)
        : size(size), adj(size), revAdj(size), nodeWeights(size, 0), totalWeight(0), inDegree(size, 0) {
}

void Graph::addNode(uint64_t id, uint64_t weight) {
    assert(id < size && "Node ID must be smaller than graph size");
    nodeWeights[id] = weight;
    totalWeight += weight;
}

void Graph::addEdge(uint64_t from, uint64_t to, uint64_t weight) {
    assert(from < size && "Node ID `from` must be smaller than graph size");
    assert(to < size && "Node ID `to` must be smaller than graph size");
    adj[from].emplace_back(to, weight);
    revAdj[to].emplace_back(from, weight);
    inDegree[to]++;
}

std::vector<std::tuple<uint64_t, uint64_t, bool>> Graph::getNeighbors(uint64_t node) const {
    std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors;
    for (const auto &[outNeighbor, edgeWeight]: adj[node]) {
        neighbors.emplace_back(outNeighbor, edgeWeight, true);
    }

    for (const auto &[inNeighbor, edgeWeight]: revAdj[node]) {
        neighbors.emplace_back(inNeighbor, edgeWeight, false);
    }

    return neighbors;
}

std::vector<std::tuple<uint64_t, uint64_t, bool>> Graph::getNeighborsSortedByEdgeWeightAsc(uint64_t node) const {
    std::vector<std::tuple<uint64_t, uint64_t, bool>> neighbors = getNeighbors(node);

    std::sort(neighbors.begin(), neighbors.end(),
              [](const std::tuple<uint64_t, uint64_t, bool> &a, const std::tuple<uint64_t, uint64_t, bool> &b) {
                  return std::get<1>(a) < std::get<1>(b);
              });

    return neighbors;
}

bool Graph::iterativeDfsHasCycle(uint64_t start) const {
    std::vector<bool> visited(size, false);
    std::vector<bool> inStack(size, false);
    std::stack<std::pair<uint64_t, size_t>> stack;

    stack.emplace(start, 0);
    visited[start] = true;
    inStack[start] = true;

    while (!stack.empty()) {
        auto [node, neighborIndex] = stack.top();

        if (neighborIndex >= adj[node].size()) {
            inStack[node] = false;
            stack.pop();
            continue;
        }

        stack.top().second++;
        uint64_t next = adj[node][neighborIndex].first;

        if (!visited[next]) {
            visited[next] = true;
            inStack[next] = true;
            stack.emplace(next, 0);
        } else if (inStack[next]) {
            return true;
        }
    }
    return false;
}

bool Graph::hasCycle() const {
    std::vector<bool> visited(adj.size(), false);

    for (uint64_t i = 0; i < adj.size(); i++) {
        if (!visited[i] && iterativeDfsHasCycle(i)) {
            return true;
        }
    }
    return false;
}

std::vector<uint64_t> Graph::topologicalSort() const {
    std::vector<uint64_t> topologicalOrder;
    topologicalOrder.reserve(size);
    std::vector<uint64_t> localInDegree = inDegree;

    std::queue<uint64_t> q;
    for (uint64_t i = 0; i < size; i++) {
        if (localInDegree[i] == 0) {
            q.push(i);
        }
    }

    while (!q.empty()) {
        uint64_t curr = q.front();
        q.pop();
        topologicalOrder.push_back(curr);

        for (const auto &[next, weight]: adj[curr]) {
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
        }
    }

    assert(topologicalOrder.size() == size && "Graph is cyclic");

    return topologicalOrder;
}

std::vector<uint64_t> Graph::computeTopLevels() const {
    std::vector<uint64_t> localInDegree = inDegree;
    std::vector<uint64_t> topLevels(size, 0);

    std::queue<uint64_t> q;
    for (uint64_t i = 0; i < size; i++) {
        if (localInDegree[i] == 0) {
            q.push(i);
        }
    }

    while (!q.empty()) {
        uint64_t curr = q.front();
        q.pop();

        for (const auto &[next, weight]: adj[curr]) {
            topLevels[next] = std::max(topLevels[next], topLevels[curr] + 1);
            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
            assert(localInDegree[next] >= 0 && "Graph is cyclic");
        }
    }

    return topLevels;
}

std::vector<uint64_t> Graph::computeTopLevels(const std::vector<uint64_t> &topologicalOrder) const {
    std::vector<uint64_t> topLevels(size, 0);

    for (uint64_t u: topologicalOrder) {
        for (const auto &[v, weight]: adj[u]) {
            topLevels[v] = std::max(topLevels[v], topLevels[u] + 1);
        }
    }

    return topLevels;
}

std::pair<std::vector<uint64_t>, std::vector<uint64_t>> Graph::topologicalSortAndTopLevels() const {
    std::vector<uint64_t> localInDegree = inDegree;
    std::vector<uint64_t> topologicalOrder;
    std::vector<uint64_t> topLevels(size, 0);
    topologicalOrder.reserve(size);

    std::queue<uint64_t> q;
    for (uint64_t i = 0; i < size; i++) {
        if (localInDegree[i] == 0) {
            q.push(i);
            topLevels[i] = 0;
        }
    }

    while (!q.empty()) {
        uint64_t curr = q.front();
        q.pop();
        topologicalOrder.push_back(curr);

        for (const auto &[next, weight]: adj[curr]) {
            topLevels[next] = std::max(topLevels[next], topLevels[curr] + 1);

            localInDegree[next]--;
            if (localInDegree[next] == 0) {
                q.push(next);
            }
            assert(localInDegree[next] >= 0 && "Graph is cyclic");
        }
    }

//    assert(topologicalOrder.size() == size && "Graph is cyclic");
    if (topologicalOrder.size() != size) {
        std::cerr << "\n\nGRAPH IS CYCLIC, current graph in err.dot\n\n";

        printToDot("/home/panagiotis/code/dag-partitioning/err.dot");

        if (hasCycle()) {
            std::cerr << "GRAPH IS INDEED CYCLIC\n\n";
        }

        assert(1 == 0);
    }

    return {topologicalOrder, topLevels};
}

void Graph::print(llvm::raw_ostream &os) const {
    assert(adj.size() == nodeWeights.size() &&
           "Adjacency-list 2D vectors' size and nodeWeights vectors' size must be equal");
    assert(size == adj.size() && "Adjacency-list 2D vectors' size must be the same as the graph size");
    for (size_t i = 0; i < size; ++i) {
        const auto &adjList = adj[i];
        uint64_t weight = nodeWeights[i];
        os << "Node: " << i << "\nWeight: " << weight << "\n";
        os << "Neighbors: ";
        for (const auto &[neighborId, edgeWeight]: adjList) {
            os << neighborId << "(" << edgeWeight << "), ";
        }
        os << "\n\n";
    }
}

void Graph::printToDot(const std::string &dotFilename) const {
    std::ofstream dotFile(dotFilename);

    dotFile << "// size=" << size << "\n";
    dotFile << "digraph cfg {\n";
    for (uint64_t i = 0; i < size; ++i) {
        dotFile << i << "[weight=" << nodeWeights[i] << "];\n";
    }
    for (uint64_t i = 0; i < size; ++i) {
        for (const auto &[neighborId, edgeWeight]: adj[i]) {
            dotFile << i << "->" << neighborId << "[weight=" << edgeWeight << "];\n";
        }
    }
    dotFile << "}\n";
}

Graph readDotFile(const std::string &dotFilename, const std::string &mappingFilename) {
    std::ifstream dotFile(dotFilename);
    std::string line;
    std::regex sizeRegex(R"(//\s*size\s*=\s*(\d+))");

    uint64_t graphSize = UINT64_MAX;

    while (std::getline(dotFile, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, sizeRegex)) {
            graphSize = std::stoi(matches[1]);
            break;
        }
    }

    assert(graphSize != UINT64_MAX && "Failed to extract graph size");

    Graph g(graphSize);

    std::map<std::string, uint64_t> nodeMap;
    std::regex nodeRegex(R"(([a-zA-Z0-9_]+)\[weight=(\d+)\];)");
    std::regex edgeRegex(R"(([a-zA-Z0-9_]+)->([a-zA-Z0-9_]+)\[weight=(\d+)\];)");
    uint64_t nodeId = 0;

    while (std::getline(dotFile, line)) {
        std::smatch matches;
        if (std::regex_search(line, matches, edgeRegex)) {
            uint64_t from = nodeMap[matches[1].str()];
            uint64_t to = nodeMap[matches[2].str()];
            uint64_t nodeWeight = std::stoull(matches[3].str());
            g.addEdge(from, to, nodeWeight);
        } else if (std::regex_search(line, matches, nodeRegex)) {
            uint64_t edgeWeight = std::stoull(matches[2].str());
            g.addNode(nodeId, edgeWeight);
            nodeMap[matches[1].str()] = nodeId;
            nodeId++;
        }
    }

    std::ofstream mapFile(mappingFilename);
    for (const auto &[name, id]: nodeMap) {
        mapFile << name << " -> " << id << "\n";
    }

    return g;
}