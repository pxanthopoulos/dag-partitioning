#include <fstream>
#include <iostream>
#include <random>
#include <set>
#include <stack>
#include <vector>

class DAGGenerator {
private:
    std::vector<std::vector<std::pair<int, int>>> adj;
    int size;
    std::random_device rd;
    std::mt19937 gen;

    bool hasCycle() {
        std::vector<int> visited(size, 0);
        std::vector<int> recursionStack(size, 0);

        for (int i = 0; i < size; i++)
            if (dfs(i, visited, recursionStack))
                return true;
        return false;
    }

    bool dfs(int v, std::vector<int> &visited,
             std::vector<int> &recursionStack) {
        if (visited[v] == 0) {
            visited[v] = 1;
            recursionStack[v] = 1;

            for (const auto &neighbor : adj[v]) {
                if (visited[neighbor.first] == 0 &&
                    dfs(neighbor.first, visited, recursionStack))
                    return true;
                else if (recursionStack[neighbor.first] == 1)
                    return true;
            }
        }
        recursionStack[v] = 0;
        return false;
    }

    bool isConnected() {
        std::vector<bool> visited(size, false);
        std::stack<int> s;
        s.push(0);
        visited[0] = true;
        int count = 1;

        while (!s.empty()) {
            int v = s.top();
            s.pop();

            for (const auto &neighbor : adj[v]) {
                if (!visited[neighbor.first]) {
                    visited[neighbor.first] = true;
                    s.push(neighbor.first);
                    count++;
                }
            }
        }

        return count == size;
    }

    void createSpanningTree() {
        std::vector<bool> used(size, false);
        std::uniform_int_distribution<> dis(0, size - 1);

        // Start with node 0
        used[0] = true;
        int usedCount = 1;

        // Connect remaining nodes
        while (usedCount < size) {
            int from = dis(gen);
            if (!used[from])
                continue;

            int to = dis(gen);
            if (!used[to] && from != to) {
                adj[from].emplace_back(to, 1);
                used[to] = true;
                usedCount++;
            }
        }
    }

public:
    explicit DAGGenerator(int n) : size(n), gen(rd()) { adj.resize(n); }

    int generate(double edgeRatio) {
        std::uniform_int_distribution<> dis(0, size - 1);
        std::set<std::pair<int, int>> edges;

        // First create a spanning tree to ensure connectivity
        createSpanningTree();

        // Copy existing edges to set
        for (int i = 0; i < size; i++) {
            for (const auto &edge : adj[i]) {
                edges.insert({i, edge.first});
            }
        }

        // Add additional random edges
        int targetEdges = (int)((double)size * edgeRatio);
        int attempts = 0;
        int maxAttempts = size * size;

        while (edges.size() < targetEdges && attempts < maxAttempts) {
            int from = dis(gen);
            int to = dis(gen);

            if (from != to && edges.find({from, to}) == edges.end()) {
                adj[from].emplace_back(to, 1);

                if (!hasCycle()) {
                    edges.insert({from, to});
                } else {
                    adj[from].pop_back();
                }
            }
            attempts++;
        }

        // Rebuild adjacency list
        adj.clear();
        adj.resize(size);
        for (const auto &edge : edges) {
            std::uniform_int_distribution<> disW(1, 20);
            int random_number = disW(gen);
            adj[edge.first].emplace_back(edge.second, random_number);
        }

        if (hasCycle()) {
            std::cout << "CYCLIC\n";
            return 1;
        }

        if (!isConnected()) {
            std::cout << "NOT CONNECTED\n";
            return 1;
        }

        return 0;
    }

    void printToDot(const std::string &dotFilename) {
        std::ofstream dotFile(dotFilename);
        dotFile << "// size=" << size << "\n";
        dotFile << "digraph cfg {\n";
        for (uint64_t i = 0; i < size; ++i) {
            std::uniform_int_distribution<> disW(1, 20);
            int random_number = disW(gen);
            dotFile << i << "[weight=" << random_number << "];\n";
        }
        for (uint64_t i = 0; i < size; ++i) {
            for (const auto &[neighborId, edgeWeight] : adj[i]) {
                dotFile << i << "->" << neighborId << "[weight=" << edgeWeight
                        << "];\n";
            }
        }
        dotFile << "}\n";
    }
};

int main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <size> <ratio * 100>\n";
        return 1;
    }

    int size = std::stoi(argv[1]);
    double ratio = (double)std::stoi(argv[2]) / 100;
    DAGGenerator generator(size);
    int result = generator.generate(ratio);

    if (result == 0) {
        generator.printToDot("./dag.dot");
    } else {
        std::cout << "FAILED\n";
    }
    return result;
}