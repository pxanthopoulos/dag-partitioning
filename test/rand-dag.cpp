#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <set>
#include <stack>
#include <vector>

class DAGGenerator {
  private:
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> adj;
    uint64_t size;
    std::random_device rd;
    std::mt19937 gen;
    uint64_t debug;

    [[nodiscard]] bool hasCycle() {
        std::vector<uint64_t> visited(size, 0);
        std::vector<uint64_t> recursionStack(size, 0);

        for (uint64_t i = 0; i < size; i++)
            if (dfs(i, visited, recursionStack))
                return true;
        return false;
    }

    [[nodiscard]] bool dfs(uint64_t v, std::vector<uint64_t> &visited,
                           std::vector<uint64_t> &recursionStack) {
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

    [[nodiscard]] bool isConnected() {
        std::vector<bool> visited(size, false);
        std::stack<uint64_t> s;
        s.push(0);
        visited[0] = true;
        uint64_t count = 1;

        while (!s.empty()) {
            uint64_t v = s.top();
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
        std::uniform_int_distribution<uint64_t> dis(0, size - 1);

        // Start with node 0
        used[0] = true;
        uint64_t usedCount = 1;

        uint64_t last10perc = 0, last5perc = 0;
        // Connect remaining nodes
        while (usedCount < size) {
            if (debug > 2) {
                uint64_t percentage = (usedCount * 100 / size) / 5 * 5;
                if (percentage != last5perc) {
                    std::cout << "At " << percentage << "% of size ..."
                              << std::endl;
                    last5perc = percentage;
                }
            } else if (debug > 1) {
                uint64_t percentage = (usedCount * 100 / size) / 10 * 10;
                if (percentage != last10perc) {
                    std::cout << "At " << percentage << "% of size ..."
                              << std::endl;
                    last10perc = percentage;
                }
            }

            uint64_t from = dis(gen);
            if (!used[from])
                continue;

            uint64_t to = dis(gen);
            if (!used[to] && from != to) {
                adj[from].emplace_back(to, 1);
                used[to] = true;
                usedCount++;
            }
        }
    }

    void createLongSpanningTree() {
        std::vector<bool> used(size, false);
        std::uniform_int_distribution<uint64_t> dis(0, size - 1);

        // Track recently added nodes
        std::deque<uint64_t> recent;
        used[0] = true;
        recent.push_back(0);
        uint64_t usedCount = 1;

        while (usedCount < size) {
            // 80% chance to use recent node as source
            uint64_t from;
            if (dis(gen) % 100 < 80 && !recent.empty()) {
                from = recent.back();
            } else {
                from = dis(gen);
                if (!used[from])
                    continue;
            }

            uint64_t to = dis(gen);
            if (!used[to] && from != to) {
                adj[from].emplace_back(to, 1);
                used[to] = true;
                usedCount++;
                recent.push_back(to);
                if (recent.size() > 10)
                    recent.pop_front();
            }
        }
    }

    [[nodiscard]] std::vector<uint64_t> calculateInDegrees() const {
        std::vector<uint64_t> inDegree(size, 0);
        for (uint64_t i = 0; i < size; i++) {
            for (const auto &[neighbor, weight] : adj[i]) {
                inDegree[neighbor]++;
            }
        }
        return inDegree;
    }

    [[nodiscard]] std::vector<uint64_t> topologicalSort() const {
        std::vector<uint64_t> order(size);
        std::vector<uint64_t> localInDegree = calculateInDegrees();
        std::queue<uint64_t> q;
        uint64_t pos = 0;

        for (uint64_t i = 0; i < size; i++) {
            if (localInDegree[i] == 0)
                q.push(i);
        }

        while (!q.empty()) {
            uint64_t curr = q.front();
            q.pop();
            order[curr] = pos++;

            for (const auto &[next, weight] : adj[curr]) {
                if (--localInDegree[next] == 0)
                    q.push(next);
            }
        }
        return order;
    }

    [[nodiscard]] std::vector<uint64_t> computeTopLevels() const {
        std::vector<uint64_t> localInDegree = calculateInDegrees();
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

            for (const auto &[next, weight] : adj[curr]) {
                topLevels[next] =
                    std::max(topLevels[next], topLevels[curr] + 1);
                localInDegree[next]--;
                if (localInDegree[next] == 0) {
                    q.push(next);
                }
            }
        }

        return topLevels;
    }

    static void printDistribution(const std::vector<uint64_t> &values) {
        std::map<uint64_t, uint64_t> freq;
        for (auto v : values)
            freq[v]++;

        uint64_t maxFreq = 0;
        for (const auto &[_, count] : freq)
            maxFreq = std::max(maxFreq, count);

        const int width = 50;

        uint64_t maxVal = freq.rbegin()->first;
        uint64_t padding = std::to_string(maxVal).length();

        for (const auto &[value, count] : freq) {
            int bars = static_cast<int>((count * width) / maxFreq);
            std::cout << std::setw((int)padding) << value << " |"
                      << std::string(bars, '#') << "\n";
        }
    }

  public:
    explicit DAGGenerator(uint64_t n, uint64_t debug)
        : size(n), gen(rd()), debug(debug) {
        adj.resize(n);
    }

    int generate(double edgeRatio) {
        std::uniform_int_distribution<uint64_t> dis(0, size - 1);
        std::set<std::pair<uint64_t, uint64_t>> edges;

        if (debug > 0)
            std::cout << "Generating spanning tree ..." << std::endl;
        // First create a spanning tree to ensure connectivity
        createLongSpanningTree();
        if (debug > 0)
            std::cout << "Finished generating spanning tree, copying edges ..."
                      << std::endl;
        // Copy existing edges to set
        for (uint64_t i = 0; i < size; i++) {
            for (const auto &edge : adj[i]) {
                edges.insert({i, edge.first});
            }
        }

        auto toporder = topologicalSort();

        // Add additional random edges
        auto targetEdges = (uint64_t)((double)size * edgeRatio);
        uint64_t attempts = 0;
        uint64_t maxAttempts = size * size;
        uint64_t last2perc = 0, last5perc = 0, last10perc = 0;
        uint64_t last2percatt = 0, last5percatt = 0, last10percatt = 0;

        if (debug > 0)
            std::cout << "Adding additional edges ..." << std::endl;
        while (edges.size() < targetEdges && attempts < maxAttempts) {
            if (debug > 3) {
                uint64_t percentage =
                    (edges.size() * 100 / targetEdges) / 2 * 2;
                if (percentage != last2perc) {
                    std::cout << "At " << percentage << "% of target edges ..."
                              << std::endl;
                    last2perc = percentage;
                }
            } else if (debug > 2) {
                uint64_t percentage =
                    (edges.size() * 100 / targetEdges) / 5 * 5;
                if (percentage != last5perc) {
                    std::cout << "At " << percentage << "% of target edges ..."
                              << std::endl;
                    last5perc = percentage;
                }
            } else if (debug > 1) {
                uint64_t percentage =
                    (edges.size() * 100 / targetEdges) / 10 * 10;
                if (percentage != last10perc) {
                    std::cout << "At " << percentage << "% of target edges ..."
                              << std::endl;
                    last10perc = percentage;
                }
            }

            if (debug > 6) {
                uint64_t percentage = (attempts * 100 / maxAttempts) / 2 * 2;
                if (percentage != last2percatt) {
                    std::cout << "At " << percentage << "% of max attempts ..."
                              << std::endl;
                    last2percatt = percentage;
                }
            } else if (debug > 5) {
                uint64_t percentage = (attempts * 100 / maxAttempts) / 5 * 5;
                if (percentage != last5percatt) {
                    std::cout << "At " << percentage << "% of max attempts ..."
                              << std::endl;
                    last5percatt = percentage;
                }
            } else if (debug > 4) {
                uint64_t percentage = (attempts * 100 / maxAttempts) / 10 * 10;
                if (percentage != last10percatt) {
                    std::cout << "At " << percentage << "% of max attempts ..."
                              << std::endl;
                    last10percatt = percentage;
                }
            }

            uint64_t from = dis(gen);
            uint64_t to = dis(gen);

            if (toporder[from] < toporder[to] &&
                edges.find({from, to}) == edges.end()) {
                edges.insert({from, to});
            }
            attempts++;
        }

        if (debug > 0)
            std::cout << "Rebuilding graph ..." << std::endl;

        // Rebuild adjacency list
        adj.clear();
        adj.resize(size);
        for (const auto &edge : edges) {
            std::uniform_int_distribution<> disW(1, 20);
            uint64_t random_number = disW(gen);
            adj[edge.first].emplace_back(edge.second, random_number);
        }

        if (debug > 0)
            std::cout << "Checking for cycles ..." << std::endl;
        if (hasCycle()) {
            std::cout << "CYCLIC\n";
            return 1;
        }

        if (debug > 0)
            std::cout << "Checking if connected ..." << std::endl;
        if (!isConnected()) {
            std::cout << "NOT CONNECTED\n";
            return 1;
        }

        if (debug > 0)
            printDistribution(computeTopLevels());

        return 0;
    }

    void printToDot(const std::string &dotFilename) {
        std::ofstream dotFile(dotFilename);
        dotFile << "// size=" << size << "\n";
        dotFile << "digraph cfg {\n";
        for (uint64_t i = 0; i < size; ++i) {
            std::uniform_int_distribution<> disW(0, 20);
            uint64_t random_number = disW(gen);
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
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0]
                  << " <size> <ratio * 100> <debug> <output_filename>\n";
        std::cerr << "  size: number of nodes in the graph\n";
        std::cerr << "  ratio * 100: edge ratio multiplied by 100 (e.g., 150 "
                     "for 1.5x ratio)\n";
        std::cerr << "  debug: debug level (0-7)\n";
        std::cerr << "  output_filename: path for the output .dot file\n";
        return 1;
    }

    uint64_t size = std::stoi(argv[1]);
    double ratio = (double)std::stoi(argv[2]) / 100;
    uint64_t debug = std::stoi(argv[3]);
    std::string outputFilename = argv[4];

    DAGGenerator generator(size, debug);
    int result = generator.generate(ratio);

    if (result == 0) {
        generator.printToDot(outputFilename);
        if (debug > 0) {
            std::cout << "Successfully generated DAG with " << size
                      << " nodes and wrote to " << outputFilename << std::endl;
        }
    } else {
        std::cout << "FAILED\n";
    }
    return result;
}