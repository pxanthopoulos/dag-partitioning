#include "Graph.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <string>

using dag_partitioning::core::Graph;
using dag_partitioning::core::readDotFile;

namespace {

/**
 * @brief A DOT file with the given text in a fresh directory
 */
class DotFile : public testing::Test {
  protected:
    std::filesystem::path dir;

    void SetUp() override {
        dir =
            std::filesystem::temp_directory_path() /
            ("dag-partitioning-" + std::string(testing::UnitTest::GetInstance()
                                                   ->current_test_info()
                                                   ->name()));
        std::filesystem::create_directories(dir);
    }

    void TearDown() override { std::filesystem::remove_all(dir); }

    Graph read(const std::string &text) {
        std::ofstream(dir / "graph.dot") << text;
        return readDotFile((dir / "graph.dot").string(),
                           (dir / "mapping.txt").string());
    }
};

} // namespace

TEST_F(DotFile, ReadsNodesAndEdges) {
    Graph graph = read("digraph cfg {\n"
                       "0[weight=1];\n"
                       "1[weight=2];\n"
                       "2[weight=3];\n"
                       "0->1[weight=4];\n"
                       "0->2[weight=5];\n"
                       "1->2[weight=6];\n"
                       "}\n");
    EXPECT_EQ(graph.size, 3u);
    EXPECT_EQ(graph.nodeWeights, (std::vector<uint64_t>{1, 2, 3}));
    EXPECT_EQ(graph.totalWeight, 6u);
    EXPECT_EQ(graph.adj[0],
              (std::vector<std::pair<uint64_t, uint64_t>>{{1, 4}, {2, 5}}));
    EXPECT_EQ(graph.revAdj[2],
              (std::vector<std::pair<uint64_t, uint64_t>>{{0, 5}, {1, 6}}));
}

TEST_F(DotFile, NumbersNamedNodesInOrderAndWritesTheMapping) {
    Graph graph = read("strict digraph  {\n"
                       "  b [weight=7];\n"
                       "  a [weight=8];\n"
                       "  b -> a [weight=9];\n"
                       "}\n");
    EXPECT_EQ(graph.nodeWeights, (std::vector<uint64_t>{7, 8}));
    EXPECT_EQ(graph.adj[0],
              (std::vector<std::pair<uint64_t, uint64_t>>{{1, 9}}));
    std::ifstream mapping(dir / "mapping.txt");
    std::string text((std::istreambuf_iterator<char>(mapping)),
                     std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("b -> 0\n"), std::string::npos);
    EXPECT_NE(text.find("a -> 1\n"), std::string::npos);
}

TEST_F(DotFile, SkipsComments) {
    Graph graph = read("// size=2\n"
                       "digraph cfg {\n"
                       "0[weight=1];\n"
                       "  // 1[weight=5];\n"
                       "1[weight=2];\n"
                       "0->1[weight=3];\n"
                       "}\n");
    EXPECT_EQ(graph.nodeWeights, (std::vector<uint64_t>{1, 2}));
}

TEST_F(DotFile, SkipsLinesThatAreNeitherNodesNorEdges) {
    // Blank lines, graph attributes and a line without a weight
    Graph graph = read("digraph cfg {\n"
                       "   \n"
                       "\t\n"
                       "graph [rankdir=LR];\n"
                       "0[weight=1];\n"
                       "x;\n"
                       "1[weight=2];\n"
                       "0->1[weight=3];\n"
                       "}\n");
    EXPECT_EQ(graph.nodeWeights, (std::vector<uint64_t>{1, 2}));
}

TEST_F(DotFile, EdgesMayComeBeforeTheirNodes) {
    Graph graph = read("digraph cfg {\n"
                       "a->b[weight=3];\n"
                       "a[weight=1];\n"
                       "b[weight=2];\n"
                       "}\n");
    EXPECT_EQ(graph.adj[0],
              (std::vector<std::pair<uint64_t, uint64_t>>{{1, 3}}));
}

TEST_F(DotFile, UndeclaredNodesThrow) {
    EXPECT_THROW(read("digraph cfg {\n"
                      "a[weight=1];\n"
                      "a->b[weight=3];\n"
                      "}\n"),
                 std::runtime_error);
}

TEST_F(DotFile, ReadsWhatPrintToDotWrites) {
    Graph graph(3);
    graph.addNode(0, 5);
    graph.addNode(1, 6);
    graph.addNode(2, 7);
    graph.addEdge(0, 2, 8);
    graph.addEdge(1, 2, 9);
    graph.printToDot((dir / "printed.dot").string());
    Graph copy = readDotFile((dir / "printed.dot").string(),
                             (dir / "mapping.txt").string());
    EXPECT_EQ(copy.nodeWeights, graph.nodeWeights);
    EXPECT_EQ(copy.adj, graph.adj);
    EXPECT_EQ(copy.revAdj, graph.revAdj);
}

TEST_F(DotFile, MissingFileThrows) {
    EXPECT_THROW(readDotFile((dir / "none.dot").string(),
                             (dir / "mapping.txt").string()),
                 std::runtime_error);
}
