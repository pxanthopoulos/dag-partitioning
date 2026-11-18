#include "CPSATMemorySolver.h"
#include "BruteForceSolver.h"
#include "TestGraphs.h"

#include <gtest/gtest.h>

using dag_partitioning::scheduling::bruteforce::BruteForceSolver;
using dag_partitioning::scheduling::cpsat::CPSATMemorySolver;
using dag_partitioning::scheduling::hypergraph::HyperGraph;
using operations_research::sat::CpSolverStatus;

TEST(CPSATMemorySolver, FindsTheLowestPeak) {
    HyperGraph graph = testgraphs::twoChains();
    CPSATMemorySolver solver(graph);
    auto [status, schedule, peak] = solver.solve(60, 1, true);
    EXPECT_EQ(status, CpSolverStatus::OPTIMAL);
    EXPECT_EQ(schedule, (std::vector<uint64_t>{1, 2, 0, 3}));
    EXPECT_EQ(peak, 11u);
    EXPECT_EQ(graph.computePeakMemory(schedule), peak);
}

TEST(CPSATMemorySolver, CountsExtraMemoryLiveOutsAndLiveIns) {
    // 0 and 1 are independent; 1 has 50 bytes of extra memory, 0 makes a
    // live-out of 30, and a live-in of 40 is read by 1
    HyperGraph graph({0, 50}, {{0, 30, {}, true}}, {{40, {1}}});
    CPSATMemorySolver solver(graph);
    auto [status, schedule, peak] = solver.solve(60, 1, true);
    EXPECT_EQ(status, CpSolverStatus::OPTIMAL);
    // Running 1 first holds the live-in and its extra memory, not the
    // live-out: 90 bytes, against 120 the other way
    EXPECT_EQ(schedule, (std::vector<uint64_t>{1, 0}));
    EXPECT_EQ(peak, 90u);
    EXPECT_EQ(graph.computePeakMemory({0, 1}), 120u);
}

TEST(CPSATMemorySolver, MatchesBruteForceForAnyNumberOfWorkers) {
    HyperGraph graph({5, 0, 7, 0, 3}, {{0, 16, {1, 3}, false},
                                       {1, 8, {4}, false},
                                       {2, 32, {3}, false},
                                       {3, 4, {4}, false},
                                       {4, 2, {}, true}});
    uint64_t best = BruteForceSolver(graph).solve().second;
    for (int32_t workers : {1, 4}) {
        CPSATMemorySolver solver(graph);
        auto [status, schedule, peak] = solver.solve(60, workers, true);
        EXPECT_EQ(status, CpSolverStatus::OPTIMAL);
        EXPECT_EQ(peak, best);
        EXPECT_EQ(graph.computePeakMemory(schedule), best);
    }
}
