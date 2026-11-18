#include "CPSATTransferMemorySolver.h"
#include "TestGraphs.h"

#include <gtest/gtest.h>

using dag_partitioning::scheduling::cpsat::CPSATTransferMemorySolver;
using operations_research::sat::CpSolverStatus;
using testgraphs::Minimal;

namespace {

struct Solution {
    std::vector<uint64_t> schedule;
    uint64_t transfers;
    uint64_t peak;
};

Solution solve(const Minimal &minimal, uint64_t limit, int32_t workers = 1) {
    CPSATTransferMemorySolver solver(minimal.graph, minimal.footprints, limit);
    auto [status, schedule, transfers] = solver.solve(60, workers, true);
    EXPECT_EQ(status, CpSolverStatus::OPTIMAL);
    EXPECT_EQ(minimal.graph.computeTransfers(schedule), transfers);
    return {schedule, transfers, minimal.footprints.peak(schedule)};
}

} // namespace

TEST(CPSATTransferMemorySolver, MovesTheLeastWithoutALimit) {
    Solution s = solve(Minimal(), UINT64_MAX);
    EXPECT_EQ(s.schedule, (std::vector<uint64_t>{0, 1, 2}));
    EXPECT_EQ(s.transfers, 8u);
    EXPECT_EQ(s.peak, 12u);
}

TEST(CPSATTransferMemorySolver, ALimitThatFitsChangesNothing) {
    for (uint64_t limit : {12u, 13u}) {
        Solution s = solve(Minimal(), limit);
        EXPECT_EQ(s.transfers, 8u);
        EXPECT_EQ(s.peak, 12u);
    }
}

TEST(CPSATTransferMemorySolver, ATighterLimitMovesMore) {
    for (uint64_t limit : {8u, 11u}) {
        Solution s = solve(Minimal(), limit);
        EXPECT_EQ(s.schedule, (std::vector<uint64_t>{1, 0, 2}));
        EXPECT_EQ(s.transfers, 16u);
        EXPECT_EQ(s.peak, 8u);
    }
}

TEST(CPSATTransferMemorySolver, BelowEveryOrderTakesTheLowestPeak) {
    for (int32_t workers : {1, 4}) {
        Solution s = solve(Minimal(), 1, workers);
        EXPECT_EQ(s.schedule, (std::vector<uint64_t>{1, 0, 2}));
        EXPECT_EQ(s.transfers, 16u);
        EXPECT_EQ(s.peak, 8u);
    }
}

TEST(CPSATTransferMemorySolver, LiveInsCountTowardsTheLimit) {
    Minimal minimal(true);
    EXPECT_EQ(solve(minimal, 12).transfers, 8u);
    Solution s = solve(minimal, 9);
    EXPECT_EQ(s.transfers, 16u);
    EXPECT_EQ(s.peak, 9u);
    // No order reaches 8: the lowest peak is 9
    s = solve(minimal, 8);
    EXPECT_EQ(s.transfers, 16u);
    EXPECT_EQ(s.peak, 9u);
}
