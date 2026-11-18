#include "BruteForceSolver.h"
#include "TestGraphs.h"

#include <gtest/gtest.h>

using dag_partitioning::scheduling::Objective;
using dag_partitioning::scheduling::bruteforce::BruteForceSolver;
using testgraphs::Minimal;

TEST(BruteForceSolver, LowestPeak) {
    auto graph = testgraphs::twoChains();
    auto [schedule, peak] = BruteForceSolver(graph).solve();
    EXPECT_EQ(schedule, (std::vector<uint64_t>{1, 2, 0, 3}));
    EXPECT_EQ(peak, 11u);
}

TEST(BruteForceSolver, FewestTransfersWithoutALimit) {
    Minimal minimal;
    auto [schedule, transfers] =
        BruteForceSolver(minimal.graph, Objective::TransferMemory,
                         &minimal.footprints)
            .solve();
    EXPECT_EQ(schedule, (std::vector<uint64_t>{0, 1, 2}));
    EXPECT_EQ(transfers, 8u);
}

TEST(BruteForceSolver, FewestTransfersUnderALimit) {
    Minimal minimal;
    auto [schedule, transfers] =
        BruteForceSolver(minimal.graph, Objective::TransferMemory,
                         &minimal.footprints, 11)
            .solve();
    EXPECT_EQ(schedule, (std::vector<uint64_t>{1, 0, 2}));
    EXPECT_EQ(transfers, 16u);
}

TEST(BruteForceSolver, BelowEveryOrderTakesTheLowestPeak) {
    Minimal minimal(true);
    auto [schedule, transfers] =
        BruteForceSolver(minimal.graph, Objective::TransferMemory,
                         &minimal.footprints, 1)
            .solve();
    EXPECT_EQ(schedule, (std::vector<uint64_t>{1, 0, 2}));
    EXPECT_EQ(transfers, 16u);
}
