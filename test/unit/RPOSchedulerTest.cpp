#include "RPOScheduler.h"
#include "TestGraphs.h"

#include <gtest/gtest.h>

using dag_partitioning::scheduling::heuristic::RPOScheduler;
using dag_partitioning::scheduling::hypergraph::HyperGraph;

TEST(RPOScheduler, ReversePostOrderOfADiamond) {
    // 0 -> 1 -> 3 and 0 -> 2 -> 3: the depth-first search finishes 3, 1, 2,
    // 0
    HyperGraph graph({1, 2, 3, 4}, {{0, 8, {1, 2}, false},
                                    {1, 8, {3}, false},
                                    {2, 8, {3}, false},
                                    {3, 8, {}, true}});
    auto [schedule, peak] = RPOScheduler(graph).solve();
    EXPECT_EQ(schedule, (std::vector<uint64_t>{0, 2, 1, 3}));
    EXPECT_EQ(peak, graph.computePeakMemory(schedule));
}

TEST(RPOScheduler, VisitsEveryRoot) {
    HyperGraph graph = testgraphs::twoChains();
    auto [schedule, peak] = RPOScheduler(graph).solve();
    ASSERT_EQ(schedule.size(), 4u);
    std::vector<uint64_t> pos(4);
    for (uint64_t step = 0; step < 4; step++)
        pos[schedule[step]] = step;
    EXPECT_LT(pos[0], pos[3]);
    EXPECT_LT(pos[1], pos[2]);
    EXPECT_LT(pos[2], pos[3]);
    EXPECT_EQ(peak, graph.computePeakMemory(schedule));
}
