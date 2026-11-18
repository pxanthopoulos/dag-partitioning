#include "Scheduler.h"
#include "Graph.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <stdexcept>
#include <tuple>

using dag_partitioning::core::Graph;
using dag_partitioning::scheduling::Objective;
using dag_partitioning::scheduling::Scheduler;
using dag_partitioning::scheduling::Solver;
using dag_partitioning::scheduling::hypergraph::HyperGraph;

namespace {

/**
 * @brief Gives the tests the coarse graph
 */
class CoarseScheduler : public Scheduler {
  public:
    using Scheduler::Scheduler;

    const HyperGraph &coarse() {
        buildCoarseGraph();
        return *coarseGraph;
    }
};

struct Edge {
    uint64_t from, to, weight;
};

Graph makeGraph(const std::vector<uint64_t> &weights,
                const std::vector<Edge> &edges) {
    Graph graph(weights.size());
    for (uint64_t node = 0; node < weights.size(); node++)
        graph.addNode(node, weights[node]);
    for (const Edge &edge : edges)
        graph.addEdge(edge.from, edge.to, edge.weight);
    return graph;
}

// The tensors of a coarse graph as (producer, size, sorted users, live-out)
using TensorTuple = std::tuple<uint64_t, uint64_t, std::vector<uint64_t>, bool>;

std::vector<TensorTuple> tensorsOf(const HyperGraph &graph) {
    std::vector<TensorTuple> tensors;
    for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
        const auto &users = graph.getTensorUsers()[t];
        std::vector<uint64_t> sorted(users.begin(), users.end());
        std::sort(sorted.begin(), sorted.end());
        tensors.emplace_back(graph.getTensorProducers()[t],
                             graph.getTensorSizes()[t], sorted,
                             graph.getLiveOut()[t]);
    }
    return tensors;
}

} // namespace

TEST(Scheduler, EscapingBuffersBecomeTensors) {
    // 0 (8 bytes) -> 1 (4 bytes), in two partitions; 1 has no consumers, so
    // it is a graph output
    Graph graph = makeGraph({8, 4}, {{0, 1, 8}});
    std::vector<uint64_t> mapping = {0, 1};
    CoarseScheduler scheduler(graph, mapping);
    const HyperGraph &coarse = scheduler.coarse();
    EXPECT_EQ(coarse.getExtraSizes(), (std::vector<uint64_t>{0, 0}));
    EXPECT_EQ(tensorsOf(coarse),
              (std::vector<TensorTuple>{{0, 8, {1}, false}, {1, 4, {}, true}}));
}

TEST(Scheduler, LocalBuffersArePackedIntoTheWeight) {
    // 0 (10) -> 1 (20) -> 2 (5) in partition 0, 2 -> 3 in partition 1: 0 and
    // 1 overlap at 1
    Graph graph =
        makeGraph({10, 20, 5, 1}, {{0, 1, 10}, {1, 2, 20}, {2, 3, 5}});
    std::vector<uint64_t> mapping = {0, 0, 0, 1};
    EXPECT_EQ(CoarseScheduler(graph, mapping).coarse().getExtraSizes(),
              (std::vector<uint64_t>{30, 0}));
    // Offsets aligned to 64 and sizes to 16: 16 at 0, 32 at 64
    EXPECT_EQ(CoarseScheduler(graph, mapping, false, false, 64, 16)
                  .coarse()
                  .getExtraSizes(),
              (std::vector<uint64_t>{96, 0}));
    // A packer of the caller's
    auto sum =
        [](std::vector<dag_partitioning::scheduling::utils::Tensor> &tensors) {
            uint64_t total = 0;
            for (const auto &tensor : tensors)
                total += tensor.size;
            return total;
        };
    EXPECT_EQ(CoarseScheduler(graph, mapping, false, false, 1, 1, sum)
                  .coarse()
                  .getExtraSizes(),
              (std::vector<uint64_t>{30, 0}));
}

TEST(Scheduler, TiedOpsShareTheirInputBuffer) {
    // 1 has weight 0 and writes into 0's 8-byte buffer, which 2 reads in
    // another partition: one tensor, from partition 0
    Graph graph = makeGraph({8, 0, 4}, {{0, 1, 8}, {1, 2, 8}});
    std::vector<uint64_t> mapping = {0, 0, 1};
    EXPECT_EQ(tensorsOf(CoarseScheduler(graph, mapping).coarse()),
              (std::vector<TensorTuple>{{0, 8, {1}, false}, {1, 4, {}, true}}));
}

TEST(Scheduler, EdgesThatNoTensorOrdersBecomeEmptyTensors) {
    // A buffer made in 0, written in place in 1 and read in 2 orders 0 before
    // 1 and 2, but not 1 before 2
    Graph graph = makeGraph({8, 0, 4}, {{0, 1, 8}, {1, 2, 8}});
    std::vector<uint64_t> mapping = {0, 1, 2};
    EXPECT_EQ(tensorsOf(CoarseScheduler(graph, mapping).coarse()),
              (std::vector<TensorTuple>{{0, 8, {1, 2}, false},
                                        {2, 4, {}, true},
                                        {1, 0, {2}, false}}));
}

TEST(Scheduler, LiveOutNodesAreTheGraphOutputs) {
    // 0 is an output although 1 consumes it; 1 is not, although nothing does
    Graph graph = makeGraph({8, 4}, {{0, 1, 8}});
    std::vector<uint64_t> mapping = {0, 1};
    CoarseScheduler scheduler(graph, mapping, false, false, 1, 1, {},
                              {true, false});
    EXPECT_EQ(tensorsOf(scheduler.coarse()),
              (std::vector<TensorTuple>{{0, 8, {1}, true}}));
    // 1's buffer stays in partition 1
    EXPECT_EQ(scheduler.coarse().getExtraSizes(),
              (std::vector<uint64_t>{0, 4}));
}

TEST(Scheduler, LiveInsAndTheirTiedUsers) {
    // A live-in of 16 bytes is read by 0 and written in place by 1 (weight
    // 0), whose result 2 reads in another partition
    Graph graph = makeGraph({4, 0, 2}, {{0, 1, 4}, {1, 2, 16}});
    std::vector<uint64_t> mapping = {0, 0, 1};
    CoarseScheduler scheduler(graph, mapping, false, false, 1, 1, {}, {},
                              {{16, {0, 1}, {1}}});
    const HyperGraph &coarse = scheduler.coarse();
    ASSERT_EQ(coarse.getLiveIns().size(), 1u);
    EXPECT_EQ(coarse.getLiveIns()[0].size, 16u);
    EXPECT_EQ(coarse.getLiveIns()[0].users, (std::vector<uint64_t>{0}));
    // 1 starts a buffer of the live-in's size
    EXPECT_EQ(tensorsOf(coarse), (std::vector<TensorTuple>{{0, 16, {1}, false},
                                                           {1, 2, {}, true}}));
}

TEST(Scheduler, NodeIdsNeedNotBeTopological) {
    // The chain 2 -> 0 -> 1, split after 0
    Graph graph = makeGraph({8, 4, 16}, {{2, 0, 16}, {0, 1, 8}});
    std::vector<uint64_t> mapping = {0, 1, 0};
    CoarseScheduler scheduler(graph, mapping);
    const HyperGraph &coarse = scheduler.coarse();
    // 2's buffer lives until 0, which runs after it in partition 0
    EXPECT_EQ(coarse.getExtraSizes(), (std::vector<uint64_t>{16, 0}));
    EXPECT_EQ(tensorsOf(coarse),
              (std::vector<TensorTuple>{{0, 8, {1}, false}, {1, 4, {}, true}}));
}

TEST(Scheduler, RejectsInvalidArguments) {
    Graph graph = makeGraph({8, 4}, {{0, 1, 8}});
    std::vector<uint64_t> mapping = {0, 1};
    EXPECT_THROW(Scheduler(graph, mapping, false, false, 1, 1, {}, {true}),
                 std::invalid_argument);
    EXPECT_THROW(
        Scheduler(graph, mapping, false, false, 1, 1, {}, {}, {{4, {2}, {}}}),
        std::invalid_argument);
}

TEST(Scheduler, RunsEverySolverAndObjective) {
    // Two chains joined at the end, one partition per node; every CP-SAT
    // result is verified against brute force
    Graph graph =
        makeGraph({10, 10, 1, 1, 2},
                  {{0, 3, 10}, {1, 2, 10}, {2, 3, 1}, {3, 4, 1}, {0, 4, 10}});
    std::vector<uint64_t> mapping = {0, 1, 2, 3, 4};
    std::vector<std::vector<uint64_t>> orders;
    std::vector<uint64_t> order = {0, 1, 2, 3, 4};
    do {
        std::vector<uint64_t> pos(5);
        for (uint64_t step = 0; step < 5; step++)
            pos[order[step]] = step;
        if (pos[0] < pos[3] && pos[1] < pos[2] && pos[2] < pos[3] &&
            pos[3] < pos[4])
            orders.push_back(order);
    } while (std::next_permutation(order.begin(), order.end()));

    for (Objective objective : {Objective::Memory, Objective::TransferMemory}) {
        for (Solver solver : {Solver::CPSAT, Solver::RPO}) {
            Scheduler scheduler(graph, mapping, false, true);
            std::vector<uint64_t> schedule;
            EXPECT_NO_THROW(
                schedule =
                    scheduler.run(60, 2, solver, objective, UINT64_MAX, true)
                        .first);
            EXPECT_NE(std::find(orders.begin(), orders.end(), schedule),
                      orders.end());
        }
    }

    // A 1-byte limit gives the lowest peak of any order
    CoarseScheduler scheduler(graph, mapping, false, true);
    auto [schedule, transfers] =
        scheduler.run(60, 2, Solver::CPSAT, Objective::TransferMemory, 1, true);
    uint64_t lowest = UINT64_MAX;
    for (const auto &candidate : orders)
        lowest = std::min(lowest, scheduler.peakFootprint(candidate));
    EXPECT_EQ(scheduler.peakFootprint(schedule), lowest);
}
