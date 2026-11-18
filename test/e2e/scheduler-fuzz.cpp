// End-to-end check of the scheduler on random partitioned graphs, with tied
// ops, live-outs, live-ins and alignments: every CP-SAT schedule is verified
// against brute force (Scheduler's verify), meets its memory limit, and
// moves no less than the unconstrained one.
//
// Usage: scheduler-fuzz <seed> <instances> [workers]
// Prints one line per schedule (instance, solver, limit, cost, peak, order),
// so that the outputs of two builds can be compared.
#include "Graph.h"
#include "Scheduling.h"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace dag_partitioning;
using scheduling::Objective;
using scheduling::Scheduler;
using scheduling::Solver;

struct Instance {
    core::Graph graph{0};
    std::vector<uint64_t> mapping;
    std::vector<bool> liveOut;
    std::vector<Scheduler::LiveIn> liveIns;
    uint64_t offsetAlignment, rangeAlignment;
};

static Instance generate(std::mt19937_64 &rng) {
    auto pick = [&](uint64_t lo, uint64_t hi) {
        return std::uniform_int_distribution<uint64_t>(lo, hi)(rng);
    };
    uint64_t n = pick(6, 36);
    const uint64_t sizes[] = {0, 16, 64, 256, 1024, 1024, 4096};
    std::vector<uint64_t> weight(n);
    std::vector<std::vector<uint64_t>> preds(n);
    for (uint64_t j = 0; j < n; j++) {
        weight[j] = sizes[pick(0, 6)];
        if (j == 0)
            continue;
        uint64_t k = pick(j < 3 ? 0 : 1, 3);
        for (uint64_t e = 0; e < k; e++) {
            // Mostly recent producers: chains with some long edges
            uint64_t back =
                pick(0, 2) ? pick(1, std::min<uint64_t>(j, 3)) : pick(1, j);
            uint64_t p = j - back;
            if (std::find(preds[j].begin(), preds[j].end(), p) ==
                preds[j].end())
                preds[j].push_back(p);
        }
    }
    // Ids are topological, or a random relabeling
    std::vector<uint64_t> id(n);
    for (uint64_t i = 0; i < n; i++)
        id[i] = i;
    if (pick(0, 2) == 0)
        std::shuffle(id.begin(), id.end(), rng);

    Instance inst;
    inst.graph = core::Graph(n);
    // A weight-0 node writes in place into its first input: its result has
    // that input's size
    std::vector<uint64_t> resultSize(n);
    for (uint64_t j = 0; j < n; j++) {
        resultSize[j] = weight[j];
        if (weight[j] == 0 && !preds[j].empty())
            resultSize[j] = resultSize[preds[j][0]];
    }
    for (uint64_t j = 0; j < n; j++)
        inst.graph.addNode(id[j], weight[j]);
    for (uint64_t j = 0; j < n; j++)
        for (uint64_t p : preds[j])
            inst.graph.addEdge(id[p], id[j], resultSize[p]);

    // Partitions: contiguous chunks of a random topological order
    std::vector<uint64_t> indeg(n, 0), order;
    std::vector<std::vector<uint64_t>> succs(n);
    for (uint64_t j = 0; j < n; j++)
        for (uint64_t p : preds[j]) {
            succs[p].push_back(j);
            indeg[j]++;
        }
    std::vector<uint64_t> ready;
    for (uint64_t j = 0; j < n; j++)
        if (!indeg[j])
            ready.push_back(j);
    while (!ready.empty()) {
        uint64_t r = pick(0, ready.size() - 1);
        uint64_t v = ready[r];
        ready.erase(ready.begin() + r);
        order.push_back(v);
        for (uint64_t s : succs[v])
            if (--indeg[s] == 0)
                ready.push_back(s);
    }
    uint64_t parts = pick(2, std::min<uint64_t>(7, n));
    std::vector<uint64_t> cuts;
    for (uint64_t i = 1; i < n; i++)
        cuts.push_back(i);
    std::shuffle(cuts.begin(), cuts.end(), rng);
    cuts.resize(parts - 1);
    std::sort(cuts.begin(), cuts.end());
    inst.mapping.assign(n, 0);
    uint64_t part = 0, c = 0;
    for (uint64_t pos = 0; pos < n; pos++) {
        while (c < cuts.size() && cuts[c] == pos) {
            part++;
            c++;
        }
        inst.mapping[id[order[pos]]] = part;
    }

    if (pick(0, 1)) {
        inst.liveOut.assign(n, false);
        for (uint64_t i = 0; i < n; i++)
            inst.liveOut[i] = pick(0, 6) == 0;
    }
    uint64_t numLiveIns = pick(0, 3);
    for (uint64_t l = 0; l < numLiveIns; l++) {
        Scheduler::LiveIn liveIn{sizes[pick(1, 6)], {}, {}};
        uint64_t k = pick(1, 3);
        for (uint64_t u = 0; u < k; u++)
            liveIn.users.push_back(pick(0, n - 1));
        std::sort(liveIn.users.begin(), liveIn.users.end());
        liveIn.users.erase(
            std::unique(liveIn.users.begin(), liveIn.users.end()),
            liveIn.users.end());
        for (uint64_t u : liveIn.users)
            if (inst.graph.nodeWeights[u] == 0 && pick(0, 1))
                liveIn.tiedUsers.push_back(u);
        inst.liveIns.push_back(liveIn);
    }
    inst.offsetAlignment = pick(0, 1) ? 1 : 64;
    inst.rangeAlignment = pick(0, 1) ? 1 : 16;
    return inst;
}

static void report(uint64_t c, const char *what, uint64_t limit,
                   const std::vector<uint64_t> &schedule, uint64_t cost,
                   uint64_t peak) {
    std::cout << c << " " << what << " limit "
              << (limit == UINT64_MAX ? std::string("none")
                                      : std::to_string(limit))
              << " cost " << cost << " peak " << peak << " order";
    for (uint64_t p : schedule)
        std::cout << " " << p;
    std::cout << "\n";
}

int main(int argc, char **argv) {
    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <seed> <instances> [workers]\n";
        return 2;
    }
    uint64_t seed = std::strtoull(argv[1], nullptr, 10);
    uint64_t count = std::strtoull(argv[2], nullptr, 10);
    int workers = argc > 3 ? std::atoi(argv[3]) : 4;
    std::mt19937_64 rng(seed);
    uint64_t failures = 0;
    auto fail = [&](uint64_t c, const std::string &why) {
        std::cerr << "FAIL instance " << c << ": " << why << "\n";
        failures++;
    };
    for (uint64_t c = 0; c < count; c++) {
        Instance inst = generate(rng);
        // Cost and peak of a schedule
        auto run = [&](const char *what, Solver solver, Objective objective,
                       uint64_t limit) -> std::pair<uint64_t, uint64_t> {
            Scheduler s(inst.graph, inst.mapping, false, true,
                        inst.offsetAlignment, inst.rangeAlignment, {},
                        inst.liveOut, inst.liveIns);
            auto [schedule, cost] =
                s.run(600, workers, solver, objective, limit, true);
            uint64_t peak = s.peakFootprint(schedule);
            report(c, what, limit, schedule, cost, peak);
            return {cost, peak};
        };
        try {
            run("memory-cpsat", Solver::CPSAT, Objective::Memory, UINT64_MAX);
            run("memory-rpo", Solver::RPO, Objective::Memory, UINT64_MAX);
            run("transfer-rpo", Solver::RPO, Objective::TransferMemory,
                UINT64_MAX);
            auto [freeCost, freePeak] =
                run("transfer-cpsat", Solver::CPSAT, Objective::TransferMemory,
                    UINT64_MAX);
            uint64_t lowPeak = run("transfer-cpsat", Solver::CPSAT,
                                   Objective::TransferMemory, 1)
                                   .second;
            if (lowPeak > freePeak)
                fail(c, "the lowest peak is above the unconstrained one");
            std::vector<uint64_t> limits = {lowPeak, (lowPeak + freePeak) / 2,
                                            freePeak - 1, freePeak};
            std::sort(limits.begin(), limits.end());
            limits.erase(std::unique(limits.begin(), limits.end()),
                         limits.end());
            for (uint64_t limit : limits) {
                if (limit < lowPeak || limit <= 1)
                    continue;
                auto [cost, peak] = run("transfer-cpsat", Solver::CPSAT,
                                        Objective::TransferMemory, limit);
                if (peak > limit)
                    fail(c, "peak " + std::to_string(peak) + " over limit " +
                                std::to_string(limit));
                if (cost < freeCost)
                    fail(c, "moves less than without a limit");
            }
        } catch (const std::exception &e) {
            fail(c, e.what());
        }
    }
    std::cerr << count << " instances, " << failures << " failures\n";
    return failures == 0 ? 0 : 1;
}
