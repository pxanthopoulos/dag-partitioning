/**
 * @file CPSATTransferMemorySolver.cpp
 * @brief CP-SAT solver that minimizes the bytes moved between device and host
 * under a limit on the peak device footprint
 */

#include "CPSATTransferMemorySolver.h"
#include "CPSAT.h"

#include "ortools/sat/cp_model_solver.h"
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>

namespace dag_partitioning {

namespace scheduling {

namespace cpsat {

namespace {

namespace sat = operations_research::sat;

/**
 * @brief The positions of the nodes, and literals that tell which node runs
 * right after which
 */
class OrderModel {
  public:
    sat::CpModelBuilder cpModel;
    std::vector<sat::IntVar> pos;

    /**
     * @brief Positions, pruned to [|anc(i)|, n-1-|des(i)|], with producers
     * before their users, warm-started from hint
     */
    OrderModel(const hypergraph::HyperGraph &graph,
               const std::vector<uint64_t> &hint)
        : n(graph.getSize()), graph(graph), hintPos(n), successors(n),
          predecessors(n) {
        pos.reserve(n);
        for (uint64_t i = 0; i < n; i++) {
            pos.push_back(cpModel.NewIntVar(
                {static_cast<int64_t>(graph.getAncestors()[i].size()),
                 static_cast<int64_t>(n - 1 -
                                      graph.getDescendants()[i].size())}));
        }
        cpModel.AddAllDifferent(pos);

        const auto &producers = graph.getTensorProducers();
        const auto &users = graph.getTensorUsers();
        for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
            for (uint64_t user : users[t])
                cpModel.AddLessThan(pos[producers[t]], pos[user]);
        }

        for (uint64_t step = 0; step < n; step++)
            hintPos[hint[step]] = step;
        for (uint64_t i = 0; i < n; i++)
            cpModel.AddHint(pos[i], static_cast<int64_t>(hintPos[i]));
    }

    /**
     * @brief Whether after can run right after before (after = n: before can
     * run last): their positions allow it, and no node has to run between
     * them, a predecessor of after that descends from before
     */
    [[nodiscard]] bool canBeAdjacent(uint64_t before, uint64_t after) const {
        const auto &ancestors = graph.getAncestors();
        uint64_t lastBefore = n - 1 - graph.getDescendants()[before].size();
        if (after == n)
            return lastBefore == n - 1;
        uint64_t firstAfter = ancestors[after].size();
        uint64_t lastAfter = n - 1 - graph.getDescendants()[after].size();
        if (std::max<uint64_t>(ancestors[before].size() + 1, firstAfter) >
                std::min(lastBefore + 1, lastAfter) ||
            ancestors[before].count(after))
            return false;
        for (uint64_t t : graph.getInputTensors()[after]) {
            uint64_t pred = graph.getTensorProducers()[t];
            if (ancestors[pred].count(before))
                return false;
        }
        return true;
    }

    /**
     * @brief Position of a node in the hint
     */
    [[nodiscard]] uint64_t hintPosition(uint64_t node) const {
        return hintPos[node];
    }

    /**
     * @brief Whether after runs right after before in the hint; after = n
     * stands for after the last step
     */
    [[nodiscard]] bool hintAdjacent(uint64_t before, uint64_t after) const {
        return (after == n ? n : hintPos[after]) == hintPos[before] + 1;
    }

    /**
     * @brief A literal that may be true only if after runs right after
     * before (after = n: before runs last), made on first use
     */
    sat::BoolVar adjacent(uint64_t before, uint64_t after) {
        auto [it, inserted] = literals.try_emplace({before, after});
        if (inserted) {
            it->second = cpModel.NewBoolVar();
            if (after == n) {
                cpModel.AddEquality(pos[before], static_cast<int64_t>(n - 1))
                    .OnlyEnforceIf(it->second);
            } else {
                cpModel.AddEquality(pos[after], pos[before] + 1)
                    .OnlyEnforceIf(it->second);
                predecessors[after].push_back(it->second);
            }
            successors[before].push_back(it->second);
            cpModel.AddHint(it->second, hintAdjacent(before, after));
        }
        return it->second;
    }

    /**
     * @brief The literal made by adjacent
     */
    [[nodiscard]] sat::BoolVar existing(uint64_t before, uint64_t after) const {
        return literals.at({before, after});
    }

    /**
     * @brief Makes every adjacency literal also true whenever after runs
     * right after before, and gives every node that can run last one for it
     */
    void makeAdjacencyExact() {
        for (uint64_t i = 0; i < n; i++) {
            if (canBeAdjacent(i, n))
                adjacent(i, n);
        }
        for (const auto &[pair, literal] : literals) {
            const auto &[before, after] = pair;
            if (after == n) {
                cpModel.AddNotEqual(pos[before], static_cast<int64_t>(n - 1))
                    .OnlyEnforceIf(literal.Not());
            } else {
                cpModel.AddNotEqual(pos[after], pos[before] + 1)
                    .OnlyEnforceIf(literal.Not());
            }
        }
    }

    /**
     * @brief Per node, the nodes it has an adjacency literal with, before and
     * after it (after may hold n)
     */
    void neighbors(std::vector<std::vector<uint64_t>> &before,
                   std::vector<std::vector<uint64_t>> &after) const {
        before.assign(n, {});
        after.assign(n, {});
        for (const auto &[pair, literal] : literals) {
            after[pair.first].push_back(pair.second);
            if (pair.second != n)
                before[pair.second].push_back(pair.first);
        }
    }

    /**
     * @brief A node has at most one successor and one predecessor
     */
    void addAtMostOneNeighbor() {
        for (uint64_t i = 0; i < n; i++) {
            cpModel.AddAtMostOne(successors[i]);
            cpModel.AddAtMostOne(predecessors[i]);
        }
    }

  private:
    uint64_t n;
    const hypergraph::HyperGraph &graph;
    std::vector<uint64_t> hintPos;
    std::map<std::pair<uint64_t, uint64_t>, sat::BoolVar> literals;
    std::vector<std::vector<sat::BoolVar>> successors, predecessors;
};

/**
 * @brief Bytes moved: a consumer of a tensor needs a fetch unless the producer
 * or another consumer runs right before it, and any fetch needs the offload
 */
sat::LinearExpr addTransfers(OrderModel &model,
                             const hypergraph::HyperGraph &graph) {
    uint64_t n = graph.getSize();
    const auto &producers = graph.getTensorProducers();
    const auto &users = graph.getTensorUsers();
    const auto &sizes = graph.getTensorSizes();
    const auto &liveOut = graph.getLiveOut();
    sat::CpModelBuilder &cpModel = model.cpModel;

    sat::LinearExpr transfers;
    for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
        if (sizes[t] == 0)
            continue;
        std::vector<uint64_t> members = {producers[t]};
        members.insert(members.end(), users[t].begin(), users[t].end());
        std::vector<uint64_t> consumers(users[t].begin(), users[t].end());
        if (liveOut[t])
            consumers.push_back(n);
        if (consumers.empty())
            continue;

        int64_t size = static_cast<int64_t>(sizes[t]);
        sat::BoolVar offload = cpModel.NewBoolVar();
        transfers += size * offload;
        bool hintOffload = false;
        for (uint64_t consumer : consumers) {
            sat::BoolVar fetch = cpModel.NewBoolVar();
            std::vector<sat::BoolVar> covered = {fetch};
            bool hintFetch = true;
            for (uint64_t member : members) {
                if (member == consumer ||
                    !model.canBeAdjacent(member, consumer))
                    continue;
                covered.push_back(model.adjacent(member, consumer));
                hintFetch &= !model.hintAdjacent(member, consumer);
            }
            cpModel.AddBoolOr(covered);
            cpModel.AddImplication(fetch, offload);
            transfers += size * fetch;
            cpModel.AddHint(fetch, hintFetch);
            hintOffload |= hintFetch;
        }
        cpModel.AddHint(offload, hintOffload);
    }
    return transfers;
}

/**
 * @brief The live-in bytes every node holds, which footprints leave out
 */
struct HeldLiveIns {
    std::vector<sat::LinearExpr> at; // Bytes held while the node runs
    std::vector<int64_t> max;        // Most bytes it can hold
};

/**
 * @brief Node i holds every live-in whose last user runs at or after it. A
 * user, or an ancestor of one, always holds it; a descendant of every user
 * never does
 */
HeldLiveIns addHeldLiveIns(OrderModel &model,
                           const hypergraph::HyperGraph &graph,
                           const std::vector<hypergraph::LiveIn> &liveIns) {
    uint64_t n = graph.getSize();
    const auto &ancestors = graph.getAncestors();
    sat::CpModelBuilder &cpModel = model.cpModel;

    HeldLiveIns held{std::vector<sat::LinearExpr>(n),
                     std::vector<int64_t>(n, 0)};
    for (const auto &liveIn : liveIns) {
        if (liveIn.users.empty())
            continue;
        int64_t size = static_cast<int64_t>(liveIn.size);
        std::vector<sat::LinearExpr> userPositions;
        uint64_t hintLast = 0;
        for (uint64_t user : liveIn.users) {
            userPositions.push_back(model.pos[user]);
            hintLast = std::max(hintLast, model.hintPosition(user));
        }
        sat::IntVar last = cpModel.NewIntVar({0, static_cast<int64_t>(n - 1)});
        cpModel.AddMaxEquality(last, userPositions);
        cpModel.AddHint(last, static_cast<int64_t>(hintLast));
        for (uint64_t i = 0; i < n; i++) {
            bool always = false, never = true;
            for (uint64_t user : liveIn.users) {
                always |= user == i || ancestors[user].count(i);
                never &= ancestors[i].count(user) > 0;
            }
            if (always) {
                held.at[i] += size;
                held.max[i] += size;
            } else if (!never) {
                // Not held only if it runs after the last user
                sat::BoolVar isHeld = cpModel.NewBoolVar();
                cpModel.AddGreaterThan(model.pos[i], last)
                    .OnlyEnforceIf(isHeld.Not());
                cpModel.AddHint(isHeld, model.hintPosition(i) <= hintLast);
                held.at[i] += size * isHeld;
                held.max[i] += size;
            }
        }
    }
    return held;
}

/**
 * @brief The combinations of neighbors that decide the footprint of each node
 *
 * Needs exact adjacency literals (OrderModel::makeAdjacencyExact).
 */
class NeighborCombinations {
  public:
    static constexpr uint64_t none = PartitionFootprints::None;

    explicit NeighborCombinations(const OrderModel &model) : model(model) {
        model.neighbors(before, after);
    }

    /**
     * @brief Calls visit(i, pred, succ) for every combination of a node's
     * relevant neighbors: pred and succ are a node it has an adjacency literal
     * with, or none; they differ unless both are none
     */
    template <typename Visit> void forEach(Visit visit) const {
        for (uint64_t i = 0; i < before.size(); i++) {
            for (uint64_t pred : withNone(before[i])) {
                for (uint64_t succ : withNone(after[i])) {
                    if (pred != none && pred == succ)
                        continue;
                    visit(i, pred, succ);
                }
            }
        }
    }

    /**
     * @brief Literals that hold exactly when pred and succ are the relevant
     * neighbors of i
     */
    [[nodiscard]] std::vector<sat::BoolVar> literals(uint64_t i, uint64_t pred,
                                                     uint64_t succ) const {
        std::vector<sat::BoolVar> result;
        if (pred != none) {
            result.push_back(model.existing(pred, i));
        } else {
            for (uint64_t node : before[i])
                result.push_back(model.existing(node, i).Not());
        }
        if (succ != none) {
            result.push_back(model.existing(i, succ));
        } else {
            for (uint64_t node : after[i])
                result.push_back(model.existing(i, node).Not());
        }
        return result;
    }

  private:
    const OrderModel &model;
    std::vector<std::vector<uint64_t>> before, after;

    static std::vector<uint64_t> withNone(std::vector<uint64_t> neighbors) {
        neighbors.insert(neighbors.begin(), none);
        return neighbors;
    }
};

/**
 * @brief A variable at least the peak, with the hint's peak as upper bound,
 * and no lower than the smallest footprint of any node, which no order goes
 * below
 */
sat::IntVar addPeak(OrderModel &model, const NeighborCombinations &combinations,
                    const HeldLiveIns &held,
                    const PartitionFootprints &footprints,
                    const std::vector<uint64_t> &hint) {
    sat::CpModelBuilder &cpModel = model.cpModel;
    std::vector<uint64_t> smallest(model.pos.size(), UINT64_MAX);
    combinations.forEach([&](uint64_t i, uint64_t pred, uint64_t succ) {
        smallest[i] =
            std::min(smallest[i], footprints.footprint(i, pred, succ));
    });
    int64_t lowest = 0;
    for (uint64_t footprint : smallest)
        lowest = std::max(lowest, static_cast<int64_t>(footprint));

    int64_t hintPeak = static_cast<int64_t>(footprints.peak(hint));
    sat::IntVar peak = cpModel.NewIntVar({lowest, hintPeak});
    cpModel.AddHint(peak, hintPeak);
    combinations.forEach([&](uint64_t i, uint64_t pred, uint64_t succ) {
        int64_t footprint =
            static_cast<int64_t>(footprints.footprint(i, pred, succ));
        if (footprint + held.max[i] > lowest)
            cpModel.AddGreaterOrEqual(peak, held.at[i] + footprint)
                .OnlyEnforceIf(combinations.literals(i, pred, succ));
    });
    return peak;
}

/**
 * @brief Forbids the combinations of neighbors whose footprint exceeds the
 * limit, and holds no more live-ins than fit under it
 */
void addLimit(OrderModel &model, const NeighborCombinations &combinations,
              const HeldLiveIns &held, const PartitionFootprints &footprints,
              uint64_t limit) {
    sat::CpModelBuilder &cpModel = model.cpModel;
    combinations.forEach([&](uint64_t i, uint64_t pred, uint64_t succ) {
        uint64_t footprint = footprints.footprint(i, pred, succ);
        if (footprint > limit) {
            std::vector<sat::BoolVar> forbidden;
            for (const auto &literal : combinations.literals(i, pred, succ))
                forbidden.push_back(literal.Not());
            cpModel.AddBoolOr(forbidden);
        } else if (footprint + static_cast<uint64_t>(held.max[i]) > limit) {
            cpModel
                .AddLessOrEqual(held.at[i],
                                static_cast<int64_t>(limit - footprint))
                .OnlyEnforceIf(combinations.literals(i, pred, succ));
        }
    });
}

} // namespace

CPSATTransferMemorySolver::CPSATTransferMemorySolver(
    const hypergraph::HyperGraph &graph, const PartitionFootprints &footprints,
    uint64_t memoryLimit, bool debug)
    : graph(graph), footprints(footprints), memoryLimit(memoryLimit),
      debug(debug) {}

bool CPSATTransferMemorySolver::Result::found() const {
    return status == operations_research::sat::CpSolverStatus::OPTIMAL ||
           status == operations_research::sat::CpSolverStatus::FEASIBLE;
}

CPSATTransferMemorySolver::Result CPSATTransferMemorySolver::solveOnce(
    Goal goal, uint64_t limit, const std::vector<uint64_t> &hint,
    uint64_t timeLimitSeconds, int32_t numWorkers, bool deterministic) const {
    auto constructionStart = std::chrono::steady_clock::now();
    uint64_t n = graph.getSize();

    OrderModel model(graph, hint);
    sat::LinearExpr transfers = addTransfers(model, graph);
    std::optional<sat::IntVar> peak;
    if (goal == Goal::Peak || limit != UINT64_MAX) {
        // The footprints need to know who runs last and the exact neighbors,
        // else a solution could dodge a limit by leaving a literal false
        model.makeAdjacencyExact();
        NeighborCombinations combinations(model);
        HeldLiveIns held =
            addHeldLiveIns(model, graph, footprints.getLiveIns());
        if (goal == Goal::Peak)
            peak = addPeak(model, combinations, held, footprints, hint);
        else
            addLimit(model, combinations, held, footprints, limit);
    }
    model.addAtMostOneNeighbor();
    if (goal == Goal::Peak)
        model.cpModel.Minimize(*peak);
    else
        model.cpModel.Minimize(transfers);

    sat::CpModelProto proto = model.cpModel.Build();
    auto solveStart = std::chrono::steady_clock::now();
    sat::CpSolverResponse response =
        solveModel(proto, timeLimitSeconds, numWorkers, deterministic, debug);
    auto solveEnd = std::chrono::steady_clock::now();

    if (debug) {
        std::cerr
            << "Solver status: " << CpSolverResponseStats(response)
            << "\nCP-SAT model construction time: "
            << std::chrono::duration<double>(solveStart - constructionStart)
                   .count()
            << " seconds\nCP-SAT solve time: "
            << std::chrono::duration<double>(solveEnd - solveStart).count()
            << " seconds" << std::endl;
    }

    Result result{response.status(), {}, 0};
    if (!result.found())
        return result;

    result.schedule.resize(n);
    for (uint64_t i = 0; i < n; i++)
        result.schedule[SolutionIntegerValue(response, model.pos[i])] = i;
    long long objective = std::llround(response.objective_value());
    if (goal == Goal::Peak) {
        result.cost = footprints.peak(result.schedule);
        assert(objective == static_cast<long long>(result.cost) &&
               "CP-SAT objective must match the peak footprint");
    } else {
        // The objective may keep an avoidable fetch or offload in a
        // solution that is only feasible, never in an optimal one
        result.cost = graph.computeTransfers(result.schedule);
        assert(objective >= static_cast<long long>(result.cost) &&
               (response.status() != sat::CpSolverStatus::OPTIMAL ||
                objective == static_cast<long long>(result.cost)) &&
               "CP-SAT objective must match computeTransfers");
        assert((limit == UINT64_MAX ||
                footprints.peak(result.schedule) <= limit) &&
               "CP-SAT solution must meet the memory limit");
    }
    return result;
}

std::tuple<operations_research::sat::CpSolverStatus, std::vector<uint64_t>,
           uint64_t>
CPSATTransferMemorySolver::solve(uint64_t timeLimitSeconds, int32_t numWorkers,
                                 bool deterministic) const {
    std::vector<uint64_t> hint = warmStart(
        graph,
        [&](const std::vector<uint64_t> &order) {
            return graph.computeTransfers(order);
        },
        debug);

    Result result = solveOnce(Goal::Transfers, memoryLimit, hint,
                              timeLimitSeconds, numWorkers, deterministic);
    sat::CpSolverStatus status = result.status;
    if (!result.found()) {
        // Nothing meets the limit (or the time limit hit first): move as little
        // as possible at the lowest peak
        Result lowest = solveOnce(Goal::Peak, UINT64_MAX, hint,
                                  timeLimitSeconds, numWorkers, deterministic);
        if (!lowest.found()) {
            throw std::runtime_error(
                "No solution found by CP-SAT solver. Status: " +
                std::to_string(static_cast<int>(lowest.status)));
        }
        result = solveOnce(Goal::Transfers, std::max(memoryLimit, lowest.cost),
                           lowest.schedule, timeLimitSeconds, numWorkers,
                           deterministic);
        if (!result.found()) {
            result.schedule = lowest.schedule;
            result.cost = graph.computeTransfers(lowest.schedule);
        }
        status = lowest.status == sat::CpSolverStatus::OPTIMAL ? result.status
                                                               : lowest.status;
    }
    if (debug) {
        std::cerr << "Transferred bytes: " << result.cost
                  << "\nPeak footprint: " << footprints.peak(result.schedule)
                  << std::endl;
    }
    return {status, result.schedule, result.cost};
}

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning
