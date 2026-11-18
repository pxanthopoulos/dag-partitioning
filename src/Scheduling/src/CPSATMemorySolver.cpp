/**
 * @file CPSATMemorySolver.cpp
 * @brief CP-SAT solver that minimizes the peak memory of a schedule
 *
 * Variables, for nodes i, tensors t (produced by p(t)), live-ins l and steps j:
 * - O[i][j]: node i runs at step j
 * - T[t][j]: tensor t is in memory at step j
 * - L[l][j]: live-in l is still in memory at step j
 * - mem: the peak
 *
 * Constraints:
 * (1) one node per step, and (2) one step per node
 * (3) O[i][j] => T[t][j] for every input t of i
 * (3b) an empty tensor holds no memory and only orders its producer before
 *      its users: pos(p(t)) < pos(u), pos(i) = sum of j O[i][j], for the
 *      users that nothing else orders after p(t)
 * (4a) T[t][j] => T[t][j-1] or O[p(t)][j]: a tensor is made, then kept
 * (4b) O[p(t)][j] => T[t][j]
 * (4c) T[t][j] => T[t][j+1] or a user of t runs at j: freed after its last
 *      user (never if live-out; at once if it has no users)
 * (5) T[t][0] => O[p(t)][0]
 * (5b) L[l][j] => L[l][j-1], and a user of l at j => L[l][j]
 * (6) per step, the tensors and live-ins held plus the extra memory of the
 *     node that runs are at most mem
 * Objective: minimize mem.
 *
 * Pruning: node i runs in [|anc(i)|, n-1-|des(i)|], and tensor t is held
 * from the earliest step of p(t) to the latest of its users (to the end if
 * live-out), surely (a constant) from the latest step of p(t) to the earliest
 * of the last user to run.
 */

#include "CPSATMemorySolver.h"
#include "CPSAT.h"

#include "ortools/sat/cp_model_solver.h"
#include <chrono>
#include <iostream>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>

namespace dag_partitioning {

namespace scheduling {

namespace cpsat {

namespace {

namespace sat = operations_research::sat;

/**
 * @brief Booleans for the steps [first, last] of something
 */
struct Window {
    uint64_t first = 0, last = 0;
    std::vector<sat::BoolVar> at;

    [[nodiscard]] bool has(uint64_t j) const { return first <= j && j <= last; }
    [[nodiscard]] sat::BoolVar operator[](uint64_t j) const {
        return at[j - first];
    }
    void create(sat::CpModelBuilder &model, uint64_t from, uint64_t to) {
        first = from;
        last = to;
        for (uint64_t j = from; j <= to; j++)
            at.push_back(model.NewBoolVar());
    }
};

/**
 * @brief The orders (producer, user) that empty tensors add: those that no
 * tensor with memory gives, and that no other predecessor of the user, a
 * descendant of the producer, implies
 */
std::vector<std::pair<uint64_t, uint64_t>>
emptyTensorOrders(const hypergraph::HyperGraph &graph) {
    const auto &producers = graph.getTensorProducers();
    const auto &users = graph.getTensorUsers();
    const auto &sizes = graph.getTensorSizes();
    std::vector<std::vector<uint64_t>> predecessors(graph.getSize());
    std::set<std::pair<uint64_t, uint64_t>> held;
    for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
        for (uint64_t user : users[t]) {
            predecessors[user].push_back(producers[t]);
            if (sizes[t] != 0)
                held.insert({producers[t], user});
        }
    }
    std::set<std::pair<uint64_t, uint64_t>> orders;
    for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
        if (sizes[t] != 0)
            continue;
        uint64_t producer = producers[t];
        for (uint64_t user : users[t]) {
            bool implied = held.count({producer, user}) > 0;
            for (uint64_t pred : predecessors[user])
                implied |= graph.getAncestors()[pred].count(producer) > 0;
            if (!implied)
                orders.insert({producer, user});
        }
    }
    return {orders.begin(), orders.end()};
}

} // namespace

CPSATMemorySolver::CPSATMemorySolver(const hypergraph::HyperGraph &graph,
                                     bool debug)
    : graph(graph), debug(debug) {}

std::tuple<operations_research::sat::CpSolverStatus, std::vector<uint64_t>,
           uint64_t>
CPSATMemorySolver::solve(uint64_t timeLimitSeconds, int32_t numWorkers,
                         bool deterministic) {
    auto constructionStart = std::chrono::steady_clock::now();
    uint64_t n = graph.getSize();
    uint64_t m = graph.getNumTensors();
    const auto &producers = graph.getTensorProducers();
    const auto &users = graph.getTensorUsers();
    const auto &liveOut = graph.getLiveOut();
    const auto &liveIns = graph.getLiveIns();

    std::vector<uint64_t> hint = warmStart(
        graph,
        [&](const std::vector<uint64_t> &order) {
            return graph.computePeakMemory(order);
        },
        debug);
    uint64_t hintPeak = graph.computePeakMemory(hint);
    std::vector<uint64_t> hintStep(n);
    for (uint64_t step = 0; step < n; step++)
        hintStep[hint[step]] = step;

    sat::CpModelBuilder model;
    std::vector<Window> O(n), T(m), L(liveIns.size());
    for (uint64_t i = 0; i < n; i++)
        O[i].create(model, graph.getAncestors()[i].size(),
                    n - 1 - graph.getDescendants()[i].size());
    for (uint64_t t = 0; t < m; t++) {
        if (graph.getTensorSizes()[t] == 0) {
            T[t].first = 1; // No steps
            continue;
        }
        uint64_t last = liveOut[t] ? n - 1 : O[producers[t]].last;
        for (uint64_t user : users[t])
            last = std::max(last, O[user].last);
        // Held for sure once the producer has run and while a user (or the
        // end, if live-out) still has to
        uint64_t heldUntil = liveOut[t] ? n - 1 : 0;
        for (uint64_t user : users[t])
            heldUntil = std::max(heldUntil, O[user].first);
        T[t].first = O[producers[t]].first;
        T[t].last = last;
        for (uint64_t j = T[t].first; j <= last; j++)
            T[t].at.push_back(O[producers[t]].last <= j && j <= heldUntil
                                  ? model.TrueVar()
                                  : model.NewBoolVar());
    }
    // (5b)
    for (uint64_t l = 0; l < liveIns.size(); l++) {
        if (liveIns[l].users.empty())
            continue;
        uint64_t last = 0;
        for (uint64_t user : liveIns[l].users)
            last = std::max(last, O[user].last);
        for (uint64_t j = 0; j <= last; j++) {
            L[l].at.push_back(model.NewBoolVar());
            if (j > 0)
                model.AddImplication(L[l][j], L[l][j - 1]);
        }
        L[l].last = last;
        for (uint64_t user : liveIns[l].users) {
            for (uint64_t j = O[user].first; j <= O[user].last; j++)
                model.AddImplication(O[user][j], L[l][j]);
        }
    }
    sat::IntVar mem = model.NewIntVar({0, static_cast<int64_t>(hintPeak)});

    // (1) and (2)
    for (uint64_t j = 0; j < n; j++) {
        std::vector<sat::BoolVar> nodes;
        for (uint64_t i = 0; i < n; i++) {
            if (O[i].has(j))
                nodes.push_back(O[i][j]);
        }
        model.AddExactlyOne(nodes);
    }
    for (uint64_t i = 0; i < n; i++)
        model.AddExactlyOne(O[i].at);

    // (3)
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = O[i].first; j <= O[i].last; j++) {
            for (uint64_t t : graph.getInputTensors()[i]) {
                if (graph.getTensorSizes()[t] == 0)
                    continue;
                if (T[t].has(j))
                    model.AddImplication(O[i][j], T[t][j]);
                else
                    model.AddEquality(O[i][j], 0);
            }
        }
    }

    // (3b)
    std::vector<std::optional<sat::IntVar>> positions(n);
    auto position = [&](uint64_t i) {
        if (!positions[i]) {
            positions[i] = model.NewIntVar({static_cast<int64_t>(O[i].first),
                                            static_cast<int64_t>(O[i].last)});
            sat::LinearExpr step;
            for (uint64_t j = O[i].first; j <= O[i].last; j++)
                step += O[i][j] * static_cast<int64_t>(j);
            model.AddEquality(*positions[i], step);
            model.AddHint(*positions[i], static_cast<int64_t>(hintStep[i]));
        }
        return *positions[i];
    };
    for (const auto &[producer, user] : emptyTensorOrders(graph))
        model.AddLessThan(position(producer), position(user));

    // (4a)
    for (uint64_t t = 0; t < m; t++) {
        const Window &producer = O[producers[t]];
        for (uint64_t j = std::max<uint64_t>(T[t].first, 1); j <= T[t].last;
             j++) {
            std::vector<sat::BoolVar> clause = {T[t][j].Not()};
            if (T[t].has(j - 1))
                clause.push_back(T[t][j - 1]);
            if (producer.has(j))
                clause.push_back(producer[j]);
            model.AddBoolOr(clause);
        }
    }

    // (4b)
    for (uint64_t t = 0; t < m; t++) {
        const Window &producer = O[producers[t]];
        for (uint64_t j = producer.first; j <= producer.last; j++) {
            if (T[t].has(j))
                model.AddImplication(producer[j], T[t][j]);
        }
    }

    // (4c)
    for (uint64_t t = 0; t < m; t++) {
        for (uint64_t j = T[t].first; j <= T[t].last; j++) {
            std::vector<sat::BoolVar> clause = {T[t][j].Not()};
            if (T[t].has(j + 1))
                clause.push_back(T[t][j + 1]);
            if (liveOut[t]) {
                if (clause.size() > 1)
                    model.AddBoolOr(clause);
                continue;
            }
            for (uint64_t user : users[t]) {
                if (O[user].has(j))
                    clause.push_back(O[user][j]);
            }
            if (users[t].empty() && O[producers[t]].has(j))
                clause.push_back(O[producers[t]][j]);
            model.AddBoolOr(clause);
        }
    }

    // (5)
    for (uint64_t t = 0; t < m; t++) {
        if (!T[t].has(0))
            continue;
        if (O[producers[t]].has(0))
            model.AddImplication(T[t][0], O[producers[t]][0]);
        else
            model.AddEquality(T[t][0], 0);
    }

    // (6)
    for (uint64_t j = 0; j < n; j++) {
        sat::LinearExpr held;
        for (uint64_t t = 0; t < m; t++) {
            if (T[t].has(j))
                held +=
                    T[t][j] * static_cast<int64_t>(graph.getTensorSizes()[t]);
        }
        for (uint64_t i = 0; i < n; i++) {
            if (O[i].has(j))
                held +=
                    O[i][j] * static_cast<int64_t>(graph.getExtraSizes()[i]);
        }
        for (uint64_t l = 0; l < L.size(); l++) {
            if (j < L[l].at.size())
                held += L[l][j] * static_cast<int64_t>(liveIns[l].size);
        }
        model.AddLessOrEqual(held, mem);
    }
    model.Minimize(mem);

    // Warm start: a tensor is held from its producer to
    // its last user (to the end if live-out)
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = O[i].first; j <= O[i].last; j++)
            model.AddHint(O[i][j], j == hintStep[i]);
    }
    for (uint64_t t = 0; t < m; t++) {
        uint64_t birth = hintStep[producers[t]];
        uint64_t death = liveOut[t] ? n - 1 : birth;
        for (uint64_t user : users[t])
            death = std::max(death, hintStep[user]);
        for (uint64_t j = T[t].first; j <= T[t].last; j++) {
            if (!(T[t][j] == model.TrueVar()))
                model.AddHint(T[t][j], birth <= j && j <= death);
        }
    }
    for (uint64_t l = 0; l < L.size(); l++) {
        uint64_t last = 0;
        for (uint64_t user : liveIns[l].users)
            last = std::max(last, hintStep[user]);
        for (uint64_t j = 0; j < L[l].at.size(); j++)
            model.AddHint(L[l][j], j <= last);
    }
    model.AddHint(mem, static_cast<int64_t>(hintPeak));

    sat::CpModelProto proto = model.Build();
    auto solveStart = std::chrono::steady_clock::now();
    sat::CpSolverResponse response =
        solveModel(proto, timeLimitSeconds, numWorkers, deterministic, debug);
    auto solveEnd = std::chrono::steady_clock::now();
    lastConstructionSeconds =
        std::chrono::duration<double>(solveStart - constructionStart).count();
    lastSolveSeconds =
        std::chrono::duration<double>(solveEnd - solveStart).count();
    if (debug) {
        std::cerr << "Solver status: " << CpSolverResponseStats(response)
                  << "\nCP-SAT model construction time: "
                  << lastConstructionSeconds << " seconds"
                  << "\nCP-SAT solve time: " << lastSolveSeconds << " seconds"
                  << std::endl;
    }

    if (response.status() != sat::CpSolverStatus::OPTIMAL &&
        response.status() != sat::CpSolverStatus::FEASIBLE) {
        throw std::runtime_error(
            "No solution found by CP-SAT solver. Status: " +
            std::to_string(static_cast<int>(response.status())));
    }
    uint64_t peak = static_cast<uint64_t>(SolutionIntegerValue(response, mem));
    std::vector<uint64_t> schedule(n);
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = O[i].first; j <= O[i].last; j++) {
            if (SolutionBooleanValue(response, O[i][j]))
                schedule[j] = i;
        }
    }
    if (debug)
        std::cerr << "Peak memory: " << peak << std::endl;
    return {response.status(), schedule, peak};
}

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning
