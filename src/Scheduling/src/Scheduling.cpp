/**
 * @file Scheduling.cpp
 * @brief Implementation of memory-aware scheduling algorithms
 */

#include "Scheduling.h"

#include <cassert>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <stdexcept>

namespace dag_partitioning {

namespace scheduling {

PartitionFootprints::PartitionFootprints(
    std::vector<std::vector<utils::Tensor>> localTensors,
    std::vector<Escaping> tensors, std::vector<std::vector<Use>> uses,
    utils::SlabPacker packer, std::vector<hypergraph::LiveIn> liveIns)
    : localTensors(std::move(localTensors)), tensors(std::move(tensors)),
      uses(std::move(uses)), packer(std::move(packer)),
      liveIns(std::move(liveIns)), slabs(this->localTensors.size()) {}

uint64_t PartitionFootprints::footprint(uint64_t partition, uint64_t pred,
                                        uint64_t succ) const {
    bool runsLast = succ == localTensors.size();
    auto isMember = [&](uint64_t tensor, uint64_t node) {
        const auto &members = tensors[tensor].members;
        return std::binary_search(members.begin(), members.end(), node);
    };

    // Kept tensors count for the whole partition
    const auto &partitionUses = uses[partition];
    std::vector<bool> kept(partitionUses.size());
    uint64_t held = 0;
    for (uint64_t u = 0; u < partitionUses.size(); u++) {
        uint64_t tensor = partitionUses[u].tensor;
        kept[u] = isMember(tensor, pred) || isMember(tensor, succ) ||
                  (runsLast && tensors[tensor].liveOut);
        if (kept[u])
            held += tensors[tensor].size;
    }
    // The last partition also holds the live-outs it doesn't use: it fetches
    // them before it ends
    if (runsLast) {
        for (uint64_t t = 0; t < tensors.size(); t++) {
            if (tensors[t].liveOut && !isMember(t, partition))
                held += tensors[t].size;
        }
    }

    // The others share the slab with the local tensors
    auto [it, inserted] = slabs[partition].try_emplace(kept, 0);
    if (inserted) {
        std::vector<utils::Tensor> slab = localTensors[partition];
        for (uint64_t u = 0; u < partitionUses.size(); u++) {
            if (!kept[u])
                slab.push_back({tensors[partitionUses[u].tensor].size,
                                partitionUses[u].first, partitionUses[u].last,
                                0});
        }
        std::sort(slab.begin(), slab.end(),
                  [](const utils::Tensor &a, const utils::Tensor &b) {
                      if (a.birth != b.birth)
                          return a.birth < b.birth;
                      return a.death < b.death;
                  });
        it->second = slab.empty() ? 0 : packer(slab);
    }
    return held + it->second;
}

uint64_t PartitionFootprints::peak(const std::vector<uint64_t> &order) const {
    // The live-in bytes held at each step: every live-in up to its last user
    std::vector<uint64_t> held(order.size() + 1, 0);
    if (!liveIns.empty()) {
        std::vector<uint64_t> pos(localTensors.size());
        for (uint64_t step = 0; step < order.size(); step++)
            pos[order[step]] = step;
        std::vector<int64_t> delta(order.size() + 1, 0);
        for (const auto &liveIn : liveIns) {
            if (liveIn.users.empty())
                continue;
            uint64_t last = 0;
            for (uint64_t user : liveIn.users)
                last = std::max(last, pos[user]);
            delta[0] += static_cast<int64_t>(liveIn.size);
            delta[last + 1] -= static_cast<int64_t>(liveIn.size);
        }
        int64_t live = 0;
        for (uint64_t step = 0; step < order.size(); step++) {
            live += delta[step];
            held[step] = static_cast<uint64_t>(live);
        }
    }
    uint64_t peak = 0;
    for (uint64_t step = 0; step < order.size(); step++) {
        peak = std::max(
            peak, held[step] +
                      footprint(order[step], step == 0 ? None : order[step - 1],
                                step + 1 == order.size() ? localTensors.size()
                                                         : order[step + 1]));
    }
    return peak;
}

namespace heuristic {

RPOScheduler::RPOScheduler(const hypergraph::HyperGraph &graph, bool debug)
    : graph(graph), debug(debug) {}

std::pair<std::vector<uint64_t>, uint64_t> RPOScheduler::solve() const {
    uint64_t n = graph.getSize();
    std::vector<uint64_t> schedule;
    schedule.reserve(n);

    // Build adjacency list for DFS (node -> children/successors)
    std::vector<std::vector<uint64_t>> children(n);
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t input : graph.getInputTensors()[i]) {
            children[graph.getTensorProducers()[input]].push_back(i);
        }
    }

    // Find root nodes (nodes with no inputs)
    std::vector<uint64_t> roots;
    for (uint64_t i = 0; i < n; i++) {
        if (graph.getInputTensors()[i].empty()) {
            roots.push_back(i);
        }
    }

    // DFS to compute reverse post-order
    std::vector<bool> visited(n, false);
    std::vector<uint64_t> postOrder;
    postOrder.reserve(n);

    // Iterative DFS using explicit stack to avoid stack overflow
    std::stack<std::tuple<uint64_t, size_t, bool>> dfsStack;

    for (uint64_t root : roots) {
        if (visited[root])
            continue;

        dfsStack.push({root, 0, false});

        while (!dfsStack.empty()) {
            auto &[node, childIdx, postVisit] = dfsStack.top();

            if (postVisit) {
                postOrder.push_back(node);
                dfsStack.pop();
                continue;
            }

            if (!visited[node]) {
                visited[node] = true;
            }

            // Find next unvisited child
            while (childIdx < children[node].size() &&
                   visited[children[node][childIdx]]) {
                childIdx++;
            }

            if (childIdx < children[node].size()) {
                uint64_t child = children[node][childIdx];
                childIdx++;
                dfsStack.push({child, 0, false});
            } else {
                postVisit = true;
            }
        }
    }

    // Reverse post-order
    for (auto it = postOrder.rbegin(); it != postOrder.rend(); ++it) {
        schedule.push_back(*it);
    }

    if (debug) {
        std::cerr << "\nRPO schedule: ";
        for (uint64_t op : schedule) {
            std::cerr << op << " ";
        }
        std::cerr << std::endl;
    }

    // Compute peak memory
    uint64_t peakMemory = graph.computePeakMemory(schedule);

    if (debug) {
        std::cerr << "RPO peak memory: " << peakMemory << std::endl;
    }

    return {schedule, peakMemory};
}

} // namespace heuristic

namespace cpsat {

CPSATMemorySolver::CPSATMemorySolver(const hypergraph::HyperGraph &graph,
                                     bool debug)
    : graph(graph), debug(debug) {
    uint64_t n = graph.getSize();
    earliestOp.resize(n);
    latestOp.resize(n);
    earliestTensor.resize(graph.getNumTensors());
    latestTensor.resize(graph.getNumTensors());
}

/**
 * Pruning:
 * (7) O[i][j] exists for j >= |anc(i)|
 * (8) T[t][j] exists for j >= |anc(prod(t))|
 * (9) O[i][j] exists for j <= n-1-|des(i)|
 * (10) T[t][j] exists for j <= max over users u of n-1-|des(u)| (n-1 for
 * live-out tensors)
 */
void CPSATMemorySolver::computePruningBounds() {
    uint64_t n = graph.getSize();
    for (uint64_t i = 0; i < n; i++) {
        earliestOp[i] = graph.getAncestors()[i].size();
        latestOp[i] = n - 1 - graph.getDescendants()[i].size();
    }
    for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
        uint64_t producer = graph.getTensorProducers()[t];
        earliestTensor[t] = earliestOp[producer];
        latestTensor[t] = graph.getLiveOut()[t] ? n - 1 : latestOp[producer];
        for (uint64_t u : graph.getTensorUsers()[t]) {
            latestTensor[t] = std::max(latestTensor[t], latestOp[u]);
        }
    }

    if (debug) {
        std::cerr << "\nPruning bounds:" << std::endl;
        for (uint64_t i = 0; i < n; i++) {
            std::cerr << "  Node " << i << ":"
                      << " earliestOp=" << earliestOp[i]
                      << ", latestOp=" << latestOp[i] << std::endl;
        }
        for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
            std::cerr << "  Tensor " << t << ":"
                      << " earliestTensor=" << earliestTensor[t]
                      << ", latestTensor=" << latestTensor[t] << std::endl;
        }
    }
}

bool CPSATMemorySolver::hasO(uint64_t i, uint64_t j) const {
    return earliestOp[i] <= j && j <= latestOp[i];
}

bool CPSATMemorySolver::hasT(uint64_t t, uint64_t j) const {
    return earliestTensor[t] <= j && j <= latestTensor[t];
}

uint64_t CPSATMemorySolver::oKey(uint64_t i, uint64_t j) const {
    return i * graph.getSize() + j;
}

uint64_t CPSATMemorySolver::tKey(uint64_t t, uint64_t j) const {
    return t * graph.getSize() + j;
}

operations_research::sat::BoolVar CPSATMemorySolver::getO(uint64_t i,
                                                          uint64_t j) const {
    return O.at(oKey(i, j));
}

operations_research::sat::BoolVar CPSATMemorySolver::getT(uint64_t t,
                                                          uint64_t j) const {
    return T.at(tKey(t, j));
}

/**
 * CP-SAT Model Formulation:
 *
 * Variables:
 * - O[i][j] ∈ {0,1}: operator i scheduled at step j
 * - T[t][j] ∈ {0,1}: tensor t in memory at step j (an operator may produce
 *   any number of tensors; p(t) is the producer of t)
 * - mem ∈ [0, memoryUpperBound]: peak memory
 *
 * Constraints:
 * (1) ExactlyOne(O[i][j] for all i) for each j -- one operator per step
 * (2) ExactlyOne(O[i][j] for all j) for each i -- each operator scheduled once
 * (3) O[i][j] => T[t][j] for all t∈IN[i]   -- inputs available (implication)
 * (4a) T[t][j] => T[t][j-1] OR O[p(t)][j]  -- tensor persistence upper bound
 * (4b) O[p(t)][j] => T[t][j]               -- tensor must exist when created
 * (4c) T[t][j] => T[t][j+1] OR OR_{c∈users} O[c][j] -- tensor freed after
 * last consumer (T[t][j] => T[t][j+1] for live-out tensors)
 * (5) T[t][0] => O[p(t)][0]       -- initial condition
 * (5b) L[l][j] => L[l][j-1], O[u][j] => L[l][j] for users u of live-in l
 *      -- a live-in tensor is on the device from the start to its last user
 * (6) Σ_t T[t][j]·S[t] + Σ_l L[l][j]·S[l] + Σᵢ ES[i]·O[i][j] ≤ mem
 *      -- memory bound
 *
 * Objective: minimize mem
 */
void CPSATMemorySolver::createVariablesAndConstraints(
    uint64_t memoryUpperBound) {
    uint64_t n = graph.getSize();
    uint64_t m = graph.getNumTensors();
    const auto &producers = graph.getTensorProducers();

    // Create O and T variables
    uint64_t numO = 0, numT = 0;
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
            if (hasO(i, j)) {
                O[oKey(i, j)] = cpModel.NewBoolVar();
                numO++;
            }
        }
    }
    for (uint64_t t = 0; t < m; t++) {
        for (uint64_t j = 0; j < n; j++) {
            if (hasT(t, j)) {
                T[tKey(t, j)] = cpModel.NewBoolVar();
                numT++;
            }
        }
    }

    // Live-in tensors: L[l][j] up to the latest step of their last user; a user
    // running at j keeps it at j and before (on the device from the start)
    const auto &liveIns = graph.getLiveIns();
    L.assign(liveIns.size(), {});
    for (uint64_t l = 0; l < liveIns.size(); l++) {
        if (liveIns[l].users.empty())
            continue;
        uint64_t latest = 0;
        for (uint64_t u : liveIns[l].users)
            latest = std::max(latest, latestOp[u]);
        for (uint64_t j = 0; j <= latest; j++) {
            L[l].push_back(cpModel.NewBoolVar());
            if (j > 0)
                cpModel.AddImplication(L[l][j], L[l][j - 1]);
        }
        for (uint64_t u : liveIns[l].users) {
            for (uint64_t j = 0; j <= latest; j++) {
                if (hasO(u, j))
                    cpModel.AddImplication(getO(u, j), L[l][j]);
            }
        }
    }

    mem = cpModel.NewIntVar({0, static_cast<int64_t>(memoryUpperBound)});

    if (debug) {
        std::cerr << "\nVariables: O=" << numO << ", T=" << numT
                  << ", total=" << (numO + numT) << std::endl;
    }

    // Constraint (1): exactly one operator per step
    for (uint64_t j = 0; j < n; j++) {
        std::vector<operations_research::sat::BoolVar> opsAtStep;
        for (uint64_t i = 0; i < n; i++) {
            if (hasO(i, j)) {
                opsAtStep.push_back(getO(i, j));
            }
        }
        cpModel.AddExactlyOne(opsAtStep);
    }

    // Constraint (2): each operator scheduled exactly once
    for (uint64_t i = 0; i < n; i++) {
        std::vector<operations_research::sat::BoolVar> stepsForOp;
        for (uint64_t j = 0; j < n; j++) {
            if (hasO(i, j)) {
                stepsForOp.push_back(getO(i, j));
            }
        }
        cpModel.AddExactlyOne(stepsForOp);
    }

    // Constraint (3): inputs must be available (O[i][j] => T[k][j])
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
            if (!hasO(i, j))
                continue;

            for (uint64_t k : graph.getInputTensors()[i]) {
                if (hasT(k, j)) {
                    cpModel.AddImplication(getO(i, j), getT(k, j));
                } else {
                    cpModel.AddEquality(getO(i, j), 0);
                }
            }
        }
    }

    // Constraint (4a): T[t][j] => T[t][j-1] OR O[p(t)][j]
    for (uint64_t t = 0; t < m; t++) {
        for (uint64_t j = 1; j < n; j++) {
            if (!hasT(t, j))
                continue;

            std::vector<operations_research::sat::BoolVar> clause;
            clause.push_back(getT(t, j).Not());

            if (hasT(t, j - 1)) {
                clause.push_back(getT(t, j - 1));
            }
            if (hasO(producers[t], j)) {
                clause.push_back(getO(producers[t], j));
            }

            cpModel.AddBoolOr(clause);
        }
    }

    // Constraint (4b): O[p(t)][j] => T[t][j]
    for (uint64_t t = 0; t < m; t++) {
        for (uint64_t j = 0; j < n; j++) {
            if (hasO(producers[t], j) && hasT(t, j)) {
                cpModel.AddImplication(getO(producers[t], j), getT(t, j));
            }
        }
    }

    // Constraint (4c): T[t][j] => still live at j+1 or a consumer runs at j
    // (so a tensor live at j has a consumer at j or later); live-out tensors
    // stay live
    for (uint64_t t = 0; t < m; t++) {
        const auto &users = graph.getTensorUsers()[t];
        for (uint64_t j = 0; j < n; j++) {
            if (!hasT(t, j))
                continue;

            std::vector<operations_research::sat::BoolVar> clause;
            clause.push_back(getT(t, j).Not());
            if (hasT(t, j + 1)) {
                clause.push_back(getT(t, j + 1));
            }
            if (graph.getLiveOut()[t]) {
                if (clause.size() > 1)
                    cpModel.AddBoolOr(clause);
                continue;
            }
            for (uint64_t c : users) {
                if (hasO(c, j)) {
                    clause.push_back(getO(c, j));
                }
            }
            // Without consumers, a tensor only lives while it is created
            if (users.empty() && hasO(producers[t], j)) {
                clause.push_back(getO(producers[t], j));
            }
            cpModel.AddBoolOr(clause);
        }
    }

    // Constraint (5): T[t][0] => O[p(t)][0]
    for (uint64_t t = 0; t < m; t++) {
        if (hasT(t, 0)) {
            if (hasO(producers[t], 0)) {
                cpModel.AddImplication(getT(t, 0), getO(producers[t], 0));
            } else {
                cpModel.AddEquality(getT(t, 0), 0);
            }
        }
    }

    // Constraint (6): memory bound at each step (exactly one operator runs)
    for (uint64_t j = 0; j < n; j++) {
        operations_research::sat::LinearExpr memExpr;
        for (uint64_t t = 0; t < m; t++) {
            if (hasT(t, j)) {
                memExpr += getT(t, j) *
                           static_cast<int64_t>(graph.getTensorSizes()[t]);
            }
        }
        for (uint64_t k = 0; k < n; k++) {
            if (hasO(k, j)) {
                memExpr +=
                    getO(k, j) * static_cast<int64_t>(graph.getExtraSizes()[k]);
            }
        }
        for (uint64_t l = 0; l < L.size(); l++) {
            if (j < L[l].size())
                memExpr += L[l][j] * static_cast<int64_t>(liveIns[l].size);
        }

        cpModel.AddLessOrEqual(memExpr, mem);
    }

    cpModel.Minimize(mem);
}

void CPSATMemorySolver::setWarmStartHints(const std::vector<uint64_t> &schedule,
                                          uint64_t peakMemory) {
    uint64_t n = graph.getSize();

    // Compute opToStep mapping
    std::vector<uint64_t> opToStep(n);
    for (uint64_t step = 0; step < n; step++) {
        opToStep[schedule[step]] = step;
    }

    // Compute tensor liveness
    std::vector<robin_hood::unordered_set<uint64_t>> live(n);
    for (uint64_t tensor = 0; tensor < graph.getNumTensors(); tensor++) {
        uint64_t birthStep = opToStep[graph.getTensorProducers()[tensor]];
        const auto &users = graph.getTensorUsers()[tensor];

        uint64_t deathStep = graph.getLiveOut()[tensor] ? n - 1 : birthStep;
        for (uint64_t user : users) {
            deathStep = std::max(deathStep, opToStep[user]);
        }

        for (uint64_t step = birthStep; step <= deathStep; step++) {
            live[step].insert(tensor);
        }
    }

    // Add hints for O variables
    for (uint64_t op = 0; op < n; op++) {
        uint64_t scheduledStep = opToStep[op];
        for (uint64_t step = 0; step < n; step++) {
            if (hasO(op, step)) {
                cpModel.AddHint(getO(op, step), step == scheduledStep ? 1 : 0);
            }
        }
    }

    // Add hints for T variables
    for (uint64_t tensor = 0; tensor < graph.getNumTensors(); tensor++) {
        for (uint64_t step = 0; step < n; step++) {
            if (hasT(tensor, step)) {
                cpModel.AddHint(getT(tensor, step),
                                live[step].count(tensor) ? 1 : 0);
            }
        }
    }

    // Add hints for the live-in tensors: on the device up to their last user
    const auto &liveIns = graph.getLiveIns();
    for (uint64_t l = 0; l < L.size(); l++) {
        uint64_t last = 0;
        for (uint64_t user : liveIns[l].users)
            last = std::max(last, opToStep[user]);
        for (uint64_t j = 0; j < L[l].size(); j++)
            cpModel.AddHint(L[l][j], j <= last ? 1 : 0);
    }

    // Add hint for mem variable
    cpModel.AddHint(mem, static_cast<int64_t>(peakMemory));
}

std::tuple<operations_research::sat::CpSolverStatus, std::vector<uint64_t>,
           uint64_t>
CPSATMemorySolver::solve(uint64_t timeLimitSeconds, int32_t numWorkers,
                         bool deterministic) {
    uint64_t n = graph.getSize();

    auto constructionStart = std::chrono::steady_clock::now();

    // Compute pruning bounds
    computePruningBounds();

    // Get heuristic schedule
    heuristic::RPOScheduler rpoScheduler(graph, debug);
    auto [heuristicSchedule, heuristicPeak] = rpoScheduler.solve();

    // Build model
    createVariablesAndConstraints(heuristicPeak);

    // Set warm start hints
    setWarmStartHints(heuristicSchedule, heuristicPeak);

    // Configure solver
    operations_research::sat::SatParameters params;
    params.set_max_time_in_seconds(static_cast<double>(timeLimitSeconds));
    params.set_num_search_workers(numWorkers);
    // Interleaved search: the same result for any number of workers, unless
    // the time limit is hit
    params.set_interleave_search(deterministic);
    params.set_log_search_progress(debug);

    if (debug) {
        std::cerr << "\nSolving with CP-SAT..." << std::endl;
    }

    // Solve
    operations_research::sat::Model model;
    model.Add(NewSatParameters(params));
    operations_research::sat::CpModelProto modelProto = cpModel.Build();

    auto solveStart = std::chrono::steady_clock::now();
    operations_research::sat::CpSolverResponse response =
        SolveCpModel(modelProto, &model);
    auto solveEnd = std::chrono::steady_clock::now();
    lastConstructionSeconds =
        std::chrono::duration<double>(solveStart - constructionStart).count();
    lastSolveSeconds =
        std::chrono::duration<double>(solveEnd - solveStart).count();

    if (debug) {
        std::cerr << "Solver status: " << CpSolverResponseStats(response)
                  << std::endl;
        std::cerr << "CP-SAT model construction time: "
                  << lastConstructionSeconds << " seconds" << std::endl;
        std::cerr << "CP-SAT solve time: " << lastSolveSeconds << " seconds"
                  << std::endl;
    }

    // Extract results
    if (response.status() ==
            operations_research::sat::CpSolverStatus::OPTIMAL ||
        response.status() ==
            operations_research::sat::CpSolverStatus::FEASIBLE) {
        uint64_t peakMemory =
            static_cast<uint64_t>(SolutionIntegerValue(response, mem));

        if (debug) {
            std::cerr
                << (response.status() ==
                            operations_research::sat::CpSolverStatus::OPTIMAL
                        ? "OPTIMAL"
                        : "FEASIBLE")
                << " solution found!" << std::endl;
            std::cerr << "Peak memory: " << peakMemory << std::endl;
        }

        std::vector<uint64_t> schedule(n, UINT64_MAX);
        for (uint64_t i = 0; i < n; i++) {
            for (uint64_t j = 0; j < n; j++) {
                if (hasO(i, j) &&
                    operations_research::sat::SolutionBooleanValue(
                        response, getO(i, j))) {
                    schedule[j] = i;
                    break;
                }
            }
        }

        if (debug) {
            std::cerr << "Schedule: ";
            for (uint64_t i = 0; i < n; i++) {
                std::cerr << schedule[i] << " ";
            }
            std::cerr << std::endl;
        }

        return {response.status(), schedule, peakMemory};
    }

    throw std::runtime_error(
        "No solution found by CP-SAT solver. Status: " +
        std::to_string(static_cast<int>(response.status())));
}

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
    namespace sat = operations_research::sat;
    auto constructionStart = std::chrono::steady_clock::now();
    uint64_t n = graph.getSize();
    const auto &producers = graph.getTensorProducers();
    const auto &users = graph.getTensorUsers();
    const auto &sizes = graph.getTensorSizes();
    const auto &liveOut = graph.getLiveOut();
    bool limitsMemory = goal == Goal::Peak || limit != UINT64_MAX;

    sat::CpModelBuilder cpModel;

    // Positions, pruned to [|anc(i)|, n-1-|des(i)|]
    std::vector<sat::IntVar> pos;
    pos.reserve(n);
    for (uint64_t i = 0; i < n; i++) {
        pos.push_back(cpModel.NewIntVar(
            {static_cast<int64_t>(graph.getAncestors()[i].size()),
             static_cast<int64_t>(n - 1 - graph.getDescendants()[i].size())}));
    }
    cpModel.AddAllDifferent(pos);

    // Producers run before their users
    for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
        for (uint64_t user : users[t])
            cpModel.AddLessThan(pos[producers[t]], pos[user]);
    }

    // Warm start from the given schedule; step n stands for after the last
    // step
    std::vector<uint64_t> hintPos(n);
    for (uint64_t step = 0; step < n; step++)
        hintPos[hint[step]] = step;
    for (uint64_t i = 0; i < n; i++)
        cpModel.AddHint(pos[i], static_cast<int64_t>(hintPos[i]));
    auto hintAdjacent = [&](uint64_t before, uint64_t after) {
        return (after == n ? n : hintPos[after]) == hintPos[before] + 1;
    };

    // adjacent[{before, after}] may be true only if after runs right after
    // before (after = n: before runs last); a node has at most one
    // successor and one predecessor
    std::map<std::pair<uint64_t, uint64_t>, sat::BoolVar> adjacent;
    std::vector<std::vector<sat::BoolVar>> successors(n), predecessors(n);
    auto getAdjacent = [&](uint64_t before, uint64_t after) {
        auto [it, inserted] = adjacent.try_emplace({before, after});
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
    };

    // A consumer of a tensor needs a fetch unless the producer or another
    // consumer runs right before it, and any fetch needs the offload
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
                if (member == consumer)
                    continue;
                covered.push_back(getAdjacent(member, consumer));
                hintFetch &= !hintAdjacent(member, consumer);
            }
            cpModel.AddBoolOr(covered);
            cpModel.AddImplication(fetch, offload);
            transfers += size * fetch;
            cpModel.AddHint(fetch, hintFetch);
            hintOffload |= hintFetch;
        }
        cpModel.AddHint(offload, hintOffload);
    }

    std::optional<sat::IntVar> peak;
    if (limitsMemory) {
        // The footprints need to know who runs last and the exact neighbors:
        // here adjacent is also true whenever after runs right after before,
        // else a solution could dodge a limit by leaving it false
        for (uint64_t i = 0; i < n; i++)
            getAdjacent(i, n);
        for (const auto &[pair, literal] : adjacent) {
            const auto &[before, after] = pair;
            if (after == n) {
                cpModel.AddNotEqual(pos[before], static_cast<int64_t>(n - 1))
                    .OnlyEnforceIf(literal.Not());
            } else {
                cpModel.AddNotEqual(pos[after], pos[before] + 1)
                    .OnlyEnforceIf(literal.Not());
            }
        }

        // The neighbors of every node that can change its footprint
        std::vector<std::vector<uint64_t>> before(n), after(n);
        for (const auto &[pair, literal] : adjacent) {
            after[pair.first].push_back(pair.second);
            if (pair.second != n)
                before[pair.second].push_back(pair.first);
        }
        const uint64_t none = PartitionFootprints::None;
        auto withNone = [&](std::vector<uint64_t> neighbors) {
            neighbors.insert(neighbors.begin(), none);
            return neighbors;
        };

        // Literals that hold exactly when pred and succ are the relevant
        // neighbors of i
        auto literalsOf = [&](uint64_t i, uint64_t pred, uint64_t succ) {
            std::vector<sat::BoolVar> literals;
            if (pred != none) {
                literals.push_back(adjacent.at({pred, i}));
            } else {
                for (uint64_t node : before[i])
                    literals.push_back(adjacent.at({node, i}).Not());
            }
            if (succ != none) {
                literals.push_back(adjacent.at({i, succ}));
            } else {
                for (uint64_t node : after[i])
                    literals.push_back(adjacent.at({i, node}).Not());
            }
            return literals;
        };

        // The live-in tensors, which footprint leaves out: node i also holds
        // every live-in whose last user runs at or after it, heldAt[i] bytes,
        // at most heldMax[i]. A user, or an ancestor of one, always holds it;
        // a descendant of every user never does
        const auto &ancestors = graph.getAncestors();
        std::vector<sat::LinearExpr> heldAt(n);
        std::vector<int64_t> heldMax(n, 0);
        for (const auto &liveIn : footprints.getLiveIns()) {
            if (liveIn.users.empty())
                continue;
            int64_t size = static_cast<int64_t>(liveIn.size);
            std::vector<sat::LinearExpr> userPositions;
            uint64_t hintLast = 0;
            for (uint64_t user : liveIn.users) {
                userPositions.push_back(pos[user]);
                hintLast = std::max(hintLast, hintPos[user]);
            }
            sat::IntVar last =
                cpModel.NewIntVar({0, static_cast<int64_t>(n - 1)});
            cpModel.AddMaxEquality(last, userPositions);
            cpModel.AddHint(last, static_cast<int64_t>(hintLast));
            for (uint64_t i = 0; i < n; i++) {
                bool always = false, never = true;
                for (uint64_t user : liveIn.users) {
                    always |= user == i || ancestors[user].count(i);
                    never &= ancestors[i].count(user) > 0;
                }
                if (always) {
                    heldAt[i] += size;
                    heldMax[i] += size;
                } else if (!never) {
                    // Not held only if it runs after the last user
                    sat::BoolVar held = cpModel.NewBoolVar();
                    cpModel.AddGreaterThan(pos[i], last)
                        .OnlyEnforceIf(held.Not());
                    cpModel.AddHint(held, hintPos[i] <= hintLast);
                    heldAt[i] += size * held;
                    heldMax[i] += size;
                }
            }
        }

        if (goal == Goal::Peak) {
            // No order goes below the smallest footprint of any node, and the
            // hint reaches its own peak
            int64_t lowest = 0;
            for (uint64_t i = 0; i < n; i++) {
                uint64_t smallest = UINT64_MAX;
                for (uint64_t pred : withNone(before[i])) {
                    for (uint64_t succ : withNone(after[i])) {
                        if (pred == none || pred != succ)
                            smallest = std::min(
                                smallest, footprints.footprint(i, pred, succ));
                    }
                }
                lowest = std::max(lowest, static_cast<int64_t>(smallest));
            }
            int64_t hintPeak = static_cast<int64_t>(footprints.peak(hint));
            peak = cpModel.NewIntVar({lowest, hintPeak});
            cpModel.AddHint(*peak, hintPeak);
            for (uint64_t i = 0; i < n; i++) {
                for (uint64_t pred : withNone(before[i])) {
                    for (uint64_t succ : withNone(after[i])) {
                        if (pred != none && pred == succ)
                            continue;
                        int64_t footprint = static_cast<int64_t>(
                            footprints.footprint(i, pred, succ));
                        if (footprint + heldMax[i] > lowest)
                            cpModel
                                .AddGreaterOrEqual(*peak, heldAt[i] + footprint)
                                .OnlyEnforceIf(literalsOf(i, pred, succ));
                    }
                }
            }
        } else {
            for (uint64_t i = 0; i < n; i++) {
                for (uint64_t pred : withNone(before[i])) {
                    for (uint64_t succ : withNone(after[i])) {
                        if (pred != none && pred == succ)
                            continue;
                        uint64_t footprint =
                            footprints.footprint(i, pred, succ);
                        if (footprint > limit) {
                            std::vector<sat::BoolVar> forbidden;
                            for (const auto &literal :
                                 literalsOf(i, pred, succ))
                                forbidden.push_back(literal.Not());
                            cpModel.AddBoolOr(forbidden);
                        } else if (footprint +
                                       static_cast<uint64_t>(heldMax[i]) >
                                   limit) {
                            // Then only as many live-ins as fit
                            cpModel
                                .AddLessOrEqual(
                                    heldAt[i],
                                    static_cast<int64_t>(limit - footprint))
                                .OnlyEnforceIf(literalsOf(i, pred, succ));
                        }
                    }
                }
            }
        }
    }

    for (uint64_t i = 0; i < n; i++) {
        cpModel.AddAtMostOne(successors[i]);
        cpModel.AddAtMostOne(predecessors[i]);
    }
    if (goal == Goal::Peak)
        cpModel.Minimize(*peak);
    else
        cpModel.Minimize(transfers);

    sat::SatParameters params;
    params.set_max_time_in_seconds(static_cast<double>(timeLimitSeconds));
    params.set_num_search_workers(numWorkers);
    // Interleaved search: the same result for any number of workers, unless
    // the time limit is hit
    params.set_interleave_search(deterministic);
    params.set_log_search_progress(debug);
    sat::Model model;
    model.Add(NewSatParameters(params));

    auto solveStart = std::chrono::steady_clock::now();
    sat::CpSolverResponse response = SolveCpModel(cpModel.Build(), &model);
    auto solveEnd = std::chrono::steady_clock::now();

    if (debug) {
        std::cerr << "Solver status: " << CpSolverResponseStats(response)
                  << std::endl;
        std::cerr << "CP-SAT model construction time: "
                  << std::chrono::duration<double>(solveStart -
                                                   constructionStart)
                         .count()
                  << " seconds" << std::endl;
        std::cerr
            << "CP-SAT solve time: "
            << std::chrono::duration<double>(solveEnd - solveStart).count()
            << " seconds" << std::endl;
    }

    Result result{response.status(), {}, 0};
    if (!result.found())
        return result;

    result.schedule.resize(n);
    for (uint64_t i = 0; i < n; i++)
        result.schedule[SolutionIntegerValue(response, pos[i])] = i;
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
    namespace sat = operations_research::sat;

    // Get heuristic schedule
    heuristic::RPOScheduler rpoScheduler(graph, debug);
    std::vector<uint64_t> hint = rpoScheduler.solve().first;

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
        std::cerr << "Transferred bytes: " << result.cost << std::endl;
        std::cerr << "Peak footprint: " << footprints.peak(result.schedule)
                  << std::endl;
    }
    return {status, result.schedule, result.cost};
}

} // namespace cpsat

namespace bruteforce {

BruteForceSolver::BruteForceSolver(const hypergraph::HyperGraph &graph,
                                   bool debug, Objective objective,
                                   const PartitionFootprints *footprints,
                                   uint64_t memoryLimit)
    : graph(graph), debug(debug), objective(objective), footprints(footprints),
      memoryLimit(memoryLimit) {}

bool BruteForceSolver::canPlace(uint64_t node,
                                const std::vector<bool> &placed) const {
    const auto &inputTensors = graph.getInputTensors();
    for (uint64_t t : inputTensors[node]) {
        if (!placed[graph.getTensorProducers()[t]])
            return false;
    }
    return true;
}

void BruteForceSolver::enumerate(std::vector<uint64_t> &current,
                                 std::vector<bool> &placed, uint64_t &bestCost,
                                 std::vector<uint64_t> &bestOrder,
                                 uint64_t &count,
                                 std::vector<std::string> &outputs,
                                 const Cost &cost) const {
    uint64_t n = graph.getSize();
    if (current.size() == n) {
        count++;
        uint64_t orderCost = cost(current);

        if (debug) {
            std::ostringstream oss;
            oss << "Order " << count - 1 << " (cost " << orderCost << "): ";
            for (const uint64_t &node : current) {
                oss << node << " ";
            }
            outputs.emplace_back(oss.str());
        }

        if (orderCost < bestCost) {
            bestCost = orderCost;
            bestOrder = current;
        }
        return;
    }

    for (uint64_t node = 0; node < n; node++) {
        if (!placed[node] && canPlace(node, placed)) {
            placed[node] = true;
            current.push_back(node);

            enumerate(current, placed, bestCost, bestOrder, count, outputs,
                      cost);

            current.pop_back();
            placed[node] = false;
        }
    }
}

std::pair<std::vector<uint64_t>, uint64_t> BruteForceSolver::solve() const {
    if (objective == Objective::Memory) {
        return best([&](const std::vector<uint64_t> &order) {
            return graph.computePeakMemory(order);
        });
    }
    if (objective != Objective::TransferMemory)
        throw std::invalid_argument("Unknown scheduling objective");

    // The lowest peak any order reaches can be above the limit
    uint64_t limit = memoryLimit;
    if (limit != UINT64_MAX) {
        uint64_t lowest = best([&](const std::vector<uint64_t> &order) {
                              return footprints->peak(order);
                          }).second;
        limit = std::max(limit, lowest);
    }
    return best([&](const std::vector<uint64_t> &order) {
        return limit == UINT64_MAX || footprints->peak(order) <= limit
                   ? graph.computeTransfers(order)
                   : UINT64_MAX;
    });
}

std::pair<std::vector<uint64_t>, uint64_t>
BruteForceSolver::best(const Cost &cost) const {
    uint64_t n = graph.getSize();
    std::vector<std::string> outputs;

    std::vector<uint64_t> current, bestOrder;
    std::vector<bool> placed(n, false);
    uint64_t bestCost = UINT64_MAX;
    uint64_t count = 0;

    enumerate(current, placed, bestCost, bestOrder, count, outputs, cost);

    if (debug) {
        std::cerr << "\nBrute-force results:\n";
        for (const auto &out : outputs) {
            std::cerr << out << std::endl;
        }
        std::cerr << "Best result: cost " << bestCost << " in order: ";
        for (const uint64_t &node : bestOrder) {
            std::cerr << node << " ";
        }
        std::cerr << std::endl;
    }

    return {bestOrder, bestCost};
}

} // namespace bruteforce

Scheduler::Scheduler(const core::Graph &originalGraph,
                     const std::vector<uint64_t> &partitionMapping, bool debug,
                     bool verify, uint64_t offsetAlignment,
                     uint64_t rangeAlignment, utils::SlabPacker packer,
                     std::vector<bool> liveOutNodes,
                     std::vector<LiveIn> liveIns)
    : originalGraph(originalGraph), partitionMapping(partitionMapping),
      debug(debug), verify(verify), offsetAlignment(offsetAlignment),
      rangeAlignment(rangeAlignment), packer(std::move(packer)),
      liveOutNodes(std::move(liveOutNodes)), liveIns(std::move(liveIns)) {
    if (!this->liveOutNodes.empty() &&
        this->liveOutNodes.size() != originalGraph.size)
        throw std::invalid_argument("liveOutNodes needs one entry per node");
    for (const LiveIn &liveIn : this->liveIns) {
        for (const auto *nodes : {&liveIn.users, &liveIn.tiedUsers}) {
            for (uint64_t node : *nodes) {
                if (node >= originalGraph.size)
                    throw std::invalid_argument("live-in user out of range");
            }
        }
    }
    if (!this->packer) {
        this->packer = [offsetAlignment,
                        rangeAlignment](std::vector<utils::Tensor> &tensors) {
            return utils::packBestFit(tensors, offsetAlignment, rangeAlignment);
        };
    }
}

std::vector<uint64_t> Scheduler::calculatePartitionWeights(
    std::vector<hypergraph::Tensor> &escapingTensors) {
    // Node ids are the order when they are topological (as a compiler emits
    // them, in block order); otherwise use a topological sort
    bool idsAreTopological = true;
    for (uint64_t from = 0; from < originalGraph.size; ++from) {
        for (const auto &[to, edgeWeight] : originalGraph.adj[from])
            idsAreTopological &= from < to;
    }
    std::vector<uint64_t> topologicalOrder(originalGraph.size);
    if (idsAreTopological)
        std::iota(topologicalOrder.begin(), topologicalOrder.end(), 0);
    else
        topologicalOrder = originalGraph.topologicalSort();

    // Create global position mapping
    std::vector<uint64_t> globalNodeToPosition(originalGraph.size);
    for (uint64_t pos = 0; pos < topologicalOrder.size(); ++pos) {
        globalNodeToPosition[topologicalOrder[pos]] = pos;
    }

    // Find number of partitions
    uint64_t numPartitions = 0;
    for (uint64_t partId : partitionMapping) {
        numPartitions = std::max(numPartitions, partId + 1);
    }

    // Group nodes by partition and sort by global topological position
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> partitionNodes(
        numPartitions);
    for (uint64_t node = 0; node < originalGraph.size; ++node) {
        uint64_t partition = partitionMapping[node];
        partitionNodes[partition].emplace_back(globalNodeToPosition[node],
                                               node);
    }

    // Sort nodes within each partition by their global topological position
    for (auto &nodes : partitionNodes) {
        std::sort(nodes.begin(), nodes.end());
    }

    // Create local position mapping: nodeToLocalPosition[node] = local position
    // within its partition
    std::vector<uint64_t> nodeToLocalPosition(originalGraph.size);
    for (uint64_t partition = 0; partition < numPartitions; ++partition) {
        for (uint64_t localPos = 0; localPos < partitionNodes[partition].size();
             ++localPos) {
            uint64_t node = partitionNodes[partition][localPos].second;
            nodeToLocalPosition[node] = localPos;
        }
    }

    // A tied op writing into a live-in starts a buffer of the live-in's size
    // (the live-in is not a node, so it would have no buffer at all)
    std::vector<uint64_t> bufferSize(originalGraph.nodeWeights.begin(),
                                     originalGraph.nodeWeights.end());
    std::vector<bool> startsBuffer(originalGraph.size, false);
    for (const LiveIn &liveIn : liveIns) {
        for (uint64_t node : liveIn.tiedUsers) {
            if (originalGraph.nodeWeights[node] == 0) {
                startsBuffer[node] = true;
                bufferSize[node] = liveIn.size;
            }
        }
    }

    // The live-ins, with the partitions of their users
    partitionLiveIns.clear();
    for (const LiveIn &liveIn : liveIns) {
        std::set<uint64_t> partitions;
        for (uint64_t user : liveIn.users)
            partitions.insert(partitionMapping[user]);
        partitionLiveIns.push_back(
            {liveIn.size,
             std::vector<uint64_t>(partitions.begin(), partitions.end())});
    }

    // Buffer of every value: a tied op (weight 0, consumed results) writes
    // into the input whose edge weight equals its result size, preferring
    // the last-defined one (the init operand)
    std::vector<uint64_t> root(originalGraph.size, UINT64_MAX);
    for (uint64_t node : topologicalOrder) {
        if (originalGraph.nodeWeights[node] > 0 || startsBuffer[node]) {
            root[node] = node;
            continue;
        }
        uint64_t resultSize = 0;
        for (const auto &[consumer, edgeWeight] : originalGraph.adj[node])
            resultSize = std::max(resultSize, edgeWeight);
        uint64_t tiedInput = UINT64_MAX;
        for (const auto &[input, edgeWeight] : originalGraph.revAdj[node]) {
            if (resultSize > 0 && edgeWeight == resultSize &&
                (tiedInput == UINT64_MAX ||
                 globalNodeToPosition[input] > globalNodeToPosition[tiedInput]))
                tiedInput = input;
        }
        if (tiedInput != UINT64_MAX)
            root[node] = root[tiedInput];
    }
    std::vector<std::vector<uint64_t>> aliases(originalGraph.size);
    for (uint64_t node = 0; node < originalGraph.size; ++node) {
        if (root[node] != UINT64_MAX)
            aliases[root[node]].push_back(node);
    }

    // Collect tensors per partition
    std::vector<std::vector<utils::Tensor>> partitionTensors(numPartitions);
    escapingTensors.clear();
    std::vector<PartitionFootprints::Escaping> footprintTensors;
    std::vector<std::vector<PartitionFootprints::Use>> footprintUses(
        numPartitions);

    for (uint64_t producer = 0; producer < originalGraph.size; ++producer) {
        if (aliases[producer].empty())
            continue;
        uint64_t producerPartition = partitionMapping[producer];

        // Find all users of the buffer's values within the same partition,
        // and the partitions it escapes to
        uint64_t lastConsumerPosition =
            nodeToLocalPosition[producer]; // Birth time
        hypergraph::Tensor escaping{
            producerPartition, bufferSize[producer], {}, false};
        // First and last use in every partition that uses it
        std::map<uint64_t, std::pair<uint64_t, uint64_t>> usedAt;
        auto use = [&](uint64_t node) {
            uint64_t local = nodeToLocalPosition[node];
            auto [it, inserted] =
                usedAt.try_emplace(partitionMapping[node], local, local);
            it->second.first = std::min(it->second.first, local);
            it->second.second = std::max(it->second.second, local);
            if (partitionMapping[node] == producerPartition) {
                lastConsumerPosition =
                    std::max(lastConsumerPosition, nodeToLocalPosition[node]);
            } else {
                escaping.users.push_back(partitionMapping[node]);
            }
        };
        for (uint64_t value : aliases[producer]) {
            use(value);
            // Graph outputs stay live until the end
            if (liveOutNodes.empty() ? originalGraph.adj[value].empty()
                                     : liveOutNodes[value])
                escaping.liveOut = true;
            for (const auto &[consumer, edgeWeight] : originalGraph.adj[value])
                use(consumer);
        }

        // Tensor is packed only if it has internal users only
        if (escaping.users.empty() && !escaping.liveOut) {
            utils::Tensor tensor;
            tensor.size = bufferSize[producer];
            tensor.birth = nodeToLocalPosition[producer];
            tensor.death = lastConsumerPosition;
            tensor.offset = 0; // Will be assigned during packing
            partitionTensors[producerPartition].push_back(tensor);
        } else {
            std::sort(escaping.users.begin(), escaping.users.end());
            escaping.users.erase(
                std::unique(escaping.users.begin(), escaping.users.end()),
                escaping.users.end());
            if (escaping.size != 0) {
                std::vector<uint64_t> members;
                for (const auto &[partition, range] : usedAt) {
                    footprintUses[partition].push_back(
                        {footprintTensors.size(), range.first, range.second});
                    members.push_back(partition);
                }
                footprintTensors.push_back(
                    {escaping.size, std::move(members), escaping.liveOut});
            }
            escapingTensors.push_back(std::move(escaping));
        }
    }

    footprints = std::make_unique<PartitionFootprints>(
        partitionTensors, std::move(footprintTensors), std::move(footprintUses),
        packer, partitionLiveIns);

    // Perform best-fit packing for each partition
    std::vector<uint64_t> partitionWeights(numPartitions, 0);

    for (uint64_t partition = 0; partition < numPartitions; ++partition) {
        auto &tensors = partitionTensors[partition];

        if (tensors.empty()) {
            continue;
        }

        // Sort tensors by birth time, then by death time (for deterministic
        // packing)
        std::sort(tensors.begin(), tensors.end(),
                  [](const utils::Tensor &a, const utils::Tensor &b) {
                      if (a.birth != b.birth)
                          return a.birth < b.birth;
                      return a.death < b.death;
                  });

        uint64_t peakMemory = packer(tensors);

        partitionWeights[partition] = peakMemory;
    }

    return partitionWeights;
}

void Scheduler::buildCoarseGraph() {
    // Find the number of partitions (assuming contiguous IDs from 0 to k-1)
    uint64_t numPartitions =
        *std::max_element(partitionMapping.begin(), partitionMapping.end()) + 1;
    assert(numPartitions != 0 && "Number of partitions must be > 0");

    std::vector<hypergraph::Tensor> tensors;
    std::vector<uint64_t> partitionWeights = calculatePartitionWeights(tensors);

    // Partition edges that no tensor orders become empty tensors (e.g. a
    // buffer allocated in A, written in place in B and read in C orders A
    // before B and C, but not B before C)
    std::set<std::pair<uint64_t, uint64_t>> ordered;
    for (const auto &tensor : tensors) {
        for (uint64_t user : tensor.users) {
            ordered.insert({tensor.producer, user});
        }
    }
    for (uint64_t from = 0; from < originalGraph.size; ++from) {
        for (const auto &[to, edgeWeight] : originalGraph.adj[from]) {
            uint64_t fromPartition = partitionMapping[from];
            uint64_t toPartition = partitionMapping[to];

            if (fromPartition != toPartition &&
                ordered.insert({fromPartition, toPartition}).second) {
                tensors.push_back({fromPartition, 0, {toPartition}, false});
            }
        }
    }

    coarseGraph = std::make_unique<hypergraph::HyperGraph>(
        std::move(partitionWeights), tensors, partitionLiveIns);
}

std::pair<std::vector<uint64_t>, uint64_t>
Scheduler::run(uint64_t timeLimitSeconds, int32_t numWorkers, Solver solver,
               Objective objective, uint64_t memoryLimit, bool deterministic) {
    auto buildStart = std::chrono::steady_clock::now();
    buildCoarseGraph();
    auto buildEnd = std::chrono::steady_clock::now();
    if (debug) {
        std::cerr
            << "Coarse graph construction time: "
            << std::chrono::duration<double>(buildEnd - buildStart).count()
            << " seconds" << std::endl;
    }
    const hypergraph::HyperGraph &hyperGraph = *coarseGraph;

    if (solver == Solver::RPO) {
        auto solveStart = std::chrono::steady_clock::now();
        heuristic::RPOScheduler rpoScheduler(hyperGraph, debug);
        auto [schedule, peakMemory] = rpoScheduler.solve();
        uint64_t cost;
        const char *costLabel;
        if (objective == Objective::Memory) {
            cost = peakMemory;
            costLabel = "Peak memory: ";
        } else if (objective == Objective::TransferMemory) {
            cost = hyperGraph.computeTransfers(schedule);
            costLabel = "Transferred bytes: ";
        } else {
            throw std::invalid_argument("Unknown scheduling objective");
        }
        auto solveEnd = std::chrono::steady_clock::now();
        if (debug) {
            std::cerr
                << "RPO solve time: "
                << std::chrono::duration<double>(solveEnd - solveStart).count()
                << " seconds" << std::endl;
            std::cerr << costLabel << cost << std::endl;
        }
        return {schedule, cost};
    }

    operations_research::sat::CpSolverStatus status;
    std::vector<uint64_t> schedule;
    uint64_t cost;
    if (objective == Objective::Memory) {
        cpsat::CPSATMemorySolver memorySolver(hyperGraph, debug);
        std::tie(status, schedule, cost) =
            memorySolver.solve(timeLimitSeconds, numWorkers, deterministic);
    } else if (objective == Objective::TransferMemory) {
        cpsat::CPSATTransferMemorySolver transferSolver(hyperGraph, *footprints,
                                                        memoryLimit, debug);
        std::tie(status, schedule, cost) =
            transferSolver.solve(timeLimitSeconds, numWorkers, deterministic);
    } else {
        throw std::invalid_argument("Unknown scheduling objective");
    }

    if (verify) {
        bruteforce::BruteForceSolver bfSolver(hyperGraph, debug, objective,
                                              footprints.get(), memoryLimit);
        auto [bfSchedule, bfCost] = bfSolver.solve();

        // Both must also stay within the limit, or the lowest peak above it
        bool peakOk = objective != Objective::TransferMemory ||
                      memoryLimit == UINT64_MAX ||
                      footprints->peak(schedule) <=
                          std::max(memoryLimit, footprints->peak(bfSchedule));
        if (cost != bfCost || !peakOk) {
            std::cerr << "\n=== VERIFICATION FAILED: DUMP START ===\n"
                      << std::endl;

            std::cerr << "=== HyperGraph ===\n";
            std::cerr << hyperGraph;
            std::cerr << std::endl;

            std::cerr << "=== CP-SAT Solution ===\n";
            std::cerr << "Cost: " << cost << std::endl;
            std::cerr << "Schedule: ";
            for (size_t i = 0; i < schedule.size(); ++i) {
                std::cerr << schedule[i] << " ";
            }
            std::cerr << "\n" << std::endl;

            std::cerr << "=== Brute-Force Solution ===\n";
            std::cerr << "Cost: " << bfCost << std::endl;
            std::cerr << "Schedule: ";
            for (size_t i = 0; i < bfSchedule.size(); ++i) {
                std::cerr << bfSchedule[i] << " ";
            }
            std::cerr << "\n" << std::endl;

            std::cerr << "\nStatus: ";
            switch (status) {
            case operations_research::sat::CpSolverStatus::OPTIMAL:
                std::cerr << "OPTIMAL" << std::endl;
                break;
            case operations_research::sat::CpSolverStatus::FEASIBLE:
                std::cerr << "FEASIBLE" << std::endl;
                break;
            case operations_research::sat::CpSolverStatus::INFEASIBLE:
                std::cerr << "INFEASIBLE" << std::endl;
                break;
            case operations_research::sat::CpSolverStatus::MODEL_INVALID:
                std::cerr << "MODEL_INVALID" << std::endl;
                break;
            case operations_research::sat::CpSolverStatus::UNKNOWN:
                std::cerr << "UNKNOWN" << std::endl;
                break;
            }

            std::cerr << "\n=== VERIFICATION FAILED: DUMP END ===\n"
                      << std::endl;

            throw std::runtime_error(
                "Verification failed: CP-SAT cost (" + std::to_string(cost) +
                ") != brute-force cost (" + std::to_string(bfCost) + ")");
        }
    }

    return {schedule, cost};
}

uint64_t Scheduler::peakFootprint(const std::vector<uint64_t> &schedule) const {
    return footprints->peak(schedule);
}

} // namespace scheduling

} // namespace dag_partitioning