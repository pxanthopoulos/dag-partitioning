/**
 * @file CPSATMemorySolver.cpp
 * @brief CP-SAT solver that minimizes the peak memory of a schedule
 */

#include "CPSATMemorySolver.h"
#include "CPSAT.h"
#include "RPOScheduler.h"

#include "ortools/sat/cp_model_solver.h"
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>

namespace dag_partitioning {

namespace scheduling {

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

    if (debug) {
        std::cerr << "\nSolving with CP-SAT..." << std::endl;
    }

    operations_research::sat::CpModelProto modelProto = cpModel.Build();
    auto solveStart = std::chrono::steady_clock::now();
    operations_research::sat::CpSolverResponse response = solveModel(
        modelProto, timeLimitSeconds, numWorkers, deterministic, debug);
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

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning
