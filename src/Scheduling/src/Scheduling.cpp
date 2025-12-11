/**
 * @file Scheduling.cpp
 * @brief Implementation of memory-aware scheduling algorithms
 */

#include "Scheduling.h"

#include <cassert>

namespace dag_partitioning {

namespace scheduling {

namespace packing {

uint64_t packBestFit(std::vector<Tensor> &tensors) {
    auto overlaps = [](const Tensor &a, const Tensor &b) {
        return !(a.death < b.birth || b.death < a.birth);
    };

    auto verify = [&overlaps](const std::vector<Tensor> &tensors) {
        for (size_t i = 0; i < tensors.size(); ++i) {
            for (size_t j = i + 1; j < tensors.size(); ++j) {
                if (overlaps(tensors[i], tensors[j])) {
                    uint64_t end1 = tensors[i].offset + tensors[i].size;
                    uint64_t end2 = tensors[j].offset + tensors[j].size;
                    if (!(end1 <= tensors[j].offset ||
                          end2 <= tensors[i].offset)) {
                        std::cerr << "CONFLICT: tensor " << i << " and " << j
                                  << "\n";
                        return false;
                    }
                }
            }
        }
        return true;
    };

    if (tensors.empty())
        return 0;

    for (size_t i = 0; i < tensors.size(); ++i) {
        Tensor &t = tensors[i];

        std::vector<std::pair<uint64_t, uint64_t>> occupied;
        for (size_t j = 0; j < i; ++j) { // all previously placed
            if (overlaps(t, tensors[j])) {
                occupied.push_back(
                    {tensors[j].offset, tensors[j].offset + tensors[j].size});
            }
        }

        // Sort and merge overlapping intervals
        std::sort(occupied.begin(), occupied.end());
        std::vector<std::pair<uint64_t, uint64_t>> merged;
        for (auto &iv : occupied) {
            if (merged.empty() || merged.back().second < iv.first) {
                merged.push_back(iv);
            } else {
                merged.back().second =
                    std::max(merged.back().second, iv.second);
            }
        }

        // Best-fit: find smallest gap >= t.size
        uint64_t bestOffset = 0;
        uint64_t bestWaste = UINT64_MAX;

        if (merged.empty()) {
            bestOffset = 0;
        } else {
            // Gap before first interval [0, merged[0].first)
            if (merged[0].first >= t.size &&
                merged[0].first - t.size < bestWaste) {
                bestWaste = merged[0].first - t.size;
                bestOffset = 0;
            }

            // Gaps between intervals
            for (size_t i = 0; i + 1 < merged.size(); ++i) {
                uint64_t gapSize = merged[i + 1].first - merged[i].second;
                if (gapSize >= t.size && gapSize - t.size < bestWaste) {
                    bestWaste = gapSize - t.size;
                    bestOffset = merged[i].second;
                }
            }

            // If no gap fits, place after last interval
            if (bestWaste == UINT64_MAX) {
                bestOffset = merged.back().second;
            }
        }

        t.offset = bestOffset;
    }

    uint64_t peak = 0;
    for (const Tensor &t : tensors) {
        peak = std::max(peak, t.offset + t.size);
    }

    assert(verify(tensors) && "Memory packing verification failed");

    return peak;
}

} // namespace packing

namespace ilp {

ILPGraph::ILPGraph(const core::Graph &graph) {
    size = graph.size;
    tensorSizes.resize(graph.size);
    extraSizes.resize(graph.size);
    inputTensors.resize(graph.size);

    for (uint64_t node = 0; node < graph.size; ++node) {
        tensorSizes[node] = 0;
        extraSizes[node] = graph.nodeWeights[node];

        for (const auto &[neighbor, edgeWeight] : graph.adj[node]) {
            inputTensors[neighbor].insert(node);
            tensorSizes[node] = std::max(tensorSizes[node], edgeWeight);
        }
    }

    computeTopology();
}

void ILPGraph::computeTopology() {
    computeAncestors();
    computeDescendants();
}

void ILPGraph::computeAncestors() {
    ancestors.clear();
    ancestors.resize(size);
    tensorUsers.clear();
    tensorUsers.resize(size);

    // Compute in-degree for topological sort
    std::vector<uint64_t> inDegree(size, 0);
    for (uint64_t i = 0; i < size; i++) {
        inDegree[i] = inputTensors[i].size();
    }

    // Topological sort using Kahn's algorithm
    std::queue<uint64_t> queue;
    for (uint64_t i = 0; i < size; i++) {
        if (inDegree[i] == 0) {
            queue.push(i);
        }
    }

    while (!queue.empty()) {
        uint64_t node = queue.front();
        queue.pop();

        // For this node, its ancestors are the union of:
        // 1. Its direct inputs
        // 2. All ancestors of its inputs
        for (uint64_t input : inputTensors[node]) {
            ancestors[node].insert(input);
            ancestors[node].insert(ancestors[input].begin(),
                                   ancestors[input].end());
            tensorUsers[input].insert(node);
        }

        // Update in-degrees and enqueue nodes
        for (uint64_t i = 0; i < size; i++) {
            if (inputTensors[i].count(node)) {
                inDegree[i]--;
                if (inDegree[i] == 0) {
                    queue.push(i);
                }
            }
        }
    }
}

void ILPGraph::computeDescendants() {
    descendants.clear();
    descendants.resize(size);

    // Build descendants from ancestors: if j is ancestor of i, then i is
    // descendant of j
    for (uint64_t i = 0; i < size; i++) {
        for (uint64_t ancestor : ancestors[i]) {
            descendants[ancestor].insert(i);
        }
    }
}

void ILPGraph::print(std::ostream &os) const {
    os << "Nodes: " << size << std::endl;

    for (uint64_t i = 0; i < size; i++) {
        os << "\nNode: " << i << std::endl;
        os << "    tensorSizes=" << tensorSizes[i]
           << ", extraSizes=" << extraSizes[i] << std::endl;

        os << "    inputs={";
        for (auto it = inputTensors[i].begin(); it != inputTensors[i].end();
             ++it) {
            if (it != inputTensors[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}, users={";
        for (auto it = tensorUsers[i].begin(); it != tensorUsers[i].end();
             ++it) {
            if (it != tensorUsers[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}" << std::endl;

        os << "    anc={";
        for (auto it = ancestors[i].begin(); it != ancestors[i].end(); ++it) {
            if (it != ancestors[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}, des={";
        for (auto it = descendants[i].begin(); it != descendants[i].end();
             ++it) {
            if (it != descendants[i].begin())
                os << ", ";
            os << *it;
        }
        os << "}" << std::endl;
    }
}

uint64_t ILPGraph::getSize() const { return size; }

const std::vector<robin_hood::unordered_set<uint64_t>> &
ILPGraph::getAncestors() const {
    return ancestors;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
ILPGraph::getDescendants() const {
    return descendants;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
ILPGraph::getTensorUsers() const {
    return tensorUsers;
}

const std::vector<uint64_t> &ILPGraph::getTensorSizes() const {
    return tensorSizes;
}

const std::vector<uint64_t> &ILPGraph::getExtraSizes() const {
    return extraSizes;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
ILPGraph::getInputTensors() const {
    return inputTensors;
}

// Stream output operator
std::ostream &operator<<(std::ostream &os, const ILPGraph &graph) {
    graph.print(os);
    return os;
}

ILPSolver::ILPSolver(const ILPGraph &graph, bool debug = false)
    : graph(graph), debug(debug) {
    solver.reset(operations_research::MPSolver::CreateSolver("SCIP"));
    if (solver)
        return;

    solver.reset(operations_research::MPSolver::CreateSolver("CBC"));
    if (solver)
        return;

    throw std::runtime_error(
        "Failed to create ILP solver: neither SCIP nor CBC available");
}

void ILPSolver::computePruningBound(std::vector<uint64_t> &earliestOp,
                                    std::vector<uint64_t> &latestOp,
                                    std::vector<uint64_t> &earliestTensor,
                                    std::vector<uint64_t> &latestTensor) const {
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        earliestOp[i] = graph.getAncestors()[i].size();
        latestOp[i] = graph.getSize() - 1 - graph.getDescendants()[i].size();
        earliestTensor[i] = graph.getAncestors()[i].size();

        // Constraint 10: tensor not needed after last consumer could run
        const auto &users = graph.getTensorUsers()[i];
        if (users.empty()) {
            latestTensor[i] = graph.getSize() - 1; // Output tensor
        } else {
            // Find user with minimum descendants (can run latest)
            uint64_t min_des = graph.getSize();
            for (uint64_t u : users) {
                min_des = std::min(min_des, graph.getDescendants()[u].size());
            }
            latestTensor[i] = graph.getSize() - 1 - min_des;
        }
    }

    if (debug) {
        std::cerr << "\nPruning bounds:" << std::endl;
        for (uint64_t i = 0; i < graph.getSize(); i++) {
            std::cerr << "  Node " << i << ":"
                      << " earliestOp=" << earliestOp[i]
                      << ", latestOp=" << latestOp[i]
                      << ", earliestTensor=" << earliestTensor[i]
                      << ", latestTensor=" << latestTensor[i] << std::endl;
        }
    }
}

/**
 * Variables:
 * - O[i][j] ∈ {0,1}: operator i scheduled at step j
 * - T[i][j] ∈ {0,1}: tensor i in memory at step j
 * - mem: peak memory
 *
 * Constraints:
 * (1) Σᵢ O[i][j] = 1  ∀j        -- one operator per step
 * (2) Σⱼ O[i][j] = 1  ∀i        -- each operator scheduled once
 * (3) O[i][j] ≤ T[k][j]  ∀k∈IN[i]  -- inputs available
 * (4a) T[i][j] ≤ T[i][j-1] + O[i][j]  -- tensor persistence upper bound
 * (4b) T[i][j] ≥ O[i][j]              -- tensor must exist when created
 * (4c) T[i][j] ≤ Σ_{c∈users[i]} Σ_{k≥j} O[c][k]  -- tensor freed after last
 * consumer
 * (5) T[i][0] ≤ O[i][0]          -- initial condition
 * (6) Σᵢ T[i][j]·S[i] + ES[k]·O[k][j] ≤ mem  -- memory bound
 *
 * Pruning:
 * (7) O[i][j] = 0 for j < |anc(i)|
 * (8) T[i][j] = 0 for j < |anc(i)|
 * (9) O[i][j] = 0 for j > n-1-|des(i)|
 */

void ILPSolver::createVariables(const std::vector<uint64_t> &earliestOp,
                                const std::vector<uint64_t> &latestOp,
                                const std::vector<uint64_t> &earliestTensor,
                                const std::vector<uint64_t> &latestTensor) {
    O.assign(graph.getSize(), std::vector<operations_research::MPVariable *>(
                                  graph.getSize(), nullptr));
    T.assign(graph.getSize(), std::vector<operations_research::MPVariable *>(
                                  graph.getSize(), nullptr));

    uint64_t num_O = 0, num_T = 0;
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        for (uint64_t j = 0; j < graph.getSize(); j++) {
            if (earliestOp[i] <= j && j <= latestOp[i]) {
                O[i][j] = solver->MakeBoolVar("O_" + std::to_string(i) + "_" +
                                              std::to_string(j));
                num_O++;
            }
            if (earliestTensor[i] <= j && j <= latestTensor[i]) {
                T[i][j] = solver->MakeBoolVar("T_" + std::to_string(i) + "_" +
                                              std::to_string(j));
                num_T++;
            }
        }
    }

    mem = solver->MakeNumVar(0, solver->infinity(), "mem");

    if (debug) {
        std::cerr << "\nVariables: O=" << num_O << ", T=" << num_T
                  << ", total=" << (num_O + num_T) << ", out of maximum "
                  << graph.getSize() * graph.getSize() * 2 << std::endl;
    }

    // Constraint (1): exactly one operator per step
    for (uint64_t j = 0; j < graph.getSize(); j++) {
        operations_research::MPConstraint *ct = solver->MakeRowConstraint(
            1, 1, "one_per_step_" + std::to_string(j));
        for (uint64_t i = 0; i < graph.getSize(); i++) {
            if (O[i][j])
                ct->SetCoefficient(O[i][j], 1);
        }
    }

    // Constraint (2): each operator scheduled exactly once
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        operations_research::MPConstraint *ct = solver->MakeRowConstraint(
            1, 1, "scheduled_once_" + std::to_string(i));
        for (uint64_t j = 0; j < graph.getSize(); j++) {
            if (O[i][j])
                ct->SetCoefficient(O[i][j], 1);
        }
    }

    // Constraint (3): inputs must be available
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        for (uint64_t j = 0; j < graph.getSize(); j++) {
            if (!O[i][j])
                continue;
            for (uint64_t k : graph.getInputTensors()[i]) {
                if (T[k][j]) {
                    operations_research::MPConstraint *ct =
                        solver->MakeRowConstraint(-solver->infinity(), 0,
                                                  "input_" + std::to_string(i) +
                                                      "_" + std::to_string(j) +
                                                      "_" + std::to_string(k));
                    ct->SetCoefficient(O[i][j], 1);
                    ct->SetCoefficient(T[k][j], -1);
                } else {
                    operations_research::MPConstraint *ct =
                        solver->MakeRowConstraint(
                            0, 0,
                            "input_pruned_" + std::to_string(i) + "_" +
                                std::to_string(j) + "_" + std::to_string(k));
                    ct->SetCoefficient(O[i][j], 1);
                }
            }
        }
    }

    // Constraint (4a): tensor persistence upper bound
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        for (uint64_t j = 1; j < graph.getSize(); j++) {
            if (!T[i][j])
                continue;
            operations_research::MPConstraint *ct = solver->MakeRowConstraint(
                -solver->infinity(), 0,
                "persist_upper_" + std::to_string(i) + "_" + std::to_string(j));
            ct->SetCoefficient(T[i][j], 1);
            if (T[i][j - 1])
                ct->SetCoefficient(T[i][j - 1], -1);
            if (O[i][j])
                ct->SetCoefficient(O[i][j], -1);
        }
    }

    // Constraint (4b): tensor must exist when created
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        for (uint64_t j = 0; j < graph.getSize(); j++) {
            if (O[i][j] && T[i][j]) {
                operations_research::MPConstraint *ct =
                    solver->MakeRowConstraint(-solver->infinity(), 0,
                                              "must_exist_" +
                                                  std::to_string(i) + "_" +
                                                  std::to_string(j));
                ct->SetCoefficient(O[i][j], 1);
                ct->SetCoefficient(T[i][j], -1);
            }
        }
    }

    // Constraint (4c): tensor can only exist if at least one consumer hasn't
    // run yet
    // T[i][j] ≤ Σ_{c ∈ users} Σ_{k≥j} O[c][k]
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        const auto &users = graph.getTensorUsers()[i];
        if (users.empty())
            continue; // Output tensor, no constraint needed

        for (uint64_t j = 0; j < graph.getSize(); j++) {
            if (!T[i][j])
                continue;

            // T[i][j] - Σ_{c} Σ_{k≥j} O[c][k] ≤ 0
            operations_research::MPConstraint *ct = solver->MakeRowConstraint(
                -solver->infinity(), 0,
                "lifetime_" + std::to_string(i) + "_" + std::to_string(j));

            ct->SetCoefficient(T[i][j], 1);

            for (uint64_t c : users) {
                for (uint64_t k = j; k < graph.getSize(); k++) {
                    if (O[c][k]) {
                        ct->SetCoefficient(O[c][k], -1);
                    }
                }
            }
        }
    }

    // Constraint (5): initial condition
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        if (T[i][0]) {
            if (O[i][0]) {
                operations_research::MPConstraint *ct =
                    solver->MakeRowConstraint(-solver->infinity(), 0,
                                              "init_" + std::to_string(i));
                ct->SetCoefficient(T[i][0], 1);
                ct->SetCoefficient(O[i][0], -1);
            } else {
                operations_research::MPConstraint *ct =
                    solver->MakeRowConstraint(0, 0,
                                              "init_zero_" + std::to_string(i));
                ct->SetCoefficient(T[i][0], 1);
            }
        }
    }

    // Constraint (6): memory bound
    for (uint64_t j = 0; j < graph.getSize(); j++) {
        for (uint64_t k = 0; k < graph.getSize(); k++) {
            if (!O[k][j])
                continue;

            // ∑ T[i][j]·S[i] + ES[k]·O[k][j] ≤ mem
            // Rearranged: ∑ T[i][j]·S[i] + ES[k]·O[k][j] - mem ≤ 0
            operations_research::MPConstraint *ct = solver->MakeRowConstraint(
                -solver->infinity(), 0,
                "mem_" + std::to_string(j) + "_" + std::to_string(k));

            // Sum of live tensor sizes
            for (uint64_t i = 0; i < graph.getSize(); i++) {
                if (T[i][j]) {
                    ct->SetCoefficient(T[i][j], graph.getTensorSizes()[i]);
                }
            }

            // Extra workspace for operator k
            ct->SetCoefficient(O[k][j], graph.getExtraSizes()[k]);

            // -mem (moved to LHS)
            ct->SetCoefficient(mem, -1);
        }
    }
}

std::tuple<operations_research::MPSolver::ResultStatus, std::vector<uint64_t>,
           uint64_t>
ILPSolver::solve(uint64_t timeLimitSeconds) {
    std::vector<uint64_t> earliestOp(graph.getSize()),
        latestOp(graph.getSize());
    std::vector<uint64_t> earliestTensor(graph.getSize()),
        latestTensor(graph.getSize());

    computePruningBound(earliestOp, latestOp, earliestTensor, latestTensor);
    createVariables(earliestOp, latestOp, earliestTensor, latestTensor);

    // Objective: minimize peak memory
    operations_research::MPObjective *objective = solver->MutableObjective();
    objective->SetCoefficient(mem, 1);
    objective->SetMinimization();

    if (debug) {
        std::cerr << "Constraints: " << solver->NumConstraints() << std::endl;
        std::cerr << "\nSolving..." << std::endl;
    }

    solver->SetTimeLimit(absl::Seconds(timeLimitSeconds));
    operations_research::MPSolver::ResultStatus status = solver->Solve();

    if (status == operations_research::MPSolver::OPTIMAL ||
        status == operations_research::MPSolver::FEASIBLE) {
        uint64_t peakMemory = static_cast<uint64_t>(mem->solution_value());

        if (debug) {
            std::cerr << (status == operations_research::MPSolver::OPTIMAL
                              ? "OPTIMAL"
                              : "FEASIBLE")
                      << " solution found!" << std::endl;
            std::cerr << "Peak memory: " << peakMemory << std::endl;
        }

        std::vector<uint64_t> schedule(graph.getSize(), -1);
        for (uint64_t i = 0; i < graph.getSize(); i++) {
            for (uint64_t j = 0; j < graph.getSize(); j++) {
                if (O[i][j] && O[i][j]->solution_value() > 0.5) {
                    schedule[i] = j;
                    break;
                }
            }
        }

        if (debug) {
            std::cerr << "Schedule: ";
            for (uint64_t i = 0; i < graph.getSize(); i++) {
                std::cerr << schedule[i] << " ";
            }
            std::cerr << std::endl;
        }

        return {status, schedule, peakMemory};
    }

    if (debug) {
        std::cerr << "No solution found. Status: " << status << std::endl;
    }

    return {status, {}, 0};
}

void ILPSolver::printVariables(std::ostream &os) const {
    for (uint64_t i = 0; i < graph.getSize(); i++) {
        for (uint64_t j = 0; j < graph.getSize(); j++) {
            if (O[i][j]) {
                os << "O[" << i << "][" << j
                   << "] = " << O[i][j]->solution_value() << std::endl;
            } else {
                os << "O[" << i << "][" << j << "] = NULL" << std::endl;
            }
            if (T[i][j]) {
                os << "T[" << i << "][" << j
                   << "] = " << T[i][j]->solution_value() << std::endl;
            } else {
                os << "T[" << i << "][" << j << "] = NULL" << std::endl;
            }
            os << std::endl;
        }
    }
}

} // namespace ilp

namespace bruteforce {

BruteForceSolver::BruteForceSolver(const ilp::ILPGraph &graph, bool debug)
    : graph(graph), debug(debug) {}

uint64_t
BruteForceSolver::computePeakMemory(const std::vector<uint64_t> &order) const {
    uint64_t n = order.size();
    robin_hood::unordered_set<uint64_t> live;
    uint64_t peak = 0;

    const auto &inputTensors = graph.getInputTensors();
    const auto &tensorSizes = graph.getTensorSizes();
    const auto &extraSizes = graph.getExtraSizes();

    for (uint64_t step = 0; step < n; step++) {
        uint64_t op = order[step];

        // Memory during execution: live tensors + new output + workspace
        uint64_t mem = extraSizes[op] + tensorSizes[op];
        for (uint64_t t : live) {
            mem += tensorSizes[t];
        }
        peak = std::max(peak, mem);

        // Add output tensor to live set
        live.insert(op);

        // Remove tensors whose last consumer just ran
        // (A tensor can be freed after all its consumers have executed)
        std::vector<uint64_t> toRemove;
        for (uint64_t t : live) {
            bool canFree = true;
            for (uint64_t i = 0; i < n; i++) {
                if (inputTensors[i].count(t)) {
                    // Operator i needs tensor t - check if i has run
                    bool hasRun = false;
                    for (uint64_t s = 0; s <= step; s++) {
                        if (order[s] == i) {
                            hasRun = true;
                            break;
                        }
                    }
                    if (!hasRun) {
                        canFree = false;
                        break;
                    }
                }
            }
            if (canFree)
                toRemove.push_back(t);
        }
        for (uint64_t t : toRemove)
            live.erase(t);
    }

    return peak;
}

bool BruteForceSolver::canPlace(uint64_t node,
                                const std::vector<bool> &placed) const {
    const auto &inputTensors = graph.getInputTensors();
    for (uint64_t dep : inputTensors[node]) {
        if (!placed[dep])
            return false;
    }
    return true;
}

void BruteForceSolver::enumerate(std::vector<uint64_t> &current,
                                 std::vector<bool> &placed, uint64_t &bestPeak,
                                 std::vector<uint64_t> &bestOrder,
                                 uint64_t &count,
                                 std::vector<std::string> &outputs) const {
    uint64_t n = graph.getSize();
    if (current.size() == n) {
        count++;
        uint64_t peak = computePeakMemory(current);

        if (debug) {
            std::ostringstream oss;
            oss << "Order " << count - 1 << " (peak " << peak << "): ";
            for (const uint64_t &node : current) {
                oss << node << " ";
            }
            outputs.emplace_back(oss.str());
        }

        if (peak < bestPeak) {
            bestPeak = peak;
            bestOrder = current;
        }
        return;
    }

    for (uint64_t node = 0; node < n; node++) {
        if (!placed[node] && canPlace(node, placed)) {
            placed[node] = true;
            current.push_back(node);

            enumerate(current, placed, bestPeak, bestOrder, count, outputs);

            current.pop_back();
            placed[node] = false;
        }
    }
}

std::pair<std::vector<uint64_t>, uint64_t> BruteForceSolver::solve() const {
    uint64_t n = graph.getSize();
    std::vector<std::string> outputs;

    std::vector<uint64_t> current, bestOrder;
    std::vector<bool> placed(n, false);
    uint64_t bestPeak = UINT64_MAX;
    uint64_t count = 0;

    enumerate(current, placed, bestPeak, bestOrder, count, outputs);

    if (debug) {
        std::cerr << "\nBrute-force results:\n";
        for (const auto &out : outputs) {
            std::cerr << out << std::endl;
        }
        std::cerr << "Best result: peak memory " << bestPeak << " in order: ";
        for (const uint64_t &node : bestOrder) {
            std::cerr << node << " ";
        }
        std::cerr << std::endl;
    }

    return {bestOrder, bestPeak};
}

} // namespace bruteforce

Scheduler::Scheduler(const core::Graph &originalGraph,
                     std::vector<uint64_t> &partitionMapping, bool debug,
                     bool verify)
    : originalGraph(originalGraph), partitionMapping(partitionMapping),
      debug(debug), verify(verify) {}

std::vector<uint64_t> Scheduler::calculatePartitionWeights() const {
    std::vector<uint64_t> topologicalOrder = originalGraph.topologicalSort();

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

    // Collect tensors per partition
    std::vector<std::vector<packing::Tensor>> partitionTensors(numPartitions);

    for (uint64_t producer = 0; producer < originalGraph.size; ++producer) {
        uint64_t producerPartition = partitionMapping[producer];

        // Find all consumers of this producer within the same partition
        uint64_t lastConsumerPosition =
            nodeToLocalPosition[producer]; // Birth time
        bool hasInternalConsumers = false;

        for (const auto &[consumer, edgeWeight] : originalGraph.adj[producer]) {
            // Only consider consumers in the same partition
            if (partitionMapping[consumer] == producerPartition) {
                hasInternalConsumers = true;
                lastConsumerPosition = std::max(lastConsumerPosition,
                                                nodeToLocalPosition[consumer]);
            }
        }

        // If there are internal consumers, add tensor to partition
        if (hasInternalConsumers) {
            packing::Tensor tensor;
            tensor.producer = producer;
            tensor.size = originalGraph.nodeWeights[producer];
            tensor.birth = nodeToLocalPosition[producer];
            tensor.death = lastConsumerPosition;
            tensor.offset = 0; // Will be assigned during packing
            partitionTensors[producerPartition].push_back(tensor);
        }
    }

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
                  [](const packing::Tensor &a, const packing::Tensor &b) {
                      if (a.birth != b.birth)
                          return a.birth < b.birth;
                      return a.death < b.death;
                  });

        // Best-fit packing algorithm
        uint64_t peakMemory = packing::packBestFit(tensors);

        partitionWeights[partition] = peakMemory;
    }

    return partitionWeights;
}

void Scheduler::buildCoarseGraph() {
    // Find the number of partitions (assuming contiguous IDs from 0 to k-1)
    uint64_t numPartitions =
        *std::max_element(partitionMapping.begin(), partitionMapping.end()) + 1;
    assert(numPartitions != 0 && "Number of partitions must be > 0");

    coarseGraph = std::make_unique<core::Graph>(numPartitions);

    std::vector<uint64_t> partitionWeights = calculatePartitionWeights();

    // Add nodes to coarse graph - partition ID is the node ID
    for (uint64_t partitionId = 0; partitionId < numPartitions; ++partitionId) {
        coarseGraph->addNode(partitionId, partitionWeights[partitionId]);
    }

    // For each edge in the original graph, create an edge between partitions
    // if needed, with total size of edges between partitions
    std::vector<std::unordered_map<uint64_t, uint64_t>> newEdges(numPartitions);
    for (uint64_t from = 0; from < originalGraph.size; ++from) {
        for (const auto &[to, edgeWeight] : originalGraph.adj[from]) {
            uint64_t fromPartition = partitionMapping[from];
            uint64_t toPartition = partitionMapping[to];

            if (fromPartition != toPartition) {
                if (newEdges[fromPartition].find(toPartition) !=
                    newEdges[fromPartition].end()) {
                    newEdges[fromPartition][toPartition] += edgeWeight;
                } else {
                    newEdges[fromPartition].emplace(toPartition, edgeWeight);
                }
            }
        }
    }

    for (uint64_t partitionId = 0; partitionId < newEdges.size();
         ++partitionId) {
        for (const auto &[neighborId, edgeWeight] : newEdges[partitionId]) {
            coarseGraph->addEdge(partitionId, neighborId, edgeWeight);
        }
    }
}

std::pair<std::vector<uint64_t>, uint64_t> Scheduler::run() {
    buildCoarseGraph();
    ilp::ILPGraph ilpGraph(*coarseGraph);

    ilp::ILPSolver solver(ilpGraph, debug);
    auto [status, schedule, peakMemory] = solver.solve();

    if (verify) {
        bruteforce::BruteForceSolver bfSolver(ilpGraph, debug);
        auto [bfSchedule, bfPeakMemory] = bfSolver.solve();

        if (peakMemory != bfPeakMemory) {
            std::cerr << "\n=== VERIFICATION FAILED: DUMP START ===\n"
                      << std::endl;

            std::cerr << "=== Coarse Graph ===\n";
            std::cerr << *coarseGraph;
            std::cerr << std::endl;

            std::cerr << "=== ILPGraph ===\n";
            std::cerr << ilpGraph;
            std::cerr << std::endl;

            std::cerr << "=== ILP Solution ===\n";
            std::cerr << "Peak Memory: " << peakMemory << std::endl;
            std::cerr << "Schedule: ";
            for (size_t i = 0; i < schedule.size(); ++i) {
                std::cerr << schedule[i] << " ";
            }
            std::cerr << "\n" << std::endl;

            std::cerr << "=== Brute-Force Solution ===\n";
            std::cerr << "Peak Memory: " << bfPeakMemory << std::endl;
            std::cerr << "Schedule: ";
            for (size_t i = 0; i < bfSchedule.size(); ++i) {
                std::cerr << bfSchedule[i] << " ";
            }
            std::cerr << "\n" << std::endl;

            solver.printVariables(std::cerr);

            std::cerr << "\nStatus of the ILP solver: ";
            switch (status) {
            case operations_research::MPSolver::ResultStatus::OPTIMAL:
                std::cerr << "OPTIMAL" << std::endl;
                break;
            case operations_research::MPSolver::ResultStatus::FEASIBLE:
                std::cerr << "FEASIBLE" << std::endl;
                break;
            case operations_research::MPSolver::ResultStatus::INFEASIBLE:
                std::cerr << "INFEASIBLE" << std::endl;
                break;
            case operations_research::MPSolver::ResultStatus::UNBOUNDED:
                std::cerr << "UNBOUNDED" << std::endl;
                break;
            case operations_research::MPSolver::ResultStatus::ABNORMAL:
                std::cerr << "ABNORMAL" << std::endl;
                break;
            case operations_research::MPSolver::ResultStatus::NOT_SOLVED:
                std::cerr << "NOT_SOLVED" << std::endl;
                break;
            case operations_research::MPSolver::ResultStatus::MODEL_INVALID:
                std::cerr << "MODEL_INVALID" << std::endl;
                break;
            }

            std::cerr << "\n=== VERIFICATION FAILED: DUMP END ===\n"
                      << std::endl;

            throw std::runtime_error("Verification failed: ILP peak memory (" +
                                     std::to_string(peakMemory) +
                                     ") != brute-force peak memory (" +
                                     std::to_string(bfPeakMemory) + ")");
        }
    }

    return {schedule, peakMemory};
}

} // namespace scheduling

} // namespace dag_partitioning