/**
 * @file Scheduling.cpp
 * @brief Implementation of memory-aware scheduling algorithms
 */

#include "Scheduling.h"

#include <cassert>
#include <sstream>
#include <stack>

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

HyperGraph::HyperGraph(const core::Graph &graph) {
    size = graph.size;
    tensorSizes.resize(size);
    extraSizes.resize(size);
    inputTensors.resize(size);

    for (uint64_t node = 0; node < size; ++node) {
        tensorSizes[node] = 0;
        extraSizes[node] = graph.nodeWeights[node];

        for (const auto &[neighbor, edgeWeight] : graph.adj[node]) {
            inputTensors[neighbor].insert(node);
            tensorSizes[node] = std::max(tensorSizes[node], edgeWeight);
        }
    }

    computeTopology();
}

void HyperGraph::computeTopology() {
    computeAncestors();
    computeDescendants();
}

void HyperGraph::computeAncestors() {
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

void HyperGraph::computeDescendants() {
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

void HyperGraph::print(std::ostream &os) const {
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

uint64_t HyperGraph::getSize() const { return size; }

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getAncestors() const {
    return ancestors;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getDescendants() const {
    return descendants;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getTensorUsers() const {
    return tensorUsers;
}

const std::vector<uint64_t> &HyperGraph::getTensorSizes() const {
    return tensorSizes;
}

const std::vector<uint64_t> &HyperGraph::getExtraSizes() const {
    return extraSizes;
}

const std::vector<robin_hood::unordered_set<uint64_t>> &
HyperGraph::getInputTensors() const {
    return inputTensors;
}

std::ostream &operator<<(std::ostream &os, const HyperGraph &graph) {
    graph.print(os);
    return os;
}

namespace rpo {

RPOScheduler::RPOScheduler(const HyperGraph &graph, bool debug)
    : graph(graph), debug(debug) {}

std::pair<std::vector<uint64_t>, uint64_t> RPOScheduler::solve() const {
    uint64_t n = graph.getSize();
    std::vector<uint64_t> schedule;
    schedule.reserve(n);

    // Build adjacency list for DFS (node -> children/successors)
    std::vector<std::vector<uint64_t>> children(n);
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t input : graph.getInputTensors()[i]) {
            children[input].push_back(i);
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

    // Compute opToStep mapping
    std::vector<uint64_t> opToStep(n);
    for (uint64_t step = 0; step < n; step++) {
        opToStep[schedule[step]] = step;
    }

    // Compute peak memory
    uint64_t peakMemory = 0;
    for (uint64_t step = 0; step < n; step++) {
        uint64_t memAtStep = 0;

        for (uint64_t tensor = 0; tensor < n; tensor++) {
            uint64_t birthStep = opToStep[tensor];
            const auto &users = graph.getTensorUsers()[tensor];

            uint64_t deathStep = birthStep;
            for (uint64_t user : users) {
                deathStep = std::max(deathStep, opToStep[user]);
            }

            if (birthStep <= step && step <= deathStep) {
                memAtStep += graph.getTensorSizes()[tensor];
            }
        }

        uint64_t op = schedule[step];
        memAtStep += graph.getExtraSizes()[op];

        peakMemory = std::max(peakMemory, memAtStep);
    }

    if (debug) {
        std::cerr << "RPO peak memory: " << peakMemory << std::endl;
    }

    return {schedule, peakMemory};
}

} // namespace rpo

namespace cpsat {

CPSATSolver::CPSATSolver(const HyperGraph &graph, bool debug)
    : graph(graph), debug(debug) {
    uint64_t n = graph.getSize();
    earliestOp.resize(n);
    latestOp.resize(n);
    earliestTensor.resize(n);
    latestTensor.resize(n);

    cpModel = operations_research::sat::CpModelBuilder();
}

/**
 * Pruning:
 * (7) O[i][j] exists for j >= |anc(i)|
 * (8) T[i][j] exists for j >= |anc(i)|
 * (9) O[i][j] exists for j <= n-1-|des(i)|
 */
void CPSATSolver::computePruningBounds() {
    uint64_t n = graph.getSize();
    for (uint64_t i = 0; i < n; i++) {
        earliestOp[i] = graph.getAncestors()[i].size();
        latestOp[i] = n - 1 - graph.getDescendants()[i].size();
        earliestTensor[i] = graph.getAncestors()[i].size();

        const auto &users = graph.getTensorUsers()[i];
        if (users.empty()) {
            latestTensor[i] = n - 1;
        } else {
            uint64_t minDes = n;
            for (uint64_t u : users) {
                minDes = std::min(minDes, graph.getDescendants()[u].size());
            }
            latestTensor[i] = n - 1 - minDes;
        }
    }

    if (debug) {
        std::cerr << "\nPruning bounds:" << std::endl;
        for (uint64_t i = 0; i < n; i++) {
            std::cerr << "  Node " << i << ":"
                      << " earliestOp=" << earliestOp[i]
                      << ", latestOp=" << latestOp[i]
                      << ", earliestTensor=" << earliestTensor[i]
                      << ", latestTensor=" << latestTensor[i] << std::endl;
        }
    }
}

bool CPSATSolver::hasO(uint64_t i, uint64_t j) const {
    return earliestOp[i] <= j && j <= latestOp[i];
}

bool CPSATSolver::hasT(uint64_t i, uint64_t j) const {
    return earliestTensor[i] <= j && j <= latestTensor[i];
}

uint64_t CPSATSolver::oKey(uint64_t i, uint64_t j) const {
    return i * graph.getSize() + j;
}

uint64_t CPSATSolver::tKey(uint64_t i, uint64_t j) const {
    return i * graph.getSize() + j;
}

operations_research::sat::BoolVar CPSATSolver::getO(uint64_t i,
                                                    uint64_t j) const {
    return O.at(oKey(i, j));
}

operations_research::sat::BoolVar CPSATSolver::getT(uint64_t i,
                                                    uint64_t j) const {
    return T.at(tKey(i, j));
}

/**
 * CP-SAT Model Formulation:
 *
 * Variables:
 * - O[i][j] ∈ {0,1}: operator i scheduled at step j
 * - T[i][j] ∈ {0,1}: tensor i in memory at step j
 * - mem ∈ [0, memoryUpperBound]: peak memory
 *
 * Constraints:
 * (1) ExactlyOne(O[i][j] for all i) for each j -- one operator per step
 * (2) ExactlyOne(O[i][j] for all j) for each i -- each operator scheduled once
 * (3) O[i][j] => T[k][j] for all k∈IN[i]   -- inputs available (implication)
 * (4a) T[i][j] => T[i][j-1] OR O[i][j]     -- tensor persistence upper bound
 * (4b) O[i][j] => T[i][j]                  -- tensor must exist when created
 * (4c) T[i][j] => OR_{c∈users} OR_{k>=j} O[c][k] -- tensor freed after last
 * consumer
 * (5) T[i][0] => O[i][0]          -- initial condition
 * (6) Σᵢ T[i][j]·S[i] + ES[k]·O[k][j] ≤ mem  -- memory bound
 *
 * Objective: minimize mem
 */
void CPSATSolver::createVariablesAndConstraints(uint64_t memoryUpperBound) {
    uint64_t n = graph.getSize();

    // Create O and T variables
    uint64_t numO = 0, numT = 0;
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
            if (hasO(i, j)) {
                O[oKey(i, j)] = cpModel.NewBoolVar();
                numO++;
            }
            if (hasT(i, j)) {
                T[tKey(i, j)] = cpModel.NewBoolVar();
                numT++;
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

    // Constraint (4a): T[i][j] => T[i][j-1] OR O[i][j]
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 1; j < n; j++) {
            if (!hasT(i, j))
                continue;

            std::vector<operations_research::sat::BoolVar> clause;
            clause.push_back(getT(i, j).Not());

            if (hasT(i, j - 1)) {
                clause.push_back(getT(i, j - 1));
            }
            if (hasO(i, j)) {
                clause.push_back(getO(i, j));
            }

            cpModel.AddBoolOr(clause);
        }
    }

    // Constraint (4b): O[i][j] => T[i][j]
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < n; j++) {
            if (hasO(i, j) && hasT(i, j)) {
                cpModel.AddImplication(getO(i, j), getT(i, j));
            }
        }
    }

    // Constraint (4c): T[i][j] => at least one consumer hasn't run yet
    for (uint64_t i = 0; i < n; i++) {
        const auto &users = graph.getTensorUsers()[i];
        if (users.empty())
            continue;

        for (uint64_t j = 0; j < n; j++) {
            if (!hasT(i, j))
                continue;

            std::vector<operations_research::sat::BoolVar> futureConsumers;
            for (uint64_t c : users) {
                for (uint64_t k = j; k < n; k++) {
                    if (hasO(c, k)) {
                        futureConsumers.push_back(getO(c, k));
                    }
                }
            }

            if (futureConsumers.empty()) {
                cpModel.AddEquality(getT(i, j), 0);
            } else {
                std::vector<operations_research::sat::BoolVar> clause;
                clause.push_back(getT(i, j).Not());
                clause.insert(clause.end(), futureConsumers.begin(),
                              futureConsumers.end());
                cpModel.AddBoolOr(clause);
            }
        }
    }

    // Constraint (5): T[i][0] => O[i][0]
    for (uint64_t i = 0; i < n; i++) {
        if (hasT(i, 0)) {
            if (hasO(i, 0)) {
                cpModel.AddImplication(getT(i, 0), getO(i, 0));
            } else {
                cpModel.AddEquality(getT(i, 0), 0);
            }
        }
    }

    // Constraint (6): memory bound at each step
    for (uint64_t j = 0; j < n; j++) {
        for (uint64_t k = 0; k < n; k++) {
            if (!hasO(k, j))
                continue;

            operations_research::sat::LinearExpr memExpr;
            for (uint64_t i = 0; i < n; i++) {
                if (hasT(i, j)) {
                    memExpr += getT(i, j) *
                               static_cast<int64_t>(graph.getTensorSizes()[i]);
                }
            }
            memExpr +=
                getO(k, j) * static_cast<int64_t>(graph.getExtraSizes()[k]);

            cpModel.AddLessOrEqual(memExpr, mem);
        }
    }

    cpModel.Minimize(mem);
}

void CPSATSolver::setWarmStartHints(const std::vector<uint64_t> &schedule,
                                    uint64_t peakMemory) {
    uint64_t n = graph.getSize();

    // Compute opToStep mapping
    std::vector<uint64_t> opToStep(n);
    for (uint64_t step = 0; step < n; step++) {
        opToStep[schedule[step]] = step;
    }

    // Compute tensor liveness
    std::vector<robin_hood::unordered_set<uint64_t>> live(n);
    for (uint64_t tensor = 0; tensor < n; tensor++) {
        uint64_t birthStep = opToStep[tensor];
        const auto &users = graph.getTensorUsers()[tensor];

        uint64_t deathStep = birthStep;
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
    for (uint64_t tensor = 0; tensor < n; tensor++) {
        for (uint64_t step = 0; step < n; step++) {
            if (hasT(tensor, step)) {
                cpModel.AddHint(getT(tensor, step),
                                live[step].count(tensor) ? 1 : 0);
            }
        }
    }

    // Add hint for mem variable
    cpModel.AddHint(mem, static_cast<int64_t>(peakMemory));
}

std::tuple<operations_research::sat::CpSolverStatus, std::vector<uint64_t>,
           uint64_t>
CPSATSolver::solve(uint64_t timeLimitSeconds, uint64_t numWorkers) {
    uint64_t n = graph.getSize();

    // Compute pruning bounds
    computePruningBounds();

    // Get heuristic schedule
    rpo::RPOScheduler rpoScheduler(graph, debug);
    auto [heuristicSchedule, heuristicPeak] = rpoScheduler.solve();

    // Build model
    createVariablesAndConstraints(heuristicPeak);

    // Set warm start hints
    setWarmStartHints(heuristicSchedule, heuristicPeak);

    // Configure solver
    operations_research::sat::SatParameters params;
    params.set_max_time_in_seconds(static_cast<double>(timeLimitSeconds));
    params.set_num_search_workers(numWorkers);
    params.set_log_search_progress(debug);

    if (debug) {
        std::cerr << "\nSolving with CP-SAT..." << std::endl;
    }

    // Solve
    operations_research::sat::Model model;
    model.Add(NewSatParameters(params));

    auto start = std::chrono::steady_clock::now();
    operations_research::sat::CpSolverResponse response =
        SolveCpModel(cpModel.Build(), &model);
    auto end = std::chrono::steady_clock::now();
    auto elapsed =
        std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

    if (debug) {
        std::cerr << "Solver status: " << CpSolverResponseStats(response)
                  << std::endl;
        std::cerr << "Measured elapsed: " << elapsed << " seconds" << std::endl;
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

namespace bruteforce {

BruteForceSolver::BruteForceSolver(const HyperGraph &graph, bool debug)
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
        std::vector<uint64_t> toRemove;
        for (uint64_t t : live) {
            bool canFree = true;
            for (uint64_t i = 0; i < n; i++) {
                if (inputTensors[i].count(t)) {
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
                     const std::vector<uint64_t> &partitionMapping, bool debug,
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

std::pair<std::vector<uint64_t>, uint64_t>
Scheduler::run(uint64_t timeLimitSeconds, uint64_t numWorkers) {
    buildCoarseGraph();
    HyperGraph hyperGraph(*coarseGraph);

    cpsat::CPSATSolver solver(hyperGraph, debug);
    auto [status, schedule, peakMemory] =
        solver.solve(timeLimitSeconds, numWorkers);

    if (verify) {
        bruteforce::BruteForceSolver bfSolver(hyperGraph, debug);
        auto [bfSchedule, bfPeakMemory] = bfSolver.solve();

        if (peakMemory != bfPeakMemory) {
            std::cerr << "\n=== VERIFICATION FAILED: DUMP START ===\n"
                      << std::endl;

            std::cerr << "=== Coarse Graph ===\n";
            std::cerr << *coarseGraph;
            std::cerr << std::endl;

            std::cerr << "=== HyperGraph ===\n";
            std::cerr << hyperGraph;
            std::cerr << std::endl;

            std::cerr << "=== CP-SAT Solution ===\n";
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
                "Verification failed: CP-SAT peak memory (" +
                std::to_string(peakMemory) + ") != brute-force peak memory (" +
                std::to_string(bfPeakMemory) + ")");
        }
    }

    return {schedule, peakMemory};
}

} // namespace scheduling

} // namespace dag_partitioning