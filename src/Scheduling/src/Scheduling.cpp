/**
 * @file Scheduling.cpp
 * @brief Implementation of memory-aware scheduling algorithms
 */

#include "Scheduling.h"

#include <cassert>
#include <queue>
#include <set>
#include <sstream>
#include <stack>

namespace dag_partitioning {

namespace scheduling {

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

CPSATSolver::CPSATSolver(const hypergraph::HyperGraph &graph, bool debug)
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
void CPSATSolver::computePruningBounds() {
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

bool CPSATSolver::hasO(uint64_t i, uint64_t j) const {
    return earliestOp[i] <= j && j <= latestOp[i];
}

bool CPSATSolver::hasT(uint64_t t, uint64_t j) const {
    return earliestTensor[t] <= j && j <= latestTensor[t];
}

uint64_t CPSATSolver::oKey(uint64_t i, uint64_t j) const {
    return i * graph.getSize() + j;
}

uint64_t CPSATSolver::tKey(uint64_t t, uint64_t j) const {
    return t * graph.getSize() + j;
}

operations_research::sat::BoolVar CPSATSolver::getO(uint64_t i,
                                                    uint64_t j) const {
    return O.at(oKey(i, j));
}

operations_research::sat::BoolVar CPSATSolver::getT(uint64_t t,
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
 * (6) Σ_t T[t][j]·S[t] + Σᵢ ES[i]·O[i][j] ≤ mem  -- memory bound
 *
 * Objective: minimize mem
 */
void CPSATSolver::createVariablesAndConstraints(uint64_t memoryUpperBound) {
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

        cpModel.AddLessOrEqual(memExpr, mem);
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

    // Add hint for mem variable
    cpModel.AddHint(mem, static_cast<int64_t>(peakMemory));
}

std::tuple<operations_research::sat::CpSolverStatus, std::vector<uint64_t>,
           uint64_t>
CPSATSolver::solve(uint64_t timeLimitSeconds, int32_t numWorkers) {
    uint64_t n = graph.getSize();

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

BruteForceSolver::BruteForceSolver(const hypergraph::HyperGraph &graph,
                                   bool debug)
    : graph(graph), debug(debug) {}

uint64_t
BruteForceSolver::computePeakMemory(const std::vector<uint64_t> &order) const {
    return graph.computePeakMemory(order);
}

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

std::vector<uint64_t> Scheduler::calculatePartitionWeights(
    std::vector<hypergraph::Tensor> &escapingTensors) const {
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

    // Buffer of every value: a tied op (weight 0, consumed results) writes
    // into the input whose edge weight equals its result size, preferring
    // the last-defined one (the init operand)
    std::vector<uint64_t> root(originalGraph.size, UINT64_MAX);
    for (uint64_t node : topologicalOrder) {
        if (originalGraph.nodeWeights[node] > 0) {
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

    for (uint64_t producer = 0; producer < originalGraph.size; ++producer) {
        if (aliases[producer].empty())
            continue;
        uint64_t producerPartition = partitionMapping[producer];

        // Find all users of the buffer's values within the same partition,
        // and the partitions it escapes to
        uint64_t lastConsumerPosition =
            nodeToLocalPosition[producer]; // Birth time
        hypergraph::Tensor escaping{
            producerPartition, originalGraph.nodeWeights[producer], {}, false};
        auto use = [&](uint64_t node) {
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
            if (originalGraph.adj[value].empty())
                escaping.liveOut = true;
            for (const auto &[consumer, edgeWeight] : originalGraph.adj[value])
                use(consumer);
        }

        // Tensor is packed only if it has internal users only
        if (escaping.users.empty() && !escaping.liveOut) {
            utils::Tensor tensor;
            tensor.size = originalGraph.nodeWeights[producer];
            tensor.birth = nodeToLocalPosition[producer];
            tensor.death = lastConsumerPosition;
            tensor.offset = 0; // Will be assigned during packing
            partitionTensors[producerPartition].push_back(tensor);
        } else {
            std::sort(escaping.users.begin(), escaping.users.end());
            escaping.users.erase(
                std::unique(escaping.users.begin(), escaping.users.end()),
                escaping.users.end());
            escapingTensors.push_back(std::move(escaping));
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
                  [](const utils::Tensor &a, const utils::Tensor &b) {
                      if (a.birth != b.birth)
                          return a.birth < b.birth;
                      return a.death < b.death;
                  });

        // Best-fit packing algorithm
        uint64_t peakMemory = utils::packBestFit(tensors);

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
        std::move(partitionWeights), tensors);
}

std::pair<std::vector<uint64_t>, uint64_t>
Scheduler::run(uint64_t timeLimitSeconds, int32_t numWorkers) {
    buildCoarseGraph();
    const hypergraph::HyperGraph &hyperGraph = *coarseGraph;

    cpsat::CPSATSolver solver(hyperGraph, debug);
    auto [status, schedule, peakMemory] =
        solver.solve(timeLimitSeconds, numWorkers);

    if (verify) {
        bruteforce::BruteForceSolver bfSolver(hyperGraph, debug);
        auto [bfSchedule, bfPeakMemory] = bfSolver.solve();

        if (peakMemory != bfPeakMemory) {
            std::cerr << "\n=== VERIFICATION FAILED: DUMP START ===\n"
                      << std::endl;

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