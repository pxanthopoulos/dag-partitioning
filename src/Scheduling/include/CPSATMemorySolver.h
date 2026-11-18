/**
 * @file CPSATMemorySolver.h
 * @brief CP-SAT solver that minimizes the peak memory of a schedule
 */

#ifndef DAG_PARTITIONING_SCHEDULING_CPSATMEMORYSOLVER_H
#define DAG_PARTITIONING_SCHEDULING_CPSATMEMORYSOLVER_H

#include "HyperGraph.h"

#include "ortools/sat/cp_model.h"
#include "robin_hood.h"
#include <cstdint>
#include <tuple>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace cpsat {

/**
 * @brief CP-SAT based solver for memory-aware scheduling
 */
class CPSATMemorySolver {
  public:
    /**
     * @brief Seconds the last solve() spent building the model (pruning
     * bounds, RPO warm start, variables, constraints, hints) and solving it
     */
    [[nodiscard]] double constructionSeconds() const {
        return lastConstructionSeconds;
    }
    [[nodiscard]] double solveSeconds() const { return lastSolveSeconds; }

  protected:
    const hypergraph::HyperGraph &graph;
    bool debug = false;
    double lastConstructionSeconds = 0;
    double lastSolveSeconds = 0;

    // Pruning bounds computed from topology
    std::vector<uint64_t> earliestOp;
    std::vector<uint64_t> latestOp;
    std::vector<uint64_t> earliestTensor;
    std::vector<uint64_t> latestTensor;

    // CP-SAT model and variables
    operations_research::sat::CpModelBuilder cpModel;
    robin_hood::unordered_map<uint64_t, operations_research::sat::BoolVar> O;
    robin_hood::unordered_map<uint64_t, operations_research::sat::BoolVar> T;
    // Per live-in tensor, per step up to the latest step of its last user:
    // whether it is still on the device
    std::vector<std::vector<operations_research::sat::BoolVar>> L;
    operations_research::sat::IntVar mem;

    /**
     * @brief Computes earliest/latest valid positions for pruning
     */
    void computePruningBounds();

    /**
     * @brief Checks if O[i][j] variable exists (not pruned)
     */
    [[nodiscard]] bool hasO(uint64_t i, uint64_t j) const;

    /**
     * @brief Checks if T[t][j] variable exists (not pruned)
     */
    [[nodiscard]] bool hasT(uint64_t t, uint64_t j) const;

    /**
     * @brief Returns key for O variable storage
     */
    [[nodiscard]] uint64_t oKey(uint64_t i, uint64_t j) const;

    /**
     * @brief Returns key for T variable storage
     */
    [[nodiscard]] uint64_t tKey(uint64_t t, uint64_t j) const;

    /**
     * @brief Retrieves O[i][j] variable
     */
    [[nodiscard]] operations_research::sat::BoolVar getO(uint64_t i,
                                                         uint64_t j) const;

    /**
     * @brief Retrieves T[t][j] variable
     */
    [[nodiscard]] operations_research::sat::BoolVar getT(uint64_t t,
                                                         uint64_t j) const;

    /**
     * @brief Creates CP-SAT variables and constraints
     *
     * Builds the complete CP-SAT model including O[i][j], T[t][j],
     * mem variables and all scheduling/memory constraints.
     *
     * @param memoryUpperBound Upper bound for peak memory variable
     */
    void createVariablesAndConstraints(uint64_t memoryUpperBound);

    /**
     * @brief Sets warm start hints from a heuristic schedule
     *
     * @param schedule Heuristic schedule where schedule[step] = operator
     * @param peakMemory Peak memory of the heuristic schedule
     */
    void setWarmStartHints(const std::vector<uint64_t> &schedule,
                           uint64_t peakMemory);

  public:
    /**
     * @brief Constructor
     * @param graph Hyper graph representation
     * @param debug Enable debug output
     */
    CPSATMemorySolver(const hypergraph::HyperGraph &graph, bool debug = false);

    virtual ~CPSATMemorySolver() = default;

    /**
     * @brief Solves the scheduling problem using CP-SAT
     *
     * @param timeLimitSeconds Time limit for solver in seconds
     * @param numWorkers Number of parallel workers for the solver
     * @param deterministic Same result for any number of workers, unless the
     * time limit is hit; slower, as the workers wait for each other
     * @return Tuple of (solver status, schedule, peak memory)
     */
    [[nodiscard]] std::tuple<operations_research::sat::CpSolverStatus,
                             std::vector<uint64_t>, uint64_t>
    solve(uint64_t timeLimitSeconds = 600, int32_t numWorkers = 1,
          bool deterministic = true);
};

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_CPSATMEMORYSOLVER_H
