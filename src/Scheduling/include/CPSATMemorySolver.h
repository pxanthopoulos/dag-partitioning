/**
 * @file CPSATMemorySolver.h
 * @brief CP-SAT solver that minimizes the peak memory of a schedule
 */

#ifndef DAG_PARTITIONING_SCHEDULING_CPSATMEMORYSOLVER_H
#define DAG_PARTITIONING_SCHEDULING_CPSATMEMORYSOLVER_H

#include "HyperGraph.h"

#include "ortools/sat/cp_model.h"
#include <cstdint>
#include <tuple>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace cpsat {

/**
 * @brief Finds the order of the nodes with the lowest peak memory
 * (hypergraph::HyperGraph::computePeakMemory), warm-started with the RPO
 * heuristic
 *
 * The model is an adaptation of the ILP of Zhong et al.
 * (arXiv:2308.13898): a boolean per node and step that runs it, and per
 * tensor and step that holds it.
 */
class CPSATMemorySolver {
  public:
    explicit CPSATMemorySolver(const hypergraph::HyperGraph &graph,
                               bool debug = false);

    /**
     * @param timeLimitSeconds Time limit of the solver
     * @param numWorkers Number of parallel workers
     * @param deterministic Same result for any number of workers, unless the
     * time limit is hit; slower, as the workers wait for each other
     * @return Solver status, order and its peak memory
     * @throws std::runtime_error If no order is found in time
     */
    [[nodiscard]] std::tuple<operations_research::sat::CpSolverStatus,
                             std::vector<uint64_t>, uint64_t>
    solve(uint64_t timeLimitSeconds = 600, int32_t numWorkers = 1,
          bool deterministic = true);

    /**
     * @brief Seconds the last solve spent building the model (RPO warm start
     * included) and solving it
     */
    [[nodiscard]] double constructionSeconds() const {
        return lastConstructionSeconds;
    }
    [[nodiscard]] double solveSeconds() const { return lastSolveSeconds; }

  private:
    const hypergraph::HyperGraph &graph;
    bool debug;
    double lastConstructionSeconds = 0;
    double lastSolveSeconds = 0;
};

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_CPSATMEMORYSOLVER_H
