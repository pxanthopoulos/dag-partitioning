/**
 * @file CPSATTransferMemorySolver.h
 * @brief CP-SAT solver that minimizes the bytes moved between device and host
 * under a limit on the peak device footprint
 */

#ifndef DAG_PARTITIONING_SCHEDULING_CPSATTRANSFERMEMORYSOLVER_H
#define DAG_PARTITIONING_SCHEDULING_CPSATTRANSFERMEMORYSOLVER_H

#include "HyperGraph.h"
#include "PartitionFootprints.h"

#include "ortools/sat/cp_model.h"
#include <cstdint>
#include <tuple>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace cpsat {

/**
 * @brief CP-SAT based solver that minimizes the bytes moved between device
 * and host (hypergraph::HyperGraph::computeTransfers), subject to a limit on
 * the peak device footprint (PartitionFootprints::peak)
 *
 * Each node gets a position, and a boolean per ordered pair of nodes that
 * share a tensor (its producer or consumers) tells if the second runs right
 * after the first. A consumer needs a fetch unless another member of the tensor
 * runs right before it, and any fetch of a tensor needs its one offload; the
 * objective is the bytes of all fetches and offloads.
 *
 * The footprint of a node depends only on its two neighbors, so the limit
 * forbids the pairs of neighbors that exceed it. If no order meets the limit,
 * the solver first minimizes the peak, then the bytes moved at that peak.
 */
class CPSATTransferMemorySolver {
  protected:
    const hypergraph::HyperGraph &graph;
    const PartitionFootprints &footprints;
    uint64_t memoryLimit = UINT64_MAX;
    bool debug = false;

    /**
     * @brief What one solve minimizes
     */
    enum class Goal {
        Transfers, // Bytes moved, with the peak at most the limit
        Peak,      // Peak, without a limit
    };

    /**
     * @brief Result of one solve; the status is UNKNOWN or INFEASIBLE if there
     * is no schedule
     */
    struct Result {
        operations_research::sat::CpSolverStatus status;
        std::vector<uint64_t> schedule;
        uint64_t cost;
        [[nodiscard]] bool found() const;
    };

    [[nodiscard]] Result solveOnce(Goal goal, uint64_t limit,
                                   const std::vector<uint64_t> &hint,
                                   uint64_t timeLimitSeconds,
                                   int32_t numWorkers,
                                   bool deterministic) const;

  public:
    /**
     * @brief Constructor
     * @param graph Hyper graph representation
     * @param footprints Footprints of the nodes of the graph
     * @param memoryLimit Limit on the peak; UINT64_MAX means no limit
     * @param debug Enable debug output
     */
    CPSATTransferMemorySolver(const hypergraph::HyperGraph &graph,
                              const PartitionFootprints &footprints,
                              uint64_t memoryLimit = UINT64_MAX,
                              bool debug = false);

    virtual ~CPSATTransferMemorySolver() = default;

    /**
     * @brief Solves the scheduling problem using CP-SAT, warm-started with
     * the RPO order or the node id order (cpsat::warmStart)
     *
     * @param timeLimitSeconds Time limit for each solve in seconds
     * @param numWorkers Number of parallel workers for the solver
     * @param deterministic Same result for any number of workers, unless the
     * time limit is hit; slower, as the workers wait for each other
     * @return Tuple of (solver status, schedule, transferred bytes); the status
     * is OPTIMAL only if every solve was
     */
    [[nodiscard]] std::tuple<operations_research::sat::CpSolverStatus,
                             std::vector<uint64_t>, uint64_t>
    solve(uint64_t timeLimitSeconds = 600, int32_t numWorkers = 1,
          bool deterministic = true) const;
};

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_CPSATTRANSFERMEMORYSOLVER_H
