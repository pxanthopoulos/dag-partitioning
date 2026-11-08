/**
 * @file Scheduling.h
 * @brief Memory-aware scheduling for DAGs
 *
 * Performs scheduling of DAG operations to minimize peak memory
 * usage or the bytes transferred between device and host. Uses a CP-SAT
 * solver (brute-force solver for small graphs for verification).
 */

#ifndef DAG_PARTITIONING_SCHEDULING_H
#define DAG_PARTITIONING_SCHEDULING_H

#include "HyperGraph.h"
#include "Utils.h"
#include <Graph.h>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model_solver.h"
#include "robin_hood.h"
#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

namespace dag_partitioning {

namespace core {
class Graph;
}

namespace scheduling {

/**
 * @brief Criterion that a schedule minimizes
 */
enum class Objective {
    Memory, // Peak memory (hypergraph::HyperGraph::computePeakMemory)
    // Bytes moved (hypergraph::HyperGraph::computeTransfers), subject to a
    // limit on the peak device footprint (PartitionFootprints::peak); if no
    // order meets the limit, subject to the lowest peak any order reaches
    TransferMemory,
};

/**
 * @brief Device memory each partition holds, as a function of its neighbors
 *
 * Partitions run one at a time, with host transfers between them (see
 * hypergraph::HyperGraph::computeTransfers). A partition holds its local
 * tensors and, for each tensor it produces or consumes:
 *  - if the partition right before or right after it also produces or consumes
 *    the tensor, or it runs last and the tensor is live-out, the tensor is kept
 *    on the device, for the whole partition;
 *  - otherwise it is on the device only from its first to its last use in the
 *    partition (the rest of the time it is on the host), so it can share memory
 *    with what isn't live at the same time.
 * The partition that runs last also holds every live-out tensor, for the whole
 * partition. So the footprint of a partition depends only on its two neighbors.
 */
class PartitionFootprints {
  public:
    /**
     * @brief A tensor used by more than one partition, or live-out
     */
    struct Escaping {
        uint64_t size;
        std::vector<uint64_t> members; // Producer and consumers, sorted
        bool liveOut;
    };

    /**
     * @brief The uses of an escaping tensor in one partition
     */
    struct Use {
        uint64_t tensor; // Index in the escaping tensors
        uint64_t first;  // Local position of the first use
        uint64_t last;   // Local position of the last use
    };

    static constexpr uint64_t None = UINT64_MAX; // No relevant neighbor

    /**
     * @param localTensors Per partition, the tensors only it uses
     * @param tensors The escaping tensors
     * @param uses Per partition, its uses of the escaping tensors
     * @param packer Packs the tensors of a partition into a slab
     */
    PartitionFootprints(std::vector<std::vector<utils::Tensor>> localTensors,
                        std::vector<Escaping> tensors,
                        std::vector<std::vector<Use>> uses,
                        utils::SlabPacker packer);

    /**
     * @brief Footprint of a partition
     * @param partition The partition
     * @param pred The partition right before it, or None
     * @param succ The partition right after it, None, or the number of
     * partitions if it runs last
     */
    [[nodiscard]] uint64_t footprint(uint64_t partition, uint64_t pred,
                                     uint64_t succ) const;

    /**
     * @brief Largest footprint of an order of all partitions
     */
    [[nodiscard]] uint64_t peak(const std::vector<uint64_t> &order) const;

  private:
    std::vector<std::vector<utils::Tensor>> localTensors;
    std::vector<Escaping> tensors;
    std::vector<std::vector<Use>> uses;
    utils::SlabPacker packer;
    // Per partition, the slab size by which of its uses are kept
    mutable std::vector<std::map<std::vector<bool>, uint64_t>> slabs;
};

namespace heuristic {

/**
 * @brief RPO (Reverse Post-Order) heuristic scheduler
 *
 * Computes a schedule using reverse post-order traversal of the DAG.
 * RPO tends to complete subtrees before moving to siblings, which
 * reduces the number of simultaneously live tensors.
 */
class RPOScheduler {
  protected:
    const hypergraph::HyperGraph &graph;
    bool debug = false;

  public:
    /**
     * @brief Constructor
     * @param graph Hyper graph representation
     * @param debug Enable debug output
     */
    RPOScheduler(const hypergraph::HyperGraph &graph, bool debug = false);

    virtual ~RPOScheduler() = default;

    /**
     * @brief Computes the RPO schedule
     * @return Pair of (schedule, peak memory)
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> solve() const;
};

} // namespace heuristic

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
     * the RPO heuristic
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

namespace bruteforce {

class BruteForceSolver {
  protected:
    const hypergraph::HyperGraph &graph;
    bool debug = false;
    Objective objective = Objective::Memory;
    const PartitionFootprints *footprints = nullptr;
    uint64_t memoryLimit = UINT64_MAX;

    /**
     * @brief Cost of an execution order; UINT64_MAX means not allowed
     */
    using Cost = std::function<uint64_t(const std::vector<uint64_t> &)>;

    /**
     * @brief Checks if a node can be placed given current placement state
     *
     * @param node Node to check
     * @param placed Vector indicating which nodes have been placed
     * @return true if all dependencies are satisfied
     */
    [[nodiscard]] bool canPlace(uint64_t node,
                                const std::vector<bool> &placed) const;

    /**
     * @brief Recursively enumerates all valid topological orders
     *
     * @param current Current partial order being built
     * @param placed Vector tracking which nodes have been placed
     * @param bestCost Best cost found so far
     * @param bestOrder Best order found so far
     * @param count Counter for total orders enumerated
     * @param outputs Storage for debug output strings
     * @param cost Cost of a complete order
     */
    void enumerate(std::vector<uint64_t> &current, std::vector<bool> &placed,
                   uint64_t &bestCost, std::vector<uint64_t> &bestOrder,
                   uint64_t &count, std::vector<std::string> &outputs,
                   const Cost &cost) const;

    /**
     * @brief Finds the best order for a cost
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t>
    best(const Cost &cost) const;

  public:
    /**
     * @brief Constructor
     *
     * @param graph Hyper graph representing the scheduling problem
     * @param debug Enable debug output
     * @param objective Criterion to minimize
     * @param footprints Footprints of the nodes (Objective::TransferMemory)
     * @param memoryLimit Limit on the peak (Objective::TransferMemory)
     */
    BruteForceSolver(const hypergraph::HyperGraph &graph, bool debug = false,
                     Objective objective = Objective::Memory,
                     const PartitionFootprints *footprints = nullptr,
                     uint64_t memoryLimit = UINT64_MAX);

    virtual ~BruteForceSolver() = default;

    /**
     * @brief Solves the scheduling problem via brute-force enumeration
     *
     * @return Pair of (schedule, cost)
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t> solve() const;
};

} // namespace bruteforce

/**
 * @brief Solver that orders the partitions
 */
enum class Solver {
    CPSAT, // CP-SAT, warm-started with the RPO heuristic
    RPO,   // The RPO heuristic alone
};

class Scheduler {
  protected:
    const core::Graph &originalGraph;              // Original graph to schedule
    const std::vector<uint64_t> &partitionMapping; // Node to partition mapping
    std::unique_ptr<hypergraph::HyperGraph> coarseGraph; // Coarse graph
    std::unique_ptr<PartitionFootprints> footprints;     // Of the partitions
    bool debug = false;
    bool verify = false;
    uint64_t offsetAlignment = 1;   // Alignment of packed tensor offsets
    uint64_t rangeAlignment = 1;    // Alignment of packed tensor sizes
    utils::SlabPacker packer;       // Packs each partition's tensors
    std::vector<bool> liveOutNodes; // Nodes whose results outlive the graph

    /**
     * @brief Calculates memory requirements for each partition
     *
     * Packs each partition's tensors with the slab packer (best fit by
     * default) to compute its peak memory usage from tensor lifetimes, and
     * collects the tensors that escape each partition. Also sets footprints.
     *
     * Lifetimes are positions within a partition in node id order when the
     * ids are a topological order (as a compiler emits them, in block order),
     * and in topological sort order otherwise.
     *
     * A node with weight 0 whose results are consumed is a tied op: its
     * result aliases the buffer of the input whose edge weight equals the
     * result size (the last-defined one if several match). All values that
     * alias one buffer share its lifetime. A buffer whose values and users
     * are all in one partition is packed into the partition's weight; any
     * other buffer escapes as a tensor of its producing partition, live
     * until its last consuming partition (until the end for graph outputs:
     * the values of liveOutNodes, or of the nodes without consumers if it
     * is empty).
     *
     * @param escapingTensors Output: escaping tensors (producer and users
     * are partition ids)
     * @return Vector of partition weights (peak memory usage)
     */
    [[nodiscard]] std::vector<uint64_t>
    calculatePartitionWeights(std::vector<hypergraph::Tensor> &escapingTensors);

    /**
     * @brief Builds the coarse graph from partitioned original graph
     *
     * Creates a hyper graph where each node represents a partition.
     * Node weights are computed using memory packing, and every partition
     * produces one tensor per buffer that escapes it, with its own size and
     * consuming partitions.
     */
    void buildCoarseGraph();

  public:
    /**
     * @brief Constructor
     * @param originalGraph Graph to schedule
     * @param partitionMapping Vector mapping each node to its partition
     * @param debug Enable debug output
     * @param verify If true, run brute-force solver and verify the cost
     * matches
     * @param offsetAlignment Alignment of packed tensor offsets
     * @param rangeAlignment Alignment of packed tensor sizes and partition
     * weights
     * @param packer Packs each partition's tensors into its weight; empty
     * means best fit (utils::packBestFit) with the two alignments
     * @param liveOutNodes For every node, whether its results are used after
     * the graph (graph outputs); empty means the nodes without consumers. A
     * node can have consumers in the graph and still be an output
     */
    Scheduler(const core::Graph &originalGraph,
              const std::vector<uint64_t> &partitionMapping, bool debug = false,
              bool verify = false, uint64_t offsetAlignment = 1,
              uint64_t rangeAlignment = 1, utils::SlabPacker packer = {},
              std::vector<bool> liveOutNodes = {});

    /**
     * @brief Virtual destructor for proper cleanup
     */
    virtual ~Scheduler() = default;

    /**
     * @brief Executes the scheduling algorithm
     *
     * Builds the coarse graph with partition weights computed via memory
     * packing and inter-partition edge weights.
     *
     * @param timeLimitSeconds Time limit for solver in seconds (CP-SAT only)
     * @param numWorkers Number of parallel workers for the solver (CP-SAT
     * only)
     * @param solver Solver that orders the partitions; RPO ignores the memory
     * limit
     * @param objective Criterion to minimize
     * @param memoryLimit Limit on the peak device footprint of
     * Objective::TransferMemory; UINT64_MAX means no limit
     * @param deterministic Same result for any number of workers, unless the
     * time limit is hit (CP-SAT only)
     * @return Pair of a vector with partition numbers ordered by their
     * scheduled execution and the cost of the schedule (peak memory usage or
     * transferred bytes)
     */
    [[nodiscard]] std::pair<std::vector<uint64_t>, uint64_t>
    run(uint64_t timeLimitSeconds = 600, int32_t numWorkers = 1,
        Solver solver = Solver::CPSAT, Objective objective = Objective::Memory,
        uint64_t memoryLimit = UINT64_MAX, bool deterministic = true);

    /**
     * @brief Peak device footprint of a schedule returned by run
     */
    [[nodiscard]] uint64_t
    peakFootprint(const std::vector<uint64_t> &schedule) const;
};

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_H