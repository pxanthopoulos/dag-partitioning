/**
 * @file Scheduler.h
 * @brief Orders the partitions of a partitioned graph
 */

#ifndef DAG_PARTITIONING_SCHEDULING_SCHEDULER_H
#define DAG_PARTITIONING_SCHEDULING_SCHEDULER_H

#include "HyperGraph.h"
#include "Objective.h"
#include "PartitionFootprints.h"
#include "Utils.h"
#include <Graph.h>

#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

/**
 * @brief Solver that orders the partitions
 */
enum class Solver {
    CPSAT, // CP-SAT, warm-started with the RPO or the node id order
    RPO,   // The RPO heuristic alone
};

class Scheduler {
  public:
    /**
     * @brief A tensor the graph reads but that is produced before it, which
     * stays on the device until its last user
     */
    struct LiveIn {
        uint64_t size;                   // Size in bytes
        std::vector<uint64_t> users;     // Nodes that read it
        std::vector<uint64_t> tiedUsers; // Users that write into it in place
    };

  protected:
    const core::Graph &originalGraph;              // Original graph to schedule
    const std::vector<uint64_t> &partitionMapping; // Node to partition mapping
    std::unique_ptr<hypergraph::HyperGraph> coarseGraph; // Coarse graph
    std::unique_ptr<PartitionFootprints> footprints;     // Of the partitions
    bool debug = false;
    bool verify = false;
    utils::SlabPacker packer;       // Packs each partition's tensors
    std::vector<bool> liveOutNodes; // Nodes whose results outlive the graph
    std::vector<LiveIn> liveIns;    // Tensors from before the graph

    /**
     * @brief Builds the coarse graph and the footprints of its nodes
     *
     * Every node of the coarse graph is a partition. Its weight is the peak
     * memory of the buffers that never leave it, packed with the slab packer
     * (best fit by default). Every buffer that escapes it is a tensor of its
     * own size, live until its last consuming partition (until the end for
     * graph outputs). Values that a tied op (weight 0) writes in place share
     * the buffer of its input. Partition edges that no tensor orders become
     * empty tensors.
     */
    void buildCoarseGraph();

  public:
    /**
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
     * @param liveIns Tensors produced before the graph that stay on the device
     * until their last user. A tied user (weight 0, writing into it in place)
     * starts a buffer of the live-in's size, which escapes like any other
     */
    Scheduler(const core::Graph &originalGraph,
              const std::vector<uint64_t> &partitionMapping, bool debug = false,
              bool verify = false, uint64_t offsetAlignment = 1,
              uint64_t rangeAlignment = 1, utils::SlabPacker packer = {},
              std::vector<bool> liveOutNodes = {},
              std::vector<LiveIn> liveIns = {});

    /**
     * @brief Builds the coarse graph and orders its nodes, the partitions
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
     * @return The partitions in execution order, and the cost of that order
     * (peak memory or transferred bytes)
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

#endif // DAG_PARTITIONING_SCHEDULING_SCHEDULER_H
