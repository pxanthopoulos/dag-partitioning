/**
 * @file HyperGraph.h
 * @brief Hyper graph of nodes and the tensors between them, for scheduling
 *
 * Every tensor (hyper edge) has one producer and any number of users; a node
 * can produce any number of tensors.
 */

#ifndef DAG_PARTITIONING_SCHEDULING_HYPERGRAPH_H
#define DAG_PARTITIONING_SCHEDULING_HYPERGRAPH_H

#include "robin_hood.h"
#include <cstdint>
#include <ostream>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace hypergraph {

/**
 * @brief Tensor (hyper edge): one producer, any number of users
 */
struct Tensor {
    uint64_t producer;           // Node that creates the tensor
    uint64_t size;               // Size in bytes
    std::vector<uint64_t> users; // Nodes that read the tensor
    bool liveOut = false;        // Live until the end of the schedule
};

/**
 * @brief Tensor from before the graph (produced outside it, and neither
 * moved to the host nor left out of the memory): on the device from the first
 * step until its last user runs
 */
struct LiveIn {
    uint64_t size;               // Size in bytes
    std::vector<uint64_t> users; // Nodes that read it
};

class HyperGraph {
  public:
    /**
     * @param extraSizes Extra memory of every node while it runs
     * @param tensors Tensors; producers and users index extraSizes
     * @param liveIns Tensors from before the graph; users index extraSizes
     */
    HyperGraph(std::vector<uint64_t> extraSizes,
               const std::vector<Tensor> &tensors,
               std::vector<LiveIn> liveIns = {});

    // Per node
    uint64_t getSize() const { return size; }
    const std::vector<uint64_t> &getExtraSizes() const { return extraSizes; }
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getInputTensors() const {
        return inputTensors;
    }
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getAncestors() const {
        return ancestors;
    }
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getDescendants() const {
        return descendants;
    }

    // Per tensor
    uint64_t getNumTensors() const { return tensorSizes.size(); }
    const std::vector<uint64_t> &getTensorSizes() const { return tensorSizes; }
    const std::vector<uint64_t> &getTensorProducers() const {
        return tensorProducers;
    }
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getTensorUsers() const {
        return tensorUsers;
    }
    const std::vector<bool> &getLiveOut() const { return liveOut; }

    const std::vector<LiveIn> &getLiveIns() const { return liveIns; }

    /**
     * @brief Peak memory of an execution order
     *
     * While order[s] runs, memory holds its extra memory and every tensor t
     * with birth(t) <= s <= death(t): birth is the step of its producer,
     * death the step of its last user (the last step if live-out, the birth
     * if it has no users). A live-in tensor is born at the first step.
     *
     * @param order Topological order of all nodes
     */
    [[nodiscard]] uint64_t
    computePeakMemory(const std::vector<uint64_t> &order) const;

    /**
     * @brief Bytes moved between device and host by an execution order
     *
     * The producer and the consumers of a tensor form runs of consecutive
     * steps (a live-out tensor counts as consumed after the last step). The
     * tensor stays on the device through a run. Each run after the first, a
     * gap, needs a fetch, and any gap needs one offload whose host copy
     * serves every fetch: a tensor with g gaps moves (1 + g) times its size,
     * or nothing if g = 0.
     *
     * @param order Topological order of all nodes
     */
    [[nodiscard]] uint64_t
    computeTransfers(const std::vector<uint64_t> &order) const;

    friend std::ostream &operator<<(std::ostream &os, const HyperGraph &graph);

  private:
    uint64_t size;
    std::vector<LiveIn> liveIns;
    std::vector<uint64_t> extraSizes;
    std::vector<robin_hood::unordered_set<uint64_t>> ancestors, descendants;
    std::vector<robin_hood::unordered_set<uint64_t>> inputTensors;
    std::vector<uint64_t> tensorSizes, tensorProducers;
    std::vector<robin_hood::unordered_set<uint64_t>> tensorUsers;
    std::vector<bool> liveOut;
};

} // namespace hypergraph

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_HYPERGRAPH_H