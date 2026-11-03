/**
 * @file HyperGraph.h
 * @brief Hyper graph representation for memory-aware scheduling
 *
 * Represents a DAG as a hyper graph where nodes are operations and
 * hyper edges represent tensor dependencies. A node may produce any
 * number of tensors, each with its own size and users. Provides methods
 * to analyze the graph topology for scheduling purposes.
 */

#ifndef DAG_PARTITIONING_SCHEDULING_HYPERGRAPH_H
#define DAG_PARTITIONING_SCHEDULING_HYPERGRAPH_H

#include "Graph.h"

#include "robin_hood.h"
#include <cstdint>
#include <iostream>
#include <vector>

namespace dag_partitioning {

namespace core {
class Graph;
}

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

class HyperGraph {
  private:
    uint64_t size;
    std::vector<uint64_t> tensorSizes;     // S_t: tensor size
    std::vector<uint64_t> tensorProducers; // Producer of tensor t
    std::vector<bool> liveOut;             // Tensor t is live until the end
    std::vector<uint64_t> extraSizes; // ES_i: extra memory during execution
    std::vector<std::vector<uint64_t>> outputTensors; // OUT_i: output tensors
    std::vector<robin_hood::unordered_set<uint64_t>>
        inputTensors; // IN_i: input tensors

    std::vector<robin_hood::unordered_set<uint64_t>> ancestors;
    std::vector<robin_hood::unordered_set<uint64_t>> descendants;
    std::vector<robin_hood::unordered_set<uint64_t>> tensorUsers;

    /**
     * @brief Computes topological relationships
     */
    void computeTopology();

    /**
     * @brief Computes transitive closure of predecessors
     */
    void computeAncestors();

    /**
     * @brief Computes transitive closure of successors
     */
    void computeDescendants();

    /**
     * @brief Outputs graph structure for debugging
     * @param os Output stream
     */
    void print(std::ostream &os) const;

  public:
    /**
     * @brief Constructs a hyper graph from its tensors
     * @param extraSizes Extra memory of every node during its execution
     * @param tensors Tensors; producers and users index extraSizes
     */
    HyperGraph(std::vector<uint64_t> extraSizes,
               const std::vector<Tensor> &tensors);

    virtual ~HyperGraph() = default;

    /**
     * @brief Returns number of nodes
     * @return Graph size
     */
    uint64_t getSize() const;

    /**
     * @brief Returns number of tensors
     */
    uint64_t getNumTensors() const;

    /**
     * @brief Returns the producer of every tensor
     */
    const std::vector<uint64_t> &getTensorProducers() const;

    /**
     * @brief Returns whether every tensor is live until the end
     */
    const std::vector<bool> &getLiveOut() const;

    /**
     * @brief Returns the tensors every node produces
     */
    const std::vector<std::vector<uint64_t>> &getOutputTensors() const;

    /**
     * @brief Returns ancestor sets for all nodes
     * @return Vector of ancestor sets
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getAncestors() const;

    /**
     * @brief Returns descendant sets for all nodes
     * @return Vector of descendant sets
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getDescendants() const;

    /**
     * @brief Returns tensor consumer information
     * @return Vector of sets indicating which operations consume each tensor
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getTensorUsers() const;

    /**
     * @brief Returns output tensor sizes
     * @return Vector of tensor sizes
     */
    const std::vector<uint64_t> &getTensorSizes() const;

    /**
     * @brief Returns extra memory required during operation execution
     * @return Vector of workspace memory sizes
     */
    const std::vector<uint64_t> &getExtraSizes() const;

    /**
     * @brief Returns input tensor dependencies
     * @return Vector of sets indicating which tensors each operation needs
     */
    const std::vector<robin_hood::unordered_set<uint64_t>> &
    getInputTensors() const;

    /**
     * @brief Peak memory of an execution order
     *
     * While order[s] runs, memory holds its extra memory and every tensor t
     * with birth(t) <= s <= death(t): birth is the step of its producer,
     * death the step of its last user (the last step if live-out, the birth
     * if it has no users).
     *
     * @param order Topological order of all nodes
     * @return Peak memory usage
     */
    [[nodiscard]] uint64_t
    computePeakMemory(const std::vector<uint64_t> &order) const;

    /**
     * @brief Bytes moved between device and host by an execution order
     *
     * A tensor stays on the device only if its one consumer runs right
     * after its producer; a live-out tensor counts as consumed after the
     * last step. Otherwise it is offloaded once and fetched once for every
     * consumer that does not run right after its producer.
     *
     * @param order Topological order of all nodes
     * @return Transferred bytes (offloads plus fetches)
     */
    [[nodiscard]] uint64_t
    computeTransfers(const std::vector<uint64_t> &order) const;

    /**
     * @brief Stream output operator
     * @param os Output stream
     * @param graph Graph to output
     * @return Output stream
     */
    friend std::ostream &operator<<(std::ostream &os, const HyperGraph &graph);
};

} // namespace hypergraph

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_HYPERGRAPH_H