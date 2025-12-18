/**
 * @file HyperGraph.h
 * @brief Hyper graph representation for memory-aware scheduling
 *
 * Represents a DAG as a hyper graph where nodes are operations and
 * hyper edges represent tensor dependencies. Provides methods to
 * analyze the graph topology for scheduling purposes.
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

class HyperGraph {
  private:
    uint64_t size;
    std::vector<uint64_t> tensorSizes; // S_i: output tensor size
    std::vector<uint64_t> extraSizes;  // ES_i: extra memory during execution
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
     * @brief Constructs a hyper graph from a DAG
     * @param graph Input graph to analyze
     */
    HyperGraph(const core::Graph &graph);

    virtual ~HyperGraph() = default;

    /**
     * @brief Returns number of nodes
     * @return Graph size
     */
    uint64_t getSize() const;

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