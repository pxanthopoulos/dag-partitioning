/**
 * @file Objective.h
 * @brief Criteria that a schedule minimizes
 */

#ifndef DAG_PARTITIONING_SCHEDULING_OBJECTIVE_H
#define DAG_PARTITIONING_SCHEDULING_OBJECTIVE_H

namespace dag_partitioning {

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

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_OBJECTIVE_H
