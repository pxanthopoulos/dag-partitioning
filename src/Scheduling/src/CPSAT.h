/**
 * @file CPSAT.h
 * @brief What the CP-SAT solvers share
 */

#ifndef DAG_PARTITIONING_SCHEDULING_CPSAT_H
#define DAG_PARTITIONING_SCHEDULING_CPSAT_H

#include "HyperGraph.h"

#include "ortools/sat/cp_model.h"
#include <cstdint>
#include <functional>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

namespace cpsat {

/**
 * @brief Solves a model
 *
 * @param model The model
 * @param timeLimitSeconds Time limit in seconds
 * @param numWorkers Number of parallel workers
 * @param deterministic Interleaved search: the same result for any number of
 * workers, unless the time limit is hit
 * @param debug Log the search
 */
[[nodiscard]] operations_research::sat::CpSolverResponse
solveModel(const operations_research::sat::CpModelProto &model,
           uint64_t timeLimitSeconds, int32_t numWorkers, bool deterministic,
           bool debug);

/**
 * @brief The order to warm-start from: the cheaper of the RPO order and the
 * node id order, if the ids are a topological order (as a compiler numbers
 * ops, in program order); RPO on ties
 */
[[nodiscard]] std::vector<uint64_t>
warmStart(const hypergraph::HyperGraph &graph,
          const std::function<uint64_t(const std::vector<uint64_t> &)> &cost,
          bool debug);

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_CPSAT_H
