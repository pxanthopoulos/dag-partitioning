/**
 * @file CPSAT.h
 * @brief What the CP-SAT solvers share
 */

#ifndef DAG_PARTITIONING_SCHEDULING_CPSAT_H
#define DAG_PARTITIONING_SCHEDULING_CPSAT_H

#include "ortools/sat/cp_model.h"
#include <cstdint>

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

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_CPSAT_H
