/**
 * @file CPSAT.cpp
 * @brief What the CP-SAT solvers share
 */

#include "CPSAT.h"

#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"

namespace dag_partitioning {

namespace scheduling {

namespace cpsat {

operations_research::sat::CpSolverResponse
solveModel(const operations_research::sat::CpModelProto &model,
           uint64_t timeLimitSeconds, int32_t numWorkers, bool deterministic,
           bool debug) {
    namespace sat = operations_research::sat;
    sat::SatParameters params;
    params.set_max_time_in_seconds(static_cast<double>(timeLimitSeconds));
    params.set_num_search_workers(numWorkers);
    params.set_interleave_search(deterministic);
    params.set_log_search_progress(debug);
    sat::Model satModel;
    satModel.Add(NewSatParameters(params));
    return SolveCpModel(model, &satModel);
}

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning
