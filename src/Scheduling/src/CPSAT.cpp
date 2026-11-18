/**
 * @file CPSAT.cpp
 * @brief What the CP-SAT solvers share
 */

#include "CPSAT.h"
#include "RPOScheduler.h"

#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include <iostream>
#include <numeric>

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

std::vector<uint64_t>
warmStart(const hypergraph::HyperGraph &graph,
          const std::function<uint64_t(const std::vector<uint64_t> &)> &cost,
          bool debug) {
    std::vector<uint64_t> rpo =
        heuristic::RPOScheduler(graph, debug).solve().first;
    std::vector<uint64_t> ids(graph.getSize());
    std::iota(ids.begin(), ids.end(), 0);
    for (uint64_t t = 0; t < graph.getNumTensors(); t++) {
        for (uint64_t user : graph.getTensorUsers()[t]) {
            if (graph.getTensorProducers()[t] > user)
                return rpo;
        }
    }
    uint64_t rpoCost = cost(rpo), idsCost = cost(ids);
    if (debug) {
        std::cerr << "Warm start: RPO " << rpoCost << ", node ids " << idsCost
                  << std::endl;
    }
    return idsCost < rpoCost ? ids : rpo;
}

} // namespace cpsat

} // namespace scheduling

} // namespace dag_partitioning
