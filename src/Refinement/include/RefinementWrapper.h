//
// Created by panagiotis on 17/1/2025.
//

#ifndef DAG_PARTITIONING_REFINEMENTWRAPPER_H
#define DAG_PARTITIONING_REFINEMENTWRAPPER_H

#include <cstdint>
#include <vector>

namespace dag_partitioning {

namespace core {
class Graph;
}

namespace refinement {

enum class RefinementMethod;

void refinementWrapper(const core::Graph &graph,
                       std::vector<uint8_t> &bisectionInfo, uint64_t &edgeCut,
                       RefinementMethod refinementMethod,
                       uint64_t refinementPasses, double upperBoundPartWeight,
                       double lowerBoundPartWeight);

} // namespace refinement

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_REFINEMENTWRAPPER_H
