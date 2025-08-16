/**
 * @file Mixed.h
 * @brief Mix of boundary FM and KL
 *
 * Runs one pass of boundary KL and one pass of boundary FM.
 */

#include "Mixed.h"

Mixed::Mixed(const Graph &graph, std::vector<uint8_t> &initialBisectionInfo,
             uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
             double upperBoundPartWeight, double lowerBoundPartWeight)
    : Refinement(graph, initialBisectionInfo, initialEdgeCut, maxNumberOfPasses,
                 upperBoundPartWeight, lowerBoundPartWeight),
      BoundaryKL(graph, initialBisectionInfo, initialEdgeCut, maxNumberOfPasses,
                 upperBoundPartWeight, lowerBoundPartWeight),
      BoundaryFM(graph, initialBisectionInfo, initialEdgeCut, maxNumberOfPasses,
                 upperBoundPartWeight, lowerBoundPartWeight) {}

bool Mixed::onePassRefinement() {
    bool changed = BoundaryKL::onePassRefinement();
    changed = changed || BoundaryFM::onePassRefinement();
    return changed;
}
