/**
 * @file Mixed.h
 * @brief Mix of boundary FM and KL
 *
 * Runs one pass of boundary KL and one pass of boundary FM.
 */

#ifndef DAG_PARTITIONING_MIXED_H
#define DAG_PARTITIONING_MIXED_H

#include "BoundaryKL.h"
#include "BoundaryFM.h"

class Mixed : public BoundaryKL, public BoundaryFM {
private:
    /**
     * @brief Runs one pass of boundary KL and one pass of boundary FM.
     *
     * Returns false (no changes) if both passes made no changes.
     *
     * @return true if improvements were made in this pass
     */
    [[nodiscard]] bool onePassRefinement() override;

public:
    /**
     * @brief Constructs the boundary KL refinement algorithm
     * @param graph Graph being refined
     * @param initialBisectionInfo Current partition assignments
     * @param initialEdgeCut Current edge cut weight
     * @param maxNumberOfPasses Maximum refinement iterations
     * @param upperBoundPartWeight Maximum allowed partition weight
     * @param lowerBoundPartWeight Minimum allowed partition weight
     */
    Mixed(const Graph &graph, std::vector<uint8_t> &initialBisectionInfo,
          uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
          double upperBoundPartWeight, double lowerBoundPartWeight);
};


#endif //DAG_PARTITIONING_MIXED_H
