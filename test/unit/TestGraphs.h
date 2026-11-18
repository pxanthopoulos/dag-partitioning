// Small scheduling problems with known answers, shared by the unit tests
#ifndef DAG_PARTITIONING_TEST_TESTGRAPHS_H
#define DAG_PARTITIONING_TEST_TESTGRAPHS_H

#include "HyperGraph.h"
#include "PartitionFootprints.h"
#include "Utils.h"

#include <vector>

namespace testgraphs {

using namespace dag_partitioning::scheduling;

inline utils::SlabPacker bestFit() {
    return [](std::vector<utils::Tensor> &tensors) {
        return utils::packBestFit(tensors);
    };
}

/**
 * @brief Partition 0 makes A, partition 1 makes B and C, partition 2 reads
 * A, then C, then B (4 bytes each; schedmem/FINDINGS.md's minimal example)
 *
 * Order 0 1 2 keeps B and C on the device for all of partition 2: it moves
 * only A (8 bytes) at a peak of 12. Order 1 0 2 keeps A and moves B and C
 * (16 bytes), at a peak of 8. With a live-in of 5 bytes read by partition 0,
 * the peaks are 12 and 9.
 */
struct Minimal {
    hypergraph::HyperGraph graph;
    PartitionFootprints footprints;

    explicit Minimal(bool withLiveIn = false)
        : graph({0, 0, 0},
                {{0, 4, {2}, false}, {1, 4, {2}, false}, {1, 4, {2}, false}},
                liveIns(withLiveIn)),
          footprints(
              {{}, {}, {}},
              {{4, {0, 2}, false}, {4, {1, 2}, false}, {4, {1, 2}, false}},
              {{{0, 0, 0}},
               {{1, 0, 0}, {2, 1, 1}},
               {{0, 0, 0}, {2, 1, 1}, {1, 2, 2}}},
              bestFit(), liveIns(withLiveIn)) {}

    static std::vector<hypergraph::LiveIn> liveIns(bool withLiveIn) {
        if (!withLiveIn)
            return {};
        return {{5, {0}}};
    }
};

/**
 * @brief Node 0 makes a (10 bytes) for node 3, node 1 makes b (10) for node
 * 2, node 2 makes c (1) for node 3. Only order 1 2 0 3 reaches the lowest
 * peak, 11; order 0 1 2 3 peaks at 21
 */
inline hypergraph::HyperGraph twoChains() {
    return hypergraph::HyperGraph(
        {0, 0, 0, 0},
        {{0, 10, {3}, false}, {1, 10, {2}, false}, {2, 1, {3}, false}});
}

} // namespace testgraphs

#endif // DAG_PARTITIONING_TEST_TESTGRAPHS_H
