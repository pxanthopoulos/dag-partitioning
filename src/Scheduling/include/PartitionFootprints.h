/**
 * @file PartitionFootprints.h
 * @brief Device memory of each partition, given the partitions next to it
 */

#ifndef DAG_PARTITIONING_SCHEDULING_PARTITIONFOOTPRINTS_H
#define DAG_PARTITIONING_SCHEDULING_PARTITIONFOOTPRINTS_H

#include "HyperGraph.h"
#include "Utils.h"

#include <cstdint>
#include <map>
#include <vector>

namespace dag_partitioning {

namespace scheduling {

/**
 * @brief Device memory each partition holds, as a function of its neighbors
 *
 * Partitions run one at a time, with host transfers between them (see
 * hypergraph::HyperGraph::computeTransfers). A partition holds its local
 * tensors and, for each tensor it produces or consumes:
 *  - if the partition right before or right after it also produces or consumes
 *    the tensor, or it runs last and the tensor is live-out, the tensor is kept
 *    on the device, for the whole partition;
 *  - otherwise it is on the device only from its first to its last use in the
 *    partition (the rest of the time it is on the host), so it can share memory
 *    with what isn't live at the same time.
 * The partition that runs last also holds every live-out tensor, for the whole
 * partition. So the footprint of a partition depends only on its two neighbors.
 * Live-in tensors (hypergraph::LiveIn, produced before the partitions and not
 * moved to the host) are held by every partition up to their last user, for
 * the whole partition; they depend on the order, so footprint leaves them
 * out and peak adds them.
 */
class PartitionFootprints {
  public:
    /**
     * @brief A tensor used by more than one partition, or live-out
     */
    struct Escaping {
        uint64_t size;
        std::vector<uint64_t> members; // Producer and consumers, sorted
        bool liveOut;
    };

    /**
     * @brief The uses of an escaping tensor in one partition
     */
    struct Use {
        uint64_t tensor; // Index in the escaping tensors
        uint64_t first;  // Local position of the first use
        uint64_t last;   // Local position of the last use
    };

    static constexpr uint64_t None = UINT64_MAX; // No relevant neighbor

    /**
     * @param localTensors Per partition, the tensors only it uses
     * @param tensors The escaping tensors
     * @param uses Per partition, its uses of the escaping tensors
     * @param packer Packs the tensors of a partition into a slab
     * @param liveIns The live-in tensors; users are partitions
     */
    PartitionFootprints(std::vector<std::vector<utils::Tensor>> localTensors,
                        std::vector<Escaping> tensors,
                        std::vector<std::vector<Use>> uses,
                        utils::SlabPacker packer,
                        std::vector<hypergraph::LiveIn> liveIns = {});

    /**
     * @brief The live-in tensors (users are partitions)
     */
    [[nodiscard]] const std::vector<hypergraph::LiveIn> &getLiveIns() const {
        return liveIns;
    }

    /**
     * @brief Footprint of a partition, without the live-in tensors
     * @param partition The partition
     * @param pred The partition right before it, or None
     * @param succ The partition right after it, None, or the number of
     * partitions if it runs last
     */
    [[nodiscard]] uint64_t footprint(uint64_t partition, uint64_t pred,
                                     uint64_t succ) const;

    /**
     * @brief Largest footprint of an order of all partitions, with the
     * live-in tensors each holds
     */
    [[nodiscard]] uint64_t peak(const std::vector<uint64_t> &order) const;

  private:
    std::vector<std::vector<utils::Tensor>> localTensors;
    std::vector<Escaping> tensors;
    std::vector<std::vector<Use>> uses;
    utils::SlabPacker packer;
    std::vector<hypergraph::LiveIn> liveIns;
    // Per partition, the slab size by which of its uses are kept
    mutable std::vector<std::map<std::vector<bool>, uint64_t>> slabs;
};

} // namespace scheduling

} // namespace dag_partitioning

#endif // DAG_PARTITIONING_SCHEDULING_PARTITIONFOOTPRINTS_H
