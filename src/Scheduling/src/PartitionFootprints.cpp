/**
 * @file PartitionFootprints.cpp
 * @brief Device memory of each partition, given the partitions next to it
 */

#include "PartitionFootprints.h"

#include <algorithm>

namespace dag_partitioning {

namespace scheduling {

PartitionFootprints::PartitionFootprints(
    std::vector<std::vector<utils::Tensor>> localTensors,
    std::vector<Escaping> tensors, std::vector<std::vector<Use>> uses,
    utils::SlabPacker packer, std::vector<hypergraph::LiveIn> liveIns)
    : localTensors(std::move(localTensors)), tensors(std::move(tensors)),
      uses(std::move(uses)), packer(std::move(packer)),
      liveIns(std::move(liveIns)), slabs(this->localTensors.size()) {}

uint64_t PartitionFootprints::footprint(uint64_t partition, uint64_t pred,
                                        uint64_t succ) const {
    bool runsLast = succ == localTensors.size();
    auto isMember = [&](uint64_t tensor, uint64_t node) {
        const auto &members = tensors[tensor].members;
        return std::binary_search(members.begin(), members.end(), node);
    };

    // Kept tensors count for the whole partition
    const auto &partitionUses = uses[partition];
    std::vector<bool> kept(partitionUses.size());
    uint64_t held = 0;
    for (uint64_t u = 0; u < partitionUses.size(); u++) {
        uint64_t tensor = partitionUses[u].tensor;
        kept[u] = isMember(tensor, pred) || isMember(tensor, succ) ||
                  (runsLast && tensors[tensor].liveOut);
        if (kept[u])
            held += tensors[tensor].size;
    }
    // The last partition also holds the live-outs it doesn't use: it fetches
    // them before it ends
    if (runsLast) {
        for (uint64_t t = 0; t < tensors.size(); t++) {
            if (tensors[t].liveOut && !isMember(t, partition))
                held += tensors[t].size;
        }
    }

    // The others share the slab with the local tensors
    auto [it, inserted] = slabs[partition].try_emplace(kept, 0);
    if (inserted) {
        std::vector<utils::Tensor> slab = localTensors[partition];
        for (uint64_t u = 0; u < partitionUses.size(); u++) {
            if (!kept[u])
                slab.push_back({tensors[partitionUses[u].tensor].size,
                                partitionUses[u].first, partitionUses[u].last,
                                0});
        }
        utils::sortByLifetime(slab);
        it->second = slab.empty() ? 0 : packer(slab);
    }
    return held + it->second;
}

uint64_t PartitionFootprints::peak(const std::vector<uint64_t> &order) const {
    // The live-in bytes held at each step: every live-in up to its last user
    std::vector<uint64_t> held(order.size() + 1, 0);
    if (!liveIns.empty()) {
        std::vector<uint64_t> pos(localTensors.size());
        for (uint64_t step = 0; step < order.size(); step++)
            pos[order[step]] = step;
        std::vector<int64_t> delta(order.size() + 1, 0);
        for (const auto &liveIn : liveIns) {
            if (liveIn.users.empty())
                continue;
            uint64_t last = 0;
            for (uint64_t user : liveIn.users)
                last = std::max(last, pos[user]);
            delta[0] += static_cast<int64_t>(liveIn.size);
            delta[last + 1] -= static_cast<int64_t>(liveIn.size);
        }
        int64_t live = 0;
        for (uint64_t step = 0; step < order.size(); step++) {
            live += delta[step];
            held[step] = static_cast<uint64_t>(live);
        }
    }
    uint64_t peak = 0;
    for (uint64_t step = 0; step < order.size(); step++) {
        peak = std::max(
            peak, held[step] +
                      footprint(order[step], step == 0 ? None : order[step - 1],
                                step + 1 == order.size() ? localTensors.size()
                                                         : order[step + 1]));
    }
    return peak;
}

} // namespace scheduling

} // namespace dag_partitioning
