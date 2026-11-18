/**
 * @file Scheduler.cpp
 * @brief Orders the partitions of a partitioned graph
 */

#include "Scheduler.h"
#include "BruteForceSolver.h"
#include "CPSATMemorySolver.h"
#include "CPSATTransferMemorySolver.h"
#include "RPOScheduler.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace dag_partitioning {

namespace scheduling {

namespace {

/**
 * @brief Where the nodes run: lifetimes are positions within a partition in
 * node id order when the ids are a topological order (as a compiler emits
 * them, in block order), and in topological sort order otherwise
 */
struct NodeOrder {
    std::vector<uint64_t> topological; // All nodes
    std::vector<uint64_t> global;      // Per node, its position in topological
    std::vector<uint64_t> local; // Per node, its position in its partition
};

NodeOrder orderNodes(const core::Graph &graph,
                     const std::vector<uint64_t> &partitionMapping,
                     uint64_t numPartitions) {
    NodeOrder order;
    bool idsAreTopological = true;
    for (uint64_t from = 0; from < graph.size; ++from) {
        for (const auto &[to, edgeWeight] : graph.adj[from])
            idsAreTopological &= from < to;
    }
    order.topological.resize(graph.size);
    if (idsAreTopological)
        std::iota(order.topological.begin(), order.topological.end(), 0);
    else
        order.topological = graph.topologicalSort();

    order.global.resize(graph.size);
    for (uint64_t pos = 0; pos < order.topological.size(); ++pos)
        order.global[order.topological[pos]] = pos;

    // The nodes of every partition, sorted by global position
    std::vector<std::vector<std::pair<uint64_t, uint64_t>>> partitionNodes(
        numPartitions);
    for (uint64_t node = 0; node < graph.size; ++node) {
        partitionNodes[partitionMapping[node]].emplace_back(order.global[node],
                                                            node);
    }
    order.local.resize(graph.size);
    for (auto &nodes : partitionNodes) {
        std::sort(nodes.begin(), nodes.end());
        for (uint64_t localPos = 0; localPos < nodes.size(); ++localPos)
            order.local[nodes[localPos].second] = localPos;
    }
    return order;
}

/**
 * @brief A buffer, named by the node that starts it, with the values that
 * alias it
 */
struct Buffer {
    uint64_t producer;  // The node that starts it
    uint64_t partition; // The partition of the producer
    uint64_t size;
    uint64_t birth; // Local position of the producer
    uint64_t death; // Local position of its last use in the partition
    // First and last use in every partition that uses it
    std::map<uint64_t, std::pair<uint64_t, uint64_t>> usedAt;
    std::vector<uint64_t> users; // The other partitions that use it, sorted
    bool liveOut;                // Used after the graph

    [[nodiscard]] bool escapes() const { return !users.empty() || liveOut; }
};

/**
 * @brief The buffers of all values
 *
 * A node with weight 0 whose results are consumed is a tied op: its result
 * aliases the buffer of the input whose edge weight equals the result size
 * (the last-defined one if several match). All values that alias one buffer
 * share its lifetime. A buffer escapes its partition if another partition
 * uses it, or if it is a graph output (the values of liveOutNodes, or of the
 * nodes without consumers if it is empty).
 *
 * A tied op writing into a live-in starts a buffer of the live-in's size (the
 * live-in is not a node, so it would have no buffer at all).
 */
std::vector<Buffer>
collectBuffers(const core::Graph &graph,
               const std::vector<uint64_t> &partitionMapping,
               uint64_t numPartitions, const std::vector<bool> &liveOutNodes,
               const std::vector<Scheduler::LiveIn> &liveIns) {
    NodeOrder order = orderNodes(graph, partitionMapping, numPartitions);

    std::vector<uint64_t> bufferSize(graph.nodeWeights.begin(),
                                     graph.nodeWeights.end());
    std::vector<bool> startsBuffer(graph.size, false);
    for (const Scheduler::LiveIn &liveIn : liveIns) {
        for (uint64_t node : liveIn.tiedUsers) {
            if (graph.nodeWeights[node] == 0) {
                startsBuffer[node] = true;
                bufferSize[node] = liveIn.size;
            }
        }
    }

    // The buffer of every value: a tied op writes into the input whose edge
    // weight equals its result size, preferring the last-defined one (the
    // init operand)
    std::vector<uint64_t> root(graph.size, UINT64_MAX);
    for (uint64_t node : order.topological) {
        if (graph.nodeWeights[node] > 0 || startsBuffer[node]) {
            root[node] = node;
            continue;
        }
        uint64_t resultSize = 0;
        for (const auto &[consumer, edgeWeight] : graph.adj[node])
            resultSize = std::max(resultSize, edgeWeight);
        uint64_t tiedInput = UINT64_MAX;
        for (const auto &[input, edgeWeight] : graph.revAdj[node]) {
            if (resultSize > 0 && edgeWeight == resultSize &&
                (tiedInput == UINT64_MAX ||
                 order.global[input] > order.global[tiedInput]))
                tiedInput = input;
        }
        if (tiedInput != UINT64_MAX)
            root[node] = root[tiedInput];
    }
    std::vector<std::vector<uint64_t>> aliases(graph.size);
    for (uint64_t node = 0; node < graph.size; ++node) {
        if (root[node] != UINT64_MAX)
            aliases[root[node]].push_back(node);
    }

    std::vector<Buffer> buffers;
    for (uint64_t producer = 0; producer < graph.size; ++producer) {
        if (aliases[producer].empty())
            continue;
        uint64_t partition = partitionMapping[producer];
        uint64_t birth = order.local[producer];
        Buffer buffer{producer, partition, bufferSize[producer],
                      birth,    birth,     {},
                      {},       false};
        auto use = [&](uint64_t node) {
            uint64_t local = order.local[node];
            auto [it, inserted] =
                buffer.usedAt.try_emplace(partitionMapping[node], local, local);
            it->second.first = std::min(it->second.first, local);
            it->second.second = std::max(it->second.second, local);
            if (partitionMapping[node] == partition)
                buffer.death = std::max(buffer.death, local);
            else
                buffer.users.push_back(partitionMapping[node]);
        };
        for (uint64_t value : aliases[producer]) {
            use(value);
            if (liveOutNodes.empty() ? graph.adj[value].empty()
                                     : liveOutNodes[value])
                buffer.liveOut = true;
            for (const auto &[consumer, edgeWeight] : graph.adj[value])
                use(consumer);
        }
        std::sort(buffer.users.begin(), buffer.users.end());
        buffer.users.erase(
            std::unique(buffer.users.begin(), buffer.users.end()),
            buffer.users.end());
        buffers.push_back(std::move(buffer));
    }
    return buffers;
}

/**
 * @brief Per partition, the buffers that never leave it
 */
std::vector<std::vector<utils::Tensor>>
localTensors(const std::vector<Buffer> &buffers, uint64_t numPartitions) {
    std::vector<std::vector<utils::Tensor>> tensors(numPartitions);
    for (const Buffer &buffer : buffers) {
        if (!buffer.escapes())
            tensors[buffer.partition].push_back(
                {buffer.size, buffer.birth, buffer.death, 0});
    }
    return tensors;
}

/**
 * @brief Peak memory of every partition: its local tensors, packed
 */
std::vector<uint64_t>
partitionWeights(std::vector<std::vector<utils::Tensor>> tensors,
                 const utils::SlabPacker &packer) {
    std::vector<uint64_t> weights(tensors.size(), 0);
    for (uint64_t partition = 0; partition < tensors.size(); ++partition) {
        if (tensors[partition].empty())
            continue;
        utils::sortByLifetime(tensors[partition]);
        weights[partition] = packer(tensors[partition]);
    }
    return weights;
}

/**
 * @brief The buffers that escape their partition, as tensors of the coarse
 * graph, live until their last consuming partition
 */
std::vector<hypergraph::Tensor>
escapingTensors(const std::vector<Buffer> &buffers) {
    std::vector<hypergraph::Tensor> tensors;
    for (const Buffer &buffer : buffers) {
        if (buffer.escapes())
            tensors.push_back(
                {buffer.partition, buffer.size, buffer.users, buffer.liveOut});
    }
    return tensors;
}

/**
 * @brief The footprints of the partitions
 */
std::unique_ptr<PartitionFootprints>
makeFootprints(const std::vector<Buffer> &buffers,
               std::vector<std::vector<utils::Tensor>> localTensors,
               const utils::SlabPacker &packer,
               std::vector<hypergraph::LiveIn> liveIns) {
    std::vector<PartitionFootprints::Escaping> escaping;
    std::vector<std::vector<PartitionFootprints::Use>> uses(
        localTensors.size());
    for (const Buffer &buffer : buffers) {
        if (!buffer.escapes() || buffer.size == 0)
            continue;
        std::vector<uint64_t> members;
        for (const auto &[partition, range] : buffer.usedAt) {
            uses[partition].push_back(
                {escaping.size(), range.first, range.second});
            members.push_back(partition);
        }
        escaping.push_back({buffer.size, std::move(members), buffer.liveOut});
    }
    return std::make_unique<PartitionFootprints>(
        std::move(localTensors), std::move(escaping), std::move(uses), packer,
        std::move(liveIns));
}

/**
 * @brief The live-ins, with the partitions of their users
 */
std::vector<hypergraph::LiveIn>
partitionLiveIns(const std::vector<Scheduler::LiveIn> &liveIns,
                 const std::vector<uint64_t> &partitionMapping) {
    std::vector<hypergraph::LiveIn> result;
    for (const Scheduler::LiveIn &liveIn : liveIns) {
        std::set<uint64_t> partitions;
        for (uint64_t user : liveIn.users)
            partitions.insert(partitionMapping[user]);
        result.push_back(
            {liveIn.size,
             std::vector<uint64_t>(partitions.begin(), partitions.end())});
    }
    return result;
}

/**
 * @brief Checks a schedule from CP-SAT against brute force: it must be an
 * order of all nodes with the reported cost, that cost must be the optimum,
 * and under Objective::TransferMemory its peak must meet the limit, or the
 * lowest peak above it. Dumps both orders and throws if not
 */
void verifySchedule(const hypergraph::HyperGraph &graph,
                    const PartitionFootprints &footprints,
                    const std::vector<uint64_t> &schedule, uint64_t cost,
                    Objective objective, uint64_t memoryLimit) {
    auto [bfSchedule, bfCost] =
        bruteforce::BruteForceSolver(graph, objective, &footprints, memoryLimit)
            .solve();

    uint64_t n = graph.getSize();
    std::vector<uint64_t> pos(n, n);
    for (uint64_t step = 0; step < schedule.size(); step++) {
        if (schedule[step] < n)
            pos[schedule[step]] = step;
    }
    bool valid = schedule.size() == n &&
                 std::find(pos.begin(), pos.end(), n) == pos.end();
    for (uint64_t t = 0; valid && t < graph.getNumTensors(); t++) {
        for (uint64_t user : graph.getTensorUsers()[t])
            valid &= pos[graph.getTensorProducers()[t]] < pos[user];
    }
    valid = valid && cost == (objective == Objective::Memory
                                  ? graph.computePeakMemory(schedule)
                                  : graph.computeTransfers(schedule));
    bool peakOk = objective == Objective::Memory || memoryLimit == UINT64_MAX ||
                  footprints.peak(schedule) <=
                      std::max(memoryLimit, footprints.peak(bfSchedule));
    if (valid && cost == bfCost && peakOk)
        return;

    auto print = [](const std::vector<uint64_t> &order) {
        std::ostringstream text;
        for (uint64_t node : order)
            text << " " << node;
        return text.str();
    };
    std::cerr << "Verification failed on:\n"
              << graph << "\nCP-SAT order (cost " << cost
              << "):" << print(schedule) << "\nBrute-force order (cost "
              << bfCost << "):" << print(bfSchedule) << std::endl;
    throw std::runtime_error(
        "Verification failed: CP-SAT cost (" + std::to_string(cost) +
        ") != brute-force cost (" + std::to_string(bfCost) + ")" +
        (valid ? "" : ", or not the cost of a valid order") +
        (peakOk ? "" : ", or over the limit"));
}

} // namespace

Scheduler::Scheduler(const core::Graph &originalGraph,
                     const std::vector<uint64_t> &partitionMapping, bool debug,
                     bool verify, uint64_t offsetAlignment,
                     uint64_t rangeAlignment, utils::SlabPacker packer,
                     std::vector<bool> liveOutNodes,
                     std::vector<LiveIn> liveIns)
    : originalGraph(originalGraph), partitionMapping(partitionMapping),
      debug(debug), verify(verify), packer(std::move(packer)),
      liveOutNodes(std::move(liveOutNodes)), liveIns(std::move(liveIns)) {
    if (!this->liveOutNodes.empty() &&
        this->liveOutNodes.size() != originalGraph.size)
        throw std::invalid_argument("liveOutNodes needs one entry per node");
    for (const LiveIn &liveIn : this->liveIns) {
        for (const auto *nodes : {&liveIn.users, &liveIn.tiedUsers}) {
            for (uint64_t node : *nodes) {
                if (node >= originalGraph.size)
                    throw std::invalid_argument("live-in user out of range");
            }
        }
    }
    if (!this->packer) {
        this->packer = [offsetAlignment,
                        rangeAlignment](std::vector<utils::Tensor> &tensors) {
            return utils::packBestFit(tensors, offsetAlignment, rangeAlignment);
        };
    }
}

void Scheduler::buildCoarseGraph() {
    // Partition ids are contiguous from 0
    uint64_t numPartitions =
        *std::max_element(partitionMapping.begin(), partitionMapping.end()) + 1;
    assert(numPartitions != 0 && "Number of partitions must be > 0");

    std::vector<Buffer> buffers = collectBuffers(
        originalGraph, partitionMapping, numPartitions, liveOutNodes, liveIns);
    std::vector<hypergraph::LiveIn> liveInsOfPartitions =
        partitionLiveIns(liveIns, partitionMapping);
    std::vector<std::vector<utils::Tensor>> local =
        localTensors(buffers, numPartitions);
    footprints = makeFootprints(buffers, local, packer, liveInsOfPartitions);
    std::vector<hypergraph::Tensor> tensors = escapingTensors(buffers);

    // Partition edges that no tensor orders become empty tensors (e.g. a
    // buffer allocated in A, written in place in B and read in C orders A
    // before B and C, but not B before C)
    std::set<std::pair<uint64_t, uint64_t>> ordered;
    for (const auto &tensor : tensors) {
        for (uint64_t user : tensor.users) {
            ordered.insert({tensor.producer, user});
        }
    }
    for (uint64_t from = 0; from < originalGraph.size; ++from) {
        for (const auto &[to, edgeWeight] : originalGraph.adj[from]) {
            uint64_t fromPartition = partitionMapping[from];
            uint64_t toPartition = partitionMapping[to];

            if (fromPartition != toPartition &&
                ordered.insert({fromPartition, toPartition}).second) {
                tensors.push_back({fromPartition, 0, {toPartition}, false});
            }
        }
    }

    coarseGraph = std::make_unique<hypergraph::HyperGraph>(
        partitionWeights(std::move(local), packer), tensors,
        std::move(liveInsOfPartitions));
}

std::pair<std::vector<uint64_t>, uint64_t>
Scheduler::run(uint64_t timeLimitSeconds, int32_t numWorkers, Solver solver,
               Objective objective, uint64_t memoryLimit, bool deterministic) {
    using Clock = std::chrono::steady_clock;
    auto secondsSince = [](Clock::time_point start) {
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    auto start = Clock::now();
    buildCoarseGraph();
    if (debug) {
        std::cerr << "Coarse graph construction time: " << secondsSince(start)
                  << " seconds" << std::endl;
    }
    const hypergraph::HyperGraph &graph = *coarseGraph;

    std::vector<uint64_t> schedule;
    uint64_t cost;
    if (solver == Solver::RPO) {
        start = Clock::now();
        schedule = heuristic::RPOScheduler(graph, debug).solve().first;
        cost = objective == Objective::Memory
                   ? graph.computePeakMemory(schedule)
                   : graph.computeTransfers(schedule);
        if (debug) {
            std::cerr << "RPO solve time: " << secondsSince(start)
                      << " seconds\n"
                      << (objective == Objective::Memory
                              ? "Peak memory: "
                              : "Transferred bytes: ")
                      << cost << std::endl;
        }
        return {schedule, cost};
    }

    if (objective == Objective::Memory) {
        std::tie(std::ignore, schedule, cost) =
            cpsat::CPSATMemorySolver(graph, debug)
                .solve(timeLimitSeconds, numWorkers, deterministic);
    } else {
        std::tie(std::ignore, schedule, cost) =
            cpsat::CPSATTransferMemorySolver(graph, *footprints, memoryLimit,
                                             debug)
                .solve(timeLimitSeconds, numWorkers, deterministic);
    }
    if (verify)
        verifySchedule(graph, *footprints, schedule, cost, objective,
                       memoryLimit);
    return {schedule, cost};
}

uint64_t Scheduler::peakFootprint(const std::vector<uint64_t> &schedule) const {
    return footprints->peak(schedule);
}

} // namespace scheduling

} // namespace dag_partitioning
