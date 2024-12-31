//
// Created by panagiotis on 23/12/2024.
//

#include "BoundaryFM.h"

BoundaryFM::BoundaryFM(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t &initialEdgeCut,
                       uint64_t maxNumberOfPasses,
                       double upperBoundPartWeight, double lowerBoundPartWeight)
        : Refinement(graph, initialBisectionInfo, initialEdgeCut, maxNumberOfPasses, upperBoundPartWeight,
                     lowerBoundPartWeight) {}

void BoundaryFM::insertMovableNodesIntoHeaps(std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                             std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                             std::vector<bool> &inHeap) const {
    for (uint64_t nodeId = 0; nodeId < workingGraph.size; ++nodeId) {
        if (inHeap[nodeId]) continue;
        if (!initialBisectionInfo[nodeId]) {
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (!initialBisectionInfo[neighborId]) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                heapV0.emplace(gain, nodeId);
                inHeap[nodeId] = true;
            }
        } else {
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.revAdj[nodeId]) {
                if (initialBisectionInfo[neighborId]) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                heapV1.emplace(gain, nodeId);
                inHeap[nodeId] = true;
            }
        }
    }
}

void BoundaryFM::insertMovableNeighborsIntoHeaps(const std::vector<bool> &initialBisectionInfo,
                                                 std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
                                                 std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
                                                 std::vector<bool> &inHeap, uint64_t movedNodeId,
                                                 bool checkOutNeighbors) const {
    if (checkOutNeighbors) {
        for (const auto &[nodeId, _]: workingGraph.adj[movedNodeId]) {
            if (inHeap[nodeId] || !initialBisectionInfo[nodeId]) continue;

            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.revAdj[nodeId]) {
                if (initialBisectionInfo[neighborId]) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                heapV1.emplace(gain, nodeId);
                inHeap[nodeId] = true;
            }
        }
    } else {
        for (const auto &[nodeId, _]: workingGraph.revAdj[movedNodeId]) {
            if (inHeap[nodeId] || initialBisectionInfo[nodeId]) continue;

            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (!initialBisectionInfo[neighborId]) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                heapV0.emplace(gain, nodeId);
                inHeap[nodeId] = true;
            }
        }
    }
}

bool BoundaryFM::isBalanceImprovedOrMaintained(bool moveFromV0, uint64_t movedNodeId, uint64_t maxNodeWeight,
                                               uint64_t sizeV0, uint64_t sizeV1, bool isBalanced) {
    // If the partition is already balanced, the move must maintain the balance
    // If not, the balance will be improved since we are moving from the maximum loaded part
    // But it must not become reversely unbalanced
    uint64_t newSizeV0, newSizeV1;
    if (moveFromV0) {
        newSizeV0 = sizeV0 - workingGraph.nodeWeights[movedNodeId];
        newSizeV1 = sizeV1 + workingGraph.nodeWeights[movedNodeId];
    } else {
        newSizeV0 = sizeV0 + workingGraph.nodeWeights[movedNodeId];
        newSizeV1 = sizeV1 - workingGraph.nodeWeights[movedNodeId];
    }
    if (isBalanced) {
        if ((double) newSizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
            (double) newSizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
            (double) newSizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
            (double) newSizeV1 > upperBoundPartWeight + (double) maxNodeWeight)
            return false;
        return true;
    } else {
        uint64_t newSizeOfPreviousBigger = sizeV0 > sizeV1 ? newSizeV0 : newSizeV1;
        uint64_t newSizeOfPreviousSmaller = sizeV0 > sizeV1 ? newSizeV1 : newSizeV0;
        if ((double) newSizeOfPreviousBigger < lowerBoundPartWeight - (double) maxNodeWeight ||
            (double) newSizeOfPreviousSmaller > upperBoundPartWeight + (double) maxNodeWeight)
            return false;
        return true;
    }
}

bool BoundaryFM::onePassRefinement() {
    std::vector<bool> moved(workingGraph.size, false);
    std::vector<bool> inHeap(workingGraph.size, false);

    std::priority_queue<std::pair<int64_t, uint64_t>> heapV0, heapV1;

    insertMovableNodesIntoHeaps(heapV0, heapV1, inHeap);

    std::vector<std::pair<uint64_t, bool>> moveSequence;
    uint64_t currentEdgeCut = initialEdgeCut;
    uint64_t bestEdgeCut = initialEdgeCut;
    size_t bestMovePrefix = 0;
    std::vector<bool> initialBisectionInfoTemp = initialBisectionInfo;

    uint64_t maxNodeWeight = workingGraph.maxNodeWeight();
    auto [sizeV0, sizeV1] = calculatePartSizes();
    bool isBalanced = checkBalance(maxNodeWeight);
    uint64_t minMaxPartSize = std::max(sizeV0, sizeV1);

    while (!heapV0.empty() || !heapV1.empty()) {
        bool moveFromV0 = false;
        uint64_t nodeId;
        int64_t gain;
        bool isV0Larger = sizeV0 > sizeV1;

        if (isBalanced) {
            if (heapV0.empty()) {
                std::tie(gain, nodeId) = heapV1.top();
                heapV1.pop();
            } else if (heapV1.empty()) {
                std::tie(gain, nodeId) = heapV0.top();
                heapV0.pop();
                moveFromV0 = true;
            } else {
                auto bestV0 = heapV0.top();
                auto bestV1 = heapV1.top();
                if (bestV0.first >= bestV1.first) {
                    std::tie(gain, nodeId) = bestV0;
                    heapV0.pop();
                    moveFromV0 = true;
                } else {
                    std::tie(gain, nodeId) = bestV1;
                    heapV1.pop();
                }
            }
        } else {
            if (isV0Larger) {
                assert(!heapV0.empty() && "Partition is unbalanced but heavy part V0 has no movable nodes");
                std::tie(gain, nodeId) = heapV0.top();
                heapV0.pop();
                moveFromV0 = true;
            } else {
                assert(!heapV1.empty() && "Partition is unbalanced but heavy part V1 has no movable nodes");
                std::tie(gain, nodeId) = heapV1.top();
                heapV1.pop();
                moveFromV0 = false;
            }
        }

        if (moved[nodeId])
            continue;
        if (!isBalanceImprovedOrMaintained(moveFromV0, nodeId, maxNodeWeight, sizeV0, sizeV1, isBalanced))
            continue;

        if (moveFromV0) {
            sizeV0 -= workingGraph.nodeWeights[nodeId];
            sizeV1 += workingGraph.nodeWeights[nodeId];
        } else {
            sizeV0 += workingGraph.nodeWeights[nodeId];
            sizeV1 -= workingGraph.nodeWeights[nodeId];
        }
        moved[nodeId] = true;
        moveSequence.emplace_back(nodeId, moveFromV0);
        currentEdgeCut -= gain;
        initialBisectionInfoTemp[nodeId] = moveFromV0;
        uint64_t maxPartSize = std::max(sizeV0, sizeV1);

        if (!isBalanced || currentEdgeCut < bestEdgeCut || maxPartSize < minMaxPartSize) {
            bestEdgeCut = currentEdgeCut;
            bestMovePrefix = moveSequence.size();
            minMaxPartSize = maxPartSize;
        }

        isBalanced = true;
        if ((double) sizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
            (double) sizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
            (double) sizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
            (double) sizeV1 > upperBoundPartWeight + (double) maxNodeWeight)
            isBalanced = false;

        if (moveFromV0) {
            // Moving from V0 to V1 can make V1 nodes unmovable
            // If an out-neighbor of the moved nodeId is in V1 and is not yet moved,
            // (if it has already been moved, we don't care if it is movable or not)
            // Then the neighbor becomes unmovable since it has one in-neighbor (the moved nodeId) that is in V1
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (!moved[neighborId] && initialBisectionInfo[neighborId] && inHeap[neighborId]) {
                    moved[neighborId] = true;
                }
            }
            // Moving from V0 to V1 can make V0 nodes movable
            // If an in-neighbor of the moved nodeId is in V0 and is not already movable,
            // (if it has already been marked movable by the initial check, we don't recheck)
            // Then the neighbor might become movable if its only out-neighbor in V0 was the moved node
            insertMovableNeighborsIntoHeaps(initialBisectionInfoTemp, heapV0, heapV1, inHeap, nodeId, false);
        } else {
            // Moving from V1 to V0 can make V0 nodes unmovable
            // If an in-neighbor of the moved nodeId is in V0 and is not yet moved,
            // (if it has already been moved, we don't care if it is movable or not)
            // Then the neighbor becomes unmovable since it has one out-neighbor (the moved nodeId) that is in V0
            for (const auto &[neighborId, _]: workingGraph.revAdj[nodeId]) {
                if (!moved[neighborId] && !initialBisectionInfo[neighborId] && inHeap[neighborId]) {
                    moved[neighborId] = true;
                }
            }
            // Moving from V1 to V0 can make V1 nodes movable
            // If an out-neighbor of the moved nodeId is in V1 and is not already movable,
            // (if it has already been marked movable by the initial check, we don't recheck)
            // Then the neighbor might become movable if its only in-neighbor in V1 was the moved node
            insertMovableNeighborsIntoHeaps(initialBisectionInfoTemp, heapV0, heapV1, inHeap, nodeId, true);
        }
    }

    if (bestMovePrefix == 0) return false;
    auto [initialSizeV0, initialSizeV1] = calculatePartSizes();
    if (initialEdgeCut < bestEdgeCut && std::max(sizeV0, sizeV1) > std::max(initialSizeV0, initialSizeV1)) return false;

    for (size_t i = 0; i < bestMovePrefix; ++i) {
        auto [nodeId, fromV0] = moveSequence[i];
        initialBisectionInfo[nodeId] = fromV0;
    }
    initialEdgeCut = bestEdgeCut;

    return true;
}
