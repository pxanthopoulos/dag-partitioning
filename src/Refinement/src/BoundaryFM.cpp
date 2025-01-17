/**
 * @file BoundaryFM.cpp
 * @brief Implementation of boundary FM refinement for directed graphs
 */

#include "BoundaryFM.h"

BoundaryFM::BoundaryFM(const Graph &graph, std::vector<uint8_t> &initialBisectionInfo,
                       uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
                       double upperBoundPartWeight, double lowerBoundPartWeight)
        : Refinement(graph, initialBisectionInfo, initialEdgeCut, maxNumberOfPasses,
                     upperBoundPartWeight, lowerBoundPartWeight) {}

void BoundaryFM::insertMovableNodesIntoHeaps(
        std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
        std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
        std::vector<uint8_t> &inHeap) const {

    // Check each vertex for movability
    for (uint64_t nodeId = 0; nodeId < workingGraph.size; ++nodeId) {
        // If the node has already been in the heap previously, don't re-add it
        if (inHeap[nodeId] == 1) continue;

        if (initialBisectionInfo[nodeId] == 0) {  // Node in V0
            // Check if movable to V1: all out-neighbors must be in V1
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (initialBisectionInfo[neighborId] == 0) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain according to equation 4.1
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                heapV0.emplace(gain, nodeId);
                inHeap[nodeId] = 1;
            }
        } else {  // Node in V1
            // Check if movable to V0: all in-neighbors must be in V0
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.revAdj[nodeId]) {
                if (initialBisectionInfo[neighborId] == 1) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain (negative of V0->V1 gain)
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                heapV1.emplace(gain, nodeId);
                inHeap[nodeId] = 1;
            }
        }
    }
}

// Called after moving a node to update movable status of its neighbors (possibly insert new movable nodes)
void BoundaryFM::insertMovableNeighborsIntoHeaps(
        const std::vector<uint8_t> &currentBisectionInfo,
        std::priority_queue<std::pair<int64_t, uint64_t>> &heapV0,
        std::priority_queue<std::pair<int64_t, uint64_t>> &heapV1,
        std::vector<uint8_t> &inHeap, uint64_t movedNodeId,
        bool checkOutNeighbors) const {

    if (checkOutNeighbors) {  // Node moved V1->V0, check out-neighbors
        for (const auto &[nodeId, _]: workingGraph.adj[movedNodeId]) {
            // If the node has already been in the heap previously, don't re-add it
            if (inHeap[nodeId] == 1 || currentBisectionInfo[nodeId] == 0) continue;

            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.revAdj[nodeId]) {
                if (currentBisectionInfo[neighborId] == 1) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain for potential move
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                heapV1.emplace(gain, nodeId);
                inHeap[nodeId] = 1;
            }
        }
    } else {  // Node moved V0->V1, check in-neighbors
        for (const auto &[nodeId, _]: workingGraph.revAdj[movedNodeId]) {
            // If the node has already been in the heap previously, don't re-add it
            if (inHeap[nodeId] == 1 || currentBisectionInfo[nodeId] == 1) continue;

            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (currentBisectionInfo[neighborId] == 0) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain for potential move
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[nodeId]) {
                    gain += (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[nodeId]) {
                    gain -= (int64_t) edgeWeight;
                }
                heapV0.emplace(gain, nodeId);
                inHeap[nodeId] = 1;
            }
        }
    }
}

bool BoundaryFM::confirmMove(
        uint8_t moveFromV0, uint64_t movedNodeId, uint64_t maxNodeWeight,
        uint64_t sizeV0, uint64_t sizeV1) const {

    // Calculate new partition sizes after proposed move
    uint64_t newSizeV0 = moveFromV0 == 1 ?
                         sizeV0 - workingGraph.nodeWeights[movedNodeId] :
                         sizeV0 + workingGraph.nodeWeights[movedNodeId];
    uint64_t newSizeV1 = moveFromV0 == 1 ?
                         sizeV1 + workingGraph.nodeWeights[movedNodeId] :
                         sizeV1 - workingGraph.nodeWeights[movedNodeId];

    // Moving from larger to smaller partition must not cause reverse imbalance
    uint64_t newSizeOfPreviousBigger = sizeV0 > sizeV1 ? newSizeV0 : newSizeV1;
    uint64_t newSizeOfPreviousSmaller = sizeV0 > sizeV1 ? newSizeV1 : newSizeV0;
    if ((double) newSizeOfPreviousBigger < lowerBoundPartWeight - (double) maxNodeWeight ||
        (double) newSizeOfPreviousSmaller > upperBoundPartWeight + (double) maxNodeWeight)
        return false;
    return true;

}

bool BoundaryFM::onePassRefinement() {
    std::vector<uint8_t> moved(workingGraph.size, 0);
    std::vector<uint8_t> inHeap(workingGraph.size, 0);

    std::priority_queue<std::pair<int64_t, uint64_t>> heapV0, heapV1;

    // Find initially movable vertices
    insertMovableNodesIntoHeaps(heapV0, heapV1, inHeap);

    std::vector<std::pair<uint64_t, bool>> moveSequence;
    uint64_t currentEdgeCut = initialEdgeCut;
    uint64_t bestEdgeCut = initialEdgeCut;
    uint64_t bestMovePrefix = 0;
    std::vector<uint8_t> initialBisectionInfoTemp = initialBisectionInfo;

    uint64_t maxNodeWeight = workingGraph.maxNodeWeight;
    auto [sizeV0, sizeV1] = calculatePartSizes();
    bool isBalanced = checkBalance(maxNodeWeight);
    uint64_t minMaxPartSize = std::max(sizeV0, sizeV1);

    // Main refinement loop - make moves until no more vertices can move
    while (!heapV0.empty() || !heapV1.empty()) {
        uint8_t moveFromV0 = 0;
        uint64_t nodeId;
        int64_t gain;
        bool isV0Larger = sizeV0 > sizeV1;

        if (isV0Larger) {
            if (heapV0.empty()) break;
            std::tie(gain, nodeId) = heapV0.top();
            heapV0.pop();
            moveFromV0 = 1;
        } else {
            if (heapV1.empty()) break;
            std::tie(gain, nodeId) = heapV1.top();
            heapV1.pop();
        }

        // If node already moved (or marked move to signify unmovable node), skip
        if (moved[nodeId] == 1)
            continue;
        // If move causes unbalance or does not improve balance, skip
        if (!confirmMove(moveFromV0, nodeId, maxNodeWeight, sizeV0, sizeV1))
            continue;

        // Tentative move
        if (moveFromV0 == 1) {
            sizeV0 -= workingGraph.nodeWeights[nodeId];
            sizeV1 += workingGraph.nodeWeights[nodeId];
        } else {
            sizeV0 += workingGraph.nodeWeights[nodeId];
            sizeV1 -= workingGraph.nodeWeights[nodeId];
        }
        moved[nodeId] = 1;
        moveSequence.emplace_back(nodeId, moveFromV0);
        currentEdgeCut -= gain;
        initialBisectionInfoTemp[nodeId] = moveFromV0;
        uint64_t maxPartSize = std::max(sizeV0, sizeV1);

        // If partition is unbalanced, always include this move to the sequence
        // If balanced, include move if it improves edge cut or improves balance quality (without worsening edge cut)
        if (!isBalanced || currentEdgeCut < bestEdgeCut ||
            (currentEdgeCut == bestEdgeCut && maxPartSize < minMaxPartSize)) {
            bestEdgeCut = currentEdgeCut;
            bestMovePrefix = moveSequence.size();
            minMaxPartSize = maxPartSize;
        }

        isBalanced = true;
        if ((double) sizeV0 < lowerBoundPartWeight ||
            (double) sizeV0 > upperBoundPartWeight ||
            (double) sizeV1 < lowerBoundPartWeight ||
            (double) sizeV1 > upperBoundPartWeight)
            isBalanced = false;

        if (moveFromV0 == 1) {
            // Moving from V0 to V1 can make V1 nodes unmovable
            // If an out-neighbor of the moved nodeId is in V1 and is not yet moved,
            // (if it has already been moved, we don't care if it is movable or not)
            // Then the neighbor becomes unmovable since it has one in-neighbor (the moved nodeId) that is in V1
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (moved[neighborId] == 0 && initialBisectionInfoTemp[neighborId] == 1 && inHeap[neighborId] == 1) {
                    moved[neighborId] = 1;
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
                if (moved[neighborId] == 0 && initialBisectionInfoTemp[neighborId] == 0 && inHeap[neighborId] == 1) {
                    moved[neighborId] = 1;
                }
            }
            // Moving from V1 to V0 can make V1 nodes movable
            // If an out-neighbor of the moved nodeId is in V1 and is not already movable,
            // (if it has already been marked movable by the initial check, we don't recheck)
            // Then the neighbor might become movable if its only in-neighbor in V1 was the moved node
            insertMovableNeighborsIntoHeaps(initialBisectionInfoTemp, heapV0, heapV1, inHeap, nodeId, true);
        }
    }

    // If no improvement was found, return false
    if (bestMovePrefix == 0) return false;

    // Check if should keep current solution
    auto [initialSizeV0, initialSizeV1] = calculatePartSizes();
    if (initialEdgeCut <= bestEdgeCut && std::max(sizeV0, sizeV1) >= std::max(initialSizeV0, initialSizeV1))
        return false;

    // Apply best move sequence
    for (uint64_t i = 0; i < bestMovePrefix; ++i) {
        auto [nodeId, fromV0] = moveSequence[i];
        initialBisectionInfo[nodeId] = fromV0;
    }
    initialEdgeCut = bestEdgeCut;

    return true;
}
