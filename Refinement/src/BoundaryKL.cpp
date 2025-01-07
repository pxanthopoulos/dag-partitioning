/**
 * @file BoundaryKL.cpp
 * @brief Implementation of boundary KL refinement for directed graphs
 */

#include "BoundaryKL.h"

BoundaryKL::BoundaryKL(const Graph &graph, std::vector<bool> &initialBisectionInfo,
                       uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
                       double upperBoundPartWeight, double lowerBoundPartWeight)
        : Refinement(graph, initialBisectionInfo, initialEdgeCut, maxNumberOfPasses,
                     upperBoundPartWeight, lowerBoundPartWeight) {}

void BoundaryKL::insertMovableNodesIntoLists(std::list<std::pair<int64_t, uint64_t>> &listV0,
                                             std::list<std::pair<int64_t, uint64_t>> &listV1,
                                             std::vector<bool> &inList) const {

    // Check each vertex for movability
    for (uint64_t nodeId = 0; nodeId < workingGraph.size; ++nodeId) {
        // If the node has already been in a vec previously, don't re-add it
        if (inList[nodeId]) continue;

        if (!initialBisectionInfo[nodeId]) {  // Node in V0
            // Check if movable to V1: all out-neighbors must be in V1
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (!initialBisectionInfo[neighborId]) {
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
                listV0.emplace_back(gain, nodeId);
                inList[nodeId] = true;
            }
        } else {  // Node in V1
            // Check if movable to V0: all in-neighbors must be in V0
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.revAdj[nodeId]) {
                if (initialBisectionInfo[neighborId]) {
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
                listV1.emplace_back(gain, nodeId);
                inList[nodeId] = true;
            }
        }
    }
}

// Called after moving a node to update movable status of its neighbors (possibly insert new movable nodes)
void BoundaryKL::insertMovableNeighborsIntoLists(const std::vector<bool> &currentBisectionInfo,
                                                 std::list<std::pair<int64_t, uint64_t>> &listV0,
                                                 std::list<std::pair<int64_t, uint64_t>> &listV1,
                                                 std::vector<bool> &inList, uint64_t movedNodeId,
                                                 bool checkOutNeighbors) const {
    if (checkOutNeighbors) {  // Node moved V0->V1, check out-neighbors
        for (const auto &[nodeId, _]: workingGraph.adj[movedNodeId]) {
            // If the node has already been in the heap previously, don't re-add it
            if (inList[nodeId] || !currentBisectionInfo[nodeId]) continue;

            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.revAdj[nodeId]) {
                if (currentBisectionInfo[neighborId]) {
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
                auto it = std::upper_bound(listV0.begin(), listV0.end(),
                                           std::make_pair(gain, nodeId),
                                           [](const auto &a, const auto &b) { return a.first > b.first; });

                listV0.insert(it, std::make_pair(gain, nodeId));
                inList[nodeId] = true;
            }
        }
    } else {  // Node moved V1->V0, check in-neighbors
        for (const auto &[nodeId, _]: workingGraph.revAdj[movedNodeId]) {
            // If the node has already been in the heap previously, don't re-add it
            if (inList[nodeId] || currentBisectionInfo[nodeId]) continue;

            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.adj[nodeId]) {
                if (!currentBisectionInfo[neighborId]) {
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
                auto it = std::upper_bound(listV1.begin(), listV1.end(),
                                           std::make_pair(gain, nodeId),
                                           [](const auto &a, const auto &b) { return a.first > b.first; });

                listV1.insert(it, std::make_pair(gain, nodeId));
                inList[nodeId] = true;
            }
        }
    }
}

std::tuple<bool, uint64_t, uint64_t, uint64_t>
BoundaryKL::findBestMovablePairBalanced(std::list<std::pair<int64_t, uint64_t>> &listV0,
                                        std::list<std::pair<int64_t, uint64_t>> &listV1,
                                        std::vector<bool> &moved, uint64_t &sizeV0,
                                        uint64_t &sizeV1) const {
    bool found = false;
    int64_t maxTotalGain = INT64_MIN;
    uint64_t bestNodeV0 = -1;
    uint64_t bestNodeV1 = -1;
    std::list<std::pair<int64_t, uint64_t>>::iterator bestItV0;
    std::list<std::pair<int64_t, uint64_t>>::iterator bestItV1;

    for (auto itV0 = listV0.begin(); itV0 != listV0.end(); ++itV0) {
        const auto &[gainV0, nodeIdV0] = *itV0;
        // Check if the node has been marked as unmovable
        if (moved[nodeIdV0]) continue;

        for (auto itV1 = listV1.begin(); itV1 != listV1.end(); ++itV1) {
            const auto &[gainV1, nodeIdV1] = *itV1;
            uint64_t nodeIdV1tmp = nodeIdV1;

            // Early exit. Before checking if the pair is swappable, we can check if the gain is better than the best so far.
            // If not, we can continue to the next pair.
            // Even better, we can break the inner loop as the gain will only decrease and go to the next outer loop iteration
            int64_t totalGain = gainV0 + gainV1;
            if (totalGain <= maxTotalGain) break;

            // Check if the node has been marked as unmovable
            if (moved[nodeIdV1]) continue;

            // If this move breaks the balance, skip it
            uint64_t newSizeV0 = sizeV0 - workingGraph.nodeWeights[nodeIdV0] + workingGraph.nodeWeights[nodeIdV1];
            uint64_t newSizeV1 = sizeV1 - workingGraph.nodeWeights[nodeIdV1] + workingGraph.nodeWeights[nodeIdV0];
            uint64_t maxNodeWeight = workingGraph.maxNodeWeight();
            if ((double) newSizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
                (double) newSizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
                (double) newSizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
                (double) newSizeV1 > upperBoundPartWeight + (double) maxNodeWeight)
                continue;

            // Check if nodeIdV0 has edge to nodeIdV1. If so, they cannot be swapped.
            auto it = std::find_if(workingGraph.adj[nodeIdV0].begin(), workingGraph.adj[nodeIdV0].end(),
                                   [nodeIdV1tmp](const auto &pair) { return pair.first == nodeIdV1tmp; });
            if (it != workingGraph.adj[nodeIdV0].end()) {
                continue;
            }

            // Pair is swappable AND improves gain
            found = true;
            bestNodeV0 = nodeIdV0;
            bestNodeV1 = nodeIdV1;
            maxTotalGain = totalGain;
            bestItV0 = itV0;
            bestItV1 = itV1;
        }
    }

    if (found) {
        listV0.erase(bestItV0);
        listV1.erase(bestItV1);
    }
    return {found, bestNodeV0, bestNodeV1, maxTotalGain};
}

std::tuple<bool, uint64_t, uint64_t, uint64_t>
BoundaryKL::findBestMovablePairUnbalanced(std::list<std::pair<int64_t, uint64_t>> &listV0,
                                          std::list<std::pair<int64_t, uint64_t>> &listV1,
                                          std::vector<bool> &moved, uint64_t &sizeV0,
                                          uint64_t &sizeV1) const {
    bool found = false;
    // Perfect balance is each size is 50% of the total. Here we calculate the difference from the perfect balance. Lower is better.
    double diffFromPerfectBalance = std::abs(50.0 - (((double) sizeV0 * 100.0) / ((double) (sizeV0 + sizeV1))));
    uint64_t bestNodeV0 = -1;
    uint64_t bestNodeV1 = -1;
    uint64_t gain = 0;
    std::list<std::pair<int64_t, uint64_t>>::iterator bestItV0;
    std::list<std::pair<int64_t, uint64_t>>::iterator bestItV1;

    for (auto itV0 = listV0.begin(); itV0 != listV0.end(); ++itV0) {
        const auto &[gainV0, nodeIdV0] = *itV0;
        // Check if the node has been marked as unmovable
        if (moved[nodeIdV0]) continue;

        for (auto itV1 = listV1.begin(); itV1 != listV1.end(); ++itV1) {
            const auto &[gainV1, nodeIdV1] = *itV1;
            uint64_t nodeIdV1tmp = nodeIdV1;

            // Check if the node has been marked as unmovable
            if (moved[nodeIdV1]) continue;

            // Check if nodeIdV0 has edge to nodeIdV1. If so, they cannot be swapped.
            auto it = std::find_if(workingGraph.adj[nodeIdV0].begin(), workingGraph.adj[nodeIdV0].end(),
                                   [nodeIdV1tmp](const auto &pair) { return pair.first == nodeIdV1tmp; });
            if (it != workingGraph.adj[nodeIdV0].end()) {
                continue;
            }

            // Check if the move causes reverse imbalance
            uint64_t newSizeV0 = sizeV0 - workingGraph.nodeWeights[nodeIdV0] + workingGraph.nodeWeights[nodeIdV1];
            uint64_t newSizeV1 = sizeV1 - workingGraph.nodeWeights[nodeIdV1] + workingGraph.nodeWeights[nodeIdV0];
            uint64_t newSizeOfPreviousBigger = sizeV0 > sizeV1 ? newSizeV0 : newSizeV1;
            uint64_t newSizeOfPreviousSmaller = sizeV0 > sizeV1 ? newSizeV1 : newSizeV0;
            uint64_t maxNodeWeight = workingGraph.maxNodeWeight();
            if ((double) newSizeOfPreviousBigger < lowerBoundPartWeight - (double) maxNodeWeight ||
                (double) newSizeOfPreviousSmaller > upperBoundPartWeight + (double) maxNodeWeight)
                continue;

            // Pair is swappable AND improves balance more than the rest
            double diffFromPerfectBalanceTemp = std::abs(
                    50.0 - (((double) newSizeV0 * 100.0) / ((double) (newSizeV0 + newSizeV1))));
            if (diffFromPerfectBalanceTemp < diffFromPerfectBalance) {
                found = true;
                bestNodeV0 = nodeIdV0;
                bestNodeV1 = nodeIdV1;
                gain = gainV0 + gainV1;
                bestItV0 = itV0;
                bestItV1 = itV1;
            }
        }
    }

    if (found) {
        listV0.erase(bestItV0);
        listV1.erase(bestItV1);
    }
    return {found, bestNodeV0, bestNodeV1, gain};
}

bool BoundaryKL::onePassRefinement() {
    std::vector<bool> moved(workingGraph.size, false);
    std::vector<bool> inList(workingGraph.size, false);

    std::list<std::pair<int64_t, uint64_t>> listV0, listV1;

    // Find initially movable vertices
    insertMovableNodesIntoLists(listV0, listV1, inList);
    listV0.sort([](const auto &a, const auto &b) { return a.first > b.first; });
    listV1.sort([](const auto &a, const auto &b) { return a.first > b.first; });

    std::vector<std::pair<uint64_t, uint64_t>> moveSequence;
    uint64_t currentEdgeCut = initialEdgeCut;
    uint64_t bestEdgeCut = initialEdgeCut;
    size_t bestMovePrefix = 0;
    std::vector<bool> initialBisectionInfoTemp = initialBisectionInfo;

    uint64_t maxNodeWeight = workingGraph.maxNodeWeight();
    auto [sizeV0, sizeV1] = calculatePartSizes();
    bool isBalanced = checkBalance(maxNodeWeight);
    uint64_t minMaxPartSize = std::max(sizeV0, sizeV1);

    // Main refinement loop - make moves until no more vertices can move
    while (!listV0.empty() && !listV1.empty()) {
        auto [found, nodeIdV0, nodeIdV1, gain] = isBalanced
                                                 ? findBestMovablePairBalanced(listV0, listV1, moved, sizeV0, sizeV1)
                                                 : findBestMovablePairUnbalanced(listV0, listV1, moved, sizeV0, sizeV1);


        // Tentative move
        sizeV0 = sizeV0 - workingGraph.nodeWeights[nodeIdV0] + workingGraph.nodeWeights[nodeIdV1];
        sizeV1 = sizeV1 - workingGraph.nodeWeights[nodeIdV1] + workingGraph.nodeWeights[nodeIdV0];
        moved[nodeIdV0] = moved[nodeIdV1] = true;
        moveSequence.emplace_back(nodeIdV0, nodeIdV1);
        currentEdgeCut -= gain;
        initialBisectionInfoTemp[nodeIdV0] = true;
        initialBisectionInfoTemp[nodeIdV1] = false;
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
        if ((double) sizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
            (double) sizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
            (double) sizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
            (double) sizeV1 > upperBoundPartWeight + (double) maxNodeWeight)
            isBalanced = false;


        // Moving from V0 to V1 can make V1 nodes unmovable
        // If an out-neighbor of the moved nodeIdV0 is in V1 and is not yet moved,
        // (if it has already been moved, we don't care if it is movable or not)
        // Then the neighbor becomes unmovable since it has one in-neighbor (the moved nodeIdV0) that is in V1
        for (const auto &[neighborId, _]: workingGraph.adj[nodeIdV0]) {
            if (!moved[neighborId] && initialBisectionInfo[neighborId] && inList[neighborId]) {
                moved[neighborId] = true;
            }
        }
        // Moving from V0 to V1 can make V0 nodes movable
        // If an in-neighbor of the moved nodeIdV0 is in V0 and is not already movable,
        // (if it has already been marked movable by the initial check, we don't recheck)
        // Then the neighbor might become movable if its only out-neighbor in V0 was the moved nodeIdV0
        insertMovableNeighborsIntoLists(initialBisectionInfoTemp, listV0, listV1, inList, nodeIdV0,
                                        false);

        // Moving from V1 to V0 can make V0 nodes unmovable
        // If an in-neighbor of the moved nodeIdV1 is in V0 and is not yet moved,
        // (if it has already been moved, we don't care if it is movable or not)
        // Then the neighbor becomes unmovable since it has one out-neighbor (the moved nodeIdV1) that is in V0
        for (const auto &[neighborId, _]: workingGraph.revAdj[nodeIdV1]) {
            if (!moved[neighborId] && !initialBisectionInfo[neighborId] && inList[neighborId]) {
                moved[neighborId] = true;
            }
        }
        // Moving from V1 to V0 can make V1 nodes movable
        // If an out-neighbor of the moved nodeIdV1 is in V1 and is not already movable,
        // (if it has already been marked movable by the initial check, we don't recheck)
        // Then the neighbor might become movable if its only in-neighbor in V1 was the moved nodeV1
        insertMovableNeighborsIntoLists(initialBisectionInfoTemp, listV0, listV1, inList, nodeIdV1, true);
    }

    // If no improvement was found, return false
    if (bestMovePrefix == 0) return false;

    // Check if should keep current solution
    auto [initialSizeV0, initialSizeV1] = calculatePartSizes();
    if (initialEdgeCut <= bestEdgeCut && std::max(sizeV0, sizeV1) >= std::max(initialSizeV0, initialSizeV1))
        return false;

    // Apply best move sequence
    for (size_t i = 0; i < bestMovePrefix; ++i) {
        auto [nodeId, fromV0] = moveSequence[i];
        initialBisectionInfo[nodeId] = fromV0;
    }
    initialEdgeCut = bestEdgeCut;

    return true;
}
