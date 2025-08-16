/**
 * @file BoundaryKL.cpp
 * @brief Implementation of boundary KL refinement for directed graphs
 */

#include "BoundaryKL.h"
#include "Graph.h"

#include <algorithm>

namespace dag_partitioning {

namespace refinement {

BoundaryKL::BoundaryKL(const core::Graph &graph,
                       std::vector<uint8_t> &initialBisectionInfo,
                       uint64_t &initialEdgeCut, uint64_t maxNumberOfPasses,
                       double upperBoundPartWeight, double lowerBoundPartWeight)
    : Refinement(graph, initialBisectionInfo, initialEdgeCut, maxNumberOfPasses,
                 upperBoundPartWeight, lowerBoundPartWeight) {}

void BoundaryKL::insertMovableNodesIntoLists(
    std::vector<std::pair<int64_t, uint64_t>> &vecV0,
    std::vector<std::pair<int64_t, uint64_t>> &vecV1,
    std::vector<uint8_t> &inList) const {

    // Check each vertex for movability
    for (uint64_t nodeId = 0; nodeId < workingGraph.size; ++nodeId) {
        // If the node has already been in a vec previously, don't re-add it
        if (inList[nodeId] == 1)
            continue;

        if (initialBisectionInfo[nodeId] == 0) { // Node in V0
            // Check if movable to V1: all out-neighbors must be in V1
            bool movable = true;
            for (const auto &[neighborId, _] : workingGraph.adj[nodeId]) {
                if (initialBisectionInfo[neighborId] == 0) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain according to equation 4.1
                int64_t gain = 0;
                for (const auto &[_, edgeWeight] : workingGraph.adj[nodeId]) {
                    gain += (int64_t)edgeWeight;
                }
                for (const auto &[_, edgeWeight] :
                     workingGraph.revAdj[nodeId]) {
                    gain -= (int64_t)edgeWeight;
                }
                vecV0.emplace_back(gain, nodeId);
                inList[nodeId] = 1;
            }
        } else { // Node in V1
            // Check if movable to V0: all in-neighbors must be in V0
            bool movable = true;
            for (const auto &[neighborId, _] : workingGraph.revAdj[nodeId]) {
                if (initialBisectionInfo[neighborId] == 1) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain (negative of V0->V1 gain)
                int64_t gain = 0;
                for (const auto &[_, edgeWeight] : workingGraph.adj[nodeId]) {
                    gain -= (int64_t)edgeWeight;
                }
                for (const auto &[_, edgeWeight] :
                     workingGraph.revAdj[nodeId]) {
                    gain += (int64_t)edgeWeight;
                }
                vecV1.emplace_back(gain, nodeId);
                inList[nodeId] = 1;
            }
        }
    }
}

// Called after moving a node to update movable status of its neighbors
// (possibly insert new movable nodes)
void BoundaryKL::insertMovableNeighborsIntoLists(
    const std::vector<uint8_t> &currentBisectionInfo,
    std::vector<std::pair<int64_t, uint64_t>> &vecV0,
    std::vector<std::pair<int64_t, uint64_t>> &vecV1,
    std::vector<uint8_t> &inList, uint64_t movedNodeId,
    bool checkOutNeighbors) const {
    if (checkOutNeighbors) { // Node moved V1->V0, check out-neighbors
        for (const auto &[nodeId, _] : workingGraph.adj[movedNodeId]) {
            // If the node has already been in the heap previously, don't re-add
            // it
            if (inList[nodeId] == 1 || currentBisectionInfo[nodeId] == 0)
                continue;

            bool movable = true;
            for (const auto &[neighborId, _] : workingGraph.revAdj[nodeId]) {
                if (currentBisectionInfo[neighborId] == 1) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain for potential move
                int64_t gain = 0;
                for (const auto &[_, edgeWeight] : workingGraph.adj[nodeId]) {
                    gain -= (int64_t)edgeWeight;
                }
                for (const auto &[_, edgeWeight] :
                     workingGraph.revAdj[nodeId]) {
                    gain += (int64_t)edgeWeight;
                }
                auto it = std::upper_bound(vecV1.begin(), vecV1.end(),
                                           std::make_pair(gain, nodeId),
                                           [](const auto &a, const auto &b) {
                                               return a.first > b.first;
                                           });

                vecV1.insert(it, std::make_pair(gain, nodeId));
                inList[nodeId] = 1;
            }
        }
    } else { // Node moved V0->V1, check in-neighbors
        for (const auto &[nodeId, _] : workingGraph.revAdj[movedNodeId]) {
            // If the node has already been in the heap previously, don't re-add
            // it
            if (inList[nodeId] == 1 || currentBisectionInfo[nodeId] == 1)
                continue;

            bool movable = true;
            for (const auto &[neighborId, _] : workingGraph.adj[nodeId]) {
                if (currentBisectionInfo[neighborId] == 0) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                // Calculate gain for potential move
                int64_t gain = 0;
                for (const auto &[_, edgeWeight] : workingGraph.adj[nodeId]) {
                    gain += (int64_t)edgeWeight;
                }
                for (const auto &[_, edgeWeight] :
                     workingGraph.revAdj[nodeId]) {
                    gain -= (int64_t)edgeWeight;
                }
                auto it = std::upper_bound(vecV0.begin(), vecV0.end(),
                                           std::make_pair(gain, nodeId),
                                           [](const auto &a, const auto &b) {
                                               return a.first > b.first;
                                           });

                vecV0.insert(it, std::make_pair(gain, nodeId));
                inList[nodeId] = 1;
            }
        }
    }
}

std::tuple<bool, uint64_t, uint64_t, int64_t>
BoundaryKL::findBestMovablePairBalanced(
    std::vector<std::pair<int64_t, uint64_t>> &vecV0,
    std::vector<std::pair<int64_t, uint64_t>> &vecV1,
    std::vector<uint8_t> &moved, uint64_t &sizeV0, uint64_t &sizeV1) const {
    bool found = false;
    int64_t maxTotalGain = INT64_MIN;
    uint64_t bestNodeV0 = -1;
    uint64_t bestNodeV1 = -1;
    uint64_t bestIdxV0 = 0;
    uint64_t bestIdxV1 = 0;

    for (uint64_t i = 0; i < vecV0.size(); ++i) {
        const auto &[gainV0, nodeIdV0] = vecV0[i];
        // Check if the node has been marked as unmovable
        if (moved[nodeIdV0] == 1)
            continue;

        for (uint64_t j = 0; j < vecV1.size(); ++j) {
            const auto &[gainV1, nodeIdV1] = vecV1[j];
            uint64_t nodeIdV1tmp = nodeIdV1;

            // Early exit. Before checking if the pair is swappable, we can
            // check if the gain is better than the best so far. If not, we can
            // continue to the next pair. Even better, we can break the inner
            // loop as the gain will only decrease and go to the next outer loop
            // iteration
            int64_t totalGain = gainV0 + gainV1;
            if (totalGain <= maxTotalGain)
                break;

            // Check if the node has been marked as unmovable
            if (moved[nodeIdV1] == 1)
                continue;

            // If this move breaks the balance, skip it
            uint64_t newSizeV0 = sizeV0 - workingGraph.nodeWeights[nodeIdV0] +
                                 workingGraph.nodeWeights[nodeIdV1];
            uint64_t newSizeV1 = sizeV1 - workingGraph.nodeWeights[nodeIdV1] +
                                 workingGraph.nodeWeights[nodeIdV0];
            if (!checkBalance(newSizeV0, newSizeV1, 0, upperBoundPartWeight,
                              lowerBoundPartWeight))
                continue;

            // Check if nodeIdV0 has edge to nodeIdV1. If so, they cannot be
            // swapped.
            auto it = std::find_if(workingGraph.adj[nodeIdV0].begin(),
                                   workingGraph.adj[nodeIdV0].end(),
                                   [nodeIdV1tmp](const auto &pair) {
                                       return pair.first == nodeIdV1tmp;
                                   });
            if (it != workingGraph.adj[nodeIdV0].end()) {
                continue;
            }

            // Pair is swappable AND improves gain
            found = true;
            bestNodeV0 = nodeIdV0;
            bestNodeV1 = nodeIdV1;
            maxTotalGain = totalGain;
            bestIdxV0 = i;
            bestIdxV1 = j;
        }
    }

    if (found) {
        if (bestIdxV0 != vecV0.size() - 1) {
            std::swap(vecV0[bestIdxV0], vecV0.back());
        }
        vecV0.pop_back();

        if (bestIdxV1 != vecV1.size() - 1) {
            std::swap(vecV1[bestIdxV1], vecV1.back());
        }
        vecV1.pop_back();
    }
    return {found, bestNodeV0, bestNodeV1, maxTotalGain};
}

std::tuple<bool, uint64_t, uint64_t, int64_t>
BoundaryKL::findBestMovablePairUnbalanced(
    std::vector<std::pair<int64_t, uint64_t>> &vecV0,
    std::vector<std::pair<int64_t, uint64_t>> &vecV1,
    std::vector<uint8_t> &moved, uint64_t &sizeV0, uint64_t &sizeV1) const {
    bool found = false;
    // Perfect balance is each size is 50% of the total. Here we calculate the
    // difference from the perfect balance. Lower is better.
    double diffFromPerfectBalance = std::abs(
        50.0 - (((double)sizeV0 * 100.0) / ((double)(sizeV0 + sizeV1))));
    uint64_t bestNodeV0 = -1;
    uint64_t bestNodeV1 = -1;
    uint64_t gain = 0;
    uint64_t bestIdxV0 = 0;
    uint64_t bestIdxV1 = 0;

    for (uint64_t i = 0; i < vecV0.size(); ++i) {
        const auto &[gainV0, nodeIdV0] = vecV0[i];
        // Check if the node has been marked as unmovable
        if (moved[nodeIdV0] == 1)
            continue;

        for (uint64_t j = 0; j < vecV1.size(); ++j) {
            const auto &[gainV1, nodeIdV1] = vecV1[j];
            uint64_t nodeIdV1tmp = nodeIdV1;

            // Check if the node has been marked as unmovable
            if (moved[nodeIdV1] == 1)
                continue;

            // Check if nodeIdV0 has edge to nodeIdV1. If so, they cannot be
            // swapped.
            auto it = std::find_if(workingGraph.adj[nodeIdV0].begin(),
                                   workingGraph.adj[nodeIdV0].end(),
                                   [nodeIdV1tmp](const auto &pair) {
                                       return pair.first == nodeIdV1tmp;
                                   });
            if (it != workingGraph.adj[nodeIdV0].end()) {
                continue;
            }

            // Check if the move causes reverse imbalance
            uint64_t newSizeV0 = sizeV0 - workingGraph.nodeWeights[nodeIdV0];
            uint64_t newSizeV1 = sizeV1 - workingGraph.nodeWeights[nodeIdV1];
            uint64_t newSizeOfPreviousBigger =
                sizeV0 > sizeV1 ? newSizeV0 : newSizeV1;
            uint64_t newSizeOfPreviousSmaller =
                sizeV0 > sizeV1 ? newSizeV1 : newSizeV0;
            if ((double)newSizeOfPreviousBigger < lowerBoundPartWeight ||
                (double)newSizeOfPreviousSmaller > upperBoundPartWeight)
                continue;

            double diffFromPerfectBalanceTemp =
                std::abs(50.0 - (((double)newSizeV0 * 100.0) /
                                 ((double)(newSizeV0 + newSizeV1))));
            if (diffFromPerfectBalanceTemp < diffFromPerfectBalance) {
                // Pair is swappable AND improves balance more than the rest
                found = true;
                bestNodeV0 = nodeIdV0;
                bestNodeV1 = nodeIdV1;
                gain = gainV0 + gainV1;
                bestIdxV0 = i;
                bestIdxV1 = j;
            }
        }
    }

    if (found) {
        if (bestIdxV0 != vecV0.size() - 1) {
            std::swap(vecV0[bestIdxV0], vecV0.back());
        }
        vecV0.pop_back();

        if (bestIdxV1 != vecV1.size() - 1) {
            std::swap(vecV1[bestIdxV1], vecV1.back());
        }
        vecV1.pop_back();
    }
    return {found, bestNodeV0, bestNodeV1, gain};
}

bool BoundaryKL::onePassRefinement() {
    std::vector<uint8_t> moved(workingGraph.size, 0);
    std::vector<uint8_t> inList(workingGraph.size, 0);

    std::vector<std::pair<int64_t, uint64_t>> vecV0, vecV1;

    // Find initially movable vertices
    insertMovableNodesIntoLists(vecV0, vecV1, inList);
    // Vectors are maintained in descending gain order to enable early exit
    // optimizations
    std::sort(vecV0.begin(), vecV0.end(), [](const auto &a, const auto &b) {
        if (a.first != b.first) {
            return a.first > b.first;
        }
        return a.second <
               b.second; // Secondary sort criteria for equal first values
    });
    std::sort(vecV1.begin(), vecV1.end(), [](const auto &a, const auto &b) {
        if (a.first != b.first) {
            return a.first > b.first;
        }
        return a.second <
               b.second; // Secondary sort criteria for equal first values
    });

    std::vector<std::pair<uint64_t, uint64_t>> moveSequence;
    uint64_t currentEdgeCut = initialEdgeCut;
    uint64_t bestEdgeCut = initialEdgeCut;
    uint64_t bestMovePrefix = 0, noImprovement = 0;
    std::vector<uint8_t> initialBisectionInfoTemp = initialBisectionInfo;

    auto [sizeV0, sizeV1] =
        calculatePartSizes(initialBisectionInfo, workingGraph);
    bool isBalanced = checkBalance(sizeV0, sizeV1, 0, upperBoundPartWeight,
                                   lowerBoundPartWeight);
    uint64_t minMaxPartSize = std::max(sizeV0, sizeV1);

    // Main refinement loop - make moves until no more vertices can move
    while (!vecV0.empty() && !vecV1.empty() &&
           noImprovement <= workingGraph.size / 4 &&
           moveSequence.size() < workingGraph.size) {
        auto [found, nodeIdV0, nodeIdV1, gain] =
            isBalanced ? findBestMovablePairBalanced(vecV0, vecV1, moved,
                                                     sizeV0, sizeV1)
                       : findBestMovablePairUnbalanced(vecV0, vecV1, moved,
                                                       sizeV0, sizeV1);
        if (!found)
            break;

        // Tentative move
        sizeV0 = sizeV0 - workingGraph.nodeWeights[nodeIdV0] +
                 workingGraph.nodeWeights[nodeIdV1];
        sizeV1 = sizeV1 - workingGraph.nodeWeights[nodeIdV1] +
                 workingGraph.nodeWeights[nodeIdV0];
        moved[nodeIdV0] = moved[nodeIdV1] = 1;
        moveSequence.emplace_back(nodeIdV0, nodeIdV1);
        currentEdgeCut -= gain;
        initialBisectionInfoTemp[nodeIdV0] = 1;
        initialBisectionInfoTemp[nodeIdV1] = 0;
        uint64_t maxPartSize = std::max(sizeV0, sizeV1);

        // If partition is unbalanced, always include this move to the sequence
        // If balanced, include move if it improves edge cut or improves balance
        // quality (without worsening edge cut)
        if (currentEdgeCut < bestEdgeCut ||
            (currentEdgeCut == bestEdgeCut && maxPartSize < minMaxPartSize)) {
            bestEdgeCut = currentEdgeCut;
            bestMovePrefix = moveSequence.size();
            minMaxPartSize = maxPartSize;
            noImprovement = 0;
        }

        isBalanced = checkBalance(sizeV0, sizeV1, 0, upperBoundPartWeight,
                                  lowerBoundPartWeight);
        noImprovement++;

        // Moving from V0 to V1 can make V1 nodes unmovable
        // If an out-neighbor of the moved nodeIdV0 is in V1 and is not yet
        // moved, (if it has already been moved, we don't care if it is movable
        // or not) Then the neighbor becomes unmovable since it has one
        // in-neighbor (the moved nodeIdV0) that is in V1
        for (const auto &[neighborId, _] : workingGraph.adj[nodeIdV0]) {
            if (moved[neighborId] == 0 &&
                initialBisectionInfoTemp[neighborId] == 1 &&
                inList[neighborId] == 1) {
                moved[neighborId] = 1;
            }
        }
        // Moving from V0 to V1 can make V0 nodes movable
        // If an in-neighbor of the moved nodeIdV0 is in V0 and is not already
        // movable, (if it has already been marked movable by the initial check,
        // we don't recheck) Then the neighbor might become movable if its only
        // out-neighbor in V0 was the moved nodeIdV0
        insertMovableNeighborsIntoLists(initialBisectionInfoTemp, vecV0, vecV1,
                                        inList, nodeIdV0, false);

        // Moving from V1 to V0 can make V0 nodes unmovable
        // If an in-neighbor of the moved nodeIdV1 is in V0 and is not yet
        // moved, (if it has already been moved, we don't care if it is movable
        // or not) Then the neighbor becomes unmovable since it has one
        // out-neighbor (the moved nodeIdV1) that is in V0
        for (const auto &[neighborId, _] : workingGraph.revAdj[nodeIdV1]) {
            if (moved[neighborId] == 0 &&
                initialBisectionInfoTemp[neighborId] == 0 &&
                inList[neighborId] == 1) {
                moved[neighborId] = 1;
            }
        }
        // Moving from V1 to V0 can make V1 nodes movable
        // If an out-neighbor of the moved nodeIdV1 is in V1 and is not already
        // movable, (if it has already been marked movable by the initial check,
        // we don't recheck) Then the neighbor might become movable if its only
        // in-neighbor in V1 was the moved nodeV1
        insertMovableNeighborsIntoLists(initialBisectionInfoTemp, vecV0, vecV1,
                                        inList, nodeIdV1, true);
    }

    // If no improvement was found, return false
    if (bestMovePrefix == 0)
        return false;

    // Apply best move sequence
    for (uint64_t i = 0; i < bestMovePrefix; ++i) {
        auto [nodeIdV0, nodeIdV1] = moveSequence[i];
        initialBisectionInfo[nodeIdV0] = 1;
        initialBisectionInfo[nodeIdV1] = 0;
    }
    initialEdgeCut = bestEdgeCut;

    return true;
}

} // namespace refinement

} // namespace dag_partitioning
