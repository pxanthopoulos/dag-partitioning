//
// Created by panagiotis on 23/12/2024.
//

#include <queue>
#include "BoundaryFM.h"

BoundaryFM::BoundaryFM(const Graph &graph, std::vector<bool> &initialBisectionInfo, uint64_t maxNumberOfPasses)
        : Refinement(graph, initialBisectionInfo, maxNumberOfPasses) {}

bool BoundaryFM::onePassRefinement() {
    std::vector<bool> moved(workingGraph.size, false);

    std::priority_queue<std::pair<int64_t, uint64_t>> heapV0, heapV1;

    for (uint64_t node = 0; node < workingGraph.size; ++node) {
        if (!initialBisectionInfo[node]) {
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.adj[node]) {
                if (!initialBisectionInfo[neighborId]) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[node]) {
                    gain += (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[node]) {
                    gain -= (int64_t) edgeWeight;
                }
                heapV0.emplace(gain, node);
            }
        } else {
            bool movable = true;
            for (const auto &[neighborId, _]: workingGraph.revAdj[node]) {
                if (initialBisectionInfo[neighborId]) {
                    movable = false;
                    break;
                }
            }

            if (movable) {
                int64_t gain = 0;
                for (const auto &[_, edgeWeight]: workingGraph.adj[node]) {
                    gain -= (int64_t) edgeWeight;
                }
                for (const auto &[_, edgeWeight]: workingGraph.revAdj[node]) {
                    gain += (int64_t) edgeWeight;
                }
                heapV1.emplace(gain, node);
            }
        }
    }

    std::vector<std::pair<uint64_t, bool>> moveSequence;
    int64_t currentGain = 0;
    int64_t bestGain = 0;
    size_t bestMovePrefix = 0;

    while (!heapV0.empty() || !heapV1.empty()) {
        bool moveFromV0 = false;
        uint64_t vertex;
        int64_t gain;

        if (heapV0.empty()) {
            std::tie(gain, vertex) = heapV1.top();
            heapV1.pop();
        } else if (heapV1.empty()) {
            std::tie(gain, vertex) = heapV0.top();
            heapV0.pop();
            moveFromV0 = true;
        } else {
            auto v0_best = heapV0.top();
            auto v1_best = heapV1.top();
            if (v0_best.first >= v1_best.first) {
                std::tie(gain, vertex) = v0_best;
                heapV0.pop();
                moveFromV0 = true;
            } else {
                std::tie(gain, vertex) = v1_best;
                heapV1.pop();
            }
        }

        if (moved[vertex]) continue;

        // Perform tentative move
        moved[vertex] = true;
        moveSequence.emplace_back(vertex, moveFromV0);
        currentGain += gain;

        // Update best solution seen
        if (currentGain > bestGain) {
            bestGain = currentGain;
            bestMovePrefix = moveSequence.size();
        }

        // After move, some vertices might become unmovable
        // Remove them from heaps (they'll be skipped when popped)
        // No need to update gains as they're static
    }

    if (bestGain <= 0) return false;

    for (size_t i = 0; i < bestMovePrefix; ++i) {
        auto [node, fromV0] = moveSequence[i];
        initialBisectionInfo[node] = !fromV0;
    }

    return true;
}
