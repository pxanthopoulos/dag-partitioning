/**
 * @file Utils.cpp
 * @brief Implementation of utility functions for memory-aware scheduling
 */

#include "Utils.h"

#include <algorithm>
#include <cassert>
#include <iostream>

namespace dag_partitioning {

namespace scheduling {

namespace utils {

uint64_t packBestFit(std::vector<Tensor> &tensors) {
    auto overlaps = [](const Tensor &a, const Tensor &b) {
        return !(a.death < b.birth || b.death < a.birth);
    };

    auto verify = [&overlaps](const std::vector<Tensor> &tensors) {
        for (size_t i = 0; i < tensors.size(); ++i) {
            for (size_t j = i + 1; j < tensors.size(); ++j) {
                if (overlaps(tensors[i], tensors[j])) {
                    uint64_t end1 = tensors[i].offset + tensors[i].size;
                    uint64_t end2 = tensors[j].offset + tensors[j].size;
                    if (!(end1 <= tensors[j].offset ||
                          end2 <= tensors[i].offset)) {
                        std::cerr << "CONFLICT: tensor " << i << " and " << j
                                  << "\n";
                        return false;
                    }
                }
            }
        }
        return true;
    };

    if (tensors.empty())
        return 0;

    for (size_t i = 0; i < tensors.size(); ++i) {
        Tensor &t = tensors[i];

        std::vector<std::pair<uint64_t, uint64_t>> occupied;
        for (size_t j = 0; j < i; ++j) { // all previously placed
            if (overlaps(t, tensors[j])) {
                occupied.push_back(
                    {tensors[j].offset, tensors[j].offset + tensors[j].size});
            }
        }

        // Sort and merge overlapping intervals
        std::sort(occupied.begin(), occupied.end());
        std::vector<std::pair<uint64_t, uint64_t>> merged;
        for (auto &iv : occupied) {
            if (merged.empty() || merged.back().second < iv.first) {
                merged.push_back(iv);
            } else {
                merged.back().second =
                    std::max(merged.back().second, iv.second);
            }
        }

        // Best-fit: find smallest gap >= t.size
        uint64_t bestOffset = 0;
        uint64_t bestWaste = UINT64_MAX;

        if (merged.empty()) {
            bestOffset = 0;
        } else {
            // Gap before first interval [0, merged[0].first)
            if (merged[0].first >= t.size &&
                merged[0].first - t.size < bestWaste) {
                bestWaste = merged[0].first - t.size;
                bestOffset = 0;
            }

            // Gaps between intervals
            for (size_t j = 0; j + 1 < merged.size(); ++j) {
                uint64_t gapSize = merged[j + 1].first - merged[j].second;
                if (gapSize >= t.size && gapSize - t.size < bestWaste) {
                    bestWaste = gapSize - t.size;
                    bestOffset = merged[j].second;
                }
            }

            // If no gap fits, place after last interval
            if (bestWaste == UINT64_MAX) {
                bestOffset = merged.back().second;
            }
        }

        t.offset = bestOffset;
    }

    uint64_t peak = 0;
    for (const Tensor &t : tensors) {
        peak = std::max(peak, t.offset + t.size);
    }

    assert(verify(tensors) && "Memory packing verification failed");

    return peak;
}

} // namespace utils

} // namespace scheduling

} // namespace dag_partitioning