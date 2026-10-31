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

uint64_t packBestFit(std::vector<Tensor> &tensors, uint64_t offsetAlignment,
                     uint64_t rangeAlignment) {
    auto align = [](uint64_t value, uint64_t alignment) {
        return (value + alignment - 1) / alignment * alignment;
    };

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
        t.size = align(t.size, rangeAlignment);

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

        // Best-fit: the smallest gap that holds t at an aligned offset, else
        // the aligned top of the last interval
        uint64_t bestOffset = UINT64_MAX;
        uint64_t bestGap = UINT64_MAX;
        uint64_t gapStart = 0;
        for (auto &[start, end] : merged) {
            uint64_t offset = align(gapStart, offsetAlignment);
            if (offset + t.size <= start && start - gapStart < bestGap) {
                bestGap = start - gapStart;
                bestOffset = offset;
            }
            gapStart = end;
        }
        if (bestOffset == UINT64_MAX)
            bestOffset = align(gapStart, offsetAlignment);

        t.offset = bestOffset;
    }

    uint64_t peak = 0;
    for (const Tensor &t : tensors) {
        peak = std::max(peak, t.offset + t.size);
    }

    assert(verify(tensors) && "Memory packing verification failed");

    return align(peak, rangeAlignment);
}

} // namespace utils

} // namespace scheduling

} // namespace dag_partitioning