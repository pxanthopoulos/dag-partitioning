#include "Utils.h"

#include <gtest/gtest.h>

using dag_partitioning::scheduling::utils::packBestFit;
using dag_partitioning::scheduling::utils::sortByLifetime;
using dag_partitioning::scheduling::utils::Tensor;

TEST(SortByLifetime, SortsByBirthThenDeath) {
    std::vector<Tensor> tensors = {
        {1, 2, 3, 0}, {2, 0, 5, 0}, {3, 2, 2, 0}, {4, 0, 1, 0}};
    sortByLifetime(tensors);
    std::vector<uint64_t> sizes;
    for (const Tensor &t : tensors)
        sizes.push_back(t.size);
    EXPECT_EQ(sizes, (std::vector<uint64_t>{4, 2, 3, 1}));
}

TEST(PackBestFit, EmptyIsZero) {
    std::vector<Tensor> tensors;
    EXPECT_EQ(packBestFit(tensors), 0u);
}

TEST(PackBestFit, OverlappingLifetimesStack) {
    std::vector<Tensor> tensors = {{10, 0, 1, 0}, {20, 1, 2, 0}};
    EXPECT_EQ(packBestFit(tensors), 30u);
    EXPECT_EQ(tensors[0].offset, 0u);
    EXPECT_EQ(tensors[1].offset, 10u);
}

TEST(PackBestFit, DisjointLifetimesShare) {
    std::vector<Tensor> tensors = {{10, 0, 0, 0}, {20, 1, 1, 0}};
    EXPECT_EQ(packBestFit(tensors), 20u);
    EXPECT_EQ(tensors[1].offset, 0u);
}

TEST(PackBestFit, BestFitTakesTheSmallestGap) {
    // At step 1 the first and third have died, leaving gaps [0, 8) and
    // [12, 16): the last goes into the smaller one
    std::vector<Tensor> tensors = {
        {8, 0, 0, 0}, {4, 0, 2, 0}, {4, 0, 0, 0}, {4, 0, 2, 0}, {4, 1, 2, 0}};
    EXPECT_EQ(packBestFit(tensors), 20u);
    EXPECT_EQ(tensors[4].offset, 12u);
}

TEST(PackBestFit, AlignsOffsetsAndSizes) {
    std::vector<Tensor> tensors = {{10, 0, 1, 0}, {20, 1, 2, 0}};
    EXPECT_EQ(packBestFit(tensors, 64, 16), 96u);
    EXPECT_EQ(tensors[0].size, 16u);
    EXPECT_EQ(tensors[1].size, 32u);
    EXPECT_EQ(tensors[1].offset, 64u);
}
