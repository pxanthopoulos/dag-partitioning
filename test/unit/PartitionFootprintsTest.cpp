#include "PartitionFootprints.h"
#include "TestGraphs.h"

#include <gtest/gtest.h>

using dag_partitioning::scheduling::PartitionFootprints;
using testgraphs::Minimal;

constexpr uint64_t None = PartitionFootprints::None;

TEST(PartitionFootprints, KeepsTensorsSharedWithANeighbor) {
    Minimal minimal;
    // After partition 1, B and C stay on the device for all of partition 2,
    // while A is fetched for its one use
    EXPECT_EQ(minimal.footprints.footprint(2, 1, 3), 8u + 4u);
    // After partition 0, A stays; C and B are fetched one after the other
    EXPECT_EQ(minimal.footprints.footprint(2, 0, 3), 4u + 4u);
    // Alone, all three share 4 bytes
    EXPECT_EQ(minimal.footprints.footprint(2, None, 3), 4u);
    EXPECT_EQ(minimal.footprints.footprint(1, None, None), 4u);
    EXPECT_EQ(minimal.footprints.footprint(1, None, 2), 8u);
}

TEST(PartitionFootprints, PeakOfAnOrder) {
    Minimal minimal;
    EXPECT_EQ(minimal.footprints.peak({0, 1, 2}), 12u);
    EXPECT_EQ(minimal.footprints.peak({1, 0, 2}), 8u);
}

TEST(PartitionFootprints, LiveInsAreHeldUntilTheirLastUser) {
    Minimal minimal(true);
    // The live-in is held by partition 0 and every partition before it
    EXPECT_EQ(minimal.footprints.peak({0, 1, 2}), 12u);
    EXPECT_EQ(minimal.footprints.peak({1, 0, 2}), 9u);
    // footprint leaves it out
    EXPECT_EQ(minimal.footprints.footprint(0, 1, 2), 4u);
}

TEST(PartitionFootprints, LastPartitionHoldsTheLiveOuts) {
    // Partition 0 makes a live-out tensor of 6 bytes; partition 1 has a
    // local tensor of 3
    PartitionFootprints footprints({{}, {{3, 0, 0, 0}}}, {{6, {0}, true}},
                                   {{{0, 0, 0}}, {}}, testgraphs::bestFit());
    EXPECT_EQ(footprints.footprint(0, None, 1), 6u);
    EXPECT_EQ(footprints.footprint(0, 1, 2), 6u);
    EXPECT_EQ(footprints.footprint(1, 0, 2), 6u + 3u);
    EXPECT_EQ(footprints.footprint(1, None, 0), 3u);
    EXPECT_EQ(footprints.peak({0, 1}), 9u);
}

TEST(PartitionFootprints, FetchedTensorsShareTheSlab) {
    // A fetched tensor is packed with the local tensors by its uses
    PartitionFootprints footprints({{}, {{5, 1, 1, 0}}}, {{7, {0, 1}, false}},
                                   {{{0, 0, 0}}, {{0, 0, 0}}},
                                   testgraphs::bestFit());
    EXPECT_EQ(footprints.footprint(1, None, None), 7u);
    EXPECT_EQ(footprints.footprint(1, 0, 2), 7u + 5u);
    // The slab size is cached by which uses are kept
    EXPECT_EQ(footprints.footprint(1, None, None), 7u);
}
