#include "Graph.h"
#include "ClusteringForbiddenEdges.h"
#include "ClusteringCycleDetection.h"
#include "ClusteringHybrid.h"
#include <iostream>
#include <chrono>

int main() {
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/2mm_10_20_30_40.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/test/dag.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/before.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/partition-graph-1.dot",
                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");

    ClusteringHybrid clusteringHyb(graph);
    clusteringHyb.runClustering();

    return 0;
}
