#include "Graph.h"
#include "ClusteringForbiddenEdges.h"
#include "ClusteringCycleDetection.h"
#include "ClusteringHybrid.h"
#include "GreedyDirectedGraphGrowing.h"
#include "UndirectedFix.h"
#include <iostream>
#include <chrono>

int main() {
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/2mm_10_20_30_40.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/test/dag.dot",
                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/before.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/partition-graph-1.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");

//    ClusteringForbiddenEdges clusteringForb(graph);
//    clusteringForb.runClustering();
//
//    ClusteringCycleDetection clusteringCyc(graph);
//    clusteringCyc.runClustering();
//
//    ClusteringHybrid clusteringHyb(graph);
//    clusteringHyb.runClustering();

    GreedyDirectedGraphGrowing GDGG(graph, 1.2 * ((double) graph.totalWeight / 2.0), 1.0);

    auto pair = GDGG.run();
//    for (int i = 0; i < pair.first.size(); i++) {
//        std::cout << i << " " << pair.first[i] << "\n";
//    }
//    std::cout << pair.second << "\n";

    UndirectedFix UF(graph, 1.2 * ((double) graph.totalWeight / 2.0), 1.0);
    UF.getUndirectedPartition();

    return 0;
}
