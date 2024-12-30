#include "Graph.h"
#include "MultilevelBisectioner.h"
#include <iostream>
#include <chrono>

int main() {
//    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/2mm_10_20_30_40.dot",
//                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");
    Graph graph = readDotFile("/home/panagiotis/code/dag-partitioning/test/dag.dot",
                              "/home/panagiotis/code/dag-partitioning/node-mappings.txt");

    MultilevelBisectioner mb_forb_undir(graph, ClusteringMethod::FORB, 20, 1, BisectionMethod::UNDIRBOTH, 1.2,
                                        RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_cyc_undir(graph, ClusteringMethod::CYC, 20, 1, BisectionMethod::UNDIRBOTH, 1.2,
                                       RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_hyb_undir(graph, ClusteringMethod::HYB, 20, 1, BisectionMethod::UNDIRBOTH, 1.2,
                                       RefinementMethod::BOUNDARYFM, 10);

    MultilevelBisectioner mb_forb_ggg(graph, ClusteringMethod::FORB, 20, 1, BisectionMethod::GGG, 1.2,
                                      RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_cyc_ggg(graph, ClusteringMethod::CYC, 20, 1, BisectionMethod::GGG, 1.2,
                                     RefinementMethod::BOUNDARYFM, 10);
    MultilevelBisectioner mb_hyb_ggg(graph, ClusteringMethod::HYB, 20, 1, BisectionMethod::GGG, 1.2,
                                     RefinementMethod::BOUNDARYFM, 10);

    auto vec = mb_forb_undir.run();
    vec = mb_cyc_undir.run();
    vec = mb_hyb_undir.run();

    vec = mb_forb_ggg.run();
    vec = mb_cyc_ggg.run();
    vec = mb_hyb_ggg.run();
    return 0;
}
