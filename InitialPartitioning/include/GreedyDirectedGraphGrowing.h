//
// Created by panagiotis on 10/12/2024.
//

#ifndef DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H
#define DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H


#include "Graph.h"

class GreedyDirectedGraphGrowing {
public:
    Graph workingGraph;

    explicit GreedyDirectedGraphGrowing(Graph &graph) : workingGraph(graph) {}
};


#endif //DAG_PARTITIONING_GREEDYDIRECTEDGRAPHGROWING_H
