//
// Created by panagiotis on 17/12/2024.
//

#include <iostream>
#include "UndirectedFix.h"
#include "scotch.h"

UndirectedFix::UndirectedFix(Graph &graph, double upperBoundPartWeight,
                             double lowerBoundPartWeight) : Bisection(graph,
                                                                      upperBoundPartWeight,
                                                                      lowerBoundPartWeight) {}

void UndirectedFix::getUndirectedPartition() {
    std::cout << SCOTCH_stratSizeof() << "\n";
}

std::pair<std::vector<bool>, uint64_t> UndirectedFix::run() {
    return {};
}
