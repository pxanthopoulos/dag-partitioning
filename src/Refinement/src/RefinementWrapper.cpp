//
// Created by panagiotis on 17/1/2025.
//

#include "RefinementWrapper.h"
#include "BoundaryFM.h"
#include "BoundaryKL.h"
#include "Mixed.h"

void refinementWrapper(const Graph &graph, std::vector<uint8_t> &bisectionInfo, uint64_t &edgeCut,
                       RefinementMethod refinementMethod, uint64_t refinementPasses, double upperBoundPartWeight,
                       double lowerBoundPartWeight) {
    if (!Refinement::checkBalance(bisectionInfo, graph, 0, upperBoundPartWeight, lowerBoundPartWeight)) {
        BoundaryFM boundaryFM = BoundaryFM(graph, bisectionInfo, edgeCut, 1,
                                           upperBoundPartWeight,
                                           lowerBoundPartWeight);
        boundaryFM.run();
    }

    std::unique_ptr<Refinement> refinement;
    switch (refinementMethod) {
        case RefinementMethod::BOUNDARYFM:
            refinement = std::make_unique<BoundaryFM>(graph, bisectionInfo, edgeCut, refinementPasses,
                                                      upperBoundPartWeight, lowerBoundPartWeight);
            break;
        case RefinementMethod::BOUNDARYKL:
            refinement = std::make_unique<BoundaryKL>(graph, bisectionInfo, edgeCut, refinementPasses,
                                                      upperBoundPartWeight, lowerBoundPartWeight);
            break;
        case RefinementMethod::MIXED:
            refinement = std::make_unique<Mixed>(graph, bisectionInfo, edgeCut, refinementPasses, upperBoundPartWeight,
                                                 lowerBoundPartWeight);
            break;
        default:
            throw std::invalid_argument("Unknown refinement type");
    }
    refinement->run();
}
