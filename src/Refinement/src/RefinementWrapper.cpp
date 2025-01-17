//
// Created by panagiotis on 17/1/2025.
//

#include "RefinementWrapper.h"
#include "BoundaryFM.h"
#include "BoundaryKL.h"
#include "Mixed.h"

std::pair<uint64_t, uint64_t>
calculatePartSizes(const Graph &graph, std::vector<uint8_t> &bisectionInfo) {
    uint64_t sizeV0 = 0, sizeV1 = 0;

    // Sum weights for each partition
    for (uint64_t i = 0; i < bisectionInfo.size(); ++i) {
        if (bisectionInfo[i] == 0)
            sizeV0 += graph.nodeWeights[i];
        else
            sizeV1 += graph.nodeWeights[i];
    }
    return {sizeV0, sizeV1};
}

bool checkBalance(const Graph &graph, std::vector<uint8_t> &bisectionInfo, double upperBoundPartWeight,
                  double lowerBoundPartWeight, uint64_t maxNodeWeight) {
    // Get current partition weights
    auto [sizeV0, sizeV1] = calculatePartSizes(graph, bisectionInfo);

    // Check balance constraints with allowance for heaviest node
    if ((double) sizeV0 < lowerBoundPartWeight - (double) maxNodeWeight ||
        (double) sizeV0 > upperBoundPartWeight + (double) maxNodeWeight ||
        (double) sizeV1 < lowerBoundPartWeight - (double) maxNodeWeight ||
        (double) sizeV1 > upperBoundPartWeight + (double) maxNodeWeight)
        return false;

    return true;
}

void refinementWrapper(const Graph &graph, std::vector<uint8_t> &bisectionInfo, uint64_t &edgeCut,
                       RefinementMethod refinementMethod, uint64_t refinementPasses, double upperBoundPartWeight,
                       double lowerBoundPartWeight) {
    if (!checkBalance(graph, bisectionInfo, upperBoundPartWeight, lowerBoundPartWeight, 0)) {
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
