#include <vector>
#include <limits>

#include "TVector3.h"

#include "Structures.hpp"

// Inputs: clusters and truePhotons as before
std::vector<int> matchClustersToTruth(
    const std::vector<Cluster>& clusters,
    const std::vector<TruePhoton>& truePhotons,
    double angleMaxRad 
) {
    const size_t Nc = clusters.size();
    const size_t Nt = truePhotons.size();

    std::vector<int> clusterToTrueIdx(Nc, -1);
    std::vector<bool> trueUsed(Nt, false);

    for (size_t ic = 0; ic < Nc; ++ic) {
        double bestAngle = std::numeric_limits<double>::infinity();
        int bestIdx = -1;
        for (size_t it = 0; it < Nt; ++it) {
            if (trueUsed[it]) continue;             // enforce unique use
            double ang = clusters[ic].centroid.Angle(truePhotons[it].dir);
            if (ang < bestAngle) {
                bestAngle = ang;
                bestIdx = static_cast<int>(it);
            }
        }
        if (bestIdx >= 0 && bestAngle <= angleMaxRad) {
            clusterToTrueIdx[ic] = bestIdx;
            trueUsed[bestIdx] = true;
        }
    }
    return clusterToTrueIdx; // size Nc, -1 means unmatched
}
