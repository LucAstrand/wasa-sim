#ifndef PHOTONMATCH_H
#define PHOTONMATCH_H

#include <vector>

#include "Structures.hpp"

std::vector<int> matchClustersToTruth(
    const std::vector<Cluster>& clusters,
    const std::vector<TruePhoton>& truePhotons,
    double angleMaxRad 
);

#endif