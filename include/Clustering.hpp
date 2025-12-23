#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <algorithm>

#include "TVector3.h"

#include "Structures.hpp"
#include "Utils.hpp"
#include "ParticleID.hpp"


std::vector<Cluster> SlidingWindowClusterHits(
    const std::vector<Hit> &hits,
    const TVector3 &vertex,
    double dEta,
    double dPhi,
    double E_seed, 
    double E_neighbor,
    int winSize = 7
);

std::vector<Cluster> SplitMergedClusterEM(
    const Cluster &merged,
    const std::vector<Hit> &hits,
    const TVector3 &vertex,
    double dEta,
    double dPhi,
    int maxIter = 50,
    double tol = 1e-4,
    double minFrac = 0.02,
    double initSigma = -1
);

std::vector<ChargedCluster> MatchHitsToTracks(
    const std::vector<ChargedTrack>& tracks,
    const std::vector<Hit>& hits,
    double thetaMax
);


#endif
