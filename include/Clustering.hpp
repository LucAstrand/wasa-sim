#ifndef CLUSTERING_H
#define CLUSTERING_H

#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>

#include "TVector3.h"

#include "Structures.hpp"
#include "Utils.hpp"
#include "ParticleID.hpp"
#include "DEDXTable.hpp"



void finalizeNeutralCluster(Cluster& cl, const TVector3 &vtx);

std::vector<Cluster> clusterNeutralHits(std::vector<Hit>& hits, const TVector3& vtx, double theta_max);


std::vector<Cluster> SlidingWindowClusterHits(
    std::vector<Hit> &hits, // Not const anymore -> We can assign ownership!
    const TVector3 &vertex,
    double dEta,
    double dPhi,
    double E_seed, 
    double E_neighbor,
    int winSize = 7
);

void MatchHitsElectron(ChargedCluster& cluster,
    std::vector<Hit>& hits,
    double moliereRadius
);
                       
void MatchHitsHadron(ChargedCluster& cluster,
    std::vector<Hit>& hits,
    double thetaMax    
);

void MatchHitsMuon(ChargedCluster& cluster,
    std::vector<Hit>& hits,
    double thetaMax
);


std::vector<ChargedCluster> MatchHitsToTracks(
    const std::vector<ChargedTrack>& tracks,
    std::vector<Hit>& hits, // Not const anymore -> We can assign ownership!
    double thetaMax, 
    const DEDXTable& dedxTable
);

std::vector<ConversionCandidate> FindConversions(
    const std::vector<ChargedCluster>& clusters,
    const TVector3& primaryVertex
);

#endif
