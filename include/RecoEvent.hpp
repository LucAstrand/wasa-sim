#ifndef RECOEVENT_H
#define RECOEVENT_H

#include "Structures.hpp"
#include "Clustering.hpp"


struct RecoEvent {
    std::vector<Cluster> clusters;
    std::vector<ChargedCluster> chargedClusters;
    std::vector<ChargedObject> chargedObjects;
    TVector3 vertex;
    double EM_energy;
};

RecoEvent ReconstructEvent(
    std::vector<Hit>& hits,
    const std::vector<ChargedTrack>& chargedTracks,
    const TVector3& vertex
);

#endif