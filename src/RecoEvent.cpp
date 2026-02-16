#include "RecoEvent.hpp"


RecoEvent ReconstructEvent(
    std::vector<Hit>& hits,
    const std::vector<ChargedTrack>& chargedTracks,
    const TVector3& vertex
) {
    RecoEvent evt;
    evt.vertex = vertex;

    evt.chargedClusters = MatchHitsToTracks(chargedTracks, hits, 25*TMath::DegToRad());

    // std::cout << "number of charged Clusters: " << evt.chargedClusters.size() << std::endl;

    evt.EM_energy = 0.0;
    // evt.clusters = clusterNeutralHits(hits, vertex, 25*TMath::DegToRad()); 
    evt.clusters = clusterNeutralHits(hits, vertex, 0.12 /* rad */); //og 0.2
    evt.clusters.erase(std::remove_if(evt.clusters.begin(), evt.clusters.end(),
                                [](const Cluster &c){ return c.p4.E() < 20.0; }),
                evt.clusters.end());    
    
    // This would be only photons
    // for (const auto& c : evt.clusters) {
    //     evt.EM_energy += c.p4.E(); // only the neutral cluster energy for now???
    // }
    for (const auto& h : hits) {
        evt.EM_energy += h.e; 
    }
    return evt;
}