#include "RecoEvent.hpp"


RecoEvent ReconstructEvent(
    std::vector<Hit>& hits,
    const std::vector<ChargedTrack>& chargedTracks,
    const TVector3& vertex, 
    const DEDXTable& dedxTable
) {
    RecoEvent evt;
    evt.vertex = vertex;

    // evt.chargedClusters = MatchHitsToTracks(chargedTracks, hits, 25*TMath::DegToRad());
    evt.chargedClusters = MatchHitsToTracks(chargedTracks, hits, 35*TMath::DegToRad(), dedxTable);

    // std::cout << "[ReconstructEvent] calling neutral clustering..." << std::endl;

    // std::cout << "number of charged Clusters: " << evt.chargedClusters.size() << std::endl;

    evt.EM_energy = 0.0;
    // evt.clusters = clusterNeutralHits(hits, vertex, 25*TMath::DegToRad()); 
    // evt.clusters = clusterNeutralHits(hits, vertex, 0.12 /* rad */); //og 0.2

    double dEta = 0.10/2;
    double dPhi = 0.10/2;
    double E_seed = 15.00;
    double E_neighbor = 0.03;
    int winSize = 7;

    evt.clusters = SlidingWindowClusterHits(hits, vertex, dEta, dPhi, E_seed, E_neighbor, winSize);

//     evt.clusters.erase(std::remove_if(evt.clusters.begin(), evt.clusters.end(),
//                                 [](const Cluster &c){ return c.p4.E() < 20.0; }),
//                 evt.clusters.end());    
    
//     // This would be only photons
//     // for (const auto& c : evt.clusters) {
//     //     evt.EM_energy += c.p4.E(); // only the neutral cluster energy for now???
//     // }
//     for (const auto& h : hits) {
//         evt.EM_energy += h.e; 
//     }
//     return evt;
// }

evt.clusters.erase(std::remove_if(evt.clusters.begin(), evt.clusters.end(),
                            [](const Cluster &c){ return c.p4.E() < 20.0; }),
            evt.clusters.end());


for (const auto& h : hits) {
    evt.EM_energy += h.e;
}

return evt;

}