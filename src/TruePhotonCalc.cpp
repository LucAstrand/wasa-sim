#include "TruePhotonCalc.hpp"

std::vector<TruePhoton> TruePhotonBuilder(
    const std::vector<TruePhotonHit> &hits,
    const TVector3 &vertex)
{
    std::vector<TruePhoton> outPhotons;

    // Sanity check
    // if (hits.size() != 2) {
    //     // std::cerr << "[TruePhotonBuilder] Expected exactly 2 hits, got "
    //     //           << hits.size() << std::endl;
    //     std::cout << "[TruePhotonBuilder] Expected exactly 2 hits, got "
    //               << hits.size() << std::endl;
    //     return outPhotons;
    // }

    for (const auto &hit : hits) {
        // Direction from vertex to hit position
        TVector3 dir(hit.x - vertex.X(),
                     hit.y - vertex.Y(),
                     hit.z - vertex.Z());
        dir = dir.Unit();

        // Photon 4-vector (massless): |p| = E
        TLorentzVector p4;
        p4.SetVect(dir * hit.e);
        p4.SetE(hit.e);

        outPhotons.push_back({dir, p4, hit.truePhotonParentID});
    }

    return outPhotons;
}

std::vector<TruePi0> TruePi0Builder(const std::vector<TruePhoton>& photons)
{
    std::unordered_map<int, TruePi0> pi0map;

    for (const auto& g : photons) {
        // if (g.parentPdg != 111) continue; // by construction the photons are from pi0s (from GEANT4)

        std::cout << "[TruePi0Builder] photon parentID: " << g.parentID << std::endl;

        auto& pi0 = pi0map[g.parentID];
        pi0.trackID = g.parentID;
        pi0.photons.push_back(&g);
        pi0.p4 += g.p4;
    }
    std::vector<TruePi0> out;
    for (auto& [id, pi0] : pi0map) {
        if (pi0.photons.size() == 2) { 
            out.push_back(std::move(pi0));
        }
    }
    return out;
}
