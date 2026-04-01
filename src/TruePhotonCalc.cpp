#include "TruePhotonCalc.hpp"

std::vector<TruePhoton> TruePhotonBuilder(
    const std::vector<TruePhotonHit> &hits,
    const TVector3 &vertex)
{
    std::vector<TruePhoton> outPhotons;

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
