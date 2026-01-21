#include "TruePhotonCalc.hpp"

std::vector<TruePhoton> TruePhotonBuilder(
    const std::vector<Hit> &hits,
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

        outPhotons.push_back({dir, p4});
    }

    return outPhotons;
}