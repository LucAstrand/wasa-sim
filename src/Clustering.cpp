#include "Clustering.hpp"


// =========================================================
//   SLIDING-WINDOW CLUSTERING (eta-phi towers)
// =========================================================
std::vector<Cluster> SlidingWindowClusterHits(
    std::vector<Hit> &hits,
    const TVector3 &vertex,
    double dEta,
    double dPhi,
    double E_seed, 
    double E_neighbor,
    int winSize)
{
    std::map<EtaPhiTowerKey, double> towers;
    for (size_t i=0; i<hits.size(); ++i) {
        if (hits[i].owner != HitOwner::None) continue; // Skip already-used hits
        double eta = calcEta(hits[i].x, hits[i].y, hits[i].z);
        double phi = calcPhi(hits[i].x, hits[i].y, hits[i].z);
        int iEta = int(std::floor(eta/dEta));
        int iPhi = int(std::floor(phi/dPhi));
        towers[{iEta,iPhi}] += hits[i].e;
    }

    std::vector<EtaPhiTowerKey> seed_keys;
    for (auto &[key, E] : towers) {
        if (E < E_seed) continue;

        bool isMax = true;
        for (int dEtaIdx=-1; dEtaIdx<=1; ++dEtaIdx) {
            for (int dPhiIdx=-1; dPhiIdx<=1; ++dPhiIdx) {
                if (dEtaIdx==0 && dPhiIdx==0) continue;
                EtaPhiTowerKey neigh{key.iEta+dEtaIdx, key.iPhi+dPhiIdx};
                if (towers.count(neigh) && towers[neigh] > E) {
                    isMax = false;
                }
            }
        }
        if (isMax) seed_keys.push_back(key);
    }

    struct SeedInfo {
        EtaPhiTowerKey key;
        double E;
    };
    std::vector<SeedInfo> seeds;
    for (auto &k : seed_keys) seeds.push_back({k, towers[k]});
    std::sort(seeds.begin(), seeds.end(), [](const SeedInfo &a, const SeedInfo &b) {
        return a.E > b.E;
    });

    std::vector<Cluster> clusters;
    std::set<EtaPhiTowerKey> assigned;
    for (auto &seed : seeds) {
        Cluster c;
        double sumE = 0, cx = 0, cy = 0, cz = 0;
        TVector3 mom(0,0,0);
        int nTowersInWindow = 0;
        int nHitsInWindow = 0;
        bool added = false;

        for (int de = -winSize; de <= winSize; ++de) {
            for (int dp = -winSize; dp <= winSize; ++dp) {
                EtaPhiTowerKey key{seed.key.iEta + de, seed.key.iPhi + dp};
                if (towers.count(key) == 0 || assigned.count(key) || towers[key] < E_neighbor) continue;

                assigned.insert(key);
                added = true;
                nTowersInWindow++;

                // Assign hits in this tower to cluster
                for (size_t i = 0; i < hits.size(); ++i) {                    
                    double eta = calcEta(hits[i].x, hits[i].y, hits[i].z);
                    double phi = calcPhi(hits[i].x, hits[i].y, hits[i].z);
                    int iEta = int(std::floor(eta / dEta));
                    int iPhi = int(std::floor(phi / dPhi));
                    if (iEta == key.iEta && iPhi == key.iPhi) {
                        auto &h = hits[i];
                        if (h.owner != HitOwner::None) continue;
                        sumE += h.e;
                        cx += h.x * h.e;
                        cy += h.y * h.e;
                        cz += h.z * h.e;
                        TVector3 v(h.x - vertex.X(), h.y - vertex.Y(), h.z - vertex.Z());
                        if (v.Mag2() > 1e-12) mom += h.e * v.Unit();
                        c.hits.push_back(&hits[i]);
                        h.owner = HitOwner::Neutral; // Mark used
                        nHitsInWindow++;
                    }
                }
            }
        }

        if (added && sumE > 0) {
            c.centroid = TVector3(cx / sumE, cy / sumE, cz / sumE);
            c.p4.SetPxPyPzE(mom.X(), mom.Y(), mom.Z(), sumE);
            clusters.push_back(c);
        }
    }
    return clusters;
}

// =========================================================
//   Neutral CLUSTERING (Angular acceptance based)
// =========================================================

void finalizeNeutralCluster(Cluster& cl, const TVector3& vtx) {
    double E_sum = 0.0;
    TVector3 centroid(0,0,0);

    for (auto* h : cl.hits) {
        centroid += h->e * TVector3(h->x, h->y, h->z);
        E_sum += h->e;
    }

    if (E_sum > 0) centroid *= (1.0 / E_sum);
    cl.centroid = centroid;

    TVector3 dir = (centroid - vtx).Unit();

    cl.p4.SetPxPyPzE(E_sum * dir.X(), E_sum * dir.Y(), E_sum * dir.Z(), E_sum);
}

// std::vector<Cluster> clusterNeutralHits(std::vector<Hit>& hits, const TVector3& vtx, double theta_max) {
//     // Sort hits by descending energy
//     std::sort(hits.begin(), hits.end(), 
//               [](const Hit& a, const Hit& b){ return a.e > b.e; });

//     std::vector<Cluster> clusters;
//     std::vector<bool> used(hits.size(), false);

//     for (size_t i = 0; i < hits.size(); ++i) {
//         if (used[i]) continue;

//         // --- Seed cluster ---
//         Cluster cl;
//         Hit& seed = hits[i];
//         used[i] = true;
//         cl.hits.push_back(&seed);

//         TVector3 seedDir = hitDirection(seed, vtx).Unit();

//         // --- Grow cluster ---
//         for (size_t j = i + 1; j < hits.size(); ++j) {
//             if (used[j]) continue;

//             Hit& h = hits[j];
//             if (h.owner != HitOwner::None) continue;
//             TVector3 hDir = hitDirection(h, vtx).Unit();

//             double angle = seedDir.Angle(hDir);
//             if (angle < theta_max) {
//                 cl.hits.push_back(&h);
//                 used[j] = true;
//                 h.owner = HitOwner::Neutral;
//             }
//         }

//         // --- Finalize cluster ---
//         finalizeNeutralCluster(cl, vtx);
//         clusters.push_back(cl);
//     }

//     return clusters;
// }

std::vector<Cluster> clusterNeutralHits(std::vector<Hit>& hits,
                                        const TVector3& vtx,
                                        double theta_max)
{
    // --- Sort hit indices by descending energy (do NOT reorder hits) ---
    std::vector<size_t> order(hits.size());
    std::iota(order.begin(), order.end(), 0);

    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return hits[a].e > hits[b].e; });

    std::vector<Cluster> clusters;
    std::vector<bool> used(hits.size(), false);

    for (size_t oi = 0; oi < order.size(); ++oi) {
        const size_t i = order[oi];
        if (used[i]) continue;
        if (hits[i].owner != HitOwner::None) continue;

        // --- Seed cluster ---
        Cluster cl;
        Hit& seed = hits[i];
        used[i] = true;
        cl.hits.push_back(&seed);

        TVector3 seedDir = hitDirection(seed, vtx).Unit();

        // --- Grow cluster ---
        for (size_t oj = oi + 1; oj < order.size(); ++oj) {
            const size_t j = order[oj];
            if (used[j]) continue;

            Hit& h = hits[j];
            if (h.owner != HitOwner::None) continue;

            TVector3 hDir = hitDirection(h, vtx).Unit();
            double angle = seedDir.Angle(hDir);

            if (angle < theta_max) {
                cl.hits.push_back(&h);
                used[j] = true;
                h.owner = HitOwner::Neutral;
            }
        }

        // --- Finalize cluster ---
        finalizeNeutralCluster(cl, vtx);
        clusters.push_back(std::move(cl));
    }

    return clusters;
}


// =========================================================
//   CHARGED OBJECT CLUSTERING (Angular acceptance based)
// =========================================================

std::vector<ChargedCluster> MatchHitsToTracks(
    const std::vector<ChargedTrack>& tracks,
    std::vector<Hit>& hits,
    double thetaMax
) {
    std::vector<ChargedCluster> clusters;
    PotentialGas tpcGas = PotentialGas::eArCO2_8020;
    for (size_t i = 0; i < tracks.size(); ++i) {
        const auto& trk = tracks[i];

        if (std::abs(trk.TruePDG) == 11) std::cout << "electron Track!" << std::endl;

        ChargedCluster c;
        c.trackID = trk.id;
        c.direction = trk.direction.Unit();
        c.TPCExitPoint = trk.exitPoint;
        c.objectTrueKE = trk.TrueKE;
        c.objectTruePDG = trk.TruePDG;
        c.objectTruedEdx = BetheBloch(trk.TruePDG, trk.TrueKE, tpcGas);
        c.totalEnergy = 0.0;
        c.clusterdEdx = trk.clusterdEdx;
        c.EdepSmeared = trk.EdepSmeared;
        c.nSigmaPion = nSigmaCalc(trk.EdepSmeared, trk.pathLength, BetheBloch(211, trk.TrueKE, tpcGas), trk.resolution);
        c.nSigmaProton = nSigmaCalc(trk.EdepSmeared, trk.pathLength, BetheBloch(2212, trk.TrueKE, tpcGas), trk.resolution);
        c.nSigmaElectron = nSigmaCalc(trk.EdepSmeared, trk.pathLength, BetheBloch(11, trk.TrueKE, tpcGas), trk.resolution);
        c.pidL = ComputePIDLikelihoods(
            c.nSigmaPion,
            c.nSigmaElectron,
            c.nSigmaProton

        );
        c.pidGuess = AssignPIDFromLikelihood(c.pidL, 0.7);

        double dirMag = trk.direction.Mag();

        clusters.push_back(c);
    }

    int nMatchedHits = 0;

    for (size_t h = 0; h < hits.size(); ++h) {
        auto& hit = hits[h];

        if (hit.owner != HitOwner::None)
            continue;

        TVector3 hitPos(hit.x, hit.y, hit.z);
        double bestAngle = std::numeric_limits<double>::max();
        double bestDot   = -999.0;
        int bestTrack    = -1;

        for (size_t i = 0; i < tracks.size(); ++i) {
            TVector3 delta = hitPos - tracks[i].exitPoint; //hitPos - tracks[i].vertex;
            // double deltaMag = delta.Mag();
            if (delta.Mag2() < 1e-12) {
                bestAngle = 0.0;
                bestTrack = static_cast<int>(i);
                break; // can't do better than 0
            }
            TVector3 hitDir = delta.Unit();
            double dot = tracks[i].direction.Unit().Dot(hitDir);
            dot = std::clamp(dot, -1.0, 1.0);
            if (dot < 0) continue;
            double theta = std::acos(dot);

            if (theta < bestAngle) {
                bestAngle = theta;
                bestTrack = static_cast<int>(i);
                bestDot   = dot;
            }
        }

        if (bestTrack >= 0 && bestAngle < thetaMax) {
            clusters[bestTrack].hits.push_back(&hit);
            clusters[bestTrack].totalEnergy += hit.e;
            hit.owner = HitOwner::Charged;
            ++nMatchedHits;
        }
    }
    return clusters;
}
