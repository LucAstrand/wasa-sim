#include "Clustering.hpp"

std::vector<Cluster> SplitMergedClusterEM(
    const Cluster &merged,
    const std::vector<Hit> &hits,
    const TVector3 &vertex,
    double dEta,           // tower size in eta
    double dPhi,           // tower size in phi
    int maxIter,
    double tol,     // relative change tolerance on log-likelihood
    double minFrac, // minimum fractional energy for a component to be considered real
    double initSigma  // initial Gaussian sigma in eta-phi units; if <0 use ~1.0 * max(dEta,dPhi)
)
{
    // 1) build per-hit eta/phi list and total energy
    struct HInfo { int idx; double eta, phi, E; };
    std::vector<HInfo> local;
    double totalE = 0;
    for (auto hitIdx : merged.hitIndices) {
        const auto &h = hits[hitIdx];
        double eta = calcEta(h.x, h.y, h.z);
        double phi = calcPhi(h.x, h.y, h.z);
        local.push_back({hitIdx, eta, phi, h.e});
        totalE += h.e;
    }
    if (totalE <= 0 || local.size() < 2) return { merged };

    // 2) initial seeds: use two highest-energy towers inside merged cluster
    // Build tower map (fine grid using dEta,dPhi)
    std::map<EtaPhiTowerKey, double> towers;
    for (auto &hi : local) {
        int iEta = int(std::floor(hi.eta / dEta));
        int iPhi = int(std::floor(hi.phi / dPhi));
        towers[{iEta, iPhi}] += hi.E;
    }
    // sort towers by energy
    std::vector<std::pair<EtaPhiTowerKey,double>> sorted;
    for (auto &kv : towers) sorted.push_back(kv);
    std::sort(sorted.begin(), sorted.end(),
              [](auto &a, auto &b){ return a.second > b.second; });
    if (sorted.empty()) return { merged };
    // pick seed1 = highest tower center, seed2 = best next tower > 1 cell away if possible
    auto towerKeyToEtaPhi = [&](const EtaPhiTowerKey &k){
        return std::pair<double,double>( (k.iEta + 0.5)*dEta, (k.iPhi + 0.5)*dPhi );
    };
    auto [eta1, phi1] = towerKeyToEtaPhi(sorted[0].first);
    double eta2=eta1, phi2=phi1;
    bool found2 = false;
    for (size_t i=1; i<sorted.size(); ++i) {
        double dEtaCells = std::abs(sorted[i].first.iEta - sorted[0].first.iEta);
        double dPhiCells = std::abs(sorted[i].first.iPhi - sorted[0].first.iPhi);
        if (dEtaCells + dPhiCells > 1) { // not immediate neighbor
            auto p = towerKeyToEtaPhi(sorted[i].first);
            eta2 = p.first; phi2 = p.second;
            found2 = true;
            break;
        }
    }
    if (!found2) {
        // fallback: pick the next-highest tower regardless
        if (sorted.size() >= 2) {
            auto p = towerKeyToEtaPhi(sorted[1].first);
            eta2 = p.first; phi2 = p.second;
        } else {
            // can't find a second seed
            return { merged };
        }
    }

    // 3) EM initialization
    double sigma = (initSigma > 0) ? initSigma : std::max(dEta, dPhi) * 1.0; // tuneable
    double alpha1 = 0.5, alpha2 = 0.5; // mixture weights (by energy fraction)
    // optionally seed alpha by tower energies if available
    double seedE1 = sorted[0].second;
    double seedE2 = (sorted.size() >= 2) ? sorted[1].second : (totalE - seedE1);
    if (seedE1 + seedE2 > 0) {
        alpha1 = seedE1 / (seedE1 + seedE2);
        alpha2 = seedE2 / (seedE1 + seedE2);
    }

    std::vector<double> w1(local.size(), 0.0), w2(local.size(), 0.0);
    double prevLL = -std::numeric_limits<double>::infinity();

    const double eps = 1e-12;

    for (int iter=0; iter<maxIter; ++iter) {
        // E-step: compute responsibilities (energy-weighted)
        double ll = 0.0; // log-likelihood (energy-weighted)
        for (size_t i=0; i<local.size(); ++i) {
            double eta = local[i].eta, phi = local[i].phi, E = local[i].E;
            double d1_eta = (eta - eta1);
            double d1_phi = dphi_wrap(phi, phi1);
            double d2_eta = (eta - eta2);
            double d2_phi = dphi_wrap(phi, phi2);
            double r1 = std::exp(-0.5 * (d1_eta*d1_eta + d1_phi*d1_phi) / (sigma*sigma));
            double r2 = std::exp(-0.5 * (d2_eta*d2_eta + d2_phi*d2_phi) / (sigma*sigma));
            double p1 = alpha1 * r1;
            double p2 = alpha2 * r2;
            double norm = p1 + p2 + eps;
            w1[i] = (p1 / norm);
            w2[i] = (p2 / norm);
            // accumulate energy-weighted log-likelihood for monitoring convergence
            ll += E * std::log(norm);
        }

        // M-step: update centroids and alphas (energy-weighted)
        double sumW1 = eps, sumW2 = eps; // energy-weighted sums
        double eta1_num = 0, phi1_num_x = 0, phi1_num_y = 0; // we'll avoid circular phi averaging by vector method
        double eta2_num = 0, phi2_num_x = 0, phi2_num_y = 0;
        for (size_t i=0; i<local.size(); ++i) {
            double E = local[i].E;
            sumW1 += w1[i] * E;
            sumW2 += w2[i] * E;
            eta1_num += w1[i] * E * local[i].eta;
            eta2_num += w2[i] * E * local[i].eta;
            // phi average using unit vectors to handle wrap
            phi1_num_x += w1[i] * E * std::cos(local[i].phi);
            phi1_num_y += w1[i] * E * std::sin(local[i].phi);
            phi2_num_x += w2[i] * E * std::cos(local[i].phi);
            phi2_num_y += w2[i] * E * std::sin(local[i].phi);
        }

        // update mixture weights
        alpha1 = sumW1 / (sumW1 + sumW2);
        alpha2 = sumW2 / (sumW1 + sumW2);

        // update centroids
        double new_eta1 = eta1_num / sumW1;
        double new_eta2 = eta2_num / sumW2;
        double new_phi1 = std::atan2(phi1_num_y, phi1_num_x);
        double new_phi2 = std::atan2(phi2_num_y, phi2_num_x);

        // small safeguard: if centroids collapse too close, we may abort
        double dEta_cent = new_eta1 - new_eta2;
        double dPhi_cent = dphi_wrap(new_phi1, new_phi2);
        double centDist2 = dEta_cent*dEta_cent + dPhi_cent*dPhi_cent;

        // update
        eta1 = new_eta1; phi1 = new_phi1;
        eta2 = new_eta2; phi2 = new_phi2;

        // optional: update sigma from weighted variance (or keep fixed)
        // compute weighted squared distances to update sigma (could help convergence)
        double var_num = 0, var_den = eps;
        for (size_t i=0; i<local.size(); ++i) {
            double E = local[i].E;
            double de1 = local[i].eta - eta1, dp1 = dphi_wrap(local[i].phi, phi1);
            double de2 = local[i].eta - eta2, dp2 = dphi_wrap(local[i].phi, phi2);
            double r1sq = (de1*de1 + dp1*dp1);
            double r2sq = (de2*de2 + dp2*dp2);
            // weight by responsibility
            var_num += E * (w1[i]*r1sq + w2[i]*r2sq);
            var_den += E;
        }
        double new_sigma = std::sqrt(std::max(var_num / var_den, 1e-6));
        // limit sigma updates to reasonable bounds:
        if (new_sigma > 5*std::max(dEta,dPhi)) new_sigma = 5*std::max(dEta,dPhi);
        if (new_sigma < 0.1*std::max(dEta,dPhi)) new_sigma = 0.1*std::max(dEta,dPhi);
        sigma = new_sigma;

        // check convergence in log-likelihood
        if (iter > 0) {
            double rel = std::abs((ll - prevLL) / (std::abs(prevLL) + 1e-12));
            if (rel < tol) break;
        }
        prevLL = ll;

        // If one component becomes too small in energy fraction -> break early
        if (alpha1 < minFrac || alpha2 < minFrac) break;

        // If centroids get extremely close compared to cell size, splitting unlikely
        if (centDist2 < 1e-6) break;
    } // end EM loop

    // Post-fit decisions: accept split only if both components have reasonable energy
    if (alpha1 < minFrac || alpha2 < minFrac) {
        // one component too small -> no reliable split
        return { merged };
    }

    // Build two clusters: compute energy share and momentum etc. using fractional weights
    std::vector<Cluster> out;
    // accumulate
    double Esum1 = 0, Esum2 = 0;
    double cx1=0, cy1=0, cz1=0, cx2=0, cy2=0, cz2=0;
    TVector3 mom1(0,0,0), mom2(0,0,0);
    std::vector<int> indices1, indices2;

    for (size_t i=0; i<local.size(); ++i) {
        const auto &hi = local[i];
        double frac1 = w1[i];
        double frac2 = w2[i];
        double e1 = hi.E * frac1;
        double e2 = hi.E * frac2;
        Esum1 += e1; Esum2 += e2;
        const auto &h = hits[hi.idx];
        cx1 += h.x * e1; cy1 += h.y * e1; cz1 += h.z * e1;
        cx2 += h.x * e2; cy2 += h.y * e2; cz2 += h.z * e2;
        TVector3 v1(h.x - vertex.X(), h.y - vertex.Y(), h.z - vertex.Z());
        TVector3 v2 = v1;
        if (v1.Mag2() > 1e-12) {
            mom1 += e1 * v1.Unit();
            mom2 += e2 * v2.Unit();
        }
        // assign hit index to the dominant fraction for bookkeeping
        if (frac1 >= frac2) indices1.push_back(hi.idx);
        if (frac2 >= frac1) indices2.push_back(hi.idx);
    }

    double tot = Esum1 + Esum2;
    if (tot <= 0) return { merged };

    // apply a final quality cut: both must have at least minFrac fraction of the merged energy
    if ( (Esum1 / tot) < minFrac || (Esum2 / tot) < minFrac ) {
        return { merged };
    }

    // create output clusters
    Cluster c1, c2;
    c1.centroid = TVector3(cx1 / Esum1, cy1 / Esum1, cz1 / Esum1);
    c1.p4.SetPxPyPzE(mom1.X(), mom1.Y(), mom1.Z(), Esum1);
    c1.hitIndices = indices1;

    c2.centroid = TVector3(cx2 / Esum2, cy2 / Esum2, cz2 / Esum2);
    c2.p4.SetPxPyPzE(mom2.X(), mom2.Y(), mom2.Z(), Esum2);
    c2.hitIndices = indices2;

    return { c1, c2 };
}




// =========================================================
//   SLIDING-WINDOW CLUSTERING (eta-phi towers)
// =========================================================
std::vector<Cluster> SlidingWindowClusterHits(
    const std::vector<Hit> &hits,
    const TVector3 &vertex,
    double dEta,
    double dPhi,
    double E_seed, 
    double E_neighbor,
    int winSize)
{
    // 1. Build towers
    std::map<EtaPhiTowerKey, double> towers;
    for (size_t i=0; i<hits.size(); ++i) {
        double eta = calcEta(hits[i].x, hits[i].y, hits[i].z);
        double phi = calcPhi(hits[i].x, hits[i].y, hits[i].z);
        int iEta = int(std::floor(eta/dEta));
        int iPhi = int(std::floor(phi/dPhi));
        towers[{iEta,iPhi}] += hits[i].e;
    }

    // 2. Find seeds (local maxima)
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

    // Sort seeds by energy descending
    struct SeedInfo {
        EtaPhiTowerKey key;
        double E;
    };
    std::vector<SeedInfo> seeds;
    for (auto &k : seed_keys) seeds.push_back({k, towers[k]});
    std::sort(seeds.begin(), seeds.end(), [](const SeedInfo &a, const SeedInfo &b) {
        return a.E > b.E;
    });

    // 3. Cluster filling around seeds (exclusive tower assignment)
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
                        const auto &h = hits[i];
                        sumE += h.e;
                        cx += h.x * h.e;
                        cy += h.y * h.e;
                        cz += h.z * h.e;
                        TVector3 v(h.x - vertex.X(), h.y - vertex.Y(), h.z - vertex.Z());
                        if (v.Mag2() > 1e-12) mom += h.e * v.Unit();
                        c.hitIndices.push_back(i);
                        nHitsInWindow++;
                    }
                }
            }
        }

        // --- DEBUG PRINT ---
        // std::cout << "Seed at (iEta,iPhi)=(" << seed.key.iEta << "," << seed.key.iPhi << ") "
        //           << "window covers " << nTowersInWindow << " towers, "
        //           << nHitsInWindow << " hits, "
        //           << "cluster E=" << sumE << std::endl;

        if (added && sumE > 0) {
            c.centroid = TVector3(cx / sumE, cy / sumE, cz / sumE);
            c.p4.SetPxPyPzE(mom.X(), mom.Y(), mom.Z(), sumE);
            clusters.push_back(c);
        }
    }
    return clusters;
}

// =========================================================
//   CHARGED OBJECT CLUSTERING (Angular acceptance based)
// =========================================================


std::vector<ChargedCluster> MatchHitsToTracks(
    const std::vector<ChargedTrack>& tracks,
    const std::vector<Hit>& hits,
    double thetaMax
) {
    std::vector<ChargedCluster> clusters;

    // Initialize one cluster per track
    for (const auto& trk : tracks) {
        ChargedCluster c;
        c.trackID   = trk.id;
        c.direction = trk.direction;
        c.clusterdEdx = trk.EdepSmeared / trk.pathLength;
        c.nSigma = nSigmaCalc(trk.EdepSmeared, trk.pathLength, trk.dEdxTheory, trk.resolution);
        clusters.push_back(c);
    }

    for (const auto& hit : hits) {
        TVector3 hitPos(hit.x, hit.y, hit.z);

        double bestAngle = std::numeric_limits<double>::max();
        int bestTrack = -1;

        for (size_t i = 0; i < tracks.size(); ++i) {
            TVector3 hitDir = (hitPos - tracks[i].vertex).Unit();
            double cosTheta = tracks[i].direction.Dot(hitDir);
            cosTheta = std::clamp(cosTheta, -1.0, 1.0);
            double theta = std::acos(cosTheta);

            if (theta < bestAngle) {
                bestAngle = theta;
                bestTrack = i;
            }
        }

        if (bestTrack >= 0 && bestAngle < thetaMax) {
            clusters[bestTrack].hits.push_back(hit);
            clusters[bestTrack].totalEnergy += hit.e;
        }
    }

    return clusters;
}