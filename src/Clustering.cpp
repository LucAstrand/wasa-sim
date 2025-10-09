#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <set>
#include <TF1.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TVector3.h>

// Structure for hits
struct Hit {
    double x, y, z, e;
    int index;
};

// Structure for clusters (neutral objects)
struct Cluster {
    std::vector<int> hitIndices;
    TVector3 centroid;
    TLorentzVector p4;
};

double AngleBetweenVectors(const TVector3& v1, const TVector3& v2) {
    double dot = v1.Dot(v2);
    double mag1 = v1.Mag();
    double mag2 = v2.Mag();
    if (mag1 * mag2 == 0) return TMath::Pi();
    return std::acos(std::min(1.0, std::max(-1.0, dot / (mag1 * mag2))));
}

// std::vector<Cluster> MomentumBasedClustering(const std::vector<Hit>& hits, const TVector3& vertex, double phiThreshold = 0.2) {
//     std::vector<Hit> hSorted = hits;
//     std::sort(hSorted.begin(), hSorted.end(), [](const Hit& a, const Hit& b) { return a.e > b.e; });

//     std::vector<Cluster> clusters;
//     std::vector<bool> used(hits.size(), false);

//     while (!hSorted.empty()) {
//         Hit hi = hSorted[0];
//         Cluster c;
//         c.hitIndices.push_back(hi.index);
//         used[hi.index] = true;

//         TVector3 v_i(hi.x - vertex.X(), hi.y - vertex.Y(), hi.z - vertex.Z());
//         if (v_i.Mag2() == 0) continue;
//         v_i = v_i.Unit();

//         for (size_t j = 1; j < hSorted.size(); ++j) {
//             if (used[hSorted[j].index]) continue;
//             TVector3 v_j(hSorted[j].x - vertex.X(), hSorted[j].y - vertex.Y(), hSorted[j].z - vertex.Z());
//             if (v_j.Mag2() == 0) continue;
//             v_j = v_j.Unit();
//             double angle = AngleBetweenVectors(v_i, v_j);
//             if (angle < phiThreshold) {
//                 c.hitIndices.push_back(hSorted[j].index);
//                 used[hSorted[j].index] = true;
//             }
//         }

//         double sumE = 0.0, cx = 0.0, cy = 0.0, cz = 0.0;
//         for (int idx : c.hitIndices) {
//             const auto& h = hits[idx];
//             sumE += h.e;
//             cx += h.x * h.e;
//             cy += h.y * h.e;
//             cz += h.z * h.e;
//         }
//         if (sumE > 0) {
//             c.centroid = TVector3(cx / sumE, cy / sumE, cz / sumE);
//             TVector3 dir = (c.centroid - vertex).Unit();
//             c.p4.SetPxPyPzE(sumE * dir.X(), sumE * dir.Y(), sumE * dir.Z(), sumE);
//             clusters.push_back(c);
//         }

//         hSorted.erase(hSorted.begin());
//         hSorted.erase(std::remove_if(hSorted.begin(), hSorted.end(),
//                                      [&used](const Hit& h) { return used[h.index]; }),
//                       hSorted.end());
//     }

//     // Filter low-energy clusters
//     clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
//                                   [](const Cluster& c) { return c.p4.E() < 10.0; }),
//                    clusters.end());

//     return clusters;
// }

std::vector<Cluster> MomentumBasedClustering(const std::vector<Hit>& hits, const TVector3& vertex, double phiThreshold = 0.05) {
    std::vector<Hit> hSorted = hits;
    // Sort by energy descending
    std::sort(hSorted.begin(), hSorted.end(), [](const Hit& a, const Hit& b) { return a.e > b.e; });

    std::vector<Cluster> clusters;
    std::vector<bool> used(hits.size(), false);

    for (const auto& seed : hSorted) {
        if (used[seed.index]) continue;
        
        // Start new cluster with seed
        Cluster c;
        c.hitIndices.push_back(seed.index);
        used[seed.index] = true;

        TVector3 seed_dir(seed.x - vertex.X(), seed.y - vertex.Y(), seed.z - vertex.Z());
        if (seed_dir.Mag2() == 0) continue;
        seed_dir = seed_dir.Unit();

        // First pass: gather all hits within cone of seed
        std::vector<int> cluster_candidates;
        for (const auto& h : hSorted) {
            if (used[h.index]) continue;
            
            TVector3 hit_dir(h.x - vertex.X(), h.y - vertex.Y(), h.z - vertex.Z());
            if (hit_dir.Mag2() == 0) continue;
            hit_dir = hit_dir.Unit();
            
            double angle = AngleBetweenVectors(seed_dir, hit_dir);
            if (angle < phiThreshold) {
                cluster_candidates.push_back(h.index);
            }
        }

        // Add candidates to cluster
        for (int idx : cluster_candidates) {
            c.hitIndices.push_back(idx);
            used[idx] = true;
        }

        // Calculate cluster properties
        double sumE = 0.0, cx = 0.0, cy = 0.0, cz = 0.0;
        for (int idx : c.hitIndices) {
            const auto& h = hits[idx];
            sumE += h.e;
            cx += h.x * h.e;
            cy += h.y * h.e;
            cz += h.z * h.e;
        }
        
        if (sumE > 0) {
            c.centroid = TVector3(cx / sumE, cy / sumE, cz / sumE);
            TVector3 dir = (c.centroid - vertex).Unit();
            c.p4.SetPxPyPzE(sumE * dir.X(), sumE * dir.Y(), sumE * dir.Z(), sumE);
            clusters.push_back(c);
        }
    }

    // Filter low-energy clusters (adjust threshold based on your physics)
    clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                  [](const Cluster& c) { return c.p4.E() < 50.0; }), // Increased threshold
                   clusters.end());

    return clusters;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    TFile *f = TFile::Open(argv[1]);
    if (!f || f->IsZombie()) return 1;

    TTree *t = (TTree *)f->Get("digitizedHits");
    if (!t) { std::cerr << "Tree digitizedHits not found" << std::endl; return 1; }

    std::vector<double> *centerXs = nullptr, *centerYs = nullptr, *centerZs = nullptr, *energies = nullptr;
    t->SetBranchAddress("centerX", &centerXs);
    t->SetBranchAddress("centerY", &centerYs);
    t->SetBranchAddress("centerZ", &centerZs);
    t->SetBranchAddress("energy", &energies);

    std::vector<double> *primaryX = nullptr, *primaryY = nullptr, *primaryZ = nullptr;
    t->SetBranchAddress("PrimaryPosX", &primaryX);
    t->SetBranchAddress("PrimaryPosY", &primaryY);
    t->SetBranchAddress("PrimaryPosZ", &primaryZ);

    Long64_t nentries = t->GetEntries();
    TH1F *hPi0Mass = new TH1F("hPi0Mass", ";M_{#gamma#gamma} [MeV];Events", 100, 0, 300);
    TH1F *hClusterE = new TH1F("hClusterE", ";Cluster E [MeV];Count", 100, 0, 500);
    TH1F *hNClusters = new TH1F("hNClusters", ";N_{clusters};Events", 20, 0, 20);
    TH1F *hOpeningAngle = new TH1F("hOpeningAngle", ";Opening Angle [rad];Entries", 100, 0, 1.0);
    TH1F *hHitE = new TH1F("hHitE", ";Hit E [MeV];Count", 100, 0, 10);

    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        t->GetEntry(ievt);

        if (!primaryX || primaryX->empty()) continue;
        TVector3 vertex((*primaryX)[0], (*primaryY)[0], (*primaryZ)[0]);

        std::vector<Hit> hits;
        if (!energies || energies->empty()) {
            std::cout << "Event " << ievt << ": No hits (skipping)" << std::endl;
            continue;
        }
        size_t nHits = energies->size();
        std::cout << "Event " << ievt << ": " << nHits << " total hits across all rings" << std::endl;
        for (size_t k = 0; k < nHits; ++k) {
            double e = (*energies)[k];
            hHitE->Fill(e);
            hits.push_back({(*centerXs)[k], (*centerYs)[k], (*centerZs)[k], e, static_cast<int>(k)});
        }

        std::set<std::tuple<double, double, double>> unique_pos;
        for (const auto& h : hits) {
            unique_pos.insert({h.x, h.y, h.z});
        }
        double min_dr = 1e9;
        for (auto it1 = unique_pos.begin(); it1 != unique_pos.end(); ++it1) {
            auto [x1, y1, z1] = *it1;
            TVector3 v1(x1 - vertex.X(), y1 - vertex.Y(), z1 - vertex.Z());
            for (auto it2 = std::next(it1); it2 != unique_pos.end(); ++it2) {
                auto [x2, y2, z2] = *it2;
                TVector3 v2(x2 - vertex.X(), y2 - vertex.Y(), z2 - vertex.Z());
                double deta = v1.PseudoRapidity() - v2.PseudoRapidity();
                double dphi = TVector2::Phi_mpi_pi(v1.Phi() - v2.Phi());
                double dr = std::sqrt(deta * deta + dphi * dphi);
                if (dr > 0 && dr < min_dr) min_dr = dr;
            }
        }
        std::cout << "Min ΔR between unique hit positions: " << min_dr << std::endl;

        // Check hit energy distribution
        std::cout << "Hit energy stats:" << std::endl;
        double totalE = 0;
        for (const auto& h : hits) totalE += h.e;
        std::cout << "  Total energy: " << totalE << " MeV" << std::endl;
        std::cout << "  Mean hit energy: " << totalE / hits.size() << " MeV" << std::endl;

        // Check angular separation between highest energy hits
        if (hits.size() >= 2) {
            auto h1 = *std::max_element(hits.begin(), hits.end(), 
                [](const Hit& a, const Hit& b) { return a.e < b.e; });
            auto hits_copy = hits;
            hits_copy.erase(std::remove_if(hits_copy.begin(), hits_copy.end(),
                [&h1](const Hit& h) { return h.index == h1.index; }), hits_copy.end());
            auto h2 = *std::max_element(hits_copy.begin(), hits_copy.end(),
                [](const Hit& a, const Hit& b) { return a.e < b.e; });
            
            TVector3 v1(h1.x - vertex.X(), h1.y - vertex.Y(), h1.z - vertex.Z());
            TVector3 v2(h2.x - vertex.X(), h2.y - vertex.Y(), h2.z - vertex.Z());
            double angle = AngleBetweenVectors(v1.Unit(), v2.Unit());
            std::cout << "  Angle between 2 highest E hits: " << angle << " rad (" 
                    << angle * 180 / TMath::Pi() << " deg)" << std::endl;
        }

        std::vector<Cluster> clusters = MomentumBasedClustering(hits, vertex, 0.3); // Adjusted phi to 0.3 rad

        double unassigned_e = 0.0;
        int unassigned_count = 0;
        for (const auto& h : hits) {
            bool assigned = false;
            for (const auto& c : clusters) {
                if (std::find(c.hitIndices.begin(), c.hitIndices.end(), h.index) != c.hitIndices.end()) {
                    assigned = true;
                    break;
                }
            }
            if (!assigned) {
                unassigned_e += h.e;
                unassigned_count++;
            }
        }
        std::cout << "Unassigned hit E: " << unassigned_e << " MeV, count: " << unassigned_count << std::endl;

        for (size_t ci = 0; ci < clusters.size(); ++ci) {
            std::cout << "  Final Cluster " << ci << " E=" << clusters[ci].p4.E()
                      << " MeV, nHits=" << clusters[ci].hitIndices.size() << std::endl;
            hClusterE->Fill(clusters[ci].p4.E());
        }
        hNClusters->Fill(clusters.size());

        std::cout << "Event " << ievt << ": " << clusters.size() << " final clusters" << std::endl;

        if (clusters.size() >= 2) {
            std::sort(clusters.begin(), clusters.end(), [](const Cluster& a, const Cluster& b) { return a.p4.E() > b.p4.E(); });
            TLorentzVector g1 = clusters[0].p4;
            TLorentzVector g2 = clusters[1].p4;
            TLorentzVector pi0 = g1 + g2;
            std::cout << "Photon 1: E=" << g1.E() << ", Px=" << g1.Px() << ", Py=" << g1.Py() << ", Pz=" << g1.Pz() << std::endl;
            std::cout << "Photon 2: E=" << g2.E() << ", Px=" << g2.Px() << ", Py=" << g2.Py() << ", Pz=" << g2.Pz() << std::endl;
            std::cout << "Pi0: M=" << pi0.M() << " MeV" << std::endl;
            hPi0Mass->Fill(pi0.M());

            TVector3 d1(g1.Px(), g1.Py(), g1.Pz());
            TVector3 d2(g2.Px(), g2.Py(), g2.Pz());
            double cosTheta = d1.Unit().Dot(d2.Unit());
            double theta = std::acos(std::min(1.0, std::max(-1.0, cosTheta)));
            std::cout << "DEBUG: E1=" << g1.E() << " E2=" << g2.E() << " theta(deg)=" << (theta * 180. / TMath::Pi())
                      << " vertex=(" << vertex.X() << "," << vertex.Y() << "," << vertex.Z() << ")" << std::endl;
            hOpeningAngle->Fill(theta);
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Results", 1600, 400);
    c1->Divide(3, 1);
    c1->cd(1); hPi0Mass->Draw();
    c1->cd(2); hClusterE->Draw();
    c1->cd(3); hNClusters->Draw();
    // c1->cd(4); hHitE->Draw();
    c1->SaveAs("pi0_mass.png");

    TCanvas *c2 = new TCanvas("c2", "Debug", 800, 400);
    c2->Divide(2, 1);
    c2->cd(1); hOpeningAngle->Draw();
    c2->cd(2); // Note: hCosTheta not defined, skipping
    c2->SaveAs("Openingangles.png");

    delete c1; delete hPi0Mass; delete hClusterE; delete hNClusters; delete hHitE;
    delete c2;
    f->Close(); delete f;
    return 0;
}