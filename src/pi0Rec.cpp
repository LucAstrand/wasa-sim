#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <numeric>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TVector3.h"  
#include "TMath.h"
#include "TString.h"   
#include "TCanvas.h"

// Hold the hits
struct Hit {
    double x, y, z;
    double e;
};

//Hold the cluster
struct Cluster {
    double energy;
    double x, y, z;
    TLorentzVector p4;
    std::vector<int> hitIndices;
};

//Helper functions

//Calculate angle between two ROOT 3-vectors
double AngleBetween(const Hit &a, const Hit &b, const Hit &vertex) {
    TVector3 v1(a.x - vertex.x, a.y - vertex.y, a.z - vertex.z);
    TVector3 v2(b.x - vertex.x, b.y - vertex.y, b.z - vertex.z);
    return v1.Angle(v2);
}

//Clustering function - Main logic
// std::vector<Cluster> ClusterHits(std::vector<Hit> &hits, const Hit &vertex, double phiMax) {
//     //sort hits in descending order
//     // seed hit as the most energetic one!
//     std::sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
//         return a.e > b.e; // sort descending by energy
//     });

//     std::vector<Cluster> clusters;
//     std::vector<bool> used(hits.size(), false);

//     for (size_t i=0; i < hits.size(); i++) {

//         if (used[i]) continue;

//         Cluster c{};
//         double sumE = 0, sumX = 0, sumY = 0, sumZ = 0;
//         std::vector<int> members; 
//         members.push_back(i);
//         used[i] = true;

//         for (size_t j=0; j < hits.size(); j++) {
//             if (used[j]) continue;
//             double ang = AngleBetween(hits[i], hits[j], vertex);
//             if (ang < phiMax) {
//                 members.push_back(j);
//                 used[j] = true;
//             }
//         }

//         for (int idx : members) {
//             sumE += hits[idx].e;
//             sumX += hits[idx].x * hits[idx].e;
//             sumY += hits[idx].y * hits[idx].e;
//             sumZ += hits[idx].z * hits[idx].e;
//         }

//         if (sumE <= 0) continue;

//         c.energy = sumE;
//         c.x = sumX / sumE;
//         c.y = sumY / sumE;
//         c.z = sumZ / sumE;
//         c.hitIndices = members;

//         //Momentum direction
//         TVector3 direction(c.x - vertex.x, c.y - vertex.y, c.z - vertex.z);
//         direction = direction.Unit();

//         double px = direction.X() * sumE;
//         double py = direction.Y() * sumE;
//         double pz = direction.Z() * sumE;
//         c.p4.SetPxPyPzE(px, py, pz, sumE);

//         clusters.push_back(c);
    
//     }

//     return clusters;

// }

// Find neighbors: cells within a certain spatial distance
std::vector<int> FindNeighbors(const std::vector<Hit> &hits, int idx, double maxDist) {
    std::vector<int> neighbors;
    const Hit &seed = hits[idx];
    for (size_t i = 0; i < hits.size(); ++i) {
        if ((int)i == idx) continue;
        double dx = hits[i].x - seed.x;
        double dy = hits[i].y - seed.y;
        double dz = hits[i].z - seed.z;
        double dist2 = dx*dx + dy*dy + dz*dz;
        if (dist2 < maxDist * maxDist) {
            neighbors.push_back(i);
        }
    }
    return neighbors;
}

// Clustering with neighbor adjacency + cluster energy cut
std::vector<Cluster> ClusterHits(std::vector<Hit> &hits,
                                 double seedThreshold,
                                 double neighborThreshold,
                                 double maxDist,
                                 double minClusterEnergy)
{
    // Sort hits by energy, descending
    std::vector<int> indices(hits.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int a, int b){ return hits[a].e > hits[b].e; });

    std::vector<Cluster> clusters;
    std::vector<bool> used(hits.size(), false);

    for (int idx : indices) {
        if (used[idx]) continue;
        if (hits[idx].e < seedThreshold) continue;

        Cluster c{};
        double sumE = 0, sumX = 0, sumY = 0, sumZ = 0;
        std::vector<int> members;

        // BFS to collect neighbors
        std::queue<int> q;
        q.push(idx);
        used[idx] = true;

        while (!q.empty()) {
            int current = q.front();
            q.pop();
            members.push_back(current);

            const Hit &h = hits[current];
            sumE += h.e;
            sumX += h.x * h.e;
            sumY += h.y * h.e;
            sumZ += h.z * h.e;

            auto neighbors = FindNeighbors(hits, current, maxDist);
            for (int n : neighbors) {
                if (used[n]) continue;
                if (hits[n].e < neighborThreshold) continue;
                used[n] = true;
                q.push(n);
            }
        }

        // Final cluster energy cut
        if (sumE < minClusterEnergy) continue;

        // Compute cluster centroid
        c.energy = sumE;
        c.x = sumX / sumE;
        c.y = sumY / sumE;
        c.z = sumZ / sumE;
        c.hitIndices = members;

        // Build TLorentzVector (photon assumption)
        TVector3 direction(c.x, c.y, c.z);
        direction = direction.Unit();
        double px = direction.X() * sumE;
        double py = direction.Y() * sumE;
        double pz = direction.Z() * sumE;
        c.p4.SetPxPyPzE(px, py, pz, sumE);

        clusters.push_back(c);
    }

    return clusters;
}



int main(int argc, char **argv) {  // Fixed main signature
    if (argc < 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;
    }
    
    const char *filename = argv[1];
    
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) return 1;

    TTree *t = (TTree*)f->Get("hibeam");
    if (!t) {
        printf("Tree 'hibeam' not found\n");
        return 1;
    }

    // For each section, we will have separate vectors
    const int NSEC = 16;
    std::vector<double> *posX[NSEC], *posY[NSEC], *posZ[NSEC], *edep[NSEC];

    // Initialize pointers to nullptr
    for (int i = 0; i < NSEC; ++i) {
        posX[i] = nullptr;
        posY[i] = nullptr;
        posZ[i] = nullptr;
        edep[i] = nullptr;
    }

    // Dynamically set branch addresses
    for (int i = 0; i < NSEC; ++i) {
        TString sec = Form("SECE%d", i);
        t->SetBranchAddress(sec + "_GlobalPosition_X", &posX[i]);
        t->SetBranchAddress(sec + "_GlobalPosition_Y", &posY[i]);
        t->SetBranchAddress(sec + "_GlobalPosition_Z", &posZ[i]);
        t->SetBranchAddress(sec + "_EDep", &edep[i]);
    }

    std::vector<double> *primaryX = nullptr, *primaryY = nullptr, *primaryZ = nullptr;
    t->SetBranchAddress("PrimaryPosX", &primaryX);
    t->SetBranchAddress("PrimaryPosY", &primaryY);
    t->SetBranchAddress("PrimaryPosZ", &primaryZ);

    Long64_t nentries = t->GetEntries();
    TH1F *hPi0Mass = new TH1F("hPi0Mass", "Reconstructed #pi^{0} mass;M_{#gamma#gamma} [GeV];Events", 100, 0, 0.3);

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);

        // Check if primary vectors are valid and non-empty
        if (!primaryX || !primaryY || !primaryZ || 
            primaryX->empty() || primaryY->empty() || primaryZ->empty()) {
            continue;
        }

        std::vector<Hit> hits;
        for (int sec = 0; sec < NSEC; ++sec) {
            if (!edep[sec] || !posX[sec] || !posY[sec] || !posZ[sec]) continue;
            
            for (size_t k = 0; k < edep[sec]->size(); ++k) {
                double E = (*edep[sec])[k];
                if (E < 1e-3) continue; // suppress noise
                hits.push_back({(*posX[sec])[k], (*posY[sec])[k], (*posZ[sec])[k], E});
            }
        }

        Hit vertex{(*primaryX)[0], (*primaryY)[0], (*primaryZ)[0], 0};

        // auto clusters = ClusterHits(hits, vertex, 0.2 /*phiMax*/);

        auto clusters = ClusterHits(hits,
                            5e-2,   // seed threshold = 50 MeV
                            2.5e-2, // neighbor threshold = 25 MeV (~20-30 MeV)
                            30.0,   // max distance between neighboring hits [mm]
                            0.05);  // minimum cluster energy [GeV]

        // Combine photon pairs into pi0 candidates
        for (size_t a = 0; a < clusters.size(); ++a) {
            for (size_t b = a+1; b < clusters.size(); ++b) {
                TLorentzVector pi0 = clusters[a].p4 + clusters[b].p4;
                hPi0Mass->Fill(pi0.M());
            }
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Pi0 Mass Distribution", 800, 600);
    hPi0Mass->Draw();
    c1->SaveAs("pi0_mass.png");
    
    // Clean up
    delete c1;
    delete hPi0Mass;
    f->Close();
    delete f;
    
    return 0;
}