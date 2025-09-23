#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <numeric>
#include <stack>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TVector3.h"  
#include "TMath.h"
#include "TString.h"   
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TLegend.h"
#include "TColor.h"
#include "TF1.h"
#include "TStyle.h"

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


std::vector<int> FindNeighbors(const std::vector<Hit> &hits, int idx, double maxDist) {
    std::vector<int> neighbors;
    const Hit &seed = hits[idx];
    double maxDist2 = maxDist * maxDist;
    
    for (size_t i = 0; i < hits.size(); ++i) {
        if (static_cast<int>(i) == idx) continue;
        
        const Hit &hit = hits[i];
        double dx = hit.x - seed.x;
        double dy = hit.y - seed.y;
        double dz = hit.z - seed.z;
        double dist2 = dx*dx + dy*dy + dz*dz;
        
        if (dist2 < maxDist2) {
            neighbors.push_back(static_cast<int>(i));
        }
    }
    return neighbors;
}

// -----------------------------------------------------
void CollectClusterHitsIterative(const std::vector<Hit> &hits,
                                 int startIdx,
                                 double neighborThreshold,
                                 double maxDist,
                                 std::vector<bool> &used,
                                 std::vector<int> &members) {
    std::stack<int> toVisit;
    toVisit.push(startIdx);
    used[startIdx] = true;

    while (!toVisit.empty()) {
        int current = toVisit.top();
        toVisit.pop();
        members.push_back(current);

        auto neighbors = FindNeighbors(hits, current, maxDist);
        for (int nIdx : neighbors) {
            if (used[nIdx]) continue;
            if (hits[nIdx].e < neighborThreshold) continue;

            used[nIdx] = true;
            toVisit.push(nIdx);
        }
    }
}


std::vector<Cluster> ClusterHits(const std::vector<Hit> &hits,
                                 const Hit &vertex,
                                 double seedThreshold,     
                                 double neighborThreshold, 
                                 double maxDist,
                                 double minClusterEnergy)  
{
    if (hits.empty()) return {};

    std::vector<int> indices(hits.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(),
              [&hits](int a, int b) { return hits[a].e > hits[b].e; });

    std::vector<Cluster> clusters;
    std::vector<bool> used(hits.size(), false);

    for (int idx : indices) {
        if (used[idx]) continue;
        if (hits[idx].e < seedThreshold) continue;

        std::vector<int> members;
        CollectClusterHitsIterative(hits, idx, neighborThreshold, maxDist, used, members);

        double sumE=0, sumX=0, sumY=0, sumZ=0;
        for (int m : members) {
            const auto &h = hits[m];
            sumE += h.e;
            sumX += h.x * h.e;
            sumY += h.y * h.e;
            sumZ += h.z * h.e;
        }

        if (sumE < minClusterEnergy) continue; 

        Cluster c;
        c.energy = sumE;
        c.x = sumX / sumE;
        c.y = sumY / sumE;
        c.z = sumZ / sumE;
        c.hitIndices = std::move(members);

        // Build 4-vector
        // Energy-weighted direction using all hits in the cluster
        TVector3 momentum(0,0,0);
        for (int idx : c.hitIndices) {
            TVector3 v(hits[idx].x - vertex.x, hits[idx].y - vertex.y, hits[idx].z - vertex.z);
            if (v.Mag2() > 1e-12) {
                v = v.Unit();           // unit vector
                momentum += hits[idx].e * v; // weight by hit energy
            }
        }

        if (momentum.Mag2() > 1e-12) {
            // momentum = momentum.Unit() * sumE;
            c.p4.SetPxPyPzE(momentum.X(), momentum.Y(), momentum.Z(), momentum.Mag());
        } else {
            // fallback
            c.p4.SetPxPyPzE(0,0,0,sumE);
        }

          clusters.push_back(std::move(c));
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
    // TH1F *hPi0Mass = new TH1F("hPi0Mass", "Reconstructed #pi^{0} mass;M_{#gamma#gamma} [GeV];Events", 100, 0, 0.3);
    TH1F *hPi0Mass = new TH1F("hPi0Mass", "Reconstructed #pi^{0} mass;M_{#gamma#gamma} [MeV];Events", 100, 0, 300);


    for (Long64_t i = 0; i < nentries; ++i) {
        // std::cout << "Event: " << i << " / " << nentries << std::endl;
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
                // I dont think this is correct, like the hits are very low energy, what we probably really wanted to do is filter on clusters not here
                // if (E < 25) continue; 
                hits.push_back({(*posX[sec])[k], (*posY[sec])[k], (*posZ[sec])[k], E});
            }
        }

        Hit vertex{(*primaryX)[0], (*primaryY)[0], (*primaryZ)[0], 0};


// // Add this diagnostic right after hit collection, before clustering:
//         if (hits.size() > 0) {
//             std::cout << "\n=== Event " << i << " Raw Hit Analysis ===" << std::endl;
//             std::cout << "Total hits: " << hits.size() << std::endl;
            
//             // Show spatial distribution of hits
//             double minX = 1e6, maxX = -1e6, minY = 1e6, maxY = -1e6, minZ = 1e6, maxZ = -1e6;
//             for (const auto& hit : hits) {
//                 minX = std::min(minX, hit.x); maxX = std::max(maxX, hit.x);
//                 minY = std::min(minY, hit.y); maxY = std::max(maxY, hit.y);  
//                 minZ = std::min(minZ, hit.z); maxZ = std::max(maxZ, hit.z);
//             }
//             std::cout << "Hit spread: X[" << minX << ", " << maxX << "] Y[" << minY << ", " << maxY << "] Z[" << minZ << ", " << maxZ << "]" << std::endl;
            
//             // Show distances between high-energy hits
//             std::vector<int> highEnergyHits;
//             for (size_t h = 0; h < hits.size(); ++h) {
//                 if (hits[h].e > 25) highEnergyHits.push_back(h);
//             }
            
//             std::cout << "High energy hits (>25 MeV): " << highEnergyHits.size() << std::endl;
//             for (int idx : highEnergyHits) {
//                 double dist = sqrt(pow(hits[idx].x - vertex.x, 2) + 
//                                   pow(hits[idx].y - vertex.y, 2) + 
//                                   pow(hits[idx].z - vertex.z, 2));
//                 std::cout << "  Hit " << idx << ": E=" << hits[idx].e 
//                           << " MeV, dist from vertex=" << dist << " cm" << std::endl;
//             }
            
//             // Calculate distances between high-energy hits
//             if (highEnergyHits.size() > 1) {
//                 std::cout << "Distances between high-energy hits:" << std::endl;
//                 for (size_t i = 0; i < highEnergyHits.size(); ++i) {
//                     for (size_t j = i+1; j < highEnergyHits.size(); ++j) {
//                         int idx1 = highEnergyHits[i], idx2 = highEnergyHits[j];
//                         double dist = sqrt(pow(hits[idx1].x - hits[idx2].x, 2) + 
//                                           pow(hits[idx1].y - hits[idx2].y, 2) + 
//                                           pow(hits[idx1].z - hits[idx2].z, 2));
//                         std::cout << "  Hit " << idx1 << " to Hit " << idx2 << ": " << dist << " cm" << std::endl;
//                     }
//             }
//         }
//     }

    const double seedThreshold   = 25;   // MeV, small enough to start a cluster
    const double neighborThreshold = 0; // MeV, include low-energy hits around seed
    const double maxDist         = 15.0;   // cm, increase if hits are spread
    const double minClusterEnergy = 30.0; // MeV, total energy in cluster

        // ------------------- Diagnostic per event -------------------
        if (!hits.empty()) {
            std::cout << "\n=== Event " << i << " Hit Energy Diagnostics ===" << std::endl;

            int nAboveSeed = 0;
            int nAboveNeighbor = 0;
            int nAboveMinCluster = 0;
            double sumE = 0;

            for (size_t h = 0; h < hits.size(); ++h) {
                const auto &hit = hits[h];
                sumE += hit.e;

                if (hit.e >= seedThreshold) nAboveSeed++;           // seedThreshold
                if (hit.e >= neighborThreshold)  nAboveNeighbor++;       // neighborThreshold
            }

            if (sumE >= minClusterEnergy) nAboveMinCluster++;         // minClusterEnergy

            std::cout << "Total hits: " << hits.size() << std::endl;
            std::cout << "Hits >= seedThreshold (25 MeV): " << nAboveSeed << std::endl;
            std::cout << "Hits >= neighborThreshold (5 MeV): " << nAboveNeighbor << std::endl;
            std::cout << "Total energy sum: " << sumE << " MeV"
                    << " (above minClusterEnergy 50 MeV? " 
                    << (sumE >= minClusterEnergy ? "yes" : "no") << ")" << std::endl;

            // // Optional: print all per-hit energies
            // std::cout << "Hit energies: ";
            // for (const auto &hit : hits) std::cout << hit.e << " ";
            // std::cout << std::endl;
        }


        auto clusters = ClusterHits(hits,
                                    vertex,
                                    seedThreshold,       // 5 MeV
                                    neighborThreshold,   // 1 MeV
                                    maxDist,             // 5 cm
                                    minClusterEnergy);   // 30 MeV

        
        // // after auto clusters = ClusterHits( ... );
        // std::cout << "Event " << i << ": Found " << clusters.size() << " clusters" << std::endl;

        // for (size_t a = 0; a < clusters.size(); ++a) {
        //     for (size_t b = a+1; b < clusters.size(); ++b) {
        //         double angle = clusters[a].p4.Vect().Angle(clusters[b].p4.Vect()) * 180.0 / TMath::Pi();
        //         std::cout << "Opening angle (deg): " << angle << std::endl;
        //     }
        // }

        // for (size_t a = 0; a < clusters.size(); ++a) {
        //     auto &A = clusters[a];
        //     TVector3 dirA = A.p4.Vect();
        //     double dirA_mag = dirA.Mag();
        //     TVector3 dirA_unit = dirA.Mag() > 0 ? dirA.Unit() : TVector3(0,0,0);
        //     std::cout << " Cluster " << a << ": E=" << A.energy
        //             << " p_mag=" << dirA_mag
        //             << " |p|/E=" << (A.p4.E() > 0 ? dirA_mag / A.p4.E() : -1)
        //             << " dir=(" << dirA_unit.X() << "," << dirA_unit.Y() << "," << dirA_unit.Z() << ")"
        //             << " #hits=" << A.hitIndices.size() << std::endl;

        //     // show the hit-level directions used to build this cluster:
        //     std::cout << "   Hit directions (unit vectors) and weights:" << std::endl;
        //     for (int hid : A.hitIndices) {
        //         TVector3 v_hit(hits[hid].x - vertex.x, hits[hid].y - vertex.y, hits[hid].z - vertex.z);
        //         double d = v_hit.Mag();
        //         if (d > 1e-12) v_hit = v_hit.Unit();
        //         std::cout << "     hit " << hid << ": E=" << hits[hid].e
        //                 << " u=(" << v_hit.X() << "," << v_hit.Y() << "," << v_hit.Z() << ")"
        //                 << std::endl;
        //     }
        // }

        // // compare TLorentz pair masses vs analytic formula
        // for (size_t a = 0; a < clusters.size(); ++a) {
        //     for (size_t b = a+1; b < clusters.size(); ++b) {
        //         double theta = clusters[a].p4.Vect().Angle(clusters[b].p4.Vect());
        //         double cosTh = std::cos(theta);
        //         double E1 = clusters[a].p4.E();
        //         double E2 = clusters[b].p4.E();
        //         double m2_analytic = 2.0 * E1 * E2 * (1.0 - cosTh);
        //         double m_analytic = (m2_analytic > 0) ? std::sqrt(m2_analytic) : 0.0;
        //         double m_tlv = (clusters[a].p4 + clusters[b].p4).M();
        //         std::cout << " Pair " << a << "-" << b 
        //                 << ": theta(deg)=" << (theta*180.0/TMath::Pi())
        //                 << " E1=" << E1 << " E2=" << E2
        //                 << " AnalyticMass=" << m_analytic 
        //                 << " TLVMass=" << m_tlv << std::endl;
        //     }
        // }


        
        // // Add this debug code in your main loop after clustering:
        std::cout << "Event " << i << ": Found " << clusters.size() << " clusters" << std::endl;
        // Count how many pairs we're checking
        int pairCount = 0;
        for (size_t a = 0; a < clusters.size(); ++a) {
            for (size_t b = a+1; b < clusters.size(); ++b) {
                TLorentzVector pi0 = clusters[a].p4 + clusters[b].p4;
                double mass = pi0.M();
                pairCount++;
                
                // Debug: print masses that seem reasonable for pi0
                if (mass > 50 && mass < 300) { // 50-300 MeV range
                    std::cout << "  Pair " << a << "-" << b 
                              << ": mass = " << mass << " MeV" << std::endl;
                    std::cout << "    Cluster A: E=" << clusters[a].p4.E() 
                              << " |p|=" << clusters[a].p4.P() << std::endl;
                    std::cout << "    Cluster B: E=" << clusters[b].p4.E() 
                              << " |p|=" << clusters[b].p4.P() << std::endl;
                }
                
                hPi0Mass->Fill(mass);
            }
        }
        std::cout << "  Checked " << pairCount << " photon pairs" << std::endl;
        
        // // // Stop after a few events for debugging
        // if (i > 10) break;

        // Combine photon pairs into pi0 candidates
        for (size_t a = 0; a < clusters.size(); ++a) {
            for (size_t b = a+1; b < clusters.size(); ++b) {
                
                double angle = clusters[a].p4.Vect().Angle(clusters[b].p4.Vect()) * 180/TMath::Pi();
                std::cout << "Opening angle: " << angle << " deg" << std::endl;                    
                TLorentzVector pi0 = clusters[a].p4 + clusters[b].p4;
                hPi0Mass->Fill(pi0.M());
            }
        }
    }



    TCanvas *c1 = new TCanvas("c1", "Pi0 Mass Distribution", 800, 600);

    // --- Optional: automatically zoom in to region around the peak ---
    int peakBin = hPi0Mass->GetMaximumBin();
    double peakCenter = hPi0Mass->GetBinCenter(peakBin);

    // Zoom to ±2σ-ish window around peak (here I just use ±50 MeV as a first guess)
    hPi0Mass->GetXaxis()->SetRangeUser(peakCenter - 50, peakCenter + 50);

    gStyle->SetOptStat(0); // turn off default stats box
    gStyle->SetOptFit(0);  // turn off auto-fit box

    // Draw & fit
    hPi0Mass->Draw();
    hPi0Mass->Fit("gaus", "Q");  // Q = quiet mode

    // Get the fit object and parameters
    TF1 *fit = hPi0Mass->GetFunction("gaus");
    if (fit) {
        double mu = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        TLegend *legend = new TLegend(0.55, 0.65, 0.85, 0.85);
        legend->AddEntry(hPi0Mass, "pi^{0} invariant mass", "l");
        legend->AddEntry((TObject*)0, Form("#mu = %.2f MeV", mu), "");
        legend->AddEntry((TObject*)0, Form("#sigma = %.2f MeV", sigma), "");
        legend->Draw();
    }

    c1->SaveAs("pi0_mass.png");


    // // DOUBLE FIT 
    // // Get the fit object and parameters
    // TF1 *fit = hPi0Mass->GetFunction("gaus");
    // if (fit) {
    //     double mu = fit->GetParameter(1);
    //     double sigma = fit->GetParameter(2);

    //     // Refine the fit range using mu ± 2σ for a more stable second fit
    //     hPi0Mass->Fit("gaus", "QR", "", mu - 2*sigma, mu + 2*sigma);

    //     // Draw legend with fit parameters
    //     TF1 *fit2 = hPi0Mass->GetFunction("gaus");
    //     if (fit2) {
    //         mu = fit2->GetParameter(1);
    //         sigma = fit2->GetParameter(2);

    //         TLegend *legend = new TLegend(0.55, 0.65, 0.85, 0.85);
    //         legend->AddEntry(hPi0Mass, "pi^{0} invariant mass", "l");
    //         legend->AddEntry((TObject*)0, Form("#mu = %.2f MeV", mu), "");
    //         legend->AddEntry((TObject*)0, Form("#sigma = %.2f MeV", sigma), "");
    //         legend->Draw();
    //     }
    // }

    // c1->SaveAs("pi0_mass.png");

    
    // Clean up
    delete c1;
    delete hPi0Mass;
    f->Close();
    delete f;
    
    return 0;
}