#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <string>
#include <set> 

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMarker.h"
#include "TLine.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TF1.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"

// =========================================================
//                 HELPER: Pretty Plot 
// =========================================================

void SetPrettyStyle() {
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.4, "Y");
    // gStyle->SetTextSize(0.045);
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
}

void PrettyPi0MassPlot(TH1F* hPi0Mass) {
    SetPrettyStyle();

    TCanvas *c = new TCanvas("cPi0", "Pi0 Mass", 800, 600);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);

    TF1 *fGaus = new TF1("fGaus", "gaus", 100, 170);
    fGaus->SetLineColor(kRed+1);
    fGaus->SetLineWidth(2);
    hPi0Mass->Fit(fGaus, "RQ");

    double mean  = fGaus->GetParameter(1);
    double sigma = fGaus->GetParameter(2);
    double errMu = fGaus->GetParError(1);
    double errSi = fGaus->GetParError(2);

    hPi0Mass->SetLineColor(kBlack);
    hPi0Mass->SetLineWidth(2);
    hPi0Mass->GetXaxis()->SetTitle("M_{#gamma#gamma} [MeV]");
    hPi0Mass->GetYaxis()->SetTitle("Events");

    hPi0Mass->Draw("HIST");
    fGaus->Draw("SAME");

    TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(hPi0Mass, "Reconstructed #pi^{0} Invariant mass", "l");
    leg->AddEntry(fGaus, "Gaussian Fit:", "l");
    leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f MeV", mean, errMu), "");
    leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f MeV", sigma, errSi), "");
    leg->SetTextSize(0.03);
    leg->Draw();

    TPaveText *info = new TPaveText(0.17, 0.70, 0.50, 0.90, "NDC");  // x1,y1,x2,y2 normalized coordinates
    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextFont(42);
    info->SetTextSize(0.04);
    // info->AddText("Hibeam Wasafull simulation");
    info->AddText("GEANT4 #pi^{0} sample");
    info->AddText("1000 simulated events");
    info->AddText("E_{kin} #in [1, 500] MeV");
    info->Draw();

    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

    c->SaveAs("Pi0Mass_pretty.png");

    // clean up
    delete leg;
    delete fGaus;
    delete c;
}

// =========================================================
//                 HIT & CLUSTER STRUCTURES
// =========================================================
struct Hit {
    double x, y, z, e;
};

struct Cluster {
    std::vector<int> hitIndices;
    TVector3 centroid;
    TLorentzVector p4;
};

// =========================================================
//          HELPER: COORDINATE CONVERSIONS
// =========================================================
struct EtaPhiTowerKey {
    int iEta;
    int iPhi;
    bool operator<(const EtaPhiTowerKey &o) const {
        if (iEta != o.iEta) return iEta < o.iEta;
        return iPhi < o.iPhi;
    }
};

inline double calcEta(double x, double y, double z) {
    TVector3 v(x, y, z);
    return v.Eta();
}

inline double calcPhi(double x, double y, double z) {
    TVector3 v(x, y, z);
    return v.Phi();
}

// =========================================================
//          HELPER: Plotting util for DEBUGGING
// =========================================================
void PlotEventDisplay(
    const std::vector<Hit>& hits,
    const std::vector<Cluster>& clusters,
    double dEta, double dPhi,
    int eventID)
{
    if (hits.empty() || clusters.empty()) return;
    // Create canvas
    TCanvas *c = new TCanvas(Form("c_evt_%d",eventID), 
                             Form("Event %d display",eventID), 800, 600);

    // Define eta-phi histogram
    TH2F *hMap = new TH2F(Form("hEvt%d",eventID),
                          Form("Event %d: Energy in #eta-#phi",eventID),
                          100, -TMath::Pi(), TMath::Pi(),   // phi bins
                          100, -5, 5);                      // eta bins

    // Fill with hits
    for (size_t i=0; i<hits.size(); ++i) {
        double eta = calcEta(hits[i].x, hits[i].y, hits[i].z);
        double phi = calcPhi(hits[i].x, hits[i].y, hits[i].z);
        hMap->Fill(phi, eta, hits[i].e);
    }

    hMap->GetXaxis()->SetTitle("#phi");
    hMap->GetYaxis()->SetTitle("#eta");
    hMap->GetZaxis()->SetTitle("E [MeV]");
    hMap->Draw("COLZ");

    // Overlay clusters with markers
    int colorIdx = 2;
    for (size_t ci=0; ci<clusters.size(); ++ci) {
        for (int idx : clusters[ci].hitIndices) {
            double eta = calcEta(hits[idx].x, hits[idx].y, hits[idx].z);
            double phi = calcPhi(hits[idx].x, hits[idx].y, hits[idx].z);
            TMarker *m = new TMarker(phi, eta, 24); 
            m->SetMarkerColor(colorIdx);
            m->SetMarkerSize(1.5);
            m->Draw("SAME");
        }
        colorIdx++;
    }

    // Optionally overlay grid lines every dEta/dPhi
    double etaMin = hMap->GetYaxis()->GetXmin();
    double etaMax = hMap->GetYaxis()->GetXmax();
    double phiMin = hMap->GetXaxis()->GetXmin();
    double phiMax = hMap->GetXaxis()->GetXmax();

    for (double e = floor(etaMin/dEta)*dEta; e <= etaMax; e+=dEta) {
        TLine *l = new TLine(phiMin, e, phiMax, e);
        l->SetLineColor(kGray);
        l->SetLineStyle(3);
        l->Draw("SAME");
    }
    for (double p = floor(phiMin/dPhi)*dPhi; p <= phiMax; p+=dPhi) {
        TLine *l = new TLine(p, etaMin, p, etaMax);
        l->SetLineColor(kGray);
        l->SetLineStyle(3);
        l->Draw("SAME");
    }

    c->SaveAs(Form("eventDisplay_evt%d.png",eventID));
}

// =========================================================
//      POST-PROCESSING CLUSTER SPLITTING
// =========================================================
std::vector<Cluster> SplitClusters(
    const std::vector<Hit> &hits,
    std::vector<Cluster> &clusters,
    const TVector3 &vertex,
    double dEta,
    double dPhi,
    double E_subseed = 10.0,  // Threshold for sub-seeds (lower than main seed)
    double r0 = 0.03,         // Angular shower scale (tune based on detector; ~Molière radius / calo radius)
    double split_E_thresh = 200.0  // Split if cluster E > this
) {
    std::vector<Cluster> final_clusters;

    for (auto &c : clusters) {
        if (c.p4.E() < split_E_thresh) {
            final_clusters.push_back(c);
            continue;
        }

        // Rebuild local towers for this cluster
        std::map<EtaPhiTowerKey, double> local_towers;
        for (int idx : c.hitIndices) {
            const auto &h = hits[idx];
            double eta = calcEta(h.x, h.y, h.z);
            double phi = calcPhi(h.x, h.y, h.z);
            int iEta = int(std::floor(eta / dEta));
            int iPhi = int(std::floor(phi / dPhi));
            local_towers[{iEta, iPhi}] += h.e;
        }

        // Find sub-seeds (local maxima)
        std::vector<EtaPhiTowerKey> sub_seeds;
        for (auto &[key, E] : local_towers) {
            if (E < E_subseed) continue;
            bool isMax = true;
            for (int de = -1; de <= 1; ++de) {
                for (int dp = -1; dp <= 1; ++dp) {
                    if (de == 0 && dp == 0) continue;
                    EtaPhiTowerKey neigh{key.iEta + de, key.iPhi + dp};
                    if (local_towers.count(neigh) && local_towers[neigh] >= E) {
                        isMax = false;
                        break;
                    }
                }
                if (!isMax) break;
            }
            if (isMax) sub_seeds.push_back(key);
        }

        size_t n_sub = sub_seeds.size();
        if (n_sub < 2) {
            final_clusters.push_back(c);
            continue;
        }

        // Limit to top 2 sub-seeds for simplicity (extend if needed for >2)
        if (n_sub > 2) {
            std::vector<std::pair<EtaPhiTowerKey, double>> sub_with_e;
            for (auto &k : sub_seeds) sub_with_e.push_back({k, local_towers[k]});
            std::sort(sub_with_e.begin(), sub_with_e.end(), [](const auto &a, const auto &b) { return a.second > b.second; });
            sub_seeds = {sub_with_e[0].first, sub_with_e[1].first};
            n_sub = 2;
        }

        // Initialize sub-clusters with centroids from seed towers
        std::vector<Cluster> sub_clusters(n_sub);
        for (size_t si = 0; si < n_sub; ++si) {
            auto &sk = sub_seeds[si];
            double sx = 0, sy = 0, sz = 0, se = 0;
            for (int idx : c.hitIndices) {
                double eta = calcEta(hits[idx].x, hits[idx].y, hits[idx].z);
                double phi = calcPhi(hits[idx].x, hits[idx].y, hits[idx].z);
                int iEta = int(std::floor(eta / dEta));
                int iPhi = int(std::floor(phi / dPhi));
                if (iEta == sk.iEta && iPhi == sk.iPhi) {
                    const auto &h = hits[idx];
                    sx += h.x * h.e;
                    sy += h.y * h.e;
                    sz += h.z * h.e;
                    se += h.e;
                }
            }
            if (se > 0) {
                sub_clusters[si].centroid = TVector3(sx / se, sy / se, sz / se);
                sub_clusters[si].p4.SetE(local_towers[sk]);  // Initial E
            }
        }

        // Iterative splitting (EM-like apportionment)
        bool converged = false;
        int max_iter = 10;
        int iter = 0;
        while (!converged && iter < max_iter) {
            iter++;
            std::vector<double> new_sumE(n_sub, 0);
            std::vector<TVector3> new_cent(n_sub, TVector3(0, 0, 0));
            std::vector<TVector3> new_mom(n_sub, TVector3(0, 0, 0));

            for (int idx : c.hitIndices) {
                const auto &h = hits[idx];
                TVector3 hpos(h.x, h.y, h.z);
                double h_eta = calcEta(h.x, h.y, h.z);
                double h_phi = calcPhi(h.x, h.y, h.z);
                double sum_f = 0;
                std::vector<double> f(n_sub, 0);
                for (size_t si = 0; si < n_sub; ++si) {
                    double s_eta = calcEta(sub_clusters[si].centroid.X(), sub_clusters[si].centroid.Y(), sub_clusters[si].centroid.Z());
                    double s_phi = calcPhi(sub_clusters[si].centroid.X(), sub_clusters[si].centroid.Y(), sub_clusters[si].centroid.Z());
                    double deta = h_eta - s_eta;
                    double dphi = fabs(h_phi - s_phi);
                    if (dphi > TMath::Pi()) dphi = 2 * TMath::Pi() - dphi;
                    double r = sqrt(deta * deta + dphi * dphi);
                    f[si] = exp(-r / r0);
                    sum_f += f[si];
                }
                if (sum_f > 0) {
                    for (size_t si = 0; si < n_sub; ++si) {
                        double w = f[si] / sum_f;
                        new_sumE[si] += w * h.e;
                        new_cent[si] += (w * h.e) * hpos;
                        TVector3 v(h.x - vertex.X(), h.y - vertex.Y(), h.z - vertex.Z());
                        if (v.Mag2() > 1e-12) new_mom[si] += (w * h.e) * v.Unit();
                    }
                }
            }

            converged = true;
            for (size_t si = 0; si < n_sub; ++si) {
                if (new_sumE[si] > 0) {
                    new_cent[si] *= ( 1.0 / new_sumE[si]);
                    if ((new_cent[si] - sub_clusters[si].centroid).Mag() > 1e-3) converged = false;
                    sub_clusters[si].centroid = new_cent[si];
                    sub_clusters[si].p4.SetPxPyPzE(new_mom[si].X(), new_mom[si].Y(), new_mom[si].Z(), new_sumE[si]);
                    // Approximate hitIndices (assign if w > 0.5 for plotting/couting)
                    // Optional: loop over hits again and push_back if w > 0.5 for this si
                }
            }
        }

        // Add valid sub-clusters
        bool split_success = false;
        for (auto &sc : sub_clusters) {
            if (sc.p4.E() >= 50.0) {  // Same threshold as main
                final_clusters.push_back(sc);
                split_success = true;
            }
        }
        if (!split_success) final_clusters.push_back(c);  // Fallback if split fails
    }
    return final_clusters;
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
    int winSize = 7)
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
        std::cout << "Seed at (iEta,iPhi)=(" << seed.key.iEta << "," << seed.key.iPhi << ") "
                  << "window covers " << nTowersInWindow << " towers, "
                  << nHitsInWindow << " hits, "
                  << "cluster E=" << sumE << std::endl;

        if (added && sumE > 0) {
            c.centroid = TVector3(cx / sumE, cy / sumE, cz / sumE);
            c.p4.SetPxPyPzE(mom.X(), mom.Y(), mom.Z(), sumE);
            clusters.push_back(c);
        }
    }
    return clusters;
}

// =========================================================
//                           MAIN
// =========================================================
int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    SetPrettyStyle();

    TFile *f = TFile::Open(argv[1]);
    if (!f || f->IsZombie()) return 1;

    TTree *t = (TTree *)f->Get("digitizedHits");
    if (!t) { std::cerr << "Tree digitizedHits not found" << std::endl; return 1; }

    // Single consolidated vectors (matches preprocessed output—no per-ring arrays!)
    std::vector<double> *centerXs = nullptr, *centerYs = nullptr, *centerZs = nullptr, *energies = nullptr;
    t->SetBranchAddress("centerX", &centerXs);
    t->SetBranchAddress("centerY", &centerYs);
    t->SetBranchAddress("centerZ", &centerZs);
    t->SetBranchAddress("energy", &energies);

    // Primary vertex (vector, but size=1 per event)
    std::vector<double> *primaryX=nullptr,*primaryY=nullptr,*primaryZ=nullptr;
    t->SetBranchAddress("PrimaryPosX",&primaryX);
    t->SetBranchAddress("PrimaryPosY",&primaryY);
    t->SetBranchAddress("PrimaryPosZ",&primaryZ);

    Long64_t nentries = t->GetEntries();
    TH1F *hPi0Mass = new TH1F("hPi0Mass",";M_{#gamma#gamma} [MeV];Events",100,0,300);
    TH1F *hClusterE = new TH1F("hClusterE",";Cluster E [MeV];Count",100,0,500);
    TH1F *hNClusters = new TH1F("hNClusters",";N_{clusters};Events",20,0,20);

    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        t->GetEntry(ievt);

        if (!primaryX||primaryX->empty()) continue;
        TVector3 vertex((*primaryX)[0],(*primaryY)[0],(*primaryZ)[0]);

        std::vector<Hit> hits;
        if (!energies || energies->empty()) {
            std::cout << "Event " << ievt << ": No hits (skipping)" << std::endl;
            continue;
        }
        size_t nHits = energies->size();
        std::cout << "Event " << ievt << ": " << nHits << " total hits across all rings" << std::endl;
        for (size_t k=0; k<nHits; ++k) {
            // Single loop over consolidated vectors (no sec/ring loop needed)
            hits.push_back({(*centerXs)[k], (*centerYs)[k], (*centerZs)[k], (*energies)[k]});
        }

        std::vector<Cluster> clusters;
        double dEta = 0.30; 
        double dPhi = 0.30;
        double E_seed = 20.00;
        double E_neighbor = 0.03;
        int winSize = 3;  
        clusters = SlidingWindowClusterHits(hits,vertex, dEta, dPhi, E_seed, E_neighbor, winSize);

        // NEW: Split large clusters
        clusters = SplitClusters(hits, clusters, vertex, dEta, dPhi);

        // Apply cluster energy threshold
        clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                    [](const Cluster &c){ return c.p4.E() < 50.0; }),
                    clusters.end());

        // If too few clusters found, relax parameters
        if (clusters.size() < 2) {
            std::cout << "Fallback: relaxing clustering params for event " << ievt << std::endl;
            double relaxed_E_seed = std::max(1.0, E_seed*0.5);
            double relaxed_E_neighbor = std::max(0.00, E_neighbor*0.5);
            int relaxed_winSize = winSize*2;
            clusters = SlidingWindowClusterHits(hits, vertex,
                                                dEta, dPhi,
                                                relaxed_E_seed,
                                                relaxed_E_neighbor,
                                                relaxed_winSize);

            // NEW: Split in fallback too
            clusters = SplitClusters(hits, clusters, vertex, dEta, dPhi);

            // Apply cut again
            clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                        [](const Cluster &c){ return c.p4.E() < 50.0; }),
                        clusters.end());
        }
        // Keep only the top 2 highest-energy clusters (for π⁰ focus)
        if (!clusters.empty()) {
            // Sort clusters by energy (descending)
            std::sort(clusters.begin(), clusters.end(), 
                    [](const Cluster &a, const Cluster &b){ return a.p4.E() > b.p4.E(); });

            if (clusters.size() > 2) clusters.resize(2);
        }

        for (size_t ci=0; ci<clusters.size(); ++ci) {
            hClusterE->Fill(clusters[ci].p4.E());
        }
        hNClusters->Fill(clusters.size());

        // =========================================================
        //                           DEBUG
        // =========================================================

        std::cout << "Event " << ievt << ": " << clusters.size() << " clusters found" << std::endl;

        for (size_t ci=0; ci<clusters.size(); ++ci) {
            std::cout << "  Cluster " << ci
                    << " E = " << clusters[ci].p4.E()
                    << " nHits = " << clusters[ci].hitIndices.size() 
                    << std::endl;
        }

        // if (ievt < 10) { // only plot first 10 events
        //     PlotEventDisplay(hits, clusters, 0.025, 0.025, ievt);
        // }

        // =========================================================
        //                           DEBUG
        // =========================================================

        for (size_t a=0;a<clusters.size();++a) {
            for (size_t b=a+1;b<clusters.size();++b) {
                TLorentzVector pi0 = clusters[a].p4 + clusters[b].p4;
                hPi0Mass->Fill(pi0.M());
            }
        }
    }

    PrettyPi0MassPlot(hPi0Mass);


    // delete c1; 
    delete hPi0Mass; delete hClusterE; delete hNClusters;
    f->Close(); delete f;
    return 0;
}