#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraph.h"

#include "CLI11.hpp"
#include "Structures.hpp"
#include "Clustering.hpp"
#include "PlotUtils.hpp"
#include "EventDisplay.hpp"
#include "Utils.hpp"
#include "TruePhotonCalc.hpp"
#include "PhotonMatch.hpp"
#include "Pi0Efficiency.hpp"
#include "Pi0Acceptance.hpp"
#include "PIDEfficiency.hpp"

// ---------------------- Helper ----------------------
template <typename T>
bool SafeSetBranch(TTree* tree, const char* branchName, T*& ptr) {
    if (tree->GetListOfBranches()->FindObject(branchName)) {
        tree->SetBranchAddress(branchName, &ptr);
        return true;
    } else {
        std::cout << "[Warning] Branch " << branchName << " not found. Skipping." << std::endl;
        return false;
    }
}
// ----------------------------------------------------

int main(int argc, char **argv) {
    if (argc < 2) {
            std::cout << "Usage: " << argv[0] << " <input.root> [options]\n"
                    << "Options:\n"
                    << "  --full-analysis       Performs everything\n"
                    << "  --pi0-analysis        Reconstruct Pi0s and invariant mass analysis\n"
                    << "  --charged-analysis    Reconstruct Charged objects and do PID studies\n"
                    << "  --truth-analysis      Performs truth level analysis for Pi0s\n";
            return 1;
        }
    
    std::string inputfile = argv[1];
    bool doPi0Analysis = false;
    bool doChargedAnalysis = false;
    bool doTruthAnalysis = false;

    for (int i = 2; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--full-analysis") {
            doPi0Analysis = true;
            doChargedAnalysis = true;
            doTruthAnalysis = true;
        }
        if (arg == "--pi0-analysis") doPi0Analysis = true;
        if (arg == "--charged-analysis") doChargedAnalysis = true;
        if (arg == "--truth-analysis") doTruthAnalysis = true;
    }

    
    SetPrettyStyle();

    TFile *f = TFile::Open(inputfile.c_str());
    if (!f || f->IsZombie()) return 1;

    TTree *t = (TTree *)f->Get("digitizedHits");
    if (!t) { std::cerr << "Tree digitizedHits not found" << std::endl; return 1; }

    //Hit Cell center coordinates
    std::vector<double> *centerXs = nullptr, *centerYs = nullptr, *centerZs = nullptr, *energies = nullptr;
    SafeSetBranch(t, "centerX", centerXs);
    SafeSetBranch(t, "centerY", centerYs);
    SafeSetBranch(t, "centerZ", centerZs);
    SafeSetBranch(t, "energy", energies);

    // Primary vertex
    std::vector<double> *primaryX = nullptr, *primaryY = nullptr, *primaryZ = nullptr, *primaryEkin = nullptr;
    std::vector<double> *primaryPx = nullptr, *primaryPy = nullptr, *primaryPz = nullptr;
    SafeSetBranch(t, "PrimaryPosX", primaryX);
    SafeSetBranch(t, "PrimaryPosY", primaryY);
    SafeSetBranch(t, "PrimaryPosZ", primaryZ);
    SafeSetBranch(t, "PrimaryEkin", primaryEkin);
    SafeSetBranch(t, "PrimaryMomX", primaryPx);
    SafeSetBranch(t, "PrimaryMomY", primaryPy);
    SafeSetBranch(t, "PrimaryMomZ", primaryPz);

    // Truth info
    std::vector<double> *truthPosX = nullptr, *truthPosY = nullptr, *truthPosZ = nullptr, *truthE = nullptr;
    SafeSetBranch(t, "truthPosX", truthPosX);
    SafeSetBranch(t, "truthPosY", truthPosY);
    SafeSetBranch(t, "truthPosZ", truthPosZ);
    SafeSetBranch(t, "truthE", truthE);

    // TPC info
    std::vector<double> *TPC_Edep = nullptr, *TPC_smearedEdep = nullptr;
    std::vector<double> *TPC_firstPosX = nullptr, *TPC_firstPosY = nullptr, *TPC_firstPosZ = nullptr;
    std::vector<double> *TPC_lastPosX = nullptr, *TPC_lastPosY = nullptr, *TPC_lastPosZ = nullptr;
    std::vector<double> *TPC_PathLength = nullptr, *TPC_dEdx = nullptr, *TPC_TrueKE = nullptr, *TPC_pdg = nullptr; // *TPC_Psm = nullptr,
    SafeSetBranch(t, "TPC_Edep", TPC_Edep);
    SafeSetBranch(t, "TPC_smearedEdep", TPC_smearedEdep);
    SafeSetBranch(t, "TPC_firstPosX", TPC_firstPosX);
    SafeSetBranch(t, "TPC_firstPosY", TPC_firstPosY);
    SafeSetBranch(t, "TPC_firstPosZ", TPC_firstPosZ);
    SafeSetBranch(t, "TPC_lastPosX", TPC_lastPosX);
    SafeSetBranch(t, "TPC_lastPosY", TPC_lastPosY);
    SafeSetBranch(t, "TPC_lastPosZ", TPC_lastPosZ);
    SafeSetBranch(t, "TPC_PathLength", TPC_PathLength);
    SafeSetBranch(t, "TPC_dEdx", TPC_dEdx);
    // SafeSetBranch(t, "TPC_Psm", TPC_Psm);
    SafeSetBranch(t, "TPC_TrueKE", TPC_TrueKE);
    SafeSetBranch(t, "TPC_pdg", TPC_pdg);

    Long64_t nentries = t->GetEntries();

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    std::vector<Hit> hits;
    std::vector<Hit> trueHits;
    std::vector<Cluster> clusters;
    std::vector<TruePhoton> truePhotons;
    std::vector<ChargedTrack> ChargedTracks;
    std::vector<ChargedCluster> chargedClusters;
    //Pi0 analysis objects
    TH1F *hPi0Mass                  = nullptr;
    TH2F *hPi0ppM_pre               = nullptr;
    TH2F *hPi0ppM_post              = nullptr;
    TH1F *hClusterE                 = nullptr;
    // TH1F *hNClusters                = nullptr;
    // TH1F *hNClusters_lowEkin        = nullptr;
    // TH1F *hNClusters_midEkin        = nullptr;
    // TH1F *hNClusters_highEkin       = nullptr;
    // TH1F *hSingleClusterE           = nullptr;
    TH1F *hPi0TrueMass              = nullptr;
    TH1F *h_mass_truthE_recoAngle   = nullptr;
    TH1F *h_mass_recoE_truthAngle   = nullptr;
    TH1F *hEffvsE                   = nullptr;

    // Charged analysis objects
    TH1F  *hNSigmaPion              = nullptr;
    TH1F  *hNSigmaProton            = nullptr;
    // TH1F  *hNSigmaElectron       = nullptr;
    TH2F  *hdEdxVsE_cluster_Pion    = nullptr;
    TH2F  *hdEdxVsE_true_Pion       = nullptr;
    TH2F  *hdEdxVsE_cluster_Proton  = nullptr;
    TH2F  *hdEdxVsE_true_Proton     = nullptr;
    TH2F  *h2_Eres                  = nullptr;
    TH1F  *hdEdxTruePion            = nullptr;
    TH1F  *hdEdxSmearPion           = nullptr;
    TH1F  *hdEdxTrueProton          = nullptr;
    TH1F  *hdEdxSmearProton         = nullptr;
    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    // Analysis helper classes (also pointers)
    Pi0Efficiency  *effPlotter           = nullptr;
    Pi0Acceptance  *accPlotter           = nullptr;
    Pi0Acceptance  *pi0AcceptanceVsEta   = nullptr;
    Pi0Acceptance  *pi0AcceptanceVsTheta = nullptr;
    PIDEfficiency  *pidEff               = nullptr;

    if (doPi0Analysis) {
        hPi0Mass                = new TH1F("hPi0Mass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
        hPi0ppM_pre             = new TH2F("hPi0ppM_pre",";M_{#gamma #gamma} [MeV];#theta (#gamma #gamma)",200, 0, 250,200, 0, 4);
        hPi0ppM_post            = new TH2F("hPi0ppM_post",";M_{#gamma #gamma} [MeV];#theta (#gamma #gamma)",200, 0, 250,200, 0, 4);
        // hClusterE               = new TH1F("hClusterE",";Cluster E [MeV];Count",100,0,500);
        // hNClusters              = new TH1F("hNClusters",";N_{clusters};Events",20,0,20);
        // hNClusters_lowEkin      = new TH1F("hNClusters",";N_{clusters};Events",10,0,10);
        // hNClusters_midEkin      = new TH1F("hNClusters",";N_{clusters};Events",10,0,10);
        // hNClusters_highEkin     = new TH1F("hNClusters",";N_{clusters};Events",10,0,10);
        // hSingleClusterE         = new TH1F("hSingleClusterE",";Cluster E [MeV];Count",100,0,500);
        hEffvsE                 = new TH1F("hEffvsE", ";#pi^0 E_{kin}; Efficiency", 100, 1, 500);
        // effPlotter              = new Pi0Efficiency(120.0, 150.0, 134.977, 20, 1, 500);
        effPlotter              = new Pi0Efficiency(120.0, 150.0, 134.977, 4, 1, 550);
    //     accPlotter              = new Pi0Acceptance(120.0, 150.0, 134.977, 4, 1, 550);
    //     pi0AcceptanceVsEta      = new Pi0Acceptance(-10, 10, 100);
    //     pi0AcceptanceVsTheta    = new Pi0Acceptance(0, TMath::Pi(), 60);
        accPlotter              = new Pi0Acceptance("E", 120, 150, 134.977, 4, 1, 550);
        pi0AcceptanceVsEta      = new Pi0Acceptance("eta", -10, 10, 100);
        pi0AcceptanceVsTheta    = new Pi0Acceptance("theta", 0, TMath::Pi(), 60);    
    }
    if (doTruthAnalysis) {
        hPi0TrueMass            = new TH1F("hPi0TrueMass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
        h_mass_truthE_recoAngle = new TH1F("h_tE_rA",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
        h_mass_recoE_truthAngle = new TH1F("h_rE_tA",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
    }
    if (doChargedAnalysis) {
        hNSigmaPion             = new TH1F("hNSigmaPion", ";n#sigma;Counts", 100, -5, 5);
        hNSigmaProton           = new TH1F("hNSigmaProton", ";n#sigma;Counts", 100, -5, 5);
        hdEdxVsE_cluster_Pion   = new TH2F("hdEdxVsE_cluster_Pion",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
        hdEdxVsE_true_Pion      = new TH2F("hdEdxVsE_true_Pion",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
        hdEdxVsE_cluster_Proton = new TH2F("hdEdxVsE_cluster_Proton",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
        hdEdxVsE_true_Proton    = new TH2F("hdEdxVsE_true_Proton",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);   
        pidEff                  = new PIDEfficiency(20, 0, 500);
        h2_Eres                 = new TH2F("h2_Eres", ";True KE [MeV];Energy Residual", 20, 0, 500, 100, -1.0, 1.0);
        hdEdxTruePion           = new TH1F("hdEdxTruePion", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        hdEdxSmearPion          = new TH1F("hdEdxSmearPion", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        hdEdxTrueProton         = new TH1F("hdEdxTrueProton", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        hdEdxSmearProton        = new TH1F("hdEdxSmearProton", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        hClusterE               = new TH1F("hClusterE",";Cluster E [MeV];Count",100,0,500);
        // TH1F *hNSigmaElectron = new TH1F("hNSigma", ";n#sigma;Counts", 100, -5, 5);
        // TH2F *hdEdxVsE = new TH2F("hdEdxVsE", ";E [MeV];dEdx", 200, 0, 1000);
        // TH2F* hdEdxVsE_cluster_Electron = new TH2F("hdEdxVsE_cluster",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
        // TH2F* hdEdxVsE_true_Electron = new TH2F("hdEdxVsE_true",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
    }


    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        t->GetEntry(ievt);

        // Get the event's vertex
        if (!primaryX||primaryX->empty()) continue;
        TVector3 vertex((*primaryX)[0],(*primaryY)[0],(*primaryZ)[0]);

        if (!primaryEkin || primaryEkin->empty()) continue; 
        double genEkin = primaryEkin->at(0); 

        hits.clear();
        if (!energies || energies->empty()) {
            std::cout << "Event " << ievt << ": No hits (skipping)" << std::endl;
            continue;
        }
        size_t nHits = energies->size();
        // std::cout << "Event " << ievt << ": " << nHits << " total hits across all rings" << std::endl;
        for (size_t k=0; k<nHits; ++k) {
            hits.push_back({(*centerXs)[k], (*centerYs)[k], (*centerZs)[k], (*energies)[k]});
        }

        if (doChargedAnalysis) {
            ChargedTracks.clear();
            size_t nChargedTracks = TPC_Edep->size();
            for (size_t k=0; k<nChargedTracks; ++k) {
                // for noe the resolution is hardcoded to be 0.15
                ChargedTracks.push_back(
                    {k, 
                    vertex, 
                    TVector3((*TPC_lastPosX)[k], (*TPC_lastPosY)[k], (*TPC_lastPosZ)[k]), 
                    TVector3((*TPC_lastPosX)[k], (*TPC_lastPosY)[k], (*TPC_lastPosZ)[k]) - TVector3((*TPC_firstPosX)[k], (*TPC_firstPosY)[k], (*TPC_firstPosZ)[k]),
                    // TVector3((*TPC_firstPosX)[k], (*TPC_firstPosY)[k], (*TPC_firstPosZ)[k]) - vertex, 
                    // (*TPC_Edep)[k], (*TPC_PathLength)[k], (*TPC_dEdx)[k], 0.15});
                    (*TPC_TrueKE)[k], (*TPC_pdg)[k], (*TPC_dEdx)[k], (*TPC_smearedEdep)[k], (*TPC_PathLength)[k], 0 /* Placeholder */, 0.15});
            }
        }

        if (doTruthAnalysis) {
            trueHits.clear();
            size_t nTrueHits = truthE->size();
            for (size_t k=0; k<nTrueHits; ++k) {
                trueHits.push_back({(*truthPosX)[k], (*truthPosY)[k], (*truthPosZ)[k], (*truthE)[k]});
            }
            
            truePhotons = TruePhotonBuilder(trueHits, vertex);
            
            if (truePhotons.size() == 2) {
                    TLorentzVector diphoton = truePhotons[0].p4 + truePhotons[1].p4;
                    if (hPi0TrueMass) hPi0TrueMass->Fill(diphoton.M());
                }
        }

        // // Charged Object Clustering 

        if (doChargedAnalysis) {
            double thetaMax = 25.0 * TMath::DegToRad();
            chargedClusters = MatchHitsToTracks(ChargedTracks, hits, thetaMax);
            for (ChargedCluster cluster : chargedClusters) {
                // std::cout << "PID Guess: " << PIDToString(cluster.pidGuess) << std::endl;
                // std::cout << "Charged Cluster Energy: " << cluster.totalEnergy << std::endl;
                // std::cout << "Charged Cluster nSigma: " << cluster.nSigma << std::endl;
                if (hNSigmaPion) hNSigmaPion->Fill(cluster.nSigmaPion);
                if (hNSigmaProton) hNSigmaProton->Fill(cluster.nSigmaProton);
                // hNSigmaElectron->Fill(cluster.nSigmaElectron);
                if (cluster.objectTruePDG == 211) {
                if (hdEdxVsE_cluster_Pion) hdEdxVsE_cluster_Pion->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y 
                if (hClusterE) hClusterE->Fill(cluster.totalEnergy);
                // if (hdEdxVsE_cluster_Pion) hdEdxVsE_cluster_Pion->Fill(cluster.objectTrueKE, cluster.clusterdEdx); // ORDER: X vs Y 
                if (hdEdxVsE_true_Pion) hdEdxVsE_true_Pion->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                if (hdEdxTruePion) hdEdxTruePion->Fill(cluster.objectTruedEdx);
                if (hdEdxSmearPion) hdEdxSmearPion->Fill(cluster.clusterdEdx);
                }
                if (cluster.objectTruePDG == 2212) {
                if (hdEdxVsE_cluster_Proton) hdEdxVsE_cluster_Proton->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
                // if (hdEdxVsE_cluster_Proton) hdEdxVsE_cluster_Proton->Fill(cluster.objectTrueKE, cluster.clusterdEdx); // ORDER: X vs Y 
                if (hdEdxVsE_true_Proton) hdEdxVsE_true_Proton->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                if (hdEdxTrueProton) hdEdxTrueProton->Fill(cluster.objectTruedEdx);
                if (hdEdxSmearProton) hdEdxSmearProton->Fill(cluster.clusterdEdx);            
                }
                // if (cluster.objectTruePDG == 11) {
                // // hdEdxVsE_cluster_Electron->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
                // hdEdxVsE_cluster_Electron->Fill(cluster.totalEnergy, cluster.objectTruedEdx); // ORDER: X vs Y 
                // hdEdxVsE_true_Electron->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                // }
                double Eres = (cluster.totalEnergy - cluster.objectTrueKE) / cluster.objectTrueKE;
                if (h2_Eres) h2_Eres->Fill(cluster.objectTrueKE, Eres);
            }            
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
        //DEBUG

        // Double check to see if ownership logic is working! 
        // Only uncomment if you want to check!
        // for (const auto& h : hits) {
        //     if (h.owner == HitOwner::Charged) {
        //         std::cout << "Charged reco claimed this hit" << std::endl;
        //     }
        //     if (h.owner == HitOwner::Neutral) {
        //         std::cout << "Neutral reco claimed this hit" << std::endl;
        //     }
        //     if (h.owner == HitOwner::None) {
        //         std::cout << "No one claimed this hit" << std::endl;
        //     }
        // }
        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Neutral Object Clustering --> To be done after the Charged Object Clustering 

        if (doPi0Analysis) {
            // double totalE_Evt = 0;
            clusters.clear();

            // double dEta = 0.10/2;
            // double dPhi = 0.10/2;
            // double E_seed = 15.00;
            // double E_neighbor = 0.03;
            // int winSize = 7;
            // clusters = SlidingWindowClusterHits(hits, vertex, dEta, dPhi, E_seed, E_neighbor, winSize);
            
            clusters = clusterNeutralHits(hits, vertex, 25 * TMath::DegToRad()); // Have to optimise the angle a bit! 
            // Apply cluster energy threshold
            clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                        [](const Cluster &c){ return c.p4.E() < 50.0; }),
                        clusters.end());
        }


        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // DOUBLE CHECK WHAT I WANT TO DO WITH THIS STUFF
        //Cluster num plots

        // for (size_t ci=0; ci<clusters.size(); ++ci) {
        //     hClusterE->Fill(clusters[ci].p4.E());
        // }
        // hNClusters->Fill(clusters.size());

        // // Fill cluster num plots based on primary Ekin
        // if (genEkin < 200) hNClusters_lowEkin->Fill(clusters.size());
        // else if (genEkin >= 200 && genEkin < 400) hNClusters_midEkin->Fill(clusters.size());
        // else hNClusters_highEkin->Fill(clusters.size());

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
       
        // // Invariant mass plot (neutral pion / double photon)

        if (doPi0Analysis) {

            // constexpr double mpi0      = 135.0;   // MeV
            constexpr double mMin      = 90.0;    // MeV
            constexpr double mMax      = 180.0;   // MeV

            constexpr double thetaMin  = 0.1;     // rad
            constexpr double thetaMax  = 3.1;     // rad

            std::cout << "\n=== New event ===" << std::endl;
            std::cout << "Number of clusters: " << clusters.size() << std::endl;

            std::vector<Pi0Candidate> candidates;

            for (size_t a = 0; a < clusters.size(); ++a) {
                for (size_t b = a + 1; b < clusters.size(); ++b) {

                    const auto& g1 = clusters[a].p4;
                    const auto& g2 = clusters[b].p4;

                    double E1 = g1.E();
                    double E2 = g2.E();


                    double mgg = (g1 + g2).M();
                    double theta = openingAngle(g1, g2);

                    if (hPi0ppM_pre) hPi0ppM_pre->Fill(mgg, theta);

                    std::cout << "Pair (" << a << "," << b << "): "
                            << "E1=" << E1
                            << " E2=" << E2
                            << " mgg=" << mgg
                            << " theta=" << theta << std::endl;

                    // if (mgg < mMin || mgg > mMax) {
                    //     std::cout << "  -> reject: mass window" << std::endl;
                    //     continue;
                    // }

                    // if (theta < thetaMin || theta > thetaMax) {
                    //     std::cout << "  -> reject: angle window" << std::endl;
                    //     continue;
                    // }

                    // if (std::pow(mgg - 1.5, 2) / (9 / 4) + std::pow(theta - 135, 2) / (45 * 45) > 1) {
                    //     std::cout << "  -> rejected by ellipse condition" << std::endl;
                    //     continue;
                    // }

                    std::cout << "  -> ACCEPTED candidate" << std::endl;

                    Pi0Candidate cand;
                    cand.i = a;
                    cand.j = b;
                    cand.mgg = mgg;
                    cand.theta = theta;
                    cand.p4 = g1 + g2;

                    candidates.push_back(cand);
                }
            }

            std::cout << "Total candidates after cuts: "
                    << candidates.size() << std::endl;

            // ------------------------------------
            // Greedy selection (no photon reuse)
            // ------------------------------------
            std::vector<bool> used(clusters.size(), false);
            std::vector<Pi0Candidate> selected;

            for (const auto& cand : candidates) {

                std::cout << "Trying candidate ("
                        << cand.i << "," << cand.j << ")" << std::endl;

                if (used[cand.i] || used[cand.j]) {
                    std::cout << "  -> rejected: photon already used" << std::endl;
                    continue;
                }

                used[cand.i] = true;
                used[cand.j] = true;

                std::cout << "  -> SELECTED π0!" << std::endl;
                selected.push_back(cand);
            }

            std::cout << "Selected π0s in event: "
                    << selected.size() << std::endl;

            // ------------------------------------
            // Fill histograms
            // ------------------------------------
            for (const auto& pi0 : selected) {
                std::cout << "Filling π0 mass: "
                        << pi0.p4.M() << std::endl;
                if (hPi0Mass && hPi0ppM_post) {
                    hPi0Mass->Fill(pi0.p4.M());
                    hPi0ppM_post->Fill(pi0.mgg, pi0.theta);
                }
            }
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Truth level + reco level plots 

        if (doTruthAnalysis) {
            std::vector<int> clusterToTrue = matchClustersToTruth(clusters, truePhotons, 10);

            // Construct hybrid TLorentzVectors
            std::vector<TLorentzVector> photons_tE_rA; // truth E + reco angle
            std::vector<TLorentzVector> photons_rE_tA; // reco E + truth angle

            for (size_t ic=0; ic<clusters.size(); ++ic) {
                int it = clusterToTrue[ic];
                if (it < 0) continue; // skip unmatched

                const Cluster &c = clusters[ic];
                const TruePhoton &t = truePhotons[it];

                // truth energy + reco angle
                photons_tE_rA.push_back(makePhotonFromEnergyAndDir(t.p4.E(), c.centroid));

                // reco energy + true angle
                photons_rE_tA.push_back(makePhotonFromEnergyAndDir(c.p4.E(), t.dir));
            }

            // Compute all pairwise invariant masses
            auto fillPairs = [](const std::vector<TLorentzVector>& phs, TH1F* hist) {
                for (size_t i=0; i<phs.size(); ++i) {
                    for (size_t j=i+1; j<phs.size(); ++j) {
                        TLorentzVector sum = phs[i] + phs[j];
                        hist->Fill(sum.M());
                    }
                }
            };

            if (h_mass_truthE_recoAngle) fillPairs(photons_tE_rA, h_mass_truthE_recoAngle);
            if (h_mass_recoE_truthAngle) fillPairs(photons_rE_tA, h_mass_recoE_truthAngle);  
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        if (doPi0Analysis) {
            if (effPlotter) effPlotter->ProcessEvent(clusters, truePhotons);
            if (accPlotter) accPlotter->ProcessEvent(clusters, truePhotons, genEkin);

            // double px = primaryPx->at(0);
            // double py = primaryPy->at(0);
            // double pz = primaryPz->at(0);
            // double p = sqrt(px*px + py*py + pz*pz);

            // double eta = 0.5 * log((p + pz) / (p - pz)); // pseudorapidity
            // double theta = acos(pz / p); // theta 

            // pi0AcceptanceVsEta.ProcessEvent(clusters, truePhotons, eta);
            // pi0AcceptanceVsTheta.ProcessEvent(clusters, truePhotons, theta);

            // pi0AcceptanceVsTheta.ProcessEventTwoHist(clusters, truePhotons, genEkin, theta);

            // double targetTheta = TMath::Pi() / 2;
            // double deltaTheta = 0.1;

            // if (genEkin == 50 && std::abs(theta - targetTheta) < deltaTheta) {
            //     // std::cout << "[Theta] " << theta << std::endl;
            //     hNClusters_lowEkin->Fill(clusters.size());
            // }
            // else if (genEkin == 500 && std::abs(theta - targetTheta) < deltaTheta) hNClusters_highEkin->Fill(clusters.size());
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        if (doChargedAnalysis) {
            pidEff->ProcessEvent(chargedClusters);
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    }
    

    if (doPi0Analysis) {
        // RECO-RECO
        PlotOptions optsPi0InvMass;
        optsPi0InvMass.doFit = true;
        optsPi0InvMass.fitMin = 100;
        optsPi0InvMass.fitMax = 170;
        optsPi0InvMass.addInfoPave = true;
        optsPi0InvMass.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        // PrettyPi0MassPlot(hPi0Mass, "Pi0Mass_Clustered.png", 100.0, 170.0);
        Plot1D({hPi0Mass}, {kBlack}, "Pi0InvMass.png", optsPi0InvMass);

        PlotOptions optsPi0ppM;
        Plot2D(hPi0ppM_pre, "Pi0ppM_pre.png", optsPi0ppM);
        Plot2D(hPi0ppM_post, "Pi0ppM_post.png", optsPi0ppM);



        // CLUSTER NUM &/or DEBUG PLOTS
        // PlotOptions optsPi0NumCluster;
        // optsPi0NumCluster.legendEntries = {"Pi0 n Clusters"};
        // PrettyPi0NumClusterPlot(hNClusters);
        // Plot1D({hNClusters}, {kBlack}, "Pi0NumClusters.png", optsPi0NumCluster);
        // Pi0ClusterNumPlotEkin(hNClusters_lowEkin, hNClusters_midEkin, hNClusters_highEkin);
        // Pi0ClusterNumPlotEkin(hNClusters_lowEkin, hNClusters_highEkin);
        // BasicHistPlot(hSingleClusterE);

        // GOTTA CHECK IF POINTER/CLASS INTERPLAY IS WORKING HERE...
        // Efficiency Plot
        // EffPlot(hEff, "Pi0_eff.png"); // OLD
        effPlotter->FinalizePlot("Pi0_efficiency_vs_Ekin.png");

        // // Acceptance Plot
        accPlotter->FinalizePlot("Pi0_acceptance_vs_Ekin.png");
        // pi0AcceptanceVsEta.FinalizePlot("plots/Pi0_acceptance_vs_Eta.png");
        // pi0AcceptanceVsTheta.FinalizePlot("plots/Pi0_acceptance_vs_Theta_50_500.png");

        //CLEANUP
        delete hPi0Mass; 
        // delete hClusterE;
        // delete hNClusters; 
        // delete hNClusters_lowEkin; 
        // delete hNClusters_midEkin; 
        // delete hNClusters_highEkin; 
        delete hEffvsE;
        delete effPlotter;
        delete accPlotter;
        delete pi0AcceptanceVsEta; 
        delete pi0AcceptanceVsTheta; 
    }

    if (doTruthAnalysis) {
        // // TRUTH-TRUTH
        PlotOptions optsPi0InvMassTT;
        optsPi0InvMassTT.addInfoPave = true;
        optsPi0InvMassTT.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        optsPi0InvMassTT.legendEntries = {"Truth-level Invariant Mass"};
        // optsPi0InvMassTT.extraLegendLines = {Form("Mean value: %.4f", hPi0TrueMass->GetXaxis()->GetBinCenter(hPi0TrueMass->GetMaximumBin() + 1))};
        optsPi0InvMassTT.extraLegendLines = {"Mean Value: ~135 MeV"};
        // TruthPi0MassPlot(hPi0TrueMass, "Pi0Mass_Truth.png");
        if (hPi0TrueMass) Plot1D({hPi0TrueMass}, {kBlack}, "Pi0InvMassTT.png", optsPi0InvMassTT);

        // Mix Plots
        PlotOptions optsPi0InvMassRecoAngle;
        optsPi0InvMassRecoAngle.doFit = true;
        optsPi0InvMassRecoAngle.fitMin = 80;
        optsPi0InvMassRecoAngle.fitMax = 180;
        optsPi0InvMassRecoAngle.addInfoPave = true;
        optsPi0InvMassRecoAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        // PrettyPi0MassPlot(h_mass_truthE_recoAngle, "Pi0Mass_truthE_recoAngle.png", 100.0, 170.0);
        if (h_mass_truthE_recoAngle) Plot1D({h_mass_truthE_recoAngle}, {kBlack}, "Pi0InvMassRecoAngle.png", optsPi0InvMassRecoAngle);

        PlotOptions optsPi0InvMassTruthAngle;
        optsPi0InvMassTruthAngle.addInfoPave = true;
        optsPi0InvMassTruthAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        optsPi0InvMassTruthAngle.legendEntries = {"Truth-level Invariant Mass"};
        optsPi0InvMassTruthAngle.extraLegendLines = {Form("Mean value: %.4f", h_mass_recoE_truthAngle->GetXaxis()->GetBinCenter(h_mass_recoE_truthAngle->GetMaximumBin() + 1))};
        // TruthPi0MassPlot(h_mass_recoE_truthAngle, "Pi0Mass_recoE_truthAngle.png");
        if (h_mass_recoE_truthAngle) Plot1D({h_mass_recoE_truthAngle}, {kBlack}, "Pi0InvMassTruthAngle.png", optsPi0InvMassTruthAngle);

        // delete hPi0TrueMass;
        delete h_mass_truthE_recoAngle; 
        delete h_mass_recoE_truthAngle; 
    }

    // PID Plots 


    if (doChargedAnalysis) {
        PlotOptions opts_nSigma;
        opts_nSigma.doFit = true;
        opts_nSigma.addLegend = true;
        // opts_nSigma.legendEntries = {
        // "#pi hypothesis",
        // "p hypothesis",
        // "e hypothesis"
        // };
        opts_nSigma.legendEntries = {
        "#pi hypothesis",
        "p hypothesis"
        };
        // std::vector<TH1*> plots1D = {hNSigmaPion, hNSigmaProton, hNSigmaElectron};
        std::vector<TH1*> plots1D = {hNSigmaPion, hNSigmaProton};
        std::vector<int> colors = {
            kRed+1,
            // kGreen+2
            kBlue+1
        };
        Plot1D(plots1D, colors, "nSigmaPlots.png", opts_nSigma);

        // PlotOptions opts_hdEdxVsE_cluster;
        // opts_hdEdxVsE_cluster.addLegend = false;
        // opts_hdEdxVsE_cluster.addInfoPave = false;
        // // opts.drawOption = "COLZ";  // If you want color instead of HIST
        // Plot2D(hdEdxVsE_cluster, "dEdxVsE_clusterE.png", opts_hdEdxVsE_cluster);

        // PlotOptions opts_hdEdxVsE_true;
        // opts_hdEdxVsE_true.addLegend = false;
        // opts_hdEdxVsE_true.addInfoPave = false;
        // // opts.drawOption = "COLZ";  // If you want color instead of HIST
        // Plot2D(hdEdxVsE_true, "dEdxVsE_trueKE.png", opts_hdEdxVsE_true);

        // PlotOptions opts_hdEdxVsE_true;
        // opts_hdEdxVsE_true.addLegend = true;
        // opts_hdEdxVsE_true.legendEntries = {"$\frac{dE}{dx}$ Pions", "$\frac{dE}{dx}$ Protons"};
        // opts_hdEdxVsE_true.addInfoPave = false;
        // std::vector<int> colors_dEdx = {kBlack+1, kRed+1};
        // // opts.drawOption = "COLZ";
        // opts_hdEdxVsE_true.drawOption = "E";
        // Plot1D({hdEdxVsE_true_Pion->ProfileX(), hdEdxVsE_true_Proton->ProfileX()}, colors_dEdx, "dEdxVsE_trueKE.png", opts_hdEdxVsE_true); 

        PlotOptions opts_hdEdxVsE_true;
        opts_hdEdxVsE_true.drawOption = "SCAT";
        opts_hdEdxVsE_true.legendEntries = {"#pi^{+}","p"};
        opts_hdEdxVsE_true.legendDrawOpt = "P";
        Plot2DOverlay({hdEdxVsE_true_Pion, hdEdxVsE_true_Proton}, {kBlack, kRed},"dedx_vs_E_overlay_true.png",opts_hdEdxVsE_true);

        PlotOptions opts_hdEdxVsE_cluster;
        opts_hdEdxVsE_cluster.drawOption = "SCAT";
        opts_hdEdxVsE_cluster.legendEntries = {"#pi^{+}","p"};
        opts_hdEdxVsE_cluster.legendDrawOpt = "P";
        // Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton, hdEdxVsE_cluster_Electron}, {kBlack, kRed, kGreen},"dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);
        Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton}, {kBlack, kRed},"dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);

        PlotOptions opts_h2_Eres;
        TProfile* pEres = h2_Eres->ProfileX();
        Plot1D({pEres}, {kBlack}, "energy_residual.png", opts_h2_Eres);

        PlotOptions opts_dEdxPlots;
        opts_dEdxPlots.addLegend = true;
        opts_dEdxPlots.legendEntries = {"True #pi", "Smeared #pi", "True p", "Smeared p"};
        std::vector<TH1*> plots1D_dEdx = {hdEdxTruePion, hdEdxSmearPion, hdEdxTrueProton, hdEdxSmearProton};
        std::vector<int> colors_dEdx = {
            kBlack,
            kRed+1,
            kGreen+2,
            kBlue+1
        };
        Plot1D(plots1D_dEdx, colors_dEdx, "dEdxPlots.png", opts_dEdxPlots);

        //Pion + 
        pidEff->FinalizePlot("plots/PIDEfficiency", 211);

        PlotOptions opts_hClusterE;
        opts_hClusterE.addLegend = true;
        opts_hClusterE.legendEntries = {"#pi"};
        Plot1D({hClusterE}, {kBlack}, "ClusterE_pion.png", opts_hClusterE);

        //CLEANUP
        delete hNSigmaPion;
        delete hNSigmaProton;
        // delete hNSigmaElectron;
        delete hdEdxVsE_cluster_Pion;
        delete hdEdxVsE_true_Pion; 
        delete hdEdxVsE_cluster_Proton;
        delete hdEdxVsE_true_Proton;
        delete hClusterE;
    }

    f->Close(); delete f;
    return 0;
}