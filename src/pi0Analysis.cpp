#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraph.h"

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
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    SetPrettyStyle();

    TFile *f = TFile::Open(argv[1]);
    if (!f || f->IsZombie()) return 1;

    TTree *t = (TTree *)f->Get("digitizedHits");
    if (!t) { std::cerr << "Tree digitizedHits not found" << std::endl; return 1; }

    // // Single consolidated vectors (matches preprocessed output—no per-ring arrays!)
    // std::vector<double> *centerXs = nullptr, *centerYs = nullptr, *centerZs = nullptr, *energies = nullptr;
    // t->SetBranchAddress("centerX", &centerXs);
    // t->SetBranchAddress("centerY", &centerYs);
    // t->SetBranchAddress("centerZ", &centerZs);
    // t->SetBranchAddress("energy", &energies);

    // // Primary vertex (vector, but size=1 per event)
    // std::vector<double> *primaryX=nullptr,*primaryY=nullptr,*primaryZ=nullptr,*primaryEkin=nullptr;
    // std::vector<double>* primaryPx = nullptr, *primaryPy = nullptr, *primaryPz = nullptr;
    // t->SetBranchAddress("PrimaryPosX",&primaryX);
    // t->SetBranchAddress("PrimaryPosY",&primaryY);
    // t->SetBranchAddress("PrimaryPosZ",&primaryZ);
    // t->SetBranchAddress("PrimaryEkin",&primaryEkin);
    // t->SetBranchAddress("PrimaryMomX", &primaryPx);
    // t->SetBranchAddress("PrimaryMomY", &primaryPy);
    // t->SetBranchAddress("PrimaryMomZ", &primaryPz);   

    // // Truth level info
    // std::vector<double> *truthPosX=nullptr, *truthPosY=nullptr, *truthPosZ=nullptr, *truthE=nullptr;
    // t->SetBranchAddress("truthPosX", &truthPosX);
    // t->SetBranchAddress("truthPosY", &truthPosY);
    // t->SetBranchAddress("truthPosZ", &truthPosZ);
    // t->SetBranchAddress("truthE", &truthE);

    // // TPC info
    // std::vector<double> *TPC_Edep=nullptr, *TPC_PosX=nullptr, *TPC_PosY=nullptr, *TPC_PosZ=nullptr;
    // std::vector<double> *TPC_PathLength=nullptr, *TPC_dEdx=nullptr, *TPC_Psm=nullptr;
    // t->SetBranchAddress("TPC_Edep", &TPC_Edep);
    // t->SetBranchAddress("TPC_PosX", &TPC_PosX);
    // t->SetBranchAddress("TPC_PosY", &TPC_PosY);
    // t->SetBranchAddress("TPC_PosZ", &TPC_PosZ);
    // t->SetBranchAddress("TPC_PathLength", &TPC_PathLength);
    // t->SetBranchAddress("TPC_dEdx_rho", &TPC_dEdx);
    // t->SetBranchAddress("TPC_Psm", &TPC_Psm);

    // Example usage in main
    // Single consolidated vectors
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
    std::vector<double> *TPC_Edep = nullptr, *TPC_PosX = nullptr, *TPC_PosY = nullptr, *TPC_PosZ = nullptr;
    std::vector<double> *TPC_PathLength = nullptr, *TPC_dEdx = nullptr, *TPC_Psm = nullptr, *TPC_TrueKE = nullptr, *TPC_pdg = nullptr;
    SafeSetBranch(t, "TPC_Edep", TPC_Edep);
    SafeSetBranch(t, "TPC_PosX", TPC_PosX);
    SafeSetBranch(t, "TPC_PosY", TPC_PosY);
    SafeSetBranch(t, "TPC_PosZ", TPC_PosZ);
    SafeSetBranch(t, "TPC_PathLength", TPC_PathLength);
    SafeSetBranch(t, "TPC_dEdx_rho", TPC_dEdx);
    SafeSetBranch(t, "TPC_Psm", TPC_Psm);
    SafeSetBranch(t, "TPC_TrueKE", TPC_TrueKE);
    SafeSetBranch(t, "TPC_pdg", TPC_pdg);

    Long64_t nentries = t->GetEntries();

    // TH1F *hPi0Mass = new TH1F("hPi0Mass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);

    // TH1F *hClusterE = new TH1F("hClusterE",";Cluster E [MeV];Count",100,0,500);
    // TH1F *hNClusters = new TH1F("hNClusters",";N_{clusters};Events",20,0,20);
    
    // TH1F *hNClusters_lowEkin = new TH1F("hNClusters",";N_{clusters};Events",10,0,10);
    // TH1F *hNClusters_midEkin = new TH1F("hNClusters",";N_{clusters};Events",10,0,10);
    // TH1F *hNClusters_highEkin = new TH1F("hNClusters",";N_{clusters};Events",10,0,10);

    
    // TH1F *hSingleClusterE = new TH1F("hSingleClusterE",";Cluster E [MeV];Count",100,0,500);
    // TH1F *hPi0TrueMass = new TH1F("hPi0TrueMass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
    // TH1F *h_mass_truthE_recoAngle = new TH1F("h_tE_rA",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
    // TH1F *h_mass_recoE_truthAngle = new TH1F("h_rE_tA",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
    // TH1F *hEffvsE = new TH1F("hEffvsE", ";#pi^0 E_{kin}; Efficiency", 100, 1, 500);

    //PID Plots

    TH1F *hNSigmaPion = new TH1F("hNSigma", ";n#sigma;Counts", 100, -5, 5);
    TH1F *hNSigmaProton = new TH1F("hNSigma", ";n#sigma;Counts", 100, -5, 5);
    // TH1F *hNSigmaElectron = new TH1F("hNSigma", ";n#sigma;Counts", 100, -5, 5);

    // TH2F *hdEdxVsE = new TH2F("hdEdxVsE", ";E [MeV];dEdx", 200, 0, 1000);
    TH2F* hdEdxVsE_cluster_Pion = new TH2F("hdEdxVsE_cluster",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
    TH2F* hdEdxVsE_true_Pion = new TH2F("hdEdxVsE_true",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
    TH2F* hdEdxVsE_cluster_Proton = new TH2F("hdEdxVsE_cluster",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
    TH2F* hdEdxVsE_true_Proton = new TH2F("hdEdxVsE_true",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
    // TH2F* hdEdxVsE_cluster_Electron = new TH2F("hdEdxVsE_cluster",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
    // TH2F* hdEdxVsE_true_Electron = new TH2F("hdEdxVsE_true",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);

    // Pi0Efficiency effPlotter(120.0, 150.0, 134.977, 20, 1, 500);
    // Pi0Efficiency effPlotter(120.0, 150.0, 134.977, 4, 1, 550);

    // Pi0Acceptance accPlotter(120.0, 150.0, 134.977, 4, 1, 550);
    // Pi0Acceptance pi0AcceptanceVsEta(-10, 10, 100);
    // Pi0Acceptance pi0AcceptanceVsTheta(0, TMath::Pi(), 60);

    PIDEfficiency pidEff(20, 0, 500);

    TH2F* h2_Eres = new TH2F("h2_Eres", ";True KE [MeV];Energy Residual", 20, 0, 500, 100, -1.0, 1.0);

    
    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        t->GetEntry(ievt);

        // Get the event's vertex
        if (!primaryX||primaryX->empty()) continue;
        TVector3 vertex((*primaryX)[0],(*primaryY)[0],(*primaryZ)[0]);

        if (!primaryEkin || primaryEkin->empty()) continue; 
        double genEkin = primaryEkin->at(0); 

        std::vector<Hit> hits;
        if (!energies || energies->empty()) {
            std::cout << "Event " << ievt << ": No hits (skipping)" << std::endl;
            continue;
        }
        size_t nHits = energies->size();
        // std::cout << "Event " << ievt << ": " << nHits << " total hits across all rings" << std::endl;
        for (size_t k=0; k<nHits; ++k) {
            hits.push_back({(*centerXs)[k], (*centerYs)[k], (*centerZs)[k], (*energies)[k]});
        }

        std::vector<ChargedTrack> ChargedTracks;
        size_t nChargedTracks = TPC_Edep->size();
        for (size_t k=0; k<nChargedTracks; ++k) {
            // for noe the resolution is hardcoded to be 0.15
            ChargedTracks.push_back(
                {k, 
                vertex, 
                TVector3((*TPC_PosX)[k], (*TPC_PosY)[k], (*TPC_PosZ)[k]), 
                TVector3((*TPC_PosX)[k], (*TPC_PosY)[k], (*TPC_PosZ)[k]) - vertex,
                // (*TPC_Edep)[k], (*TPC_PathLength)[k], (*TPC_dEdx)[k], 0.15});
                (*TPC_TrueKE)[k], (*TPC_pdg)[k], (*TPC_Psm)[k], (*TPC_PathLength)[k], (*TPC_dEdx)[k], 0.15});
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // // TRUE HITS --> True photon hits:
        // std::vector<Hit> trueHits;
        // size_t nTrueHits = truthE->size();
        // // std::cout << "[DEBUG] Event " << ievt << ": " << nTrueHits << " total true hits" << std::endl;
        // for (size_t k=0; k<nTrueHits; ++k) {
        //     trueHits.push_back({(*truthPosX)[k], (*truthPosY)[k], (*truthPosZ)[k], (*truthE)[k]});
        // }

        
        // auto truePhotons = TruePhotonBuilder(trueHits, vertex);
        
        // if (truePhotons.size() == 2) {
        //         TLorentzVector diphoton = truePhotons[0].p4 + truePhotons[1].p4;
        //         // std::cout << "[DIPHOTON] " << diphoton.M() << std::endl;
        //         hPi0TrueMass->Fill(diphoton.M());
        //     }
            
        // DEBUG print bin contents
        // hPi0TrueMass->Print("all");

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // // Charged Object Clustering 

        double thetaMax = 5.0 * TMath::DegToRad();
        auto chargedClusters = MatchHitsToTracks(ChargedTracks, hits, thetaMax);

        for (ChargedCluster cluster : chargedClusters) {
            // std::cout << "PID Guess: " << PIDToString(cluster.pidGuess) << std::endl;
            // std::cout << "Charged Cluster Energy: " << cluster.totalEnergy << std::endl;
            // std::cout << "Charged Cluster nSigma: " << cluster.nSigma << std::endl;
            hNSigmaPion->Fill(cluster.nSigmaPion);
            hNSigmaProton->Fill(cluster.nSigmaProton);
            // hNSigmaElectron->Fill(cluster.nSigmaElectron);
            if (cluster.objectTruePDG == 211) {
            // hdEdxVsE_cluster_Pion->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y 
            hdEdxVsE_cluster_Pion->Fill(cluster.totalEnergy, cluster.objectTruedEdx); // ORDER: X vs Y 
            hdEdxVsE_true_Pion->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
            }
            if (cluster.objectTruePDG == 2212) {
            // hdEdxVsE_cluster_Proton->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
            hdEdxVsE_cluster_Proton->Fill(cluster.totalEnergy, cluster.objectTruedEdx); // ORDER: X vs Y 
            hdEdxVsE_true_Proton->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
            }
            // if (cluster.objectTruePDG == 11) {
            // // hdEdxVsE_cluster_Proton->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
            // hdEdxVsE_cluster_Electron->Fill(cluster.totalEnergy, cluster.objectTruedEdx); // ORDER: X vs Y 
            // hdEdxVsE_true_Electron->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
            // }
            double Eres = (cluster.totalEnergy - cluster.objectTrueKE) / cluster.objectTrueKE;
            h2_Eres->Fill(cluster.objectTrueKE, Eres);
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
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

        // std::vector<Cluster> clusters;
        // double dEta = 0.10/2;
        // double dPhi = 0.10/2;
        // double E_seed = 15.00;
        // double E_neighbor = 0.03;
        // int winSize = 7;

        // double totalE_Evt = 0;

        // clusters = SlidingWindowClusterHits(hits, vertex, dEta, dPhi, E_seed, E_neighbor, winSize);

        // // Apply cluster energy threshold
        // clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
        //                             [](const Cluster &c){ return c.p4.E() < 50.0; }),
        //             clusters.end());


        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

        // for (size_t a=0;a<clusters.size();++a) {
        //     for (size_t b=a+1;b<clusters.size();++b) {
        //         TLorentzVector pi0 = clusters[a].p4 + clusters[b].p4;
        //         hPi0Mass->Fill(pi0.M());
        //     }
        // }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Truth level + reco level plots 

        // std::vector<int> clusterToTrue = matchClustersToTruth(clusters, truePhotons, 10);

        // // Construct hybrid TLorentzVectors
        // std::vector<TLorentzVector> photons_tE_rA; // truth E + reco angle
        // std::vector<TLorentzVector> photons_rE_tA; // reco E + truth angle

        // for (size_t ic=0; ic<clusters.size(); ++ic) {
        //     int it = clusterToTrue[ic];
        //     if (it < 0) continue; // skip unmatched

        //     const Cluster &c = clusters[ic];
        //     const TruePhoton &t = truePhotons[it];

        //     // truth energy + reco angle
        //     photons_tE_rA.push_back(makePhotonFromEnergyAndDir(t.p4.E(), c.centroid));

        //     // reco energy + true angle
        //     photons_rE_tA.push_back(makePhotonFromEnergyAndDir(c.p4.E(), t.dir));
        // }

        // // Compute all pairwise invariant masses
        // auto fillPairs = [](const std::vector<TLorentzVector>& phs, TH1F* hist) {
        //     for (size_t i=0; i<phs.size(); ++i) {
        //         for (size_t j=i+1; j<phs.size(); ++j) {
        //             TLorentzVector sum = phs[i] + phs[j];
        //             hist->Fill(sum.M());
        //         }
        //     }
        // };

        // fillPairs(photons_tE_rA, h_mass_truthE_recoAngle);
        // fillPairs(photons_rE_tA, h_mass_recoE_truthAngle);  

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Eff and Acc Plots 
        
        // effPlotter.ProcessEvent(clusters, truePhotons);
        // accPlotter.ProcessEvent(clusters, truePhotons, genEkin);

        // double px = primaryPx->at(0);
        // double py = primaryPy->at(0);
        // double pz = primaryPz->at(0);
        // double p = sqrt(px*px + py*py + pz*pz);

        // double eta = 0.5 * log((p + pz) / (p - pz)); // pseudorapidity
        // double theta = acos(pz / p); // theta 

        // // pi0AcceptanceVsEta.ProcessEvent(clusters, truePhotons, eta);
        // // pi0AcceptanceVsTheta.ProcessEvent(clusters, truePhotons, theta);

        // // pi0AcceptanceVsTheta.ProcessEventTwoHist(clusters, truePhotons, genEkin, theta);

        // // double targetTheta = TMath::Pi() / 2;
        // // double deltaTheta = 0.1;

        // // if (genEkin == 50 && std::abs(theta - targetTheta) < deltaTheta) {
        // //     // std::cout << "[Theta] " << theta << std::endl;
        // //     hNClusters_lowEkin->Fill(clusters.size());
        // // }
        // // else if (genEkin == 500 && std::abs(theta - targetTheta) < deltaTheta) hNClusters_highEkin->Fill(clusters.size());

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        pidEff.ProcessEvent(chargedClusters);

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    }
    

    // // RECO-RECO
    // PlotOptions optsPi0InvMass;
    // optsPi0InvMass.doFit = true;
    // optsPi0InvMass.fitMin = 100;
    // optsPi0InvMass.fitMax = 170;
    // optsPi0InvMass.addInfoPave = true;
    // optsPi0InvMass.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
    // // PrettyPi0MassPlot(hPi0Mass, "Pi0Mass_Clustered.png", 100.0, 170.0);
    // Plot1D({hPi0Mass}, {kBlack}, "Pi0InvMass.png", optsPi0InvMass);

    // // TRUTH-TRUTH
    // PlotOptions optsPi0InvMassTT;
    // optsPi0InvMass.addInfoPave = true;
    // optsPi0InvMassTT.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
    // optsPi0InvMassTT.legendEntries = {"Truth-level Invariant Mass"};
    // optsPi0InvMassTT.extraLegendLines = {Form("Mean value: %.4f", hPi0Mass->GetXaxis()->GetBinCenter(hPi0Mass->GetMaximumBin() + 1))};
    // // optsPi0InvMassTT.extraLegendLines = {"Mean Value: ~135 MeV"};
    // // TruthPi0MassPlot(hPi0TrueMass, "Pi0Mass_Truth.png");
    // Plot1D({hPi0TrueMass}, {kBlack}, "Pi0InvMassTT.png", optsPi0InvMassTT);

    // // Mix Plots
    // PlotOptions optsPi0InvMassRecoAngle;
    // optsPi0InvMassRecoAngle.doFit = true;
    // optsPi0InvMassRecoAngle.fitMin = 80;
    // optsPi0InvMassRecoAngle.fitMax = 180;
    // optsPi0InvMass.addInfoPave = true;
    // optsPi0InvMassRecoAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
    // // PrettyPi0MassPlot(h_mass_truthE_recoAngle, "Pi0Mass_truthE_recoAngle.png", 100.0, 170.0);
    // Plot1D({h_mass_truthE_recoAngle}, {kBlack}, "Pi0InvMassRecoAngle.png", optsPi0InvMassRecoAngle);

    // PlotOptions optsPi0InvMassTruthAngle;
    // optsPi0InvMassTruthAngle.addInfoPave = true;
    // optsPi0InvMassTruthAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
    // optsPi0InvMassTruthAngle.legendEntries = {"Truth-level Invariant Mass"};
    // optsPi0InvMassTruthAngle.extraLegendLines = {Form("Mean value: %.4f", h_mass_recoE_truthAngle->GetXaxis()->GetBinCenter(h_mass_recoE_truthAngle->GetMaximumBin() + 1))};
    // // TruthPi0MassPlot(h_mass_recoE_truthAngle, "Pi0Mass_recoE_truthAngle.png");
    // Plot1D({h_mass_recoE_truthAngle}, {kBlack}, "Pi0InvMassTruthAngle.png", optsPi0InvMassTruthAngle);

    // // CLUSTER NUM &/or DEBUG PLOTS
    // PlotOptions optsPi0NumCluster;
    // optsPi0NumCluster.legendEntries = {"Pi0 n Clusters"};
    // // PrettyPi0NumClusterPlot(hNClusters);
    // Plot1D({hNClusters}, {kBlack}, "Pi0NumClusters.png", optsPi0NumCluster);
    // // Pi0ClusterNumPlotEkin(hNClusters_lowEkin, hNClusters_midEkin, hNClusters_highEkin);
    // // Pi0ClusterNumPlotEkin(hNClusters_lowEkin, hNClusters_highEkin);
    // // BasicHistPlot(hSingleClusterE);

    // // Efficiency Plot
    // // EffPlot(hEff, "Pi0_eff.png"); // OLD
    // effPlotter.FinalizePlot("Pi0_efficiency_vs_Ekin.png");

    // // // Acceptance Plot
    // accPlotter.FinalizePlot("Pi0_acceptance_vs_Ekin.png");
    // // pi0AcceptanceVsEta.FinalizePlot("plots/Pi0_acceptance_vs_Eta.png");
    // // pi0AcceptanceVsTheta.FinalizePlot("plots/Pi0_acceptance_vs_Theta_50_500.png");

    // PID Plots 
    // nSigmaPlot(hNSigma, "nSigma.png", -3, 3); //OLD
    // dEdxVsEPlot(hdEdxVsE, "dEdxVsE.png"); //OLD

    PlotOptions opts_nSigma;
    opts_nSigma.doFit = true;
    opts_nSigma.addLegend = true;
    opts_nSigma.legendEntries = {
    "#pi hypothesis",
    "p hypothesis",
    "e hypothesis"
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

    // PlotOptions opts_hdEdxVsE_true;
    // opts_hdEdxVsE_true.drawOption = "SCAT";
    // opts_hdEdxVsE_true.legendEntries = {"#pi^{+}","p"};
    // Plot2DOverlay({hdEdxVsE_true_Pion, hdEdxVsE_true_Proton}, {kBlack, kRed},"dedx_vs_E_overlay_true.png",opts_hdEdxVsE_true);

    // PlotOptions opts_hdEdxVsE_cluster;
    // opts_hdEdxVsE_cluster.drawOption = "SCAT";
    // opts_hdEdxVsE_cluster.legendEntries = {"#pi^{+}","p", "e"};
    // // Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton, hdEdxVsE_cluster_Electron}, {kBlack, kRed, kGreen},"dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);
    // Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton}, {kBlack, kRed},"dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);

    // PlotOptions opts_h2_Eres;
    // TProfile* pEres = h2_Eres->ProfileX();
    // Plot1D({pEres}, {kBlack}, "energy_residual.png", opts_h2_Eres);

    //Pion + 
    pidEff.FinalizePlot("plots/PIDEfficiency", 211);

    // delete hPi0Mass; 
    // delete hClusterE;
    // delete hNClusters; 
    // delete hNClusters_lowEkin; 
    // delete hNClusters_midEkin; 
    // delete hNClusters_highEkin; 
    // delete hPi0TrueMass;
    // delete h_mass_truthE_recoAngle; 
    // delete h_mass_recoE_truthAngle; 
    // delete hEffvsE;
    delete hNSigmaPion;
    delete hNSigmaProton;
    // delete hNSigmaElectron;
    delete hdEdxVsE_cluster_Pion;
    delete hdEdxVsE_true_Pion; 
    delete hdEdxVsE_cluster_Proton;
    delete hdEdxVsE_true_Proton;
    f->Close(); delete f;
    return 0;
}