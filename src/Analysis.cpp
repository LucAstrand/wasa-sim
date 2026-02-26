#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TError.h"
#include "TSystem.h"

#include <mcpl.h>

#include <algorithm>

#include "CLI11.hpp"
#include "progressbar.hpp"

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
#include "RecoEvent.hpp"
#include "Calibration.hpp"

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

void MyErrorHandler(int level, Bool_t abort, const char* location, const char* msg)
{
    DefaultErrorHandler(level, abort, location, msg);

    // only for the nbins warning:
    if (msg && strstr(msg, "nbins is <=0")) {
        std::cerr << "\n=== BACKTRACE for nbins<=0 ===\n";
        gSystem->StackTrace();
        std::cerr << "==============================\n\n";
    }
}
// ----------------------------------------------------

int main(int argc, char **argv) {
    if (argc < 3) {
            std::cout << "Usage: " << argv[0] << " <input.root> <mcpl_data.mcpl> [options]\n"
                    << "Options:\n"
                    << "  --calibration         Performs calibration procedure\n"
                    << "  --full-analysis       Performs everything\n"
                    << "  --pi0-analysis        Reconstruct Pi0s and invariant mass analysis\n"
                    << "  --charged-analysis    Reconstruct Charged objects and do PID studies\n"
                    << "  --truth-analysis      Performs truth level analysis for Pi0s\n"
                    << "  --truth-mix-analysis  Performs truth-reco level mixed analysis for Pi0s\n"
                    << "  --event-variables     Computes the event level variables\n";
            return 1;
        }
    SetErrorHandler(MyErrorHandler);
    
    std::string root_inputfile = argv[1];
    std::string mcpl_inputfile = argv[2];
    bool doCalibration = false;
    bool doPi0Analysis = false;
    bool doChargedAnalysis = false;
    bool doTruthAnalysis = false;
    bool doTruthAndMixPlots = false;
    bool doEventVariables = false;

    for (int i = 3; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--calibration") doCalibration = true;
        if (arg == "--full-analysis") {
            if (!std::filesystem::exists("chargedKE.root")) {
                std::cerr << "No calibration root file detected, please run './Analysis <input.root> <mcpl_data.mcpl> --calibration' first" << std::endl;
                return 1;
            }
            doPi0Analysis = true;
            doChargedAnalysis = true;
            doTruthAnalysis = true; // NOTE for the truth-reco mix plots you also need to have the flag "--truth-mix-analysis"
            doEventVariables = true;
        }
        if (arg == "--pi0-analysis") doPi0Analysis = true;
        if (arg == "--charged-analysis") doChargedAnalysis = true;
        if (arg == "--truth-analysis") doTruthAnalysis = true;
        if (arg == "--truth-mix-analysis") {
            doTruthAndMixPlots = true; 
            doTruthAnalysis = true;
        }
            if (arg == "--event-variables") doEventVariables = true;
    }

    SetPrettyStyle();

    TFile *f = TFile::Open(root_inputfile.c_str());
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
    std::vector<int> *primaryPDG = nullptr;
    SafeSetBranch(t, "PrimaryPosX", primaryX);
    SafeSetBranch(t, "PrimaryPosY", primaryY);
    SafeSetBranch(t, "PrimaryPosZ", primaryZ);
    SafeSetBranch(t, "PrimaryEkin", primaryEkin);
    SafeSetBranch(t, "PrimaryMomX", primaryPx);
    SafeSetBranch(t, "PrimaryMomY", primaryPy);
    SafeSetBranch(t, "PrimaryMomZ", primaryPz);
    SafeSetBranch(t, "PrimaryPDG", primaryPDG);

    // Truth info
    std::vector<double> *truePhotonPosX = nullptr, *truePhotonPosY = nullptr, *truePhotonPosZ = nullptr, *truePhotonE = nullptr;
    std::vector<int> *truePhotonTrackID = nullptr, *truePhotonParentID = nullptr;
    SafeSetBranch(t, "truthPosX", truePhotonPosX);
    SafeSetBranch(t, "truthPosY", truePhotonPosY);
    SafeSetBranch(t, "truthPosZ", truePhotonPosZ);
    SafeSetBranch(t, "truthE", truePhotonE);
    SafeSetBranch(t, "TruePhotonTrackID", truePhotonTrackID);
    SafeSetBranch(t, "TruePhotonParentID", truePhotonParentID);

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
    // mcpl pre-processing to get the #Pi0s per event:

    std::vector<int> pi0_per_event;
    std::vector<int> chPi_per_event;
    std::map<int, double> Ekin_per_event;
    mcpl_file_t mcplFile = mcpl_open_file(mcpl_inputfile.c_str());
    int current_event = -1;
    const mcpl_particle_t* p;
    while ((p = mcpl_read(mcplFile))) {
        
        // we might miss one if we proceed like this... because if it starts on 0 then the second event will be the first 1...
        // NOTE here it is fine with the userflags of "cleaned_HIBEAM_WASA_rwag_signal_GBL_jbar_100k_9012_newUF.mcpl" as they start on a 1...
        if (p->userflags == 1) { 
            ++current_event;
            pi0_per_event.push_back(0);
            chPi_per_event.push_back(0);
            // pi0_EKin_per_event.push_back(0);
        }
        
        //guard
        if (current_event < 0) continue;
        
        if (abs(p->pdgcode) == 211 || p->pdgcode == 111 || p->pdgcode == 22) {
            Ekin_per_event[current_event] += p->ekin; // Only primary mesons and photons (pi0, pi+- and gammas)
        }

        if (p->pdgcode == 111) {
            pi0_per_event[current_event]++;
            // pi0_EKin_per_event[current_event] = p->ekin;
        }
        if (std::abs(p->pdgcode) == 211) {
            chPi_per_event[current_event]++;
        }
    }

    mcpl_close_file(mcplFile);

    // Reminder: This gives us pi0_per_event[ievt] --> # Pi0s in MCPL ievt & corresponding EKin in other vector with same indexing.

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    
    //Containers used throughout the analysis 
    std::vector<Hit> hits;
    std::vector<TruePhotonHit> trueHits;
    // std::vector<Cluster> clusters;
    std::vector<TruePhoton> truePhotons;
    std::vector<TruePi0> truePi0s;
    std::vector<Pi0Candidate> selected;
    std::vector<ChargedTrack> chargedTracks;
    // std::vector<ChargedCluster> chargedClusters;
    // std::vector<ChargedObject> chargedObjects;
    std::vector<TVector3> momenta;
    std::vector<double> weights;
    //Pi0 analysis objects
    TH1F *hPi0Mass                       = nullptr;
    TH2F *hPi0ppM_pre                    = nullptr;
    TH2F *hPi0ppM_post                   = nullptr;
    //--> Truth level
    TH1F *hPi0TrueMass                   = nullptr;
    //-->//--> Truth-Reco
    TH1F *h_mass_truthE_recoAngle        = nullptr;
    TH1F *h_mass_recoE_truthAngle        = nullptr;
    //--> Pi0 Acceptance and Efficiency 
    Pi0Efficiency  *effPlotter           = nullptr;
    Pi0Acceptance  *accPlotter           = nullptr;
    Pi0Acceptance  *pi0AcceptanceVsEta   = nullptr;
    Pi0Acceptance  *pi0AcceptanceVsTheta = nullptr;
    // Charged analysis objects
    //--> Pions
    TH1F  *hNSigmaPion                   = nullptr;
    TH2F  *hdEdxVsE_cluster_Pion         = nullptr;
    TH2F  *hdEdxVsE_true_Pion            = nullptr;
    TH1F  *hdEdxTruePion                 = nullptr;
    TH1F  *hdEdxSmearPion                = nullptr;
    //--> Protons
    TH1F  *hNSigmaProton                 = nullptr;
    TH2F  *hdEdxVsE_cluster_Proton       = nullptr;
    TH2F  *hdEdxVsE_true_Proton          = nullptr;
    TH1F  *hdEdxTrueProton               = nullptr;
    TH1F  *hdEdxSmearProton              = nullptr;
    //--> Electrons
    // TH1F  *hNSigmaElectron               = nullptr;
    // TH2F  *hdEdxVsE_cluster_Electron     = nullptr;
    // TH2F  *hdEdxVsE_true_Electron        = nullptr;
    // TH1F  *hdEdxTrueElectron             = nullptr;
    // TH1F  *hdEdxSmearElectron            = nullptr;
    //--> Global Charged / PID related
    TH2F  *h2_Eres                       = nullptr;
    TH1F *hClusterE                      = nullptr;
    PIDEfficiency  *pidEff               = nullptr;
    // Event level variables
    TH1F  *hEventInvariantMass           = nullptr; 
    TH1F  *hEventSphericity              = nullptr; 
    //--> Total energy related
    TH1F  *hEvis                         = nullptr;
    TH1F  *hEventCorrectedTotE           = nullptr; 
    TH2F  *hErecoVsEtrue                 = nullptr;
    TH2F  *hEvisVsEtrue                  = nullptr;
    TH1F  *hDiffErecoEtrue               = nullptr;
    //--> Pion multiplicity related
    TH2I *hNPionMultiplicity             = nullptr;
    TH2I *hChPionMultiplicity            = nullptr;

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    //Calibration procedure
    if (doCalibration) {
        DoCalibration(
            t,
            static_cast<int>(nentries * 0.7),
            pi0_per_event,
            centerXs, centerYs, centerZs, energies,
            primaryX, primaryY, primaryZ, primaryEkin,
            TPC_firstPosX, TPC_firstPosY, TPC_firstPosZ,
            TPC_lastPosX,  TPC_lastPosY,  TPC_lastPosZ,
            TPC_TrueKE, TPC_pdg, TPC_dEdx,
            TPC_smearedEdep, TPC_PathLength,
            "chargedKE.root"
        );

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


    if (doPi0Analysis) {
        hPi0Mass                = new TH1F("hPi0Mass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
        hPi0ppM_pre             = new TH2F("hPi0ppM_pre",";M_{#gamma #gamma} [MeV];#theta_{#gamma #gamma} [rad]",200, 0, 250,200, 0, 4);
        hPi0ppM_post            = new TH2F("hPi0ppM_post",";M_{#gamma #gamma} [MeV];#theta_{#gamma #gamma} [rad]",200, 0, 250,200, 0, 4);
        // Pi0 Acceptance and Efficiency
        //Reminder Order: nbins, xmin, xmax 
        effPlotter              = new Pi0Efficiency(4, 1, 500);
        accPlotter              = new Pi0Acceptance("E", 4, 1, 550);
        pi0AcceptanceVsEta      = new Pi0Acceptance("eta", 100,-10, 10);
        pi0AcceptanceVsTheta    = new Pi0Acceptance("theta", 60, 0, TMath::Pi());    
    }
    if (doTruthAnalysis) {
        hPi0TrueMass            = new TH1F("hPi0TrueMass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
        //-->//--> Truth-Reco Mix Plots 
        h_mass_truthE_recoAngle = new TH1F("h_tE_rA",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
        h_mass_recoE_truthAngle = new TH1F("h_rE_tA",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
    }
    if (doChargedAnalysis) {
        //--> Pions
        hNSigmaPion             = new TH1F("hNSigmaPion", ";n#sigma;Counts", 100, -5, 5);
        // hdEdxVsE_cluster_Pion   = new TH2F("hdEdxVsE_cluster_Pion",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
        hdEdxVsE_cluster_Pion   = new TH2F("hdEdxVsE_cluster_Pion",";E [MeV];dE/dx [MeV/cm]",100, 0, 500,100, 0, 0.1);
        hdEdxVsE_true_Pion      = new TH2F("hdEdxVsE_true_Pion",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
        hdEdxTruePion           = new TH1F("hdEdxTruePion", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        hdEdxSmearPion          = new TH1F("hdEdxSmearPion", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        //--> Protons
        hNSigmaProton           = new TH1F("hNSigmaProton", ";n#sigma;Counts", 100, -5, 5);
        // hdEdxVsE_cluster_Proton = new TH2F("hdEdxVsE_cluster_Proton",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1);
        hdEdxVsE_cluster_Proton = new TH2F("hdEdxVsE_cluster_Proton",";E [MeV];dE/dx [MeV/cm]",100, 0, 500,100, 0, 0.1);
        hdEdxVsE_true_Proton    = new TH2F("hdEdxVsE_true_Proton",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.1); 
        hdEdxTrueProton         = new TH1F("hdEdxTrueProton", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        hdEdxSmearProton        = new TH1F("hdEdxSmearProton", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        //--> Electrons
        // hNSigmaElectron         = new TH1F("hNSigmaElectron", ";n#sigma;Counts", 100, -5, 5);
        // hdEdxVsE_cluster_Electron = new TH2F("hdEdxVsE_cluster_Electron",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.2);
        // hdEdxVsE_true_Electron    = new TH2F("hdEdxVsE_true_Electron",";E [MeV];dE/dx [MeV/cm]",200, 0, 500,200, 0, 0.2);   
        // hdEdxTrueElectron       = new TH1F("hdEdxTrueElectron", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        // hdEdxSmearElectron      = new TH1F("hdEdxSmearElectron", ";dE / dx [MeV/cm];Counts", 100, 0, 0.04);
        //--> Global Charged / PID related
        pidEff                  = new PIDEfficiency(20, 0, 500);
        h2_Eres                 = new TH2F("h2_Eres", ";True KE [MeV];Energy Residual", 20, 0, 500, 100, -1.0, 1.0);
        hClusterE               = new TH1F("hClusterE",";Cluster E [MeV];Count",100,0,500);
    }
    if (doEventVariables) {
        hEventInvariantMass     = new TH1F("hEventInvariantMass",";Invariant Mass [MeV];Events",100,0,5000);
        hEventSphericity        = new TH1F("hEventSphericity",";Sphericity;Events",100,0,1);
        //--> Total energy related
        hEvis                   = new TH1F("hEventTotE",";Total Energy[MeV];Events",100,0,3000);
        hEventCorrectedTotE     = new TH1F("hEventCorrectedTotE",";Total Corrected Energy[MeV];Events",100,0,3000);
        hErecoVsEtrue           = new TH2F("hErecoVsEtrue","; E_{true} [MeV];E_{reco} [MeV]",100, 0, 3000, 100, 0, 3000);
        hEvisVsEtrue            = new TH2F("hEvisVsEtrue","; E_{true} [MeV];E_{vis} [MeV]",100, 0, 3000, 100, 0, 3000);
        hDiffErecoEtrue         = new TH1F("hERecoVsETrue", "; E_{reco} - E_{true} [MEV]; Counts", 100,-1000, 1000);
        //--> Pion multiplicity related
        hNPionMultiplicity      = new TH2I("hNPionMultiplicity", ";True #pi^0 Multiplicity; Reconstructed #pi^0 Multiplicity", 10, 0,9,10,0,9);
        hChPionMultiplicity     = new TH2I("hChPionMultiplicity", ";True #pi^{#pm} Multiplicity; Reconstructed #pi^{#pm} Multiplicity", 10, 0,9,10,0,9);
    }

    ChargedKECalibration calibration("chargedKE.root");

    progressbar bar(nentries); 

    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        bar.update(); // Visual feedback on loop progress. 
        t->GetEntry(ievt);
        
        int nPi0 = 0;
        if (ievt < (Long64_t)pi0_per_event.size()) {
            nPi0 = pi0_per_event[ievt];
        }

        // Truth-level event vertex
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
        for (size_t k=0; k<nHits; ++k) {
            hits.push_back({(*centerXs)[k], (*centerYs)[k], (*centerZs)[k], (*energies)[k]});
        }

        if (doChargedAnalysis) {
            chargedTracks.clear();
            size_t nChargedTracks = TPC_Edep->size();
            for (size_t k=0; k<nChargedTracks; ++k) {
                if (!TPC_Edep || !TPC_firstPosX || !TPC_lastPosX) continue; 
                chargedTracks.push_back(
                    {k, 
                    vertex, 
                    TVector3((*TPC_lastPosX)[k], (*TPC_lastPosY)[k], (*TPC_lastPosZ)[k]), 
                    TVector3((*TPC_lastPosX)[k], (*TPC_lastPosY)[k], (*TPC_lastPosZ)[k]) - TVector3((*TPC_firstPosX)[k], (*TPC_firstPosY)[k], (*TPC_firstPosZ)[k]),
                    (*TPC_TrueKE)[k], (*TPC_pdg)[k], (*TPC_dEdx)[k], (*TPC_smearedEdep)[k], (*TPC_PathLength)[k], 0 /* Placeholder */, 0.15});
            }
        }

        if (doTruthAnalysis) {
            trueHits.clear();
            truePi0s.clear();
            size_t nTrueHits = truePhotonE->size();
            for (size_t k=0; k<nTrueHits; ++k) {
                trueHits.push_back({(*truePhotonPosX)[k], (*truePhotonPosY)[k], (*truePhotonPosZ)[k], (*truePhotonE)[k], (*truePhotonTrackID)[k], (*truePhotonParentID)[k]});
            }
            
            truePhotons = TruePhotonBuilder(trueHits, vertex);
            truePi0s = TruePi0Builder(truePhotons);

            for (TruePi0 tpi0 : truePi0s) {
                if (hPi0TrueMass) hPi0TrueMass->Fill(tpi0.p4.M());
            }
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Reconstruct Event

        RecoEvent reco = ReconstructEvent(hits, chargedTracks, vertex);

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



        // // Charged Object Clustering

        if (doChargedAnalysis) {

            for (const ChargedCluster& cluster : reco.chargedClusters) {

                // n-Sigma plots
                if (hNSigmaPion) hNSigmaPion->Fill(cluster.nSigmaPion);
                if (hNSigmaProton) hNSigmaProton->Fill(cluster.nSigmaProton);
                // if (hNSigmaElectron) hNSigmaElectron->Fill(cluster.nSigmaElectron);

                //--> Pions
                if (cluster.objectTruePDG == 211 || cluster.objectTruePDG == -211) {
                if (hdEdxVsE_cluster_Pion) hdEdxVsE_cluster_Pion->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y 
                // if (hdEdxVsE_cluster_Pion) hdEdxVsE_cluster_Pion->Fill(cluster.objectTrueKE, cluster.clusterdEdx); // ORDER: X vs Y 
                if (hClusterE) hClusterE->Fill(cluster.totalEnergy);
                if (hdEdxVsE_true_Pion) hdEdxVsE_true_Pion->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                if (hdEdxTruePion) hdEdxTruePion->Fill(cluster.objectTruedEdx);
                if (hdEdxSmearPion) hdEdxSmearPion->Fill(cluster.clusterdEdx);
                }
                //--> Protons
                if (cluster.objectTruePDG == 2212) {
                if (hdEdxVsE_cluster_Proton) hdEdxVsE_cluster_Proton->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
                // if (hdEdxVsE_cluster_Proton) hdEdxVsE_cluster_Proton->Fill(cluster.objectTrueKE, cluster.clusterdEdx); // ORDER: X vs Y 
                if (hdEdxVsE_true_Proton) hdEdxVsE_true_Proton->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                if (hdEdxTrueProton) hdEdxTrueProton->Fill(cluster.objectTruedEdx);
                if (hdEdxSmearProton) hdEdxSmearProton->Fill(cluster.clusterdEdx);            
                }
                //--> Electrons
                // if (cluster.objectTruePDG == 11 || cluster.objectTruePDG == -11) {
                // if (hdEdxVsE_cluster_Electron) hdEdxVsE_cluster_Electron->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y 
                // if (hClusterE) hClusterE->Fill(cluster.totalEnergy);
                // if (hdEdxVsE_cluster_Pion) hdEdxVsE_cluster_Pion->Fill(cluster.objectTrueKE, cluster.clusterdEdx); // ORDER: X vs Y 
                // if (hdEdxVsE_true_Electron) hdEdxVsE_true_Electron->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                // if (hdEdxTrueElectron) hdEdxTrueElectron->Fill(cluster.objectTruedEdx);
                // if (hdEdxSmearElectron) hdEdxSmearElectron->Fill(cluster.clusterdEdx);
                // }

                // E res plots
                double Eres = (cluster.totalEnergy - cluster.objectTrueKE) / cluster.objectTrueKE;
                if (h2_Eres) h2_Eres->Fill(cluster.objectTrueKE, Eres);
            }            
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Invariant mass plot (neutral pion / double photon)

        if (doPi0Analysis) {
            selected.clear();

            double param_h = 135;
            double param_k = 1.6;
            double param_a = 45;
            double param_b = 1.55;
            double relipse = 0;

            std::vector<Pi0Candidate> candidates;

            for (size_t a = 0; a < reco.clusters.size(); ++a) {
                for (size_t b = a + 1; b < reco.clusters.size(); ++b) {

                    const auto& g1 = reco.clusters[a].p4;
                    const auto& g2 = reco.clusters[b].p4;

                    double E1 = g1.E();
                    double E2 = g2.E();

                    double mgg = (g1 + g2).M();
                    double theta = openingAngle(g1, g2);

                    if (hPi0ppM_pre) hPi0ppM_pre->Fill(mgg, theta);

                    relipse = pow((mgg-param_h),2)/pow(param_a,2) + pow( (theta - param_k), 2)/pow(param_b,2);
                    if (relipse > 1) continue;                    

                    Pi0Candidate cand;
                    // cand.i = a;
                    // cand.j = b;
                    cand.c1 = &reco.clusters[a];
                    cand.c2 = &reco.clusters[b];
                    cand.mgg = mgg;
                    cand.theta = theta;
                    cand.p4 = g1 + g2;
                    cand.score = relipse;

                    candidates.push_back(cand);
                }
            }
            std::sort(candidates.begin(), candidates.end(), [](const Pi0Candidate& a, const Pi0Candidate& b) {
                return a.score < b.score;
            });
                
            std::unordered_set<const Cluster*> used;
            for (const auto& cand : candidates) {
                if (used.count(cand.c1) || used.count(cand.c2))
                    continue;

                used.insert(cand.c1);
                used.insert(cand.c2);
                selected.push_back(cand);
            }

            for (const auto& pi0 : selected) {
                if (hPi0Mass && hPi0ppM_post) {
                    hPi0Mass->Fill(pi0.p4.M());
                    hPi0ppM_post->Fill(pi0.mgg, pi0.theta);
                }
            }

            reco.nPionMultiplicity = selected.size();

        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Truth level + reco level plots 

        if (doTruthAndMixPlots) {
            std::vector<int> clusterToTrue = matchClustersToTruth(reco.clusters, truePhotons, 10);

            // Construct hybrid TLorentzVectors
            std::vector<TLorentzVector> photons_tE_rA; // truth E + reco angle
            std::vector<TLorentzVector> photons_rE_tA; // reco E + truth angle

            for (size_t ic=0; ic<reco.clusters.size(); ++ic) {
                int it = clusterToTrue[ic];
                if (it < 0) continue; // skip unmatched

                const Cluster &c = reco.clusters[ic];
                const TruePhoton &t = truePhotons[it];

                // truth energy + reco angle
                photons_tE_rA.push_back(makePhotonFromEnergyAndDir(t.p4.E(), c.centroid));

                // reco energy + true angle
                photons_rE_tA.push_back(makePhotonFromEnergyAndDir(c.p4.E(), t.dir));
            }

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
            // std::cout << "[PI0 Eff] truePi0s: " << truePi0s.size() << " selected: " << selected.size() << std::endl;
            if (effPlotter) effPlotter->ProcessEvent(truePi0s, selected, reco.clusters);


            // if (accPlotter) accPlotter->ProcessEvent(clusters, truePhotons, genEkin);
            if (accPlotter) accPlotter->ProcessSignalEvent(truePi0s, genEkin, nPi0);


            // double px = primaryPx->at(0);
            // double py = primaryPy->at(0);
            // double pz = primaryPz->at(0);
            // double p = sqrt(px*px + py*py + pz*pz);

            // double eta = 0.5 * log((p + pz) / (p - pz)); // pseudorapidity
            // double theta = acos(pz / p); // theta 

            // pi0AcceptanceVsEta.ProcessEvent(clusters, truePhotons, eta);
            // pi0AcceptanceVsTheta.ProcessEvent(clusters, truePhotons, theta);

            // pi0AcceptanceVsTheta.ProcessEventTwoHist(clusters, truePhotons, genEkin, theta);

        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        if (doChargedAnalysis) {
            pidEff->ProcessEvent(reco.chargedClusters);
            
            // assign a charged pion multiplicity based on PID Guess --> TPC info
            for (ChargedCluster ch : reco.chargedClusters) {
                if (ch.pidGuess == PID::Pion) {
                    reco.chPionMultiplicity += 1;
                }
            }

        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        if (doEventVariables) {
            
            double eVis = reco.EM_energy; // this is the sum of all the hit energies in the calorimeter, from any source.
            double eTrue = 0.0;
            size_t n = primaryEkin->size();
            for (size_t i = 0; i < n; ++i) {
                eTrue += (*primaryEkin)[i];
            }            
            double calibratedKE = calibration.GetMeanKE(reco.chargedClusters.size(), eVis);
            double eReco = eVis + calibratedKE;
            
            if (hEvisVsEtrue) hEvisVsEtrue->Fill(eTrue, eVis);
            if (hErecoVsEtrue) hErecoVsEtrue->Fill(eTrue, eReco);
            if (hDiffErecoEtrue) hDiffErecoEtrue->Fill(eReco - eTrue);
            if (hEvis) hEvis->Fill(eVis);
            if (hEventCorrectedTotE) hEventCorrectedTotE->Fill(eReco);

            if (hNPionMultiplicity) hNPionMultiplicity->Fill(pi0_per_event[ievt], reco.nPionMultiplicity); 
            if (hChPionMultiplicity) hChPionMultiplicity->Fill(chPi_per_event[ievt], reco.chPionMultiplicity); 

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
            //OLD STUFF WHICH PROBS IS WRONG

            // momenta.clear();
            // weights.clear();

            // TLorentzVector p4_event;
            // p4_event.SetPxPyPzE(0,0,0,0);

            // for (const auto& cl : reco.clusters) {
            //     p4_event += cl.p4;
            //     TVector3 p = cl.p4.Vect();
            //     momenta.push_back(p);
            //     weights.push_back(p.Mag());
            // }

            // for (const auto& ch : reco.chargedObjects) {
            //     p4_event += ch.p4;
            //     TVector3 p = ch.p4.Vect();
            //     momenta.push_back(p);
            //     weights.push_back(p.Mag());
            // }

            // if (hEventInvariantMass) hEventInvariantMass->Fill(p4_event.M());

            // TMatrixDSym S(3);
            // double norm = 0.0;

            // for (size_t i = 0; i < momenta.size(); ++i) {
            //     const TVector3& p = momenta[i];
            //     double w = weights[i];

            //     S(0,0) += w * p.X() * p.X();
            //     S(0,1) += w * p.X() * p.Y();
            //     S(0,2) += w * p.X() * p.Z();
            //     S(1,1) += w * p.Y() * p.Y();
            //     S(1,2) += w * p.Y() * p.Z();
            //     S(2,2) += w * p.Z() * p.Z();
                
            //     norm += w * p.Mag2();
            // }

            // S(1,0) = S(0,1);
            // S(2,0) = S(0,2);
            // S(2,1) = S(1,2);

            // if (norm > 0) S *= (1/norm);

            // TVectorD eigenVals(3);
            // TMatrixD eigenVecs = S.EigenVectors(eigenVals);

            // std::vector<double> lambdas = {eigenVals[0], eigenVals[1], eigenVals[2]};
            // std::sort(lambdas.begin(), lambdas.end(), std::greater<>());

            // double sphericity = 1.5 * (lambdas[1] + lambdas[2]);

            // if (hEventSphericity) hEventSphericity->Fill(sphericity);


        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    }    

    if (doPi0Analysis) {
        // Pi0 analysis
        PlotOptions optsPi0InvMass;
        optsPi0InvMass.doFit = true;
        optsPi0InvMass.fitMin = 100;
        optsPi0InvMass.fitMax = 170;
        optsPi0InvMass.addInfoPave = true;
        optsPi0InvMass.infoLines = {"Signal Dataset", Form("%d Events", nentries)}; //.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        // PrettyPi0MassPlot(hPi0Mass, "Pi0Mass_Clustered.png", 100.0, 170.0);
        Plot1D({hPi0Mass}, {kBlack}, "Neutral/Pi0InvMass.png", optsPi0InvMass);
    
        PlotOptions optsPi0ppM;
        optsPi0ppM.overlayProfileX = false;
        optsPi0ppM.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        Plot2D(hPi0ppM_pre, "Neutral/Pi0ppM_pre.png", optsPi0ppM);
        Plot2D(hPi0ppM_post, "Neutral/Pi0ppM_post.png", optsPi0ppM);

        //--> Pi0 Acceptance and Efficiency 
        effPlotter->FinalizePlot("Neutral/Pi0_efficiency_vs_Ekin.png");
        accPlotter->FinalizePlot("Neutral/Pi0_acceptance_vs_Ekin.png");
        // pi0AcceptanceVsEta.FinalizePlot("plots/Pi0_acceptance_vs_Eta.png");
        // pi0AcceptanceVsTheta.FinalizePlot("plots/Pi0_acceptance_vs_Theta_50_500.png");

        //CLEANUP
        delete hPi0Mass; 
        delete hPi0ppM_pre;
        delete hPi0ppM_post;
        delete effPlotter;
        delete accPlotter;
        delete pi0AcceptanceVsEta; 
        delete pi0AcceptanceVsTheta; 
    }

    if (doTruthAnalysis) {
        //--> Truth level
        PlotOptions optsPi0InvMassTT;
        optsPi0InvMassTT.addInfoPave = true;
        optsPi0InvMassTT.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsPi0InvMassTT.legendEntries = {"Truth-level Invariant Mass"};
        if (hPi0TrueMass) Plot1D({hPi0TrueMass}, {kBlack}, "Truth/Pi0InvMassTT.png", optsPi0InvMassTT);

        delete hPi0TrueMass;
    }

    if (doTruthAndMixPlots) {
        //-->//--> Truth-Reco Mix Plots
        PlotOptions optsPi0InvMassRecoAngle;
        optsPi0InvMassRecoAngle.doFit = true;
        optsPi0InvMassRecoAngle.fitMin = 80;
        optsPi0InvMassRecoAngle.fitMax = 180;
        // optsPi0InvMassRecoAngle.addInfoPave = true;
        // optsPi0InvMassRecoAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        // PrettyPi0MassPlot(h_mass_truthE_recoAngle, "Pi0Mass_truthE_recoAngle.png", 100.0, 170.0);
        if (h_mass_truthE_recoAngle) Plot1D({h_mass_truthE_recoAngle}, {kBlack}, "Truth/Pi0InvMassRecoAngle.png", optsPi0InvMassRecoAngle);

        PlotOptions optsPi0InvMassTruthAngle;
        // optsPi0InvMassTruthAngle.addInfoPave = true;
        // optsPi0InvMassTruthAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        optsPi0InvMassTruthAngle.legendEntries = {"Truth-level Invariant Mass"};
        optsPi0InvMassTruthAngle.extraLegendLines = {Form("Mean value: %.4f", h_mass_recoE_truthAngle->GetXaxis()->GetBinCenter(h_mass_recoE_truthAngle->GetMaximumBin() + 1))};
        // TruthPi0MassPlot(h_mass_recoE_truthAngle, "Pi0Mass_recoE_truthAngle.png");
        if (h_mass_recoE_truthAngle) Plot1D({h_mass_recoE_truthAngle}, {kBlack}, "Truth/Pi0InvMassTruthAngle.png", optsPi0InvMassTruthAngle);

        delete h_mass_truthE_recoAngle; 
        delete h_mass_recoE_truthAngle; 
    }

    // PID Plots 


    if (doChargedAnalysis) {
        PlotOptions opts_nSigma;
        opts_nSigma.doFit = true;
        opts_nSigma.addLegend = true;
        opts_nSigma.legendEntries = {
        "#pi^{#pm} hypothesis",
        "p hypothesis"
        };
        std::vector<TH1*> plots1D = {hNSigmaPion, hNSigmaProton};
        std::vector<int> colors = {kRed+1, kBlue+1};
        // opts_nSigma.legendEntries = {
        // "#pi^{#pm} hypothesis",
        // "e^{#pm} hypothesis",
        // "p hypothesis"
        // };
        // std::vector<TH1*> plots1D = {hNSigmaPion, hNSigmaElectron, hNSigmaProton};
        // std::vector<int> colors = {kRed+1, kGreen+1, kBlue+1};
        Plot1D(plots1D, colors, "Charged/nSigmaPlots.png", opts_nSigma);

        PlotOptions opts_hdEdxVsE_true;
        opts_hdEdxVsE_true.drawOption = "SCAT";
        opts_hdEdxVsE_true.legendEntries = {"#pi^{#pm}","p"};
        opts_hdEdxVsE_true.legendDrawOpt = "P";
        Plot2DOverlay({hdEdxVsE_true_Pion, hdEdxVsE_true_Proton}, {kBlack, kRed},"Charged/dedx_vs_E_overlay_true.png",opts_hdEdxVsE_true);

        PlotOptions opts_hdEdxVsE_cluster;
        opts_hdEdxVsE_cluster.drawOption = "SCAT";
        opts_hdEdxVsE_cluster.legendEntries = {"#pi^{#pm}","p"};
        // opts_hdEdxVsE_cluster.legendEntries = {"#pi^{#pm}", "e^{#pm}","p"};
        opts_hdEdxVsE_cluster.legendDrawOpt = "P";
        // Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Electron, hdEdxVsE_cluster_Proton}, {kBlack, kGreen, kRed},"Charged/dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);
        Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton}, {kBlack, kRed},"Charged/dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);

        // PlotOptions opts_hdEdxVsE_true;
        // opts_hdEdxVsE_true.drawOption = "SCAT";
        // opts_hdEdxVsE_true.legendEntries = {"e^{#pm}"};
        // opts_hdEdxVsE_true.legendDrawOpt = "P";
        // Plot2DOverlay({hdEdxVsE_true_Electron}, {kGreen},"Charged/dedx_vs_E_overlay_true.png",opts_hdEdxVsE_true);

        // PlotOptions opts_hdEdxVsE_cluster;
        // opts_hdEdxVsE_cluster.drawOption = "SCAT";
        // opts_hdEdxVsE_cluster.legendEntries = {"e^{#pm}"};
        // opts_hdEdxVsE_cluster.legendDrawOpt = "P";
        // Plot2DOverlay({ hdEdxVsE_cluster_Electron}, {kBlack, kGreen, kRed},"Charged/dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);

        PlotOptions opts_h2_Eres;
        TProfile* pEres = h2_Eres->ProfileX(Form("pEres_%s", h2_Eres->GetName()));
        Plot1D({pEres}, {kBlack}, "Charged/energy_residual.png", opts_h2_Eres);

        PlotOptions opts_dEdxPlots;
        opts_dEdxPlots.addLegend = true;
        // opts_dEdxPlots.legendEntries = {"True #pi^{#pm}", "Smeared #pi^{#pm}", "True e^{#pm}", "Smeared e^{#pm}", "True p", "Smeared p"};
        // std::vector<TH1*> plots1D_dEdx = {hdEdxTruePion, hdEdxSmearPion, hdEdxTrueElectron, hdEdxSmearElectron, hdEdxTrueProton, hdEdxSmearProton};
        // std::vector<int> colors_dEdx = {
        //     kBlack,
        //     kBlack+1,
        //     kGreen,
        //     kGreen+1,
        //     kRed,
        //     kRed+1
        // };
        opts_dEdxPlots.legendEntries = {"True #pi^{#pm}", "Smeared #pi^{#pm}", "True p", "Smeared p"};
        std::vector<TH1*> plots1D_dEdx = {hdEdxTruePion, hdEdxSmearPion, hdEdxTrueProton, hdEdxSmearProton};
        // std::vector<int> colors_dEdx = {
        //     kBlack,
        //     kBlack+1,
        //     kRed,
        //     kRed+1
        // };
        std::vector<int> colors_dEdx = {
            kBlack,
            kBlue,
            kRed,
            kGreen
        };
        Plot1D(plots1D_dEdx, colors_dEdx, "Charged/dEdxPlots.png", opts_dEdxPlots);

        pidEff->FinalizePlot("plots/Charged/PIDEfficiency", 211);

        PlotOptions opts_hClusterE;
        opts_hClusterE.addLegend = true;
        opts_hClusterE.legendEntries = {"#pi"};
        Plot1D({hClusterE}, {kBlack}, "Charged/ClusterE_pion.png", opts_hClusterE);

        //CLEANUP
        delete hNSigmaPion;
        delete hNSigmaProton;
        delete hdEdxVsE_cluster_Pion;
        delete hdEdxVsE_true_Pion; 
        delete hdEdxVsE_cluster_Proton;
        delete hdEdxVsE_true_Proton;
        delete hClusterE;
        delete pidEff;
    }

    if (doEventVariables) {
        // PlotOptions opts_hEventInvariantMass;
        // opts_hEventInvariantMass.addLegend = true;
        // opts_hEventInvariantMass.legendEntries = {"Event invariant mass"};
        // opts_hEventInvariantMass.addInfoPave = true;
        // Plot1D({hEventInvariantMass}, {kBlack}, "EventVar/eventInvariantMass.png", opts_hEventInvariantMass);

        // PlotOptions opts_hEventSphericity;
        // opts_hEventSphericity.addLegend = true;
        // opts_hEventSphericity.legendEntries = {"Event sphericity"};
        // opts_hEventSphericity.addInfoPave = true;
        // Plot1D({hEventSphericity}, {kBlack}, "EventVar/eventSphericity.png", opts_hEventSphericity);
        
        //--> Total energy related
        PlotOptions opts_hEventTotE;
        opts_hEventTotE.addLegend = true;
        opts_hEventTotE.legendEntries = {"E_{reco}", "E_{vis}"};
        // opts_hEventTotE.addInfoPave = true;
        Plot1D({hEventCorrectedTotE, hEvis}, {kBlack, kRed}, "EventVar/EventTotE.png", opts_hEventTotE);

        PlotOptions opts_hEvisVsEtrue;
        opts_hEvisVsEtrue.addLegend = true;
        opts_hEvisVsEtrue.legendEntries = {"E_{true} vs E_{vis}"};
        opts_hEvisVsEtrue.legendDrawOpt = "p";
        // opts_hEvisVsEtrue.addInfoPave = true;
        opts_hEvisVsEtrue.overlayProfileX = true;
        opts_hEvisVsEtrue.profileColor = kRed;
        Plot2D(hEvisVsEtrue, "EventVar/EvisVsEtrue.png", opts_hEvisVsEtrue);
        
        PlotOptions opts_hErecoVsEtrue;
        opts_hErecoVsEtrue.addLegend = true;
        opts_hErecoVsEtrue.legendEntries = {"E_{true} vs E_{reco}"};
        opts_hErecoVsEtrue.legendDrawOpt = "p";
        // opts_hErecoVsEtrue.addInfoPave = true;
        opts_hErecoVsEtrue.overlayProfileX = true;
        opts_hErecoVsEtrue.profileColor = kRed;
        Plot2D(hErecoVsEtrue, "EventVar/ErecoVsEtrue.png", opts_hErecoVsEtrue);
        
        PlotOptions opts_hDiffErecoEtrue;
        opts_hDiffErecoEtrue.addLegend = true;
        opts_hDiffErecoEtrue.legendEntries = {"E_{reco} - E_{true}"};
        // opts_hDiffErecoEtrue.addInfoPave = true;
        Plot1D({hDiffErecoEtrue}, {kBlack}, "EventVar/DiffErecoEtrue.png", opts_hDiffErecoEtrue);

        //--> Pion multiplicity related
        PlotOptions opts_hNPionMultiplicity;
        opts_hNPionMultiplicity.drawOption = "COLZ";
        Plot2D(hNPionMultiplicity, "EventVar/NPionMultiplicity.png", opts_hNPionMultiplicity);

        PlotOptions opts_hChPionMultiplicity;
        opts_hChPionMultiplicity.drawOption = "COLZ";
        Plot2D(hChPionMultiplicity, "EventVar/ChPionMultiplicity.png", opts_hChPionMultiplicity);

        // delete hEventInvariantMass;
        // delete hEventSphericity;
        delete hEventCorrectedTotE;
        delete hErecoVsEtrue;
        delete hEvisVsEtrue;
        delete hDiffErecoEtrue;
        delete hChPionMultiplicity;
        delete hNPionMultiplicity;
    }

    f->Close(); delete f;
    return 0;
}