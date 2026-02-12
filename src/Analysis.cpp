#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1.h"
#include "TF1.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"

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
                    << "  --event-variables     Computes the event level variables\n";
            return 1;
        }
    
    std::string root_inputfile = argv[1];
    std::string mcpl_inputfile = argv[2];
    bool doCalibration = false;
    bool doPi0Analysis = false;
    bool doChargedAnalysis = false;
    bool doTruthAnalysis = false;
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
            doTruthAnalysis = true;
            doEventVariables = true;
        }
        if (arg == "--pi0-analysis") doPi0Analysis = true;
        if (arg == "--charged-analysis") doChargedAnalysis = true;
        if (arg == "--truth-analysis") doTruthAnalysis = true;
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
    std::map<int, double> Ekin_per_event;
    mcpl_file_t mcplFile = mcpl_open_file(mcpl_inputfile.c_str());
    int current_event = -1;
    const mcpl_particle_t* p;
    while ((p = mcpl_read(mcplFile))) {

        if (p->userflags == 1) {
            ++current_event;
            pi0_per_event.push_back(0);
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
    TH1F *hClusterE                      = nullptr;
    TH1F *hPi0TrueMass                   = nullptr;
    TH1F *h_mass_truthE_recoAngle        = nullptr;
    TH1F *h_mass_recoE_truthAngle        = nullptr;
    TH1F *hEffvsE                        = nullptr;
    // Charged analysis objects
    TH1F  *hNSigmaPion                   = nullptr;
    TH1F  *hNSigmaProton                 = nullptr;
    TH2F  *hdEdxVsE_cluster_Pion         = nullptr;
    TH2F  *hdEdxVsE_true_Pion            = nullptr;
    TH2F  *hdEdxVsE_cluster_Proton       = nullptr;
    TH2F  *hdEdxVsE_true_Proton          = nullptr;
    TH2F  *h2_Eres                       = nullptr;
    TH1F  *hdEdxTruePion                 = nullptr;
    TH1F  *hdEdxSmearPion                = nullptr;
    TH1F  *hdEdxTrueProton               = nullptr;
    TH1F  *hdEdxSmearProton              = nullptr;
    // Analysis helper classes (also pointers)
    Pi0Efficiency  *effPlotter           = nullptr;
    Pi0Acceptance  *accPlotter           = nullptr;
    Pi0Acceptance  *pi0AcceptanceVsEta   = nullptr;
    Pi0Acceptance  *pi0AcceptanceVsTheta = nullptr;
    PIDEfficiency  *pidEff               = nullptr;
    // Event level variables
    TH1F  *hEventInvariantMass           = nullptr; 
    TH1F  *hEventSphericity              = nullptr; 
    TH1F  *hEventCorrectedTotE           = nullptr; 

    TH2F  *hTrueVsVis                    = nullptr;
    TH2F  *hEventClosureTest             = nullptr;
    TH1F  *hDiffERecoVsETrue             = nullptr;

    // std::map<int, TH2F*> hClosureTestByNchs;
    // TH2F  *hEventClosureTest             = nullptr;
    // TH1F  *hDiffERecoVsETrue             = nullptr;
    // TH1F  *hDiffEVisVsETrue              = nullptr;
    // TH2F  *hDiffERecoETrueVsETRUE        = nullptr;


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
        hPi0ppM_pre             = new TH2F("hPi0ppM_pre",";M_{#gamma #gamma} [MeV];#theta (#gamma #gamma)",200, 0, 250,200, 0, 4);
        hPi0ppM_post            = new TH2F("hPi0ppM_post",";M_{#gamma #gamma} [MeV];#theta (#gamma #gamma)",200, 0, 250,200, 0, 4);
        hEffvsE                 = new TH1F("hEffvsE", ";#pi^0 E_{kin}; Efficiency", 100, 1, 500);
        effPlotter              = new Pi0Efficiency(4, 1, 500);
        accPlotter              = new Pi0Acceptance("E", 4, 1, 550);
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
    }
    if (doEventVariables) {
        hEventInvariantMass     = new TH1F("hEventInvariantMass",";Invariant Mass [MeV];Events",100,0,5000);
        hEventSphericity        = new TH1F("hEventSphericity",";Sphericity;Events",100,0,1);
        hEventCorrectedTotE     = new TH1F("hEventCorrectedTotE",";Total Corrected Energy[MeV];Events",100,0,5000);
        // hTrueVsVis              = new TH2F("hTrueVsVis", "; E_{vis} [MeV]; E_{true} [MeV]", 120, 0, 2000, 120, -1000, 1500);
        hTrueVsVis              = new TH2F("hTrueVsVis", "; E_{vis} [MeV]; E_{true} [MeV]", 100, 0, 2000, 100, 0, 3000);
        hEventClosureTest       = new TH2F("hEventClosureTest","; E_{true} [MeV];E_{vis} + linear correction [MeV]",100, 0, 3000, 100, 0, 3000);
        hDiffERecoVsETrue       = new TH1F("hERecoVsETrue", "; E_{reco} - E_{true} [MEV]; Counts", 100,-1000, 1000);


        // hEventClosureTest       = new TH2F("hEventClosureTest",";Truth primary EKin [MeV];E_{EM} + template correction [MeV]",100, 0, 3000,100, 0, 3000);
        // hEventClosureTest       = new TH2F("hEventClosureTest",";Truth primary EKin [MeV];E_{VIS} [MeV]",100, 0, 3000,100, 0, 3000);
        // hDiffERecoVsETrue       = new TH1F("hERecoVsETrue", "; E_{RECO} - E_{TRUE} [MEV]; Counts", 100,-1000, 1000);
        // hDiffEVisVsETrue        = new TH1F("hEVisVsETrue", "; E_{VIS} - E_{TRUE} [MEV]; Counts", 100,-1000, 1000);
        // hDiffERecoETrueVsETRUE  = new TH2F("hDiffERecoETrueVsETRUE",";Truth primary EKin [MeV];E_{RECO} - E_{TRUE} [MeV]",100, 0, 3000,100, -1000, 1000);
        // hDiffERecoETrueVsETRUE  = new TH2F("hDiffERecoETrueVsETRUE",";E_{VIS} [MeV];E_{TRUE} - E_{VIS}  [MeV]",100, 0, 3000,100, -1000, 1000);

        // hClosureTestByNchs[1]   = new TH2F("h2_Nch1", ";E_{EM} [MEV];True Ch KE [MeV]", 50,0,2000, 50,0,2000);
        // hClosureTestByNchs[2]   = new TH2F("h2_Nch2", ";E_{EM} [MEV];True Ch KE [MeV]", 50,0,2000, 50,0,2000);
        // hClosureTestByNchs[3]   = new TH2F("h2_Nch3", ";E_{EM} [MEV];True Ch KE [MeV]", 50,0,2000, 50,0,2000);
        // hClosureTestByNchs[4]   = new TH2F("h2_Nch4p", ";E_{EM} [MEV];True Ch KE [MeV]", 50,0,2000, 50,0,2000);
    }

    ChargedKECalibration calibration("chargedKE.root");

    progressbar bar(nentries); 

    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        bar.update(); // just to get some visual feedback on loop progress. 
        t->GetEntry(ievt);
        
        int nPi0 = 0;
        if (ievt < (Long64_t)pi0_per_event.size()) {
            nPi0 = pi0_per_event[ievt];
        }

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
            chargedTracks.clear();
            size_t nChargedTracks = TPC_Edep->size();
            // std::cout << "[CHARGED] Number of charged tracks: " << nChargedTracks << std::endl;
            for (size_t k=0; k<nChargedTracks; ++k) {
                if (!TPC_Edep || !TPC_firstPosX || !TPC_lastPosX) continue; // safety
                // for now the resolution is hardcoded to be 0.15
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
            // truePhotons.clear();
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
            // double thetaMax = 25.0 * TMath::DegToRad();
            // // chargedClusters.clear();
            // chargedObjects.clear();
            // chargedClusters = MatchHitsToTracks(chargedTracks, hits, thetaMax);

            for (const ChargedCluster& cluster : reco.chargedClusters) {
                
                // double calibratedKE = calibration.GetMeanKE(reco.chargedClusters.size(), reco.EM_energy);
                // double mass = PDGMassMeV(PIDToPDG(cluster.pidGuess));
                // double p_mag = std::sqrt(calibratedKE * (calibratedKE + 2 * mass));
                // TVector3 dir = cluster.direction.Unit();
                // //Build charged objects
                // TLorentzVector charged_p4;
                // charged_p4.SetPxPyPzE(
                //     dir.X() * p_mag,
                //     dir.Y() * p_mag,
                //     dir.Z() * p_mag,
                //     calibratedKE + mass
                // );

                // reco.chargedObjects.push_back({cluster.trackID, charged_p4, &cluster});




                if (hNSigmaPion) hNSigmaPion->Fill(cluster.nSigmaPion);
                if (hNSigmaProton) hNSigmaProton->Fill(cluster.nSigmaProton);
                if (cluster.objectTruePDG == 211 || cluster.objectTruePDG == -211) {
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
                double Eres = (cluster.totalEnergy - cluster.objectTrueKE) / cluster.objectTrueKE;
                if (h2_Eres) h2_Eres->Fill(cluster.objectTrueKE, Eres);
            }            
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Neutral Object Clustering

        // if (doPi0Analysis) {
        //     // double totalE_Evt = 0;
        //     clusters.clear();

        //     // double dEta = 0.10/2;
        //     // double dPhi = 0.10/2;
        //     // double E_seed = 15.00;
        //     // double E_neighbor = 0.03;
        //     // int winSize = 7;
        //     // clusters = SlidingWindowClusterHits(hits, vertex, dEta, dPhi, E_seed, E_neighbor, winSize);
            
        //     clusters = clusterNeutralHits(hits, vertex, 25 * TMath::DegToRad()); // Have to optimise the angle a bit! 
        //     // Apply cluster energy threshold
        //     clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
        //                                 [](const Cluster &c){ return c.p4.E() < 50.0; }),
        //                 clusters.end());
        // }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
       
        // // Invariant mass plot (neutral pion / double photon)

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
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Truth level + reco level plots 

        // should introduce the selection procedure here too! Right now this is WRONG

        if (doTruthAnalysis) {
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
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        if (doEventVariables) {
            
            // Corrected total energy
            // double totCorrectedE = 0.0; // confusing name and also possibly somethig wrong here!
            // double totTrueKE = 0.0;
            double totTpcDeposit = 0.0;

            // for (const auto& cl : reco.clusters) {
            //     totCorrectedE += cl.p4.Energy();
            // }
            for (const auto& ch : reco.chargedClusters) {
                // totCorrectedE += ch.totalEnergy;
                // totTrueKE += ch.objectTrueKE;
                totTpcDeposit += ch.EdepSmeared;
            }
            double eVis = reco.EM_energy + totTpcDeposit;
            // double calibratedKE = calibration.GetMeanKE(reco.chargedClusters.size(), eVis);
            // double eReco = eVis + calibratedKE;

            double eTrue = 0.0;
            size_t n = std::min(primaryEkin->size(), primaryPDG->size());
            // // for (size_t i = 0; i < n; ++i) {
            // //     double px = (*primaryPx)[i];
            // //     double py = (*primaryPy)[i];
            // //     double pz = (*primaryPz)[i];
            // //     double p  = std::sqrt(px*px + py*py + pz*pz);
            // //     double m  = PDGMassMeV((*primaryPDG)[i]);
            // //     double E  = std::sqrt(p*p + m*m);
            // //     eTrue += E;
            // // }
            for (size_t i = 0; i < n; ++i) {
                eTrue += (*primaryEkin)[i];
            }

            // if (hTrueVsVis) hTrueVsVis->Fill(eVis, eTrue - eVis); // bad name and axes titles! 
            if (hTrueVsVis) hTrueVsVis->Fill(eVis, eTrue);  
            // TProfile* p = hTrueVsVis->ProfileX("pTrueVsVis");
            // TF1* f = new TF1("f_lin", "[0] + [1]*x", 0, 2000);

            // p->Fit(f, "Q");
            // double A = f->GetParameter(0);
            // double B = f->GetParameter(1);
            // double Aerr = f->GetParError(0);
            // double Berr = f->GetParError(1);
            // std::cout << "Params A, B : " << A << ", " << B << std::endl;

            // Apply the fit 
            double A = 1346.54;
            double B = -1.06458;
            double correction = A + B * eVis;
            double eReco = eVis + correction;

            // // double eDiff = eReco - Ekin_per_event[ievt];
            // // if (hEventClosureTest) hEventClosureTest->Fill(Ekin_per_event[ievt], eReco);
            double eDiffReco = eReco - eTrue;
            // double eDiffVis = eVis - eTrue;
            if (hEventClosureTest) hEventClosureTest->Fill(eTrue, eReco);
            if (hDiffERecoVsETrue) hDiffERecoVsETrue->Fill(eDiffReco);
            // if (hDiffEVisVsETrue) hDiffEVisVsETrue->Fill(eDiffVis);

            // // if (hDiffERecoETrueVsETRUE) hDiffERecoETrueVsETRUE->Fill(eTrue, eReco - eTrue);
            // if (hDiffERecoETrueVsETRUE) hDiffERecoETrueVsETRUE->Fill(eVis, eTrue - eVis);



            // if (ievt < 10) {

            //     std::cout << "Event " << ievt << "\n";
            //     std::cout << "Truth mcpl: " << Ekin_per_event[ievt] << "\n";
            //     std::cout << "Truth primary: " << Etrue << "\n";
            //     std::cout << "EM:    " << reco.EM_energy << "\n";
            //     std::cout << "TPC:   " << totTpcDeposit << "\n";
            //     std::cout << "Template KE: " << calibratedKE << "\n";
            //     std::cout << "Reco:  " << eReco << "\n";
            //     std::cout << "Diff:  " << eDiff << "\n\n";
            // }

            // if (hEventCorrectedTotE) hEventCorrectedTotE->Fill(totCorrectedE + calibratedKE);

            // if (hEventClosureTest) hEventClosureTest->Fill(reco.EM_energy, totTrueKE);
            // if (hEventClosureTest) hEventClosureTest->Fill(Ekin_per_event[ievt], reco.EM_energy + calibratedKE);

            // int Nch = reco.chargedClusters.size();
            // int key = (Nch >= 4) ? 4 : Nch;
            // if (hClosureTestByNchs.count(key)) hClosureTestByNchs[key]->Fill(reco.EM_energy, totTrueKE);

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

            // TProfile* profile = hTrueVsVis->ProfileX("pTrueVsVis");
            // TF1* fitfunc = new TF1("f_lin", "[0] + [1]*x", 200, 1600);

            // profile->Fit(fitfunc, "Q");
            // double A = fitfunc->GetParameter(0);
            // double B = fitfunc->GetParameter(1);
            // double Aerr = fitfunc->GetParError(0);
            // double Berr = fitfunc->GetParError(1);
            // std::cout << "Params A, B : " << A << ", " << B << std::endl;
    

    if (doPi0Analysis) {
        // RECO-RECO
        PlotOptions optsPi0InvMass;
        optsPi0InvMass.doFit = true;
        optsPi0InvMass.fitMin = 100;
        optsPi0InvMass.fitMax = 170;
        optsPi0InvMass.addInfoPave = true;
        optsPi0InvMass.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        // PrettyPi0MassPlot(hPi0Mass, "Pi0Mass_Clustered.png", 100.0, 170.0);
        Plot1D({hPi0Mass}, {kBlack}, "Neutral/Pi0InvMass.png", optsPi0InvMass);

        PlotOptions optsPi0ppM;
        Plot2D(hPi0ppM_pre, "Neutral/Pi0ppM_pre.png", optsPi0ppM);
        Plot2D(hPi0ppM_post, "Neutral/Pi0ppM_post.png", optsPi0ppM);


        // Efficiency Plot
        effPlotter->FinalizePlot("Neutral/Pi0_efficiency_vs_Ekin.png");

        // // Acceptance Plot
        accPlotter->FinalizePlot("Neutral/Pi0_acceptance_vs_Ekin.png");
        // pi0AcceptanceVsEta.FinalizePlot("plots/Pi0_acceptance_vs_Eta.png");
        // pi0AcceptanceVsTheta.FinalizePlot("plots/Pi0_acceptance_vs_Theta_50_500.png");

        //CLEANUP
        delete hPi0Mass; 
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
        if (hPi0TrueMass) Plot1D({hPi0TrueMass}, {kBlack}, "Truth/Pi0InvMassTT.png", optsPi0InvMassTT);

        // Mix Plots
        PlotOptions optsPi0InvMassRecoAngle;
        optsPi0InvMassRecoAngle.doFit = true;
        optsPi0InvMassRecoAngle.fitMin = 80;
        optsPi0InvMassRecoAngle.fitMax = 180;
        optsPi0InvMassRecoAngle.addInfoPave = true;
        optsPi0InvMassRecoAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        // PrettyPi0MassPlot(h_mass_truthE_recoAngle, "Pi0Mass_truthE_recoAngle.png", 100.0, 170.0);
        if (h_mass_truthE_recoAngle) Plot1D({h_mass_truthE_recoAngle}, {kBlack}, "Truth/Pi0InvMassRecoAngle.png", optsPi0InvMassRecoAngle);

        PlotOptions optsPi0InvMassTruthAngle;
        optsPi0InvMassTruthAngle.addInfoPave = true;
        optsPi0InvMassTruthAngle.infoLines = {"GEANT4 pi0 sample", "5000 events", "E_{kin} #in", "[50, 150, 300,", "450, 500] MeV"};
        optsPi0InvMassTruthAngle.legendEntries = {"Truth-level Invariant Mass"};
        optsPi0InvMassTruthAngle.extraLegendLines = {Form("Mean value: %.4f", h_mass_recoE_truthAngle->GetXaxis()->GetBinCenter(h_mass_recoE_truthAngle->GetMaximumBin() + 1))};
        // TruthPi0MassPlot(h_mass_recoE_truthAngle, "Pi0Mass_recoE_truthAngle.png");
        if (h_mass_recoE_truthAngle) Plot1D({h_mass_recoE_truthAngle}, {kBlack}, "Truth/Pi0InvMassTruthAngle.png", optsPi0InvMassTruthAngle);

        delete hPi0TrueMass;
        delete h_mass_truthE_recoAngle; 
        delete h_mass_recoE_truthAngle; 
    }

    // PID Plots 


    if (doChargedAnalysis) {
        PlotOptions opts_nSigma;
        opts_nSigma.doFit = true;
        opts_nSigma.addLegend = true;
        opts_nSigma.legendEntries = {
        "#pi hypothesis",
        "p hypothesis"
        };
        std::vector<TH1*> plots1D = {hNSigmaPion, hNSigmaProton};
        std::vector<int> colors = {
            kRed+1,
            kBlue+1
        };
        Plot1D(plots1D, colors, "Charged/nSigmaPlots.png", opts_nSigma);

        PlotOptions opts_hdEdxVsE_true;
        opts_hdEdxVsE_true.drawOption = "SCAT";
        opts_hdEdxVsE_true.legendEntries = {"#pi^{+}","p"};
        opts_hdEdxVsE_true.legendDrawOpt = "P";
        Plot2DOverlay({hdEdxVsE_true_Pion, hdEdxVsE_true_Proton}, {kBlack, kRed},"Charged/dedx_vs_E_overlay_true.png",opts_hdEdxVsE_true);

        PlotOptions opts_hdEdxVsE_cluster;
        opts_hdEdxVsE_cluster.drawOption = "SCAT";
        opts_hdEdxVsE_cluster.legendEntries = {"#pi^{+}","p"};
        opts_hdEdxVsE_cluster.legendDrawOpt = "P";
        Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton}, {kBlack, kRed},"Charged/dedx_vs_E_overlay_cluster.png",opts_hdEdxVsE_cluster);

        PlotOptions opts_h2_Eres;
        TProfile* pEres = h2_Eres->ProfileX();
        Plot1D({pEres}, {kBlack}, "Charged/energy_residual.png", opts_h2_Eres);

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
        Plot1D(plots1D_dEdx, colors_dEdx, "Charged/dEdxPlots.png", opts_dEdxPlots);

        //Pion + 
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

        // PlotOptions opts_hEventCorrectedTotE;
        // opts_hEventCorrectedTotE.addLegend = true;
        // opts_hEventCorrectedTotE.legendEntries = {"Total corrected event energy"};
        // opts_hEventCorrectedTotE.addInfoPave = true;
        // Plot1D({hEventCorrectedTotE}, {kBlack}, "EventVar/EventCorrectedTotE.png", opts_hEventCorrectedTotE);

        PlotOptions opts_hEventClosureTest_full;
        opts_hEventClosureTest_full.drawOption = "COLZ";
        opts_hEventClosureTest_full.overlayProfileX = true;
        opts_hEventClosureTest_full.profileColor = kRed;
        opts_hEventClosureTest_full.addTopLatex = true;
        // opts_hEventClosureTest.addLegend = true;
        // opts_hEventClosureTest.legendEntries = {"E_{EM} vs KE^{true}"};
        Plot2D({hEventClosureTest}, "EventVar/EventClosureTest.png", opts_hEventClosureTest_full);

        // for (auto const& [nCh, hist] : hClosureTestByNchs) {
        //     if (!hist) continue;
        //     PlotOptions opts_hEventClosureTest;
        //     opts_hEventClosureTest.drawOption = "COLZ";
        //     opts_hEventClosureTest.overlayProfileX = true;
        //     opts_hEventClosureTest.profileColor = kRed;
        //     opts_hEventClosureTest.addTopLatex = true;
        //     std::string nchlabel;
        //     if (nCh < 4) nchlabel = Form("N_{ch} = %d", nCh);
        //     else nchlabel = "N_{ch} #geq 4";
        //     opts_hEventClosureTest.addLegend = true;
        //     opts_hEventClosureTest.legendEntries = {nchlabel};
        //     std::string outname = Form("EventVar/ClosureTest/EventClosureTest_%d.png", nCh);
        //     Plot2D({hist}, outname, opts_hEventClosureTest); 
        // }

        PlotOptions opts_hDiffERecoVisVsETrue;
        opts_hDiffERecoVisVsETrue.addLegend = true;
        // opts_hDiffERecoVisVsETrue.legendEntries = {"E_{RECO} - E_{TRUE}", "E_{VIS} - E_{TRUE}"};
        opts_hDiffERecoVisVsETrue.legendEntries = {"E_{RECO} - E_{TRUE}"};
        // opts_hDiffERecoVisVsETrue.extraLegendLines = {Form("1) MEAN: %.2f, RMS: %.2f", hDiffERecoVsETrue->GetMean(), hDiffERecoVsETrue->GetRMS()),
        //                                               Form("2) MEAN: %.2f, RMS: %.2f", hDiffEVisVsETrue->GetMean(), hDiffEVisVsETrue->GetRMS())};
        opts_hDiffERecoVisVsETrue.extraLegendLines = {Form("1) MEAN: %.2f, RMS: %.2f", hDiffERecoVsETrue->GetMean(), hDiffERecoVsETrue->GetRMS())};
        opts_hDiffERecoVisVsETrue.addInfoPave = true;
        opts_hDiffERecoVisVsETrue.legendX1 = 0.65;
        opts_hDiffERecoVisVsETrue.legendX2 = 0.98;
        // Plot1D({hDiffERecoVsETrue, hDiffEVisVsETrue}, {kBlack, kRed}, "EventVar/ERecoVisVsETrue.png", opts_hDiffERecoVisVsETrue);
        Plot1D({hDiffERecoVsETrue}, {kBlack}, "EventVar/DiffERecoETrue.png", opts_hDiffERecoVisVsETrue);

        PlotOptions opts_hTrueVsVis;
        // opts_hTrueVsVis.addLegend = true;
        // opts_hTrueVsVis.legendEntries = {"E_{vis}"}
        opts_hTrueVsVis.addInfoPave = true;
        opts_hTrueVsVis.overlayProfileX = true;
        opts_hTrueVsVis.overlayFitLine = true;
        Plot2D(hTrueVsVis, "EventVar/ETrueVsEVis.png", opts_hTrueVsVis);

        // PlotOptions opts_hDiffERecoETrueVsETRUE;
        // opts_hDiffERecoETrueVsETRUE.drawOption = "COLZ";
        // opts_hDiffERecoETrueVsETRUE.overlayProfileX = true;
        // opts_hDiffERecoETrueVsETRUE.profileColor = kRed;
        // opts_hDiffERecoETrueVsETRUE.addTopLatex = true;
        // // opts_hEventClosureTest.addLegend = true;
        // // opts_hEventClosureTest.legendEntries = {"E_{EM} vs KE^{true}"};
        // Plot2D({hDiffERecoETrueVsETRUE}, "EventVar/DiffERecoETrueVsETRUE.png", opts_hDiffERecoETrueVsETRUE);

        // delete hEventInvariantMass;
        // delete hEventSphericity;
        delete hEventCorrectedTotE;
        // for (auto const& [nCh, hist] : hClosureTestByNchs) delete hClosureTestByNchs[nCh];
        delete hEventClosureTest;
        delete hDiffERecoVsETrue;
        delete hTrueVsVis;
        // delete hDiffEVisVsETrue;
        // delete hDiffERecoETrueVsETRUE;
    }

    f->Close(); delete f;
    return 0;
}