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

#include <algorithm>

#include "CLI11.hpp"
#include "progressbar.hpp"

#include "BranchManager.hpp"
#include "AnalysisHistograms.hpp"
#include "Structures.hpp"
#include "Clustering.hpp"
#include "EventLoop.hpp"
#include "PlotUtils.hpp"
#include "EventDisplay.hpp"
#include "Utils.hpp"
#include "TruePhotonCalc.hpp"
#include "PhotonMatch.hpp"
#include "Pi0Efficiency.hpp"
#include "Acceptance.hpp"
#include "PIDEfficiency.hpp"
#include "RecoEvent.hpp"
#include "Calibration.hpp"
#include "EventVariables.hpp"
#include "SelectionCuts.hpp"

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

int main(int argc, char **argv) {
    if (argc < 4) {
            std::cout << "Usage: " << argv[0] << " <input.root> <background.root> <vertices.root> [options]\n"
                    << "Options:\n"
                    << "  --calibration         Performs calibration procedure\n"
                    << "  --full-analysis       Performs everything\n"
                    << "  --pi0-analysis        Reconstruct Pi0s and invariant mass analysis\n"
                    << "  --charged-analysis    Reconstruct Charged objects and do PID studies\n"
                    << "  --truth-analysis      Performs truth level analysis for Pi0s\n"
                    << "  --truth-mix-analysis  Performs truth-reco level mixed analysis for Pi0s\n"
                    << "  --event-variables     Computes the event level variables\n"
                    << "  --do-selection        Performs Selection\n";
            return 1;
        }
    SetErrorHandler(MyErrorHandler);

    std::string root_inputfile = argv[1];
    std::string background_inputfile = argv[2];
    std::string vertices_inputfile = argv[3];
    bool doCalibration = false;
    bool doSelection = false;

    AnalysisConfig cfg;

    for (int i = 4; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--calibration") doCalibration = true;
        if (arg == "--full-analysis") {
            if (!std::filesystem::exists("chargedKE.root")) {
                std::cerr << "No calibration root file detected, please run './Analysis <input.root> <background.root> <vertices.root> --calibration' first" << std::endl;
                return 1;
            }
            cfg.doPi0Analysis = true;
            cfg.doChargedAnalysis = true;
            cfg.doTruthAnalysis = true; // NOTE for the truth-reco mix plots you also need to have the flag "--truth-mix-analysis"
            cfg.doEventVariables = true;
        }
        if (arg == "--pi0-analysis") cfg.doPi0Analysis = true;
        if (arg == "--charged-analysis") cfg.doChargedAnalysis = true;
        if (arg == "--truth-analysis") cfg.doTruthAnalysis = true;
        if (arg == "--truth-mix-analysis") {
            cfg.doTruthAndMix = true;
            cfg.doTruthAnalysis = true;
        }
            if (arg == "--event-variables") cfg.doEventVariables = true;
            if (arg == "--do-selection") doSelection = true;

    }

    SetPrettyStyle();

    //Open Signal input file and branch setup
    TFile *f = TFile::Open(root_inputfile.c_str());
    if (!f || f->IsZombie()) return 1;
    TTree *t = (TTree *)f->Get("digitizedHits");
    if (!t) { std::cerr << "Tree digitizedHits not found" << std::endl; return 1;}
    BranchManagerInput br;
    br.SetBranches(t);

    // Histogram booking 
    Pi0Histograms hPi0;
    TruthHistograms hTruth;
    ChargedHistograms hCharged;
    EventVarHistograms hEvt;

    if (cfg.doPi0Analysis) hPi0.Book();
    if (cfg.doTruthAnalysis) hTruth.Book();
    if (cfg.doChargedAnalysis) hCharged.Book();
    if (cfg.doEventVariables) hEvt.Book();

    Long64_t nentries = t->GetEntries();

    //Open Vertex input file and branch setup
    TFile *vtxFile = TFile::Open(vertices_inputfile.c_str());
    if (!vtxFile || vtxFile->IsZombie()) return 1;
    TTree *vtxTree = (TTree *)vtxFile->Get("vertices");
    if (!vtxTree) {std::cerr << "Tree vertices not found" << std::endl; return 1;}
    BranchManagerVertex brVtx;   
    brVtx.SetBranches(vtxTree);
    brVtx.LoadVertices(vtxTree, nentries);

    vtxFile->Close(); delete vtxFile;

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    // Calibration & Theory dEdx tables

    DEDXTable dedxTABLE("dedx_tables_Ar80CO2.root");

    //Calibration procedure ----->>>> HAVE TO DOUBLE CHECK THIS LOGIC. A lot has changed SINCE THEN!!!!!
    // 2026-03-23: removed the unused pi0_per event. The only other change is vertex related, but since this is jsut an energy counting thing does it even matter?
    if (doCalibration) {
        DoCalibration(
            t,
            static_cast<int>(nentries * 0.7),
            // pi0_per_event,
            br.centerXs, br.centerYs, br.centerZs, br.energies,
            br.primaryX, br.primaryY, br.primaryZ, br.primaryEkin,
            br.TPC_firstPosX, br.TPC_firstPosY, br.TPC_firstPosZ,
            br.TPC_lastPosX,  br.TPC_lastPosY,  br.TPC_lastPosZ,
            br.TPC_TrueKE, br.TPC_pdg, br.TPC_nSteps, br.TPC_dEdx,
            br.TPC_smearedEdep, br.TPC_PathLength,
            dedxTABLE,
            "chargedKE.root"
        );

    }

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    ChargedKECalibration calibration("chargedKE.root");

    SelectionHistograms hSelSig, hSelBkg;

    if (doSelection) {
        hSelSig.Book("Sig");
        hSelBkg.Book("Bkg");
    }

    std::cout << "Processing signal..." << std::endl;
    RunSignalLoop(t, brVtx, dedxTABLE, calibration, cfg,
                &hPi0, &hTruth, &hCharged, &hEvt,
                doSelection ? &hSelSig : nullptr);

    if (doSelection) {
        //Open Background input file and branch setup
        TFile* fBkg = TFile::Open(background_inputfile.c_str());
        if (!fBkg || fBkg->IsZombie()) {
            std::cerr << "Cannot open background file" << std::endl; 
            return 1; 
        }
        TTree* tBkg = (TTree*)fBkg->Get("digitizedHits");
        if (!tBkg) { 
            std::cerr << "Tree digitizedHits not found in background file" << std::endl; 
            return 1; 
        }
        std::cout << "Processing background..." << std::endl;
        RunBackgroundLoop(tBkg, dedxTABLE, calibration, hSelBkg);
        fBkg->Close(); delete fBkg;

        // Selection plots Signal-Background(cosmic) overlaid
        hSelSig.PlotOverlay(hSelBkg, "Selection/");
    }

    // Other Analysis-related plotting 

    if (cfg.doPi0Analysis)     { hPi0.Plot(nentries, "preSelection/"); hPi0.Cleanup();     }
    if (cfg.doTruthAnalysis)   { hTruth.Plot("preSelection/");         hTruth.Cleanup();   }
    if (cfg.doChargedAnalysis) { hCharged.Plot("preSelection/");       hCharged.Cleanup(); }
    if (cfg.doEventVariables)  { hEvt.Plot("preSelection/");           hEvt.Cleanup();     }

    f->Close(); delete f;
    return 0;
}