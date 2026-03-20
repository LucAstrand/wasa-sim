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
                    << "  --event-variables     Computes the event level variables\n";
            return 1;
        }
    SetErrorHandler(MyErrorHandler);

    std::string root_inputfile = argv[1];
    std::string background_inputfile = argv[2];
    std::string vertices_inputfile = argv[3];
    bool doCalibration = false;
    bool doPi0Analysis = false;
    bool doChargedAnalysis = false;
    bool doTruthAnalysis = false;
    bool doTruthAndMixPlots = false;
    bool doEventVariables = false;

    for (int i = 4; i<argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--calibration") doCalibration = true;
        if (arg == "--full-analysis") {
            if (!std::filesystem::exists("chargedKE.root")) {
                std::cerr << "No calibration root file detected, please run './Analysis <input.root> <background.root> <vertices.root> --calibration' first" << std::endl;
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
    if (!t) { std::cerr << "Tree digitizedHits not found" << std::endl; return 1;}
    // Branch setup
    BranchManagerInput br;
    br.SetBranches(t);

    // Histogram booking 
    Pi0Histograms hPi0;
    TruthHistograms hTruth;
    ChargedHistograms hCharged;
    EventVarHistograms hEvt;

    if (doPi0Analysis) hPi0.Book();
    if (doTruthAnalysis) hTruth.Book();
    if (doChargedAnalysis) hCharged.Book();
    if (doEventVariables) hEvt.Book();

    Long64_t nentries = t->GetEntries();

    TFile *vtxFile = TFile::Open(vertices_inputfile.c_str());
    if (!vtxFile || vtxFile->IsZombie()) return 1;

    TTree *vtxTree = (TTree *)vtxFile->Get("vertices");
    if (!vtxTree) {std::cerr << "Tree vertices not found" << std::endl; return 1;}

    BranchManagerVertex brVtx;   
    brVtx.SetBranches(vtxTree);
    brVtx.LoadVertices(vtxTree, nentries);

    vtxFile->Close();

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    //Containers used throughout the analysis
    std::vector<Hit> hits;
    std::vector<TruePhotonHit> trueHits;
    std::vector<TruePhoton> truePhotons;
    std::vector<TruePi0> truePi0s;
    std::vector<Pi0Candidate> selected;
    std::vector<ChargedTrack> chargedTracks;
    std::vector<primaryPi0> primaryPi0s;
    std::vector<primaryChPi> primaryChPis;
    std::vector<int> pi0_per_event;
    std::vector<int> chPi_per_event;
    std::vector<TrueChPiInCal> trueChPiInCals;
    std::vector<TrueChPiDecayed> trueChPiDecayed;

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    DEDXTable dedxTABLE("dedx_tables_Ar80CO2.root");

    //Calibration procedure ----->>>> HAVE TO DOUBLE CHECK THIS LOGIC. A lot has changed SINCE THEN!!!!!
    if (doCalibration) {
        DoCalibration(
            t,
            static_cast<int>(nentries * 0.7),
            pi0_per_event,
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

    progressbar bar(nentries);

    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        bar.update(); // Visual feedback on loop progress.
        t->GetEntry(ievt);

        int nPi0 = 0;
        if (ievt < (Long64_t)pi0_per_event.size()) {
            nPi0 = pi0_per_event[ievt];
        }

        // EVENT Vertices
        if (!br.primaryX||br.primaryX->empty()) continue;
        TVector3 tVertex((*br.primaryX)[0], (*br.primaryY)[0], (*br.primaryZ)[0]);
        TVector3 vertex = brVtx.GetVertex(ievt, tVertex);  // falls back to truth if no reco vertex

        // Build true pi0s and Pi+-
        size_t nPrimaries = br.primaryPDG->size();
        primaryPi0s.clear();
        primaryChPis.clear();
        for (size_t i=0; i<nPrimaries; ++i) {
            pi0_per_event.push_back(0);
            chPi_per_event.push_back(0);
            if ((*br.primaryPDG)[i] == 111) { // pi0
                pi0_per_event[ievt]++;
                primaryPi0s.push_back({(*br.primaryTrackID)[i], TLorentzVector(TVector3((*br.primaryPx)[i], (*br.primaryPy)[i], (*br.primaryPz)[i]), (*br.primaryEkin)[i])});
            }
            if (std::abs((*br.primaryPDG)[i]) == 211) { // pi+-
                chPi_per_event[ievt]++;
                primaryChPis.push_back({(*br.primaryTrackID)[i], TLorentzVector(TVector3((*br.primaryPx)[i], (*br.primaryPy)[i], (*br.primaryPz)[i]), (*br.primaryEkin)[i])});
            }
        }
        hits.clear();
        if (!br.energies || br.energies->empty()) {
            std::cout << "Event " << ievt << ": No hits (skipping)" << std::endl;
            continue;
        }
        size_t nHits = br.energies->size();
        for (size_t k=0; k<nHits; ++k) {
            hits.push_back({(*br.centerXs)[k], (*br.centerYs)[k], (*br.centerZs)[k], (*br.energies)[k]});
        }

        if (doChargedAnalysis) {

            trueChPiInCals.clear();
            size_t nTrueChPiInCal = br.TrueChargedPionTrackID->size();
            for (size_t i=0; i<nTrueChPiInCal; ++i) {
                TrueChPiInCal chpi;
                chpi.trackID = (*br.TrueChargedPionTrackID)[i];
                chpi.throughTPC = (*br.TrueChargedPionThroughTPC)[i];
                trueChPiInCals.push_back(chpi);
            }
            trueChPiDecayed.clear();
            size_t nTrueChPiDecayed = br.TrueChargedPionDecayedTrackID->size();
            for (size_t j=0; j<nTrueChPiDecayed; ++j) {
                TrueChPiDecayed chpi;
                chpi.trackID = (*br.TrueChargedPionDecayedTrackID)[j];
                chpi.beforeCal = (*br.TrueChargedPionDecayedBeforeCal)[j];
                chpi.beforeTPC = (*br.TrueChargedPionDecayedBeforeTPC)[j];
                trueChPiDecayed.push_back(chpi);
            }
            chargedTracks.clear();
            size_t nChargedTracks = br.TPC_Edep->size();
            for (size_t k=0; k<nChargedTracks; ++k) {
                if (!br.TPC_Edep || !br.TPC_firstPosX || !br.TPC_lastPosX) continue;
                chargedTracks.push_back(
                    {(*br.TPC_trackID)[k],
                    vertex,
                    TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k]),
                    TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k]) - TVector3((*br.TPC_firstPosX)[k], (*br.TPC_firstPosY)[k], (*br.TPC_firstPosZ)[k]),
                    (*br.TPC_TrueKE)[k], (*br.TPC_pdg)[k], (*br.TPC_dEdx)[k], (*br.TPC_smearedEdep)[k], (*br.TPC_PathLength)[k], 0 /* Placeholder */, 0.15, (*br.TPC_nSteps)[k]});
            }
        }

        if (doTruthAnalysis) {
            trueHits.clear();
            truePi0s.clear();
            size_t nTrueHits = br.truePhotonE->size();
            for (size_t k=0; k<nTrueHits; ++k) {
                trueHits.push_back({(*br.truePhotonPosX)[k], (*br.truePhotonPosY)[k], (*br.truePhotonPosZ)[k], (*br.truePhotonE)[k], (*br.truePhotonTrackID)[k], (*br.truePhotonParentID)[k]});
            }
            truePhotons = TruePhotonBuilder(trueHits, vertex);
            truePi0s = TruePi0Builder(truePhotons);
            for (TruePi0 tpi0 : truePi0s) {
                if (hTruth.hTrueMass) hTruth.hTrueMass->Fill(tpi0.p4.M());
            }
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Reconstruct Event
        RecoEvent reco = ReconstructEvent(hits, chargedTracks, vertex, dedxTABLE);

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Charged cluster angle metrics
        if (doChargedAnalysis) {
            for (const auto& cluster : reco.chargedClusters) {
                if (std::abs(cluster.objectTruePDG) != 211 && 
                    std::abs(cluster.objectTruePDG) != 2212) continue;

                TVector3 exitPoint  = cluster.TPCExitPoint;
                TVector3 dir        = cluster.direction;
                if (dir.Mag2() < 1e-12) continue;

                // Check all hits - both matched and unmatched
                for (const auto& h : hits) {
                    TVector3 hitPos(h.x, h.y, h.z);
                    TVector3 delta = hitPos - exitPoint;
                    if (delta.Mag2() < 1e-12) continue;

                    double dot = dir.Unit().Dot(delta.Unit());
                    dot = std::clamp(dot, -1.0, 1.0);
                    if (dot < 0) continue;
                    double theta_deg = std::acos(dot) * TMath::RadToDeg();

                    bool isMatched = (h.owner == HitOwner::Charged);
                    if (isMatched) {
                        if (hCharged.hAngle_matched) hCharged.hAngle_matched->Fill(theta_deg);
                        if (hCharged.hAngleVsE) hCharged.hAngleVsE->Fill(h.e, theta_deg);
                    } else {
                        if (hCharged.hAngle_unmatched) hCharged.hAngle_unmatched->Fill(theta_deg);
                    }
                }
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

            if (doChargedAnalysis) {
                // add photons built from conversions!
                std::vector<ConversionCandidate> conversions = FindConversions(reco.chargedClusters, vertex);
                reco.clusters.reserve(reco.clusters.size() + conversions.size());

                for (const auto& conv : conversions) {
                    Cluster neutralCluster;
                    neutralCluster.p4       = conv.p4;
                    neutralCluster.centroid = conv.conversionVertex;
                    neutralCluster.isFromConversion = true;
                    reco.clusters.push_back(neutralCluster);
                }
                // After FindConversions:
                for (size_t i = 0; i < reco.chargedClusters.size(); ++i) {
                    auto& cluster = reco.chargedClusters[i];
                    
                    if (std::abs(cluster.objectTruePDG) != 11) continue;
                    
                    // Check if this electron was claimed by a conversion
                    bool usedInConversion = false;
                    for (const auto& conv : conversions) {
                        if (conv.track1_idx == (int)i || 
                            conv.track2_idx == (int)i) {
                            usedInConversion = true;
                            break;
                        }
                    }
                    if (usedInConversion) {
                        cluster.isUsedInConversion = true;
                    } else {
                        // Orphan electron - not part of a conversion pair
                        // Flag it so downstream code ignores it for PID/pi0 reco
                        cluster.isOrphanElectron = true;
                        // Its calorimeter energy is already counted in EM_energy
                        // via the hits loop - nothing else needed
                    }
                }
            }

            for (size_t a = 0; a < reco.clusters.size(); ++a) {
                for (size_t b = a + 1; b < reco.clusters.size(); ++b) {

                    const auto& g1 = reco.clusters[a].p4;
                    const auto& g2 = reco.clusters[b].p4;

                    double E1 = g1.E();
                    double E2 = g2.E();

                    double mgg = (g1 + g2).M();
                    double theta = openingAngle(g1, g2);

                    if (hPi0.hppM_pre) hPi0.hppM_pre->Fill(mgg, theta);

                    relipse = pow((mgg-param_h),2)/pow(param_a,2) + pow( (theta - param_k), 2)/pow(param_b,2);
                    if (relipse > 1) continue;

                    Pi0Candidate cand;
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
                if (hPi0.hMass && hPi0.hppM_post) {
                    hPi0.hMass->Fill(pi0.p4.M());
                    hPi0.hppM_post->Fill(pi0.mgg, pi0.theta);
                }
            }

            reco.nPionMultiplicity = selected.size();

        }
        
        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Charged Object Clustering
        if (doChargedAnalysis) {

            for (const ChargedCluster& cluster : reco.chargedClusters) {
                if (cluster.isOrphanElectron) continue;  // skip orphan electrons
                if (cluster.isUsedInConversion) continue; // skip electrons used in conversion reco

                // n-Sigma plots
                if (hCharged.hNSigmaPion) hCharged.hNSigmaPion->Fill(cluster.nSigmaPion);
                if (hCharged.hNSigmaProton) hCharged.hNSigmaProton->Fill(cluster.nSigmaProton);
                //--> Pions
                if (cluster.objectTruePDG == 211 || cluster.objectTruePDG == -211) {
                if (hCharged.hdEdxVsE_cluster_Pion) hCharged.hdEdxVsE_cluster_Pion->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
                // if (hCharged.hdEdxVsE_cluster_Pion) hCharged.hdEdxVsE_cluster_Pion->Fill(cluster.objectTrueKE, cluster.clusterdEdx); // ORDER: X vs Y
                if (hCharged.hClusterE) hCharged.hClusterE->Fill(cluster.totalEnergy);
                if (hCharged.hdEdxVsE_true_Pion) hCharged.hdEdxVsE_true_Pion->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                if (hCharged.hdEdxTruePion) hCharged.hdEdxTruePion->Fill(cluster.objectTruedEdx);
                if (hCharged.hdEdxSmearPion) hCharged.hdEdxSmearPion->Fill(cluster.clusterdEdx);
                }
                //--> Protons
                if (cluster.objectTruePDG == 2212) {
                if (hCharged.hdEdxVsE_cluster_Proton) hCharged.hdEdxVsE_cluster_Proton->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
                // if (hCharged.hdEdxVsE_cluster_Proton) hCharged.hdEdxVsE_cluster_Proton->Fill(cluster.objectTrueKE, cluster.clusterdEdx); // ORDER: X vs Y
                if (hCharged.hdEdxVsE_true_Proton) hCharged.hdEdxVsE_true_Proton->Fill(cluster.objectTrueKE, cluster.objectTruedEdx);
                if (hCharged.hdEdxTrueProton) hCharged.hdEdxTrueProton->Fill(cluster.objectTruedEdx);
                if (hCharged.hdEdxSmearProton) hCharged.hdEdxSmearProton->Fill(cluster.clusterdEdx);
                }
                
                // E res plots
                double totalRecoE = cluster.totalEnergy + cluster.EdepSmeared;
                // double Eres = (cluster.totalEnergy - cluster.objectTrueKE) / cluster.objectTrueKE;
                double Eres = (totalRecoE - cluster.objectTrueKE) / cluster.objectTrueKE;
                if (hCharged.h2_Eres) hCharged.h2_Eres->Fill(cluster.objectTrueKE, Eres);
            }
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

            if (hTruth.h_tE_rA) fillPairs(photons_tE_rA, hTruth.h_tE_rA);
            if (hTruth.h_rE_tA) fillPairs(photons_rE_tA, hTruth.h_rE_tA);
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Pi0 eff & acc
        if (doPi0Analysis) {
            if (hPi0.effPlotter) hPi0.effPlotter->ProcessEvent(truePi0s, selected, reco.clusters);
            if (hPi0.accPlotter) hPi0.accPlotter->Pi0ProcessSignalEvent(truePi0s, primaryPi0s);

            // THESE NEED TO BE FIXED!!!

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

        // Pi+- eff & acc
        if (doChargedAnalysis) {
            hCharged.pidEff->ProcessEvent(reco.chargedClusters);
            // chAccPlotter->ChPiProcessSignalEvent(reco.chargedClusters, primaryChPis);

            hCharged.chAccGlobal->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, 0);
            hCharged.chAccCondTPC->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, 1);
            hCharged.chAccCondNoTPC->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, 2);

            // assign a charged pion multiplicity based on PID Guess --> TPC info
            for (ChargedCluster ch : reco.chargedClusters) {
                if (ch.pidGuess == PID::Pion) {
                    //safeguards for electrons they should not happen but still...
                    if (ch.isOrphanElectron) continue;
                    if (ch.isUsedInConversion) continue;
                    reco.chPionMultiplicity += 1;
                }
            }
            for (primaryChPi chPi: primaryChPis) {
                if (hCharged.hPionTheta) hCharged.hPionTheta->Fill(chPi.p4.Theta());
                if (hCharged.hPionCosTheta) hCharged.hPionCosTheta->Fill(std::cos(chPi.p4.Theta()));
                if (hCharged.hPionPhi) hCharged.hPionPhi->Fill(chPi.p4.Phi());
            }

        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Event variables
        if (doEventVariables) {

            double eVis = reco.EM_energy; // this is the sum of all the hit energies in the calorimeter, from any source.
            double eTrue = 0.0;
            size_t n = br.primaryEkin->size();
            for (size_t i = 0; i < n; ++i) {
                eTrue += (*br.primaryEkin)[i];
            }
            double calibratedKE = calibration.GetMeanKE(reco.chargedClusters.size(), eVis);
            double eReco = eVis + calibratedKE;

            if (hEvt.hEvisVsEtrue) hEvt.hEvisVsEtrue->Fill(eTrue, eVis);
            if (hEvt.hErecoVsEtrue) hEvt.hErecoVsEtrue->Fill(eTrue, eReco);
            if (hEvt.hDiffErecoEtrue) hEvt.hDiffErecoEtrue->Fill(eReco - eTrue);
            if (hEvt.hEvis) hEvt.hEvis->Fill(eVis);
            if (hEvt.hEcorrected) hEvt.hEcorrected->Fill(eReco);

            if (hEvt.hNPionMult) hEvt.hNPionMult->Fill(pi0_per_event[ievt], reco.nPionMultiplicity);
            if (hEvt.hChPionMult) hEvt.hChPionMult->Fill(chPi_per_event[ievt], reco.chPionMultiplicity);

        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

    }

    // Plotting 

    if (doPi0Analysis)     { hPi0.Plot(nentries);   hPi0.Cleanup();     }
    if (doTruthAnalysis)   { hTruth.Plot();         hTruth.Cleanup();   }
    if (doChargedAnalysis) { hCharged.Plot();       hCharged.Cleanup(); }
    if (doEventVariables)  { hEvt.Plot();           hEvt.Cleanup();     }


    f->Close(); delete f;
    return 0;
}