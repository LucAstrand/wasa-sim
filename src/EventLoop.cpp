#include "EventLoop.hpp"

void RunSignalLoop(
    TTree* tree,
    BranchManagerVertex& brVtx,
    const DEDXTable& dedxTable,
    const ChargedKECalibration& calibration,
    const AnalysisConfig& cfg,
    Pi0Histograms* hPi0,
    TruthHistograms* hTruth,
    ChargedHistograms* hCharged,
    EventVarHistograms* hEvt,
    SelectionHistograms* hSel,
    CorrelationMatrix* hCorr)
{
    BranchManagerInput br;
    br.SetBranches(tree);

    Long64_t nentries = tree->GetEntries();
    progressbar bar(nentries);

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

    pi0_per_event.resize(nentries, 0);
    chPi_per_event.resize(nentries, 0);

    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        bar.update(); // Visual feedback on loop progress.
        tree->GetEntry(ievt);

        // int nPi0 = 0;
        // if (ievt < (Long64_t)pi0_per_event.size()) {
        //     nPi0 = pi0_per_event[ievt];
        // }

        // EVENT Vertices
        if (!br.primaryX||br.primaryX->empty()) continue;
        TVector3 tVertexVec((*br.primaryX)[0], (*br.primaryY)[0], (*br.primaryZ)[0]);
        Vtx vertex = brVtx.GetVertex(ievt, tVertexVec);  // falls back to truth if no reco vertex

        // Build true pi0s and Pi+-
        size_t nPrimaries = br.primaryPDG->size();
        primaryPi0s.clear();
        primaryChPis.clear();
        // pi0_per_event.push_back(0);
        // chPi_per_event.push_back(0);
        double pi0M  = 134.977;
        double pipmM = 139.570;
        for (size_t i=0; i<nPrimaries; ++i) {
            if ((*br.primaryPDG)[i] == 111) { // pi0
                // pi0_per_event[ievt]++;
                pi0_per_event.at(ievt)++;
                // std::cout << "[HELLO] pPi0 Ekin: " << (*br.primaryEkin)[i] << std::endl;
                primaryPi0s.push_back({(*br.primaryTrackID)[i], TLorentzVector(TVector3((*br.primaryPx)[i], (*br.primaryPy)[i], (*br.primaryPz)[i]), (*br.primaryEkin)[i] + pi0M)});
            }
            if (std::abs((*br.primaryPDG)[i]) == 211) { // pi+-
                // chPi_per_event[ievt]++;
                chPi_per_event.at(ievt)++;
                primaryChPis.push_back({(*br.primaryTrackID)[i], TLorentzVector(TVector3((*br.primaryPx)[i], (*br.primaryPy)[i], (*br.primaryPz)[i]), (*br.primaryEkin)[i] + pipmM)});
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

        if (cfg.doChargedAnalysis) {

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
            size_t nChargedTracks = br.TPC_trackID->size();
            for (size_t k=0; k<nChargedTracks; ++k) {
                if (!br.TPC_trackID || !br.TPC_firstPosX || !br.TPC_lastPosX) continue;
                // chargedTracks.push_back(
                //     {(*br.TPC_trackID)[k],
                //     vertex.vertexVec,
                //     TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k]),
                //     TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k]) - TVector3((*br.TPC_firstPosX)[k], (*br.TPC_firstPosY)[k], (*br.TPC_firstPosZ)[k]),
                //     (*br.TPC_TrueKE)[k], (*br.TPC_pdg)[k], (*br.TPC_smearedDedx)[k], (*br.TPC_theoryDedx)[k], 0.15});
                ChargedTrack trk;
                trk.id            = (*br.TPC_trackID)[k];
                trk.vertex        = vertex.vertexVec;
                trk.exitPoint     = TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k]);
                trk.direction     = trk.exitPoint
                                - TVector3((*br.TPC_firstPosX)[k], (*br.TPC_firstPosY)[k], (*br.TPC_firstPosZ)[k]);
                trk.direction     = trk.direction.Unit();
                trk.TrueKE        = (*br.TPC_TrueKE)[k];
                trk.TruePDG       = (*br.TPC_pdg)[k];
                trk.smearedDedx   = (*br.TPC_smearedDedx)[k];
                trk.theoryDedx    = (*br.TPC_theoryDedx)[k];
                trk.resolution    = 0.15;
                chargedTracks.push_back(trk);
            }
        }

        if (cfg.doTruthAnalysis) {
            trueHits.clear();
            truePi0s.clear();
            size_t nTrueHits = br.truePhotonE->size();
            for (size_t k=0; k<nTrueHits; ++k) {
                trueHits.push_back({(*br.truePhotonPosX)[k], (*br.truePhotonPosY)[k], (*br.truePhotonPosZ)[k], (*br.truePhotonE)[k], (*br.truePhotonTrackID)[k], (*br.truePhotonParentID)[k]});
            }
            truePhotons = TruePhotonBuilder(trueHits, vertex.vertexVec);
            truePi0s = TruePi0Builder(truePhotons);
            // for (TruePi0 tpi0 : truePi0s) {
            //     if (hTruth->hTrueMass) hTruth->hTrueMass->Fill(tpi0.p4.M());
            // }

            for (primaryPi0 pPi0 : primaryPi0s) {
                // std::cout << "[HELLO] pPi0 M: " << pPi0.p4.M() << std::endl;
                if (hTruth->hTrueMass) hTruth->hTrueMass->Fill(pPi0.p4.M());
            }
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Reconstruct Event
        RecoEvent reco = ReconstructEvent(hits, chargedTracks, vertex.vertexVec, dedxTable);

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Charged cluster angle metrics
        if (cfg.doChargedAnalysis) {
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
                        if (hCharged->hAngle_matched) hCharged->hAngle_matched->Fill(theta_deg);
                        if (hCharged->hAngleVsE) hCharged->hAngleVsE->Fill(h.e, theta_deg);
                    } else {
                        if (hCharged->hAngle_unmatched) hCharged->hAngle_unmatched->Fill(theta_deg);
                    }
                }
            }
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Invariant mass plot (neutral pion / double photon)
        if (cfg.doPi0Analysis) {
            selected.clear();

            double param_h = 135;
            double param_k = 1.6;
            double param_a = 45;
            double param_b = 1.55;
            double relipse = 0;

            std::vector<Pi0Candidate> candidates;

            if (cfg.doChargedAnalysis) {
                // add photons built from conversions!
                std::vector<ConversionCandidate> conversions = FindConversions(reco.chargedClusters, vertex.vertexVec);
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

                    if (hPi0->hppM_pre) hPi0->hppM_pre->Fill(mgg, theta);

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
                if (hPi0->hMass && hPi0->hppM_post) {
                    hPi0->hMass->Fill(pi0.p4.M());
                    hPi0->hppM_post->Fill(pi0.mgg, pi0.theta);
                }
            }

            reco.nPionMultiplicity = selected.size();

        }
        
        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Charged Object Clustering
        if (cfg.doChargedAnalysis) {

            for (const ChargedCluster& cluster : reco.chargedClusters) {
                if (cluster.isOrphanElectron) continue;  // skip orphan electrons
                if (cluster.isUsedInConversion) continue; // skip electrons used in conversion reco

                // n-Sigma plots
                if (hCharged->hNSigmaPion) hCharged->hNSigmaPion->Fill(cluster.nSigmaPion);
                if (hCharged->hNSigmaProton) hCharged->hNSigmaProton->Fill(cluster.nSigmaProton);
                //--> Pions
                if (cluster.objectTruePDG == 211 || cluster.objectTruePDG == -211) {
                // if (hCharged->hdEdxVsE_cluster_Pion) hCharged->hdEdxVsE_cluster_Pion->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
                if (hCharged->hdEdxVsE_cluster_Pion) hCharged->hdEdxVsE_cluster_Pion->Fill(cluster.objectTrueKE, cluster.objectSmearedDedx); // ORDER: X vs Y
                if (hCharged->hClusterE) hCharged->hClusterE->Fill(cluster.totalEnergy);
                if (hCharged->hdEdxVsE_true_Pion) hCharged->hdEdxVsE_true_Pion->Fill(cluster.objectTrueKE, cluster.objectTheoryDedx);
                if (hCharged->hdEdxTruePion) hCharged->hdEdxTruePion->Fill(cluster.objectTheoryDedx);
                if (hCharged->hdEdxSmearPion) hCharged->hdEdxSmearPion->Fill(cluster.objectSmearedDedx);
                }
                //--> Protons
                if (cluster.objectTruePDG == 2212) {
                // if (hCharged->hdEdxVsE_cluster_Proton) hCharged->hdEdxVsE_cluster_Proton->Fill(cluster.totalEnergy, cluster.clusterdEdx); // ORDER: X vs Y
                if (hCharged->hdEdxVsE_cluster_Proton) hCharged->hdEdxVsE_cluster_Proton->Fill(cluster.objectTrueKE, cluster.objectSmearedDedx); // ORDER: X vs Y
                if (hCharged->hdEdxVsE_true_Proton) hCharged->hdEdxVsE_true_Proton->Fill(cluster.objectTrueKE, cluster.objectTheoryDedx);
                if (hCharged->hdEdxTrueProton) hCharged->hdEdxTrueProton->Fill(cluster.objectTheoryDedx);
                if (hCharged->hdEdxSmearProton) hCharged->hdEdxSmearProton->Fill(cluster.objectSmearedDedx);
                }
                
                // E res plots
                double totalRecoE = cluster.totalEnergy;//+ cluster.EdepSmeared;
                // double Eres = (cluster.totalEnergy - cluster.objectTrueKE) / cluster.objectTrueKE;
                double Eres = (totalRecoE - cluster.objectTrueKE) / cluster.objectTrueKE;
                if (hCharged->h2_Eres) hCharged->h2_Eres->Fill(cluster.objectTrueKE, Eres);
            }
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Truth level + reco level plots
        if (cfg.doTruthAndMix) {
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

            if (hTruth->h_tE_rA) fillPairs(photons_tE_rA, hTruth->h_tE_rA);
            if (hTruth->h_rE_tA) fillPairs(photons_rE_tA, hTruth->h_rE_tA);
        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Pi0 eff & acc
        if (cfg.doPi0Analysis) {
            if (hPi0->effPlotter) hPi0->effPlotter->ProcessEvent(truePi0s, selected, reco.clusters);
            if (hPi0->accEkin) hPi0->accEkin->Pi0ProcessSignalEvent(truePi0s, primaryPi0s);
            if (hPi0->accTheta) hPi0->accTheta->Pi0ProcessSignalEvent(truePi0s, primaryPi0s);
            if (hPi0->acc2D) hPi0->acc2D->Pi0ProcessSignalEvent(truePi0s, primaryPi0s);

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
        if (cfg.doChargedAnalysis) {
            hCharged->pidEff->ProcessEvent(reco.chargedClusters);
            // chAccPlotter->ChPiProcessSignalEvent(reco.chargedClusters, primaryChPis);
            //Acc
            hCharged->chAccGlobal->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, 0);
            hCharged->chAccCondTPC->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, 1);
            hCharged->chAccCondNoTPC->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, 2);
            hCharged->accNChTracks->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, 0);
            //Eff
            hCharged->chEffTheta->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, reco.chargedClusters);
            hCharged->chEffNChTracks->ChPiProcessSignalEvent(trueChPiInCals, primaryChPis, reco.chargedClusters);
            // assign a charged pion multiplicity based on PID Guess --> TPC info
            for (ChargedCluster ch : reco.chargedClusters) {
                if (ch.pidGuess == PID::Pion) {
                    //safeguards for electrons they should not happen but still...
                    if (ch.isOrphanElectron) continue;
                    if (ch.isUsedInConversion) continue;
                    reco.chPionMultiplicity += 1;
                }
            }
            std::unordered_map<int, const primaryChPi*> p4Map;
            p4Map.reserve(primaryChPis.size());
            for (const auto& p : primaryChPis) p4Map.emplace(p.trackID, &p);

            for (const auto& d : trueChPiInCals) {
                if (!d.throughTPC) continue;
                auto it = p4Map.find(d.trackID);
                if (it == p4Map.end()) continue;
                double val = it->second->p4.Theta();
                if (hCharged->hPionThetaWASA) hCharged->hPionThetaWASA->Fill(val);
            }

            for (primaryChPi chPi: primaryChPis) {
                if (hCharged->hPionTheta) hCharged->hPionTheta->Fill(chPi.p4.Theta());
                if (hCharged->hPionCosTheta) hCharged->hPionCosTheta->Fill(std::cos(chPi.p4.Theta()));
                if (hCharged->hPionPhi) hCharged->hPionPhi->Fill(chPi.p4.Phi());                
            }

        }

        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

        // Event variables
        if (cfg.doEventVariables) {

            double eVis = reco.EM_energy; // this is the sum of all the hit energies in the calorimeter, from any source.
            double eTrue = 0.0;
            size_t n = br.primaryEkin->size();
            for (size_t i = 0; i < n; ++i) {
                eTrue += (*br.primaryEkin)[i];
            }
            double calibratedKE = calibration.GetMeanKE(reco.chargedClusters.size(), eVis);
            double eReco = eVis + calibratedKE;

            if (hEvt->hEvisVsEtrue) hEvt->hEvisVsEtrue->Fill(eTrue, eVis);
            if (hEvt->hErecoVsEtrue) hEvt->hErecoVsEtrue->Fill(eTrue, eReco);
            if (hEvt->hDiffErecoEtrue) hEvt->hDiffErecoEtrue->Fill(eReco - eTrue);
            if (hEvt->hEvis) hEvt->hEvis->Fill(eVis);
            if (hEvt->hEcorrected) hEvt->hEcorrected->Fill(eReco);

            if (hEvt->hNPionMult) hEvt->hNPionMult->Fill(pi0_per_event[ievt], reco.nPionMultiplicity);
            if (hEvt->hChPionMult) hEvt->hChPionMult->Fill(chPi_per_event[ievt], reco.chPionMultiplicity);

            EventVariables ev = ComputeEventVariables(ievt, reco, vertex, calibration);
            // if (ievt < 5) {  // print first 5 events
            //     std::cout << "[EventVars] nCharged=" << ev.nChargedTracks
            //             << " nNeutral=" << ev.nNeutralClusters
            //             << " Etotal=" << ev.totalRecoEnergy
            //             << " sphericity=" << ev.sphericity << std::endl;
            // }

            hSel->Fill(ev);
            hCorr->FillCorrelation(ev);

        }
        
        //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
        
    }
    std::cout << std::endl;
}

void RunBackgroundLoop(
    TTree* tree,
    const DEDXTable& dedxTable,
    const ChargedKECalibration& calibration,
    SelectionHistograms& hSel,
    CorrelationMatrix& hCorr)
{
    BranchManagerInput br;
    br.SetBranches(tree);

    Long64_t nentries = tree->GetEntries();
    progressbar bar(nentries);

    std::vector<Hit>          hits;
    std::vector<ChargedTrack> chargedTracks;

    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        bar.update();
        tree->GetEntry(ievt);
        
        // Skip empty events - these occur when MCPL ran out of 
        // particles before GEANT4 finished its requested events
        if (!br.primaryPDG || br.primaryPDG->empty()) continue;
        if (!br.energies   || br.energies->empty())   continue;

        // Have to fix vertexing for the background! 
        // TVector3 vertex(0, 0, 0);
        Vtx vertex;
        vertex.vertexVec = TVector3{0,0,0};

        // Build hits
        hits.clear();
        if (!br.energies || br.energies->empty()) continue;
        for (size_t k = 0; k < br.energies->size(); ++k)
            hits.push_back({(*br.centerXs)[k], (*br.centerYs)[k],
                            (*br.centerZs)[k], (*br.energies)[k]});

        // Build charged tracks from TPC
        chargedTracks.clear();
        if (br.TPC_trackID) {
            for (size_t k = 0; k < br.TPC_trackID->size(); ++k) {
                // chargedTracks.push_back({
                //     (*br.TPC_trackID)[k], vertex.vertexVec,
                //     TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k]),
                //     TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k])
                //   - TVector3((*br.TPC_firstPosX)[k],(*br.TPC_firstPosY)[k],(*br.TPC_firstPosZ)[k]),
                //     (*br.TPC_TrueKE)[k], (*br.TPC_pdg)[k], (*br.TPC_dEdx)[k],
                //     (*br.TPC_smearedEdep)[k], (*br.TPC_PathLength)[k],
                //     0, 0.15, (*br.TPC_nSteps)[k]
                // });
                ChargedTrack trk;
                trk.id            = k;
                trk.vertex        = vertex.vertexVec;
                trk.exitPoint     = TVector3((*br.TPC_lastPosX)[k], (*br.TPC_lastPosY)[k], (*br.TPC_lastPosZ)[k]);
                trk.direction     = trk.exitPoint
                                - TVector3((*br.TPC_firstPosX)[k], (*br.TPC_firstPosY)[k], (*br.TPC_firstPosZ)[k]);
                trk.direction     = trk.direction.Unit();
                trk.TrueKE        = (*br.TPC_TrueKE)[k];
                trk.TruePDG       = (*br.TPC_pdg)[k];
                //TEMPTEMPTEMPTEMP
                // trk.smearedDedx   = (*br.TPC_smearedDedx)[k];
                // trk.theoryDedx    = (*br.TPC_theoryDedx)[k];
                trk.resolution    = 0.15;
                chargedTracks.push_back(trk);
            }
        }

        RecoEvent reco = ReconstructEvent(hits, chargedTracks, vertex.vertexVec, dedxTable);

        EventVariables ev = ComputeEventVariables(ievt, reco, vertex, calibration);
        // if (ievt < 5) {  // print first 5 events
        //     std::cout << "[EventVars] nCharged=" << ev.nChargedTracks
        //             << " nNeutral=" << ev.nNeutralClusters
        //             << " Etotal=" << ev.totalRecoEnergy
        //             << " sphericity=" << ev.sphericity << std::endl;
        // }
        hSel.Fill(ev);
        hCorr.FillCorrelation(ev);
    }
    std::cout << std::endl;
}