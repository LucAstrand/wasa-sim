#include "Calibration.hpp"


ChargedKECalibration::ChargedKECalibration() {
    h_ke = new TH3D(
        "h_ke",
        "Charged KE template;N_{ch};E_{EM} [MeV];KE_{ch} [MeV]",
        10, 0, 10,
        50, 0, 2000,
        100, 0, 2000
    );
    h_ke->Sumw2();
}

ChargedKECalibration::ChargedKECalibration(const std::string& filename) {
    TFile f(filename.c_str(), "READ");
    h_ke = dynamic_cast<TH3D*>(f.Get("h_ke"));
    if (!h_ke) {
        throw std::runtime_error("ChargedKECalibration: h_ke not found");
    }
    h_ke = static_cast<TH3D*>(h_ke->Clone());
    h_ke->SetDirectory(nullptr);
    f.Close();
}


void ChargedKECalibration::FillCalibration(int Nch, double Eem, double KEtrue) {
    h_ke->Fill(Nch, Eem, KEtrue);
}

void ChargedKECalibration::Finalize() {
    // For now we do not do anything on the end product, could make some fits or other...
}

double ChargedKECalibration::GetMeanKE(int Nch, double Eem) const {
    int bx = h_ke->GetXaxis()->FindBin(Nch);
    int by = h_ke->GetYaxis()->FindBin(Eem);
    std::unique_ptr<TH1D> proj(
        h_ke->ProjectionZ("_pz", bx, bx, by, by)
    );
    return proj->GetMean();
}

double ChargedKECalibration::GetRMSKE(int Nch, double Eem) const {
    int bx = h_ke->GetXaxis()->FindBin(Nch);
    int by = h_ke->GetYaxis()->FindBin(Eem);
    std::unique_ptr<TH1D> proj(
        h_ke->ProjectionZ("_pz", bx, bx, by, by)
    );
    return proj->GetRMS();
}

void ChargedKECalibration::Save(const std::string& filename) const {
    TFile f(filename.c_str(), "RECREATE");
    h_ke->Write();
    f.Close();
}


void DoCalibration(
    TTree* t,
    Long64_t nentries,
    const std::vector<int>& pi0_per_event,
    const std::vector<double>* centerXs,
    const std::vector<double>* centerYs,
    const std::vector<double>* centerZs,
    const std::vector<double>* energies,
    const std::vector<double>* primaryX,
    const std::vector<double>* primaryY,
    const std::vector<double>* primaryZ,
    const std::vector<double>* primaryEkin,
    const std::vector<double>* TPC_firstPosX,
    const std::vector<double>* TPC_firstPosY,
    const std::vector<double>* TPC_firstPosZ,
    const std::vector<double>* TPC_lastPosX,
    const std::vector<double>* TPC_lastPosY,
    const std::vector<double>* TPC_lastPosZ,
    const std::vector<double>* TPC_TrueKE,
    const std::vector<double>*    TPC_pdg,
    const std::vector<double>* TPC_dEdx,
    const std::vector<double>* TPC_smearedEdep,
    const std::vector<double>* TPC_PathLength,
    const std::string& outFile
) {
    ChargedKECalibration calibration;

    std::vector<Hit> hits;
    std::vector<ChargedTrack> chargedTracks;

    for (Long64_t ievt = 0; ievt < nentries; ++ievt) {
        t->GetEntry(ievt);

        // --- Pi0 count (optional but kept) ---
        int nPi0 = (ievt < (Long64_t)pi0_per_event.size())
                 ? pi0_per_event[ievt]
                 : 0;

        // --- Vertex ---
        if (!primaryX || primaryX->empty()) continue;
        TVector3 vertex((*primaryX)[0],
                         (*primaryY)[0],
                         (*primaryZ)[0]);

        // --- Generator KE (not used directly, but sanity-checked) ---
        if (!primaryEkin || primaryEkin->empty()) continue;

        // --- Hits ---
        if (!energies || energies->empty()) continue;
        hits.clear();
        for (size_t k = 0; k < energies->size(); ++k) {
            hits.push_back({
                (*centerXs)[k],
                (*centerYs)[k],
                (*centerZs)[k],
                (*energies)[k]
            });
        }

        // --- Charged tracks ---
        if (!TPC_firstPosX || !TPC_lastPosX) continue;
        chargedTracks.clear();

        size_t nChargedTracks = TPC_firstPosX->size();
        for (size_t k=0; k<nChargedTracks; ++k) {
            ChargedTrack trk;
            trk.id            = k;
            trk.vertex        = vertex;
            trk.exitPoint     = TVector3((*TPC_lastPosX)[k], (*TPC_lastPosY)[k], (*TPC_lastPosZ)[k]);
            trk.direction     = trk.exitPoint
                            - TVector3((*TPC_firstPosX)[k], (*TPC_firstPosY)[k], (*TPC_firstPosZ)[k]);
            trk.direction     = trk.direction.Unit();
            trk.TrueKE        = (*TPC_TrueKE)[k];
            trk.TruePDG       = (*TPC_pdg)[k];
            trk.clusterdEdx   = (*TPC_dEdx)[k];
            trk.EdepSmeared   = (*TPC_smearedEdep)[k];
            trk.pathLength    = (*TPC_PathLength)[k];
            trk.dEdxTheory    = 0.0;
            trk.resolution    = 0.15;

            chargedTracks.push_back(trk);
        }

        // --- Reconstruct ---
        RecoEvent reco = ReconstructEvent(hits, chargedTracks, vertex);

        // --- True charged pion KE ---
        double KEtrue = 0.0;
        int Nch = reco.chargedClusters.size();

        for (size_t i = 0; i < TPC_pdg->size(); ++i) {
            if ((*TPC_pdg)[i] == 211)
                KEtrue += (*TPC_TrueKE)[i]; // this is a global thing, like multiple pions!!!! 
        }

        calibration.FillCalibration(
            Nch,
            reco.EM_energy,
            KEtrue
        );
    }

    calibration.Finalize();
    calibration.Save(outFile);

    std::cout << "[Calibration] Saved to " << outFile << std::endl;
}


