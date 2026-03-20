#ifndef SELECTIONHISTOGRAMS_H
#define SELECTIONHISTOGRAMS_H

#include "TH1F.h"
#include "EventVariables.hpp"


struct SelectionHistograms {
    TH1F* hNCharged       = nullptr;
    TH1F* hNNeutral       = nullptr;
    TH1F* hEtotal         = nullptr;
    TH1F* hSphericity     = nullptr;
    TH1F* hMaxTrackAngle  = nullptr;
    TH1F* hVertexRadius   = nullptr;
    TH1F* hNPi0           = nullptr;

    void Book(const std::string& prefix) {
        hNCharged      = new TH1F((prefix+"_nCharged").c_str(),
                            ";N charged tracks;Events", 10, -0.5, 9.5);
        hNNeutral      = new TH1F((prefix+"_nNeutral").c_str(),
                            ";N neutral clusters;Events", 10, -0.5, 9.5);
        hEtotal        = new TH1F((prefix+"_Etotal").c_str(),
                            ";Total reco energy [MeV];Events", 100, 0, 3000);
        hSphericity    = new TH1F((prefix+"_sphericity").c_str(),
                            ";Sphericity;Events", 50, 0, 1);
        hMaxTrackAngle = new TH1F((prefix+"_maxAngle").c_str(),
                            ";Max track opening angle [rad];Events", 60, 0, TMath::Pi());
        hVertexRadius  = new TH1F((prefix+"_vtxR").c_str(),
                            ";Vertex radius [cm];Events", 50, 0, 20);
        hNPi0          = new TH1F((prefix+"_nPi0").c_str(),
                            ";N #pi^{0} candidates;Events", 8, -0.5, 7.5);
    }

    void Fill(const EventVariables& ev) {
        hNCharged      ->Fill(ev.nChargedTracks);
        hNNeutral      ->Fill(ev.nNeutralClusters);
        hEtotal        ->Fill(ev.totalRecoEnergy);
        hSphericity    ->Fill(ev.sphericity);
        hMaxTrackAngle ->Fill(ev.maxTrackAngle);
        hVertexRadius  ->Fill(ev.vertexRadius);
        hNPi0          ->Fill(ev.nPi0Candidates);
    }
};

// In analysis main file we would call it like this! 
// SelectionHistograms hSig, hBkg;
// hSig.Book("sig");
// hBkg.Book("bkg");

#endif