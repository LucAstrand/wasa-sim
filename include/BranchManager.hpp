#ifndef BRANCHMANAGER_H
#define BRANCHMANAGER_H

#include "TTree.h"
#include <vector>
#include <iostream>

#include "Structures.hpp"

template <typename T>
bool SafeSetBranch(TTree* tree, const char* name, T*& ptr) {
    if (tree->GetListOfBranches()->FindObject(name)) {
        tree->SetBranchAddress(name, &ptr);
        return true;
    }
    std::cout << "[Warning] Branch " << name << " not found. Skipping.\n";
    return false;
}

struct BranchManagerInput {
    // Calorimeter hits
    std::vector<double> *centerXs  = nullptr;
    std::vector<double> *centerYs  = nullptr;
    std::vector<double> *centerZs  = nullptr;
    std::vector<double> *energies  = nullptr;

    // Primaries
    std::vector<double> *primaryX       = nullptr;
    std::vector<double> *primaryY       = nullptr;
    std::vector<double> *primaryZ       = nullptr;
    std::vector<double> *primaryEkin    = nullptr;
    std::vector<double> *primaryPx      = nullptr;
    std::vector<double> *primaryPy      = nullptr;
    std::vector<double> *primaryPz      = nullptr;
    std::vector<int>    *primaryPDG     = nullptr;
    std::vector<int>    *primaryTrackID = nullptr;

    // True photons
    std::vector<double> *truePhotonPosX    = nullptr;
    std::vector<double> *truePhotonPosY    = nullptr;
    std::vector<double> *truePhotonPosZ    = nullptr;
    std::vector<double> *truePhotonE       = nullptr;
    std::vector<double> *truePhotonTrackID = nullptr;
    std::vector<double> *truePhotonParentID= nullptr;

    // True charged pions
    std::vector<double> *TrueChargedPionX               = nullptr;
    std::vector<double> *TrueChargedPionY               = nullptr;
    std::vector<double> *TrueChargedPionZ               = nullptr;
    std::vector<double> *TrueChargedPionE               = nullptr;
    std::vector<double> *TrueChargedPionTrackID         = nullptr;
    std::vector<double> *TrueChargedPionThroughTPC      = nullptr;
    std::vector<double> *TrueChargedPionCreationX       = nullptr;
    std::vector<double> *TrueChargedPionCreationY       = nullptr;
    std::vector<double> *TrueChargedPionCreationZ       = nullptr;
    std::vector<double> *TrueChargedPionEndX            = nullptr;
    std::vector<double> *TrueChargedPionEndY            = nullptr;
    std::vector<double> *TrueChargedPionEndZ            = nullptr;
    std::vector<double> *TrueChargedPionDecayedBeforeCal= nullptr;
    std::vector<double> *TrueChargedPionDecayedBeforeTPC= nullptr;
    std::vector<double> *TrueChargedPionDecayedTrackID  = nullptr;

    // TPC
    std::vector<double> *TPC_trackID     = nullptr;
    std::vector<double> *TPC_Edep        = nullptr;
    std::vector<double> *TPC_smearedEdep = nullptr;
    std::vector<double> *TPC_firstPosX   = nullptr;
    std::vector<double> *TPC_firstPosY   = nullptr;
    std::vector<double> *TPC_firstPosZ   = nullptr;
    std::vector<double> *TPC_lastPosX    = nullptr;
    std::vector<double> *TPC_lastPosY    = nullptr;
    std::vector<double> *TPC_lastPosZ    = nullptr;
    std::vector<double> *TPC_PathLength  = nullptr;
    std::vector<double> *TPC_dEdx        = nullptr;
    std::vector<double> *TPC_TrueKE      = nullptr;
    std::vector<double> *TPC_pdg         = nullptr;
    std::vector<int>    *TPC_nSteps      = nullptr;

    void SetBranches(TTree* t) {
        // Calorimeter
        SafeSetBranch(t, "centerX", centerXs);
        SafeSetBranch(t, "centerY", centerYs);
        SafeSetBranch(t, "centerZ", centerZs);
        SafeSetBranch(t, "energy",  energies);
        // Primaries
        SafeSetBranch(t, "PrimaryPosX",   primaryX);
        SafeSetBranch(t, "PrimaryPosY",   primaryY);
        SafeSetBranch(t, "PrimaryPosZ",   primaryZ);
        SafeSetBranch(t, "PrimaryEkin",   primaryEkin);
        SafeSetBranch(t, "PrimaryMomX",   primaryPx);
        SafeSetBranch(t, "PrimaryMomY",   primaryPy);
        SafeSetBranch(t, "PrimaryMomZ",   primaryPz);
        SafeSetBranch(t, "PrimaryPDG",    primaryPDG);
        SafeSetBranch(t, "PrimaryTrackID",primaryTrackID);
        // True photons
        SafeSetBranch(t, "truthPosX",          truePhotonPosX);
        SafeSetBranch(t, "truthPosY",          truePhotonPosY);
        SafeSetBranch(t, "truthPosZ",          truePhotonPosZ);
        SafeSetBranch(t, "truthE",             truePhotonE);
        SafeSetBranch(t, "TruePhotonTrackID",  truePhotonTrackID);
        SafeSetBranch(t, "TruePhotonParentID", truePhotonParentID);
        // True charged pions
        SafeSetBranch(t, "TrueChargedPionX",               TrueChargedPionX);
        SafeSetBranch(t, "TrueChargedPionY",               TrueChargedPionY);
        SafeSetBranch(t, "TrueChargedPionZ",               TrueChargedPionZ);
        SafeSetBranch(t, "TrueChargedPionE",               TrueChargedPionE);
        SafeSetBranch(t, "TrueChargedPionTrackID",         TrueChargedPionTrackID);
        SafeSetBranch(t, "TrueChargedPionThroughTPC",      TrueChargedPionThroughTPC);
        SafeSetBranch(t, "TrueChargedPionCreationX",       TrueChargedPionCreationX);
        SafeSetBranch(t, "TrueChargedPionCreationY",       TrueChargedPionCreationY);
        SafeSetBranch(t, "TrueChargedPionCreationZ",       TrueChargedPionCreationZ);
        SafeSetBranch(t, "TrueChargedPionEndX",            TrueChargedPionEndX);
        SafeSetBranch(t, "TrueChargedPionEndY",            TrueChargedPionEndY);
        SafeSetBranch(t, "TrueChargedPionEndZ",            TrueChargedPionEndZ);
        SafeSetBranch(t, "TrueChargedPionDecayedBeforeCal",TrueChargedPionDecayedBeforeCal);
        SafeSetBranch(t, "TrueChargedPionDecayedBeforeTPC",TrueChargedPionDecayedBeforeTPC);
        SafeSetBranch(t, "TrueChargedPionDecayedTrackID",  TrueChargedPionDecayedTrackID);
        // TPC
        SafeSetBranch(t, "TPC_trackID",    TPC_trackID);
        SafeSetBranch(t, "TPC_Edep",       TPC_Edep);
        SafeSetBranch(t, "TPC_smearedEdep",TPC_smearedEdep);
        SafeSetBranch(t, "TPC_firstPosX",  TPC_firstPosX);
        SafeSetBranch(t, "TPC_firstPosY",  TPC_firstPosY);
        SafeSetBranch(t, "TPC_firstPosZ",  TPC_firstPosZ);
        SafeSetBranch(t, "TPC_lastPosX",   TPC_lastPosX);
        SafeSetBranch(t, "TPC_lastPosY",   TPC_lastPosY);
        SafeSetBranch(t, "TPC_lastPosZ",   TPC_lastPosZ);
        SafeSetBranch(t, "TPC_PathLength", TPC_PathLength);
        SafeSetBranch(t, "TPC_dEdx",       TPC_dEdx);
        SafeSetBranch(t, "TPC_TrueKE",     TPC_TrueKE);
        SafeSetBranch(t, "TPC_pdg",        TPC_pdg);
        SafeSetBranch(t, "TPC_nSteps",     TPC_nSteps);
    }
};

struct BranchManagerVertex {
    // Raw branch variables
    Long64_t event_id    = -1;
    Long64_t vertex_index = -1;
    Long64_t n_tracks    = 0;
    double   vtx_x       = 0;
    double   vtx_y       = 0;
    double   vtx_z       = 0;
    double   chi2ndf     = 0;
    UChar_t  accepted    = 0;

    // Processed result - best vertex per event
    std::vector<Vtx> bestVtx;

    void SetBranches(TTree* t) {
        // Note: plain variables use SetBranchAddress directly
        // not SafeSetBranch since they're not pointers-to-pointers
        t->SetBranchAddress("event_id",     &event_id);
        t->SetBranchAddress("vertex_index", &vertex_index);
        t->SetBranchAddress("n_tracks",     &n_tracks);
        t->SetBranchAddress("vtx_x",        &vtx_x);
        t->SetBranchAddress("vtx_y",        &vtx_y);
        t->SetBranchAddress("vtx_z",        &vtx_z);
        t->SetBranchAddress("chi2ndf",      &chi2ndf);
        t->SetBranchAddress("accepted",     &accepted);
    }

    // Call this after SetBranches, passing the number of signal events
    // so bestVtx is sized correctly for indexing by event_id
    void LoadVertices(TTree* vtxTree, Long64_t nentries) {
        bestVtx.resize(nentries);

        const Long64_t nv = vtxTree->GetEntries();
        for (Long64_t i = 0; i < nv; ++i) {
            vtxTree->GetEntry(i);
            if (accepted == 0) continue;
            if (event_id < 0 || event_id >= nentries) continue;

            Vtx& cur = bestVtx[event_id];
            const bool better =
                (!cur.has) ||
                (n_tracks > cur.n_tracks) ||
                (n_tracks == cur.n_tracks && chi2ndf < cur.chi2ndf);

            if (better) {
                cur.has      = true;
                // cur.x        = vtx_x;
                // cur.y        = vtx_y;
                // cur.z        = vtx_z;
                cur.vertexVec = TVector3{vtx_x, vtx_y, vtx_z};
                cur.n_tracks = n_tracks;
                cur.chi2ndf  = chi2ndf;
            }
        }

        std::cout << "Loaded accepted vertices for "
                  << std::count_if(bestVtx.begin(), bestVtx.end(),
                                   [](const Vtx& v){ return v.has; })
                  << " / " << nentries << " events\n";
    }

    // // Convenience getter for a given event
    // TVector3 GetVertex(Long64_t ievt, const TVector3& fallback) const {
    //     if (ievt >= 0 && ievt < (Long64_t)bestVtx.size() 
    //         && bestVtx[ievt].has) {
    //         return TVector3(bestVtx[ievt].x, 
    //                        bestVtx[ievt].y, 
    //                        bestVtx[ievt].z);
    //     }
    //     return fallback;
    // }
    // Convenience getter for a given event
    Vtx GetVertex(Long64_t ievt, const TVector3& fallback) const {
        if (ievt >= 0 && ievt < (Long64_t)bestVtx.size() 
            && bestVtx[ievt].has) {
            return bestVtx[ievt];
        }
        else {
            Vtx vertex;
            vertex.vertexVec = fallback;
            return vertex;
        };
    }
};

#endif