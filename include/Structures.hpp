#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

#include "ParticleID.hpp"

// struct Hit {
//     double x, y, z, e;
// };

struct Vtx {
    bool has = false;
    double x=0, y=0, z=0;
    int n_tracks = 0;
    double chi2ndf = 1e99;
};

enum class HitOwner {
    None,
    Charged,
    Neutral
};

struct Hit {
    double x, y, z;
    double e;
    HitOwner owner = HitOwner::None;
};

struct TruePhotonHit {
    double x, y, z;
    double e;
    double truePhotonTrackID;
    double truePhotonParentID;
    HitOwner owner = HitOwner::None;
};

struct TruePhoton {
    TVector3 dir;
    TLorentzVector p4;
    double parentID;
};

struct TruePi0 {
    int trackID;
    TLorentzVector p4;
    std::vector<const TruePhoton*> photons;
};

struct TrueChPiInCal {
    int trackID;
    bool throughTPC;
};

struct TrueChPiDecayed {
    int trackID;
    bool beforeCal;
    bool beforeTPC;
};

struct primaryPi0 {
    int trackID;
    TLorentzVector p4;
};

struct primaryChPi {
    int trackID;
    TLorentzVector p4; 
};

struct Cluster {
    std::vector<Hit*> hits; // Pointer now! use -> to access the quantities!
    TVector3 centroid;
    TLorentzVector p4;
    bool isFromConversion = false;
};

struct Pi0Candidate {
    const Cluster* c1;
    const Cluster* c2;
    double mgg;
    double theta;
    double score;
    TLorentzVector p4;
};


struct RecoPi0 {
    const Cluster* c1;
    const Cluster* c2;
    TLorentzVector p4;
    double mass;
};

struct ChargedTrack {
    double id;              // track ID // I might just build this dynamically not sure yet If I want to use GEANT4's TrackID --> Avoid Truth level stuff?
    TVector3 vertex;
    TVector3 exitPoint;     // TPC exit
    TVector3 direction;     // cached (exit - entry).Unit()
    double TrueKE;
    double TruePDG;
    double clusterdEdx;
    double EdepSmeared;
    double pathLength;
    double dEdxTheory;
    double resolution;
    int nSteps; 
};

struct ChargedCluster {
    int trackID;
    std::vector<Hit*> hits; // Pointer now! use -> to access the quantities!
    double totalEnergy = 0.0;
    TVector3 direction;
    TVector3 TPCExitPoint;
    int nSteps;
    double objectTrueKE; // true info from sim for control plots 
    double objectTruedEdx; // true info from sim for control plots
    int objectTruePDG;
    double clusterdEdx;
    double EdepSmeared;
    double nSigmaPion;
    double nSigmaProton;
    double nSigmaElectron;
    PIDLikelihoods pidL;
    PID pidGuess;
    bool isOrphanElectron = false;
    bool isUsedInConversion = false;
};

struct ChargedObject {
    int trackID;
    TLorentzVector p4; 
    const ChargedCluster* c; // For any reconstruction detail
};

struct ConversionCandidate {
    TVector3 conversionVertex;  // where the pair originated
    TLorentzVector p4;          // reconstructed photon 4-vector
    int track1_idx;
    int track2_idx;
    double openingAngle;
    double invMass;
};

struct EtaPhiTowerKey {
    int iEta;
    int iPhi;
    bool operator<(const EtaPhiTowerKey &o) const {
        if (iEta != o.iEta) return iEta < o.iEta;
        return iPhi < o.iPhi;
    }
};

struct EventVariables {
    // Multiplicity
    int nChargedTracks       = 0;   // TPC track count (excl electrons)
    int nNeutralClusters     = 0;   // neutral cluster count
    int nTotalObjects        = 0;   // nCharged + nNeutral
    
    // Energy
    double EM_energy         = 0.0; // total calorimeter energy
    double chargedEnergy     = 0.0; // sum of charged cluster energies
    double totalRecoEnergy   = 0.0; // EM_energy //not this + charged TPC edep
    double correctedEnergy   = 0.0; // EM_energy + calibratedKE (Calibrated missing EKin)
    
    // Topology
    double sphericity        = 0.0; // 0=pencil-like, 1=isotropic
    double maxTrackAngle     = 0.0; // largest opening angle between any two tracks
    double vertexRadius      = 0.0; // transverse distance of vertex from beam axis
    double meanDCA           = 0.0; // mean distance of closest approach of tracks to vertex
    
    // PID based
    int nPionCandidates      = 0;   // tracks with pidGuess == Pion
    int nProtonCandidates    = 0;
    int nPi0Candidates       = 0;   // reconstructed pi0s
};

#endif
