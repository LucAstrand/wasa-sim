#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

#include "ParticleID.hpp"

// struct Hit {
//     double x, y, z, e;
// };

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
    int truePhotonTrackID;
    int truePhotonParentID;
    HitOwner owner = HitOwner::None;
};

struct Cluster {
    std::vector<Hit*> hits; // Pointer now! use -> to access the quantities!
    TVector3 centroid;
    TLorentzVector p4;
};

// struct Pi0Candidate {
//     size_t i;
//     size_t j;
//     double mgg;
//     double theta;
//     double score;
//     TLorentzVector p4;
// };

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
    size_t id;              // track ID // I might just build this dynamically not sure yet If I want to use GEANT4's TrackID --> Avoid Truth level stuff?
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
};

struct ChargedCluster {
    int trackID;
    std::vector<Hit*> hits; // Pointer now! use -> to access the quantities!
    double totalEnergy = 0.0;
    TVector3 direction;
    double objectTrueKE; // true info from sim for control plots 
    double objectTruedEdx; // true info from sim for control plots
    int objectTruePDG;
    double clusterdEdx;
    double nSigmaPion;
    double nSigmaProton;
    double nSigmaElectron;
    PIDLikelihoods pidL;
    PID pidGuess;
};

struct TruePhoton {
    TVector3 dir;
    TLorentzVector p4;
    // int parentPDG;
    int parentID;
};

struct TruePi0 {
    int trackID;
    TLorentzVector p4;
    std::vector<const TruePhoton*> photons;
};

struct EtaPhiTowerKey {
    int iEta;
    int iPhi;
    bool operator<(const EtaPhiTowerKey &o) const {
        if (iEta != o.iEta) return iEta < o.iEta;
        return iPhi < o.iPhi;
    }
};

#endif
