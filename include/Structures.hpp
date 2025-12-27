#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <vector>
#include "TLorentzVector.h"
#include "TVector3.h"

struct Hit {
    double x, y, z, e;
};

struct Cluster {
    std::vector<int> hitIndices;
    TVector3 centroid;
    TLorentzVector p4;
};

struct ChargedTrack {
    size_t id;                 // track ID // I might just build this dynamically not sure yet If I want to use GEANT4's TrackID --> Avoid Truth level stuff?
    TVector3 vertex;
    TVector3 exitPoint;     // TPC exit
    TVector3 direction;     // cached (exit - vertex).Unit()
    double EdepSmeared;
    double pathLength;
    double dEdxTheory;
    double resolution;
};

struct ChargedCluster {
    int trackID;
    std::vector<Hit> hits;
    double totalEnergy = 0.0;
    TVector3 direction;
    double clusterdEdx;
    double nSigma;
};

struct TruePhoton {
    TVector3 dir;
    TLorentzVector p4;
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
