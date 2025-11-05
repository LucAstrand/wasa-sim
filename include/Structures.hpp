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
