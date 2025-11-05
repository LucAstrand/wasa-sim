#ifndef UTILS_H
#define UTILS_H

#include "TVector3.h"
#include "TLorentzVector.h"
#include <cmath>

inline double calcEta(double x, double y, double z) {
    return TVector3(x, y, z).Eta();
}

inline double calcPhi(double x, double y, double z) {
    return TVector3(x, y, z).Phi();
}

inline double dphi_wrap(double a, double b) {
    double d = a - b;
    while (d <= -M_PI) d += 2*M_PI;
    while (d >  M_PI)  d -= 2*M_PI;
    return d;
}

inline TLorentzVector makePhotonFromEnergyAndDir(double E, const TVector3 &dir) {
        TVector3 u = dir;
        if (u.Mag() == 0) return TLorentzVector(0,0,0,0);
        u.SetMag(1.0);
        TVector3 pvec = u * E;
        TLorentzVector p4;
        p4.SetPxPyPzE(pvec.X(), pvec.Y(), pvec.Z(), E);
        return p4;
}

#endif
