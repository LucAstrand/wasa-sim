#ifndef UTILS_H
#define UTILS_H

#include "Structures.hpp"

#include "TVector3.h"
#include "TLorentzVector.h"
#include <cmath>

inline TVector3 hitDirection(const Hit& h, const TVector3& vtx) {
    return TVector3(h.x, h.y, h.z) - vtx;
}

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

inline double openingAngle(const TLorentzVector& g1, const TLorentzVector& g2) {
    return g1.Vect().Angle(g2.Vect());
}

inline double expectedAngle(double E1, double E2) {
    constexpr double mpi0 = 134.9768; // GeV
    constexpr double sigmaM = 30;    // example: 10 MeV mass resolution
    constexpr double sigmaTheta = 0.5;//0.01; // example: ~10 mrad angular resolution

    double cosTheta = 1.0 - (mpi0 * mpi0) / (2.0 * E1 * E2);
    //safety
    if (cosTheta > 1.0) cosTheta = 1.0;
    if (cosTheta > -1.0) cosTheta = -1.0;

    return std::acos(cosTheta);
}

// inline double chi2Pi0(const TLorentzVector& g1, const TLorentzVector& g2) {
//     constexpr double mpi0 = 134.9768; // GeV
//     constexpr double sigmaM = 30;    // example: 10 MeV mass resolution
//     constexpr double sigmaTheta = 0.5;//0.01; // example: ~10 mrad angular resolution
//     TLorentzVector pi0 = g1 + g2;
//     double mgg = pi0.M();
//     double theta = openingAngle(g1,g2);
//     double expectedtheta = expectedAngle(g1.E(),g2.E());
//     double chi2m = (mgg - mpi0) * (mgg - mpi0) / (sigmaM * sigmaM);
//     double chi2t = (theta - expectedtheta) * (theta - expectedtheta) / (sigmaTheta * sigmaTheta);

//     std::cout << "chi2Pi0: "
//             << "E1=" << g1.E()
//             << " E2=" << g2.E()
//             << " mgg=" << (g1+g2).M()
//             << " theta=" << theta
//             << " thetaExp=" << expectedtheta
//             << " chi2=" << chi2m + chi2t
//             << std::endl;
//     return chi2m + chi2t;
// }

#endif
