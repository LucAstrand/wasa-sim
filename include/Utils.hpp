#ifndef UTILS_H
#define UTILS_H

#include "TVector3.h"
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

#endif
