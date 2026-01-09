#ifndef PARTICLEID_H
#define PARTICLEID_H

#include <iostream>
#include <math.h>

double nSigmaCalc(
    double EdepSmeared, 
    double pathLength,
    double dEdxTheory,
    double resolution
);

enum class PotentialGas { eArCO2_9010, eArCO2_8020, eAir };

inline double PDGMassMeV(int pdg)
{
    switch (std::abs(pdg)) {
        case 11:   return 0.5109989461;   // e+-
        case 211:  return 139.57039;      // pi+-
        case 2212: return 938.27208816;   // p
        default:
            std::cerr << "Unknown PDG code " << pdg << std::endl;
            return -1.0;
    }
}

double BetheBloch(
    int pdg,
    double kineticE_MeV,
    PotentialGas gas
);


#endif