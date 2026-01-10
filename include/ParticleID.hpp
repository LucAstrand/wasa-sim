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

inline double LikelihoodFromNSigma(double nSigma);

struct PIDLikelihoods {
    double Lpi;
    double Lp;
    double Le;

    double Ppi;
    double Pp;
    double Pe;
};

PIDLikelihoods ComputePIDLikelihoods(
    double nSigmaPi,
    double nSigmaP,
    double nSigmaE
);

enum class PID { Pion, Proton, Electron, Unknown };

inline const char* PIDToString(PID pid)
{
    switch(pid) {
        case PID::Pion:     return "Pion";
        case PID::Proton:   return "Proton";
        case PID::Electron: return "Electron";
        case PID::Unknown:  return "Unknown";
        default:            return "InvalidPID";
    }
}

PID AssignPIDFromLikelihood(
    const PIDLikelihoods& L,
    double minProb = 0.7
);


#endif