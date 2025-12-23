#include "ParticleID.hpp"

double nSigmaCalc(
    double EdepSmeared, 
    double pathLength,
    double dEdxTheory,
    double resolution) 
{
    double dEdx_smeared = EdepSmeared / pathLength;
    double sigma = resolution * dEdxTheory;
    double nSigma = (dEdx_smeared - dEdxTheory) / sigma;
    return nSigma; 
} 