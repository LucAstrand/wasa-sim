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
// Implemented exaclty like the sim code, It might not be 100% correct  
double BetheBloch(
    int pdg,
    double kineticE,        // MeV
    PotentialGas gas
)
{
    // --- constants 
    const double c = 1.0;
    const double z = 1.0;
    const double Me = 0.511;           // MeV
    const double cBetheBloch = 0.307075;

    // --- gas defaults (Ar/CO2 90/10) ---
    double ExPotencialGas = 1.77e-4;
    double EffZA          = 0.455;
    double GasDensity     = 1.70e-3;

    switch (gas) {
        case PotentialGas::eArCO2_9010:
            ExPotencialGas = 1.77e-4;
            EffZA          = 0.455;
            GasDensity     = 1.70e-3;
            break;

        case PotentialGas::eArCO2_8020:
            ExPotencialGas = 1.674e-4;
            EffZA          = 0.46;
            GasDensity     = 1.773e-3;
            break;

        case PotentialGas::eAir:
            ExPotencialGas = 1.0;
            EffZA          = 1.0;
            GasDensity     = 1.0;
            break;
    }

    // --- particle mass ---
    double MassParticle = PDGMassMeV(pdg);
    if (MassParticle <= 0.0) return 0.0;

    // --- kinematics ---
    double TotalE  = MassParticle + kineticE;
    double gamma   = TotalE / MassParticle;
    double beta2   = 1.0 - 1.0 / (gamma * gamma);
    if (beta2 <= 0.0) return 0.0;

    double beta = std::sqrt(beta2);

    // --- Tmax 
    double Tmax =
        (2.0 * Me * c * c * beta2 * gamma * gamma) /
        (1.0 + 2.0 * gamma * (Me / MassParticle)
             + (Me / MassParticle) * (Me / MassParticle));

    // --- Bethe-Bloch
    double dEdx =
        cBetheBloch * z * z * (EffZA / beta2) *
        (std::log(
             (2.0 * Me * c * c * beta2 * gamma * gamma * Tmax)
             / (ExPotencialGas * ExPotencialGas)
         ) - 2.0 * beta2);

    // --- density ---
    dEdx *= GasDensity;   // MeV/cm

    return dEdx;
}


inline double LikelihoodFromNSigma(double nSigma) {

    return std::exp(-0.5 * nSigma * nSigma);

}

// PIDLikelihoods ComputePIDLikelihoods(double nSigmaPi,
//                                      double nSigmaP,
//                                      double nSigmaE)

PIDLikelihoods ComputePIDLikelihoods(double nSigmaPi,
                                     double nSigmaP)
{
    PIDLikelihoods out;

    out.Lpi = LikelihoodFromNSigma(nSigmaPi);
    out.Lp  = LikelihoodFromNSigma(nSigmaP);
    // out.Le  = LikelihoodFromNSigma(nSigmaE);

    // double sum = out.Lpi + out.Lp + out.Le;
    double sum = out.Lpi + out.Lp;

    if (sum > 0) {
        out.Ppi = out.Lpi / sum;
        out.Pp  = out.Lp  / sum;
        // out.Pe  = out.Le  / sum;
    } else {
        // out.Ppi = out.Pp = out.Pe = 0.0;
        out.Ppi = out.Pp = 0.0;
    }

    return out;
}

PID AssignPIDFromLikelihood(const PIDLikelihoods& L,
                            double minProb)
{
    // if (L.Ppi > minProb && L.Ppi >= L.Pp && L.Ppi >= L.Pe)
    //     return PID::Pion;

    // if (L.Pp > minProb && L.Pp >= L.Ppi && L.Pp >= L.Pe)
    //     return PID::Proton;

    // if (L.Pe > minProb && L.Pe >= L.Ppi && L.Pe >= L.Pp)
    //     return PID::Electron;
    if (L.Ppi > minProb && L.Ppi >= L.Pp)
        return PID::Pion;

    if (L.Pp > minProb && L.Pp >= L.Ppi)
        return PID::Proton;

    return PID::Unknown;
}