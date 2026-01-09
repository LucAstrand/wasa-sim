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

double BetheBloch(
    int pdg,
    double kineticE,        // MeV
    PotentialGas gas
)
{
    // ---- Physical constants (MeV, cm) ----
    constexpr double me = 0.5109989461;   // electron mass [MeV]
    constexpr double K  = 0.307075;       // MeV cm^2 / mol
    constexpr double z  = 1.0;            // charge (assumed ±1)

    // ---- Gas parameters (defaults: Ar/CO2 90/10) ----
    double I = 1.77e-4;     // Mean excitation energy [MeV]
    double ZA = 0.455;      // Z/A
    double rho = 1.70e-3;   // g/cm^3

    switch (gas) {
        case PotentialGas::eArCO2_9010:
            I   = 1.77e-4;
            ZA  = 0.455;
            rho = 1.70e-3;
            break;

        case PotentialGas::eArCO2_8020:
            I   = 1.674e-4;
            ZA  = 0.46;
            rho = 1.773e-3;
            break;

        case PotentialGas::eAir:
            I   = 8.57e-5;
            ZA  = 0.499;
            rho = 1.225e-3;
            break;
    }

    // ---- Particle mass ----
    double m = PDGMassMeV(pdg);
    if (m <= 0.0) return 0.0;

    // ---- Kinematics ----
    double E  = kineticE + m;
    double gamma = E / m;
    double beta2 = 1.0 - 1.0 / (gamma * gamma);
    if (beta2 <= 0.0) return 0.0;

    double beta = std::sqrt(beta2);

    // ---- Tmax (max energy transfer) ----
    double massRatio = me / m;
    double Tmax =
        (2.0 * me * beta2 * gamma * gamma) /
        (1.0 + 2.0 * gamma * massRatio + massRatio * massRatio);

    // ---- Bethe-Bloch ----
    double argument =
        (2.0 * me * beta2 * gamma * gamma * Tmax) / (I * I);

    if (argument <= 1.0) return 0.0;

    double dEdx =
        K * z * z * ZA * (1.0 / beta2) *
        (0.5 * std::log(argument) - beta2);

    // ---- Apply density ----
    dEdx *= rho;   // MeV / cm

    return dEdx;
}
