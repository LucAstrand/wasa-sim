#include "ParticleID.hpp"

double nSigmaCalc(
    double dEdxSmeared, 
    double dEdxTheory,
    double resolution
)
{
    // double scaleFactor = 0.7; // this is arbitrary and should eb tuned... Its because we do the truncated mean... 
    double sigma = resolution * dEdxTheory; //* scaleFactor;
    // double nSigma = (dEdxSmeared - dEdxTheory * scaleFactor) / sigma;
    double nSigma = (dEdxSmeared - dEdxTheory) / sigma;
    return nSigma; 
} 

inline double LikelihoodFromNSigma(double nSigma) {

    return std::exp(-0.5 * nSigma * nSigma);

}

PIDLikelihoods ComputePIDLikelihoods(double nSigmaPi,
                                     double nSigmaE,
                                     double nSigmaP,
                                     double nSigmaM) 

// PIDLikelihoods ComputePIDLikelihoods(double nSigmaPi,
//                                      double nSigmaP)
{
    PIDLikelihoods out;

    out.Lpi = LikelihoodFromNSigma(nSigmaPi);
    out.Lp  = LikelihoodFromNSigma(nSigmaP);
    out.Le  = LikelihoodFromNSigma(nSigmaE);
    out.Lm  = LikelihoodFromNSigma(nSigmaM);

    double sum = out.Lpi + out.Lp + out.Le;
    // double sum = out.Lpi + out.Lp;

    if (sum > 0) {
        out.Ppi = out.Lpi / sum;
        out.Pp  = out.Lp  / sum;
        out.Pe  = out.Le  / sum;
    } else {
        out.Ppi = out.Pp = out.Pe = 0.0;
        // out.Ppi = out.Pp = 0.0;
    }

    return out;
}

PID AssignPIDFromLikelihood(const PIDLikelihoods& L,
                            double minProb)
{
    // // All species:
    // if (L.Ppi > minProb && L.Ppi >= L.Pp && L.Ppi >= L.Pe && L.Ppi >= L.Pm)
    //     return PID::Pion;

    // else if (L.Pp > minProb && L.Pp >= L.Ppi && L.Pp >= L.Pe && L.Pp >= L.Pm)
    //     return PID::Proton;

    // else if (L.Pe > minProb && L.Pe >= L.Ppi && L.Pe >= L.Pp && L.Pe >= L.Pm)
    //     return PID::Electron;

    // else if (L.Pm > minProb && L.Pm >= L.Ppi && L.Pm >= L.Pp && L.Pm >= L.Pp)
    //     return PID::Muon;

    // // Only Pions, Protons and Electrons
    if (L.Ppi > minProb && L.Ppi >= L.Pp && L.Ppi >= L.Pe)
        return PID::Pion;

    else if (L.Pp > minProb && L.Pp >= L.Ppi && L.Pp >= L.Pe)
        return PID::Proton;

    else if (L.Pe > minProb && L.Pe >= L.Ppi && L.Pe >= L.Pp)
        return PID::Electron;

    // // Only Pions and Protons
    // if (L.Ppi > minProb && L.Ppi >= L.Pp)
    //     return PID::Pion;

    // else if (L.Pp > minProb && L.Pp >= L.Ppi)
    //     return PID::Proton;

    return PID::Unknown;
}