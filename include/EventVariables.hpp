#ifndef EVENTVARIABLES_H
#define EVENTVARIABLES_H

#include "Structures.hpp"
#include "RecoEvent.hpp"
#include "Calibration.hpp"

#include "TMatrixDSym.h"
#include "TVectorD.h"

EventVariables ComputeEventVariables(
    const int& eventNumber, 
    const RecoEvent& reco,
    const TVector3& vertex, 
    const ChargedKECalibration& calibration
);

#endif