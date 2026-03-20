#ifndef EVENTVARIABLES_H
#define EVENTVARIABLES_H

#include "Structures.hpp"
#include "RecoEvent.hpp"
#include "TMatrixDSym.h"
#include "TVectorD.h"

EventVariables ComputeEventVariables(
    const RecoEvent& reco,
    const TVector3 vertex
);

#endif