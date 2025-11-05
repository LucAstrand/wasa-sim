#ifndef TRUEPHOTONCALC_H
#define TRUEPHOTONCALC_H

#include <vector>
#include <map>
#include <set>

#include "TVector3.h"
#include "TLorentzVector.h"

#include "Structures.hpp"
#include "Utils.hpp"

std::vector<TruePhoton> TruePhotonBuilder(
    const std::vector<Hit> &hits,
    const TVector3 &vertex);

#endif