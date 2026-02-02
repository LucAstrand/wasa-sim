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
    const std::vector<TruePhotonHit> &hits,
    const TVector3 &vertex);

std::vector<TruePi0> TruePi0Builder(
    const std::vector<TruePhoton>& photons
);

#endif