#ifndef SELECTIONCUTS_H
#define SELECTIONCUTS_H

#include "EventVariables.hpp"

struct SelectionCuts {
    // for now some placeholders
    int minChargedTracks      = 2.0;
    double minTotalEnergy     = 300.0; //--> MeV
    double minCorrectedEnergy = 300.0; 
    int minNeutralClusters    = 2;
    double minSphericity      = 0.3;   // we need to establish if this, since we use proxy momentum is okay! 
    double maxVertexRadius    = 5.0;   //--> cm
    int minPi0Candidates      = 1;     // could also have a min pi+- candidates etc...
};


bool PassesSelection(const EventVariables& ev,
                     const SelectionCuts& cuts)
{
    if (ev.nChargedTracks   < cuts.minChargedTracks) return false;
    if (ev.totalRecoEnergy  < cuts.minTotalEnergy) return false;
    if (ev.correctedEnergy  < cuts.minCorrectedEnergy) return false;
    if (ev.nNeutralClusters < cuts.minNeutralClusters) return false;
    if (ev.sphericity       < cuts.minSphericity) return false;
    if (ev.vertexRadius     < cuts.maxVertexRadius) return false;
    if (ev.nPi0Candidates   < cuts.minPi0Candidates) return false;
    return true;
}

#endif