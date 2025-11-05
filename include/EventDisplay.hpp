#ifndef EVENTDISPLAY_H
#define EVENTDISPLAY_H

#include <vector>
#include "Structures.hpp"

void PlotEventDisplay(
    const std::vector<Hit>& hits,
    const std::vector<Cluster>& clusters,
    double dEta, double dPhi,
    int eventID);

#endif
