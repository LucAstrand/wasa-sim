#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include "TH1F.h"

void SetPrettyStyle();
void PrettyPi0MassPlot(TH1F* hPi0Mass, TString plotname);
void TruthPi0MassPlot(TH1F* hPi0Mass, TString plotname);
void PrettyPi0NumClusterPlot(TH1F* hNCluster);
void BasicHistPlot(TH1F* histogram);

#endif
