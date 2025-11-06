#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include "TH1F.h"
#include "TGraph.h"
#include "TEfficiency.h"

void SetPrettyStyle();
void PrettyPi0MassPlot(TH1F* hPi0Mass, TString plotname, double fitMin, double fitMax);
void TruthPi0MassPlot(TH1F* hPi0Mass, TString plotname);
void PrettyPi0NumClusterPlot(TH1F* hNCluster);
void BasicHistPlot(TH1F* histogram);
void EffPlot(TH1F* hEff, TString plotname);



#endif
