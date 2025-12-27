#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TEfficiency.h"

void SetPrettyStyle();
void PrettyPi0MassPlot(TH1F* hPi0Mass, TString plotname, double fitMin, double fitMax);
void TruthPi0MassPlot(TH1F* hPi0Mass, TString plotname);
void PrettyPi0NumClusterPlot(TH1F* hNCluster);
// void Pi0ClusterNumPlotEkin(TH1F* hNClusters_lowE, TH1F* hNClusters_midE, TH1F* hNClusters_highE);
void Pi0ClusterNumPlotEkin(TH1F* hNClusters_lowE, TH1F* hNClusters_highE);
void BasicHistPlot(TH1F* histogram);
void nSigmaPlot(TH1F* hNSigma, TString plotname, double fitMin, double fitMax);
void dEdxVsEPlot(TH2F* hdEdxVsE, TString plotname);



#endif
