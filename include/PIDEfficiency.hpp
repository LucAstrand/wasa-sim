#ifndef PID_EFFICIENCY_HPP
#define PID_EFFICIENCY_HPP

#include <vector>
#include <string>
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "Structures.hpp" // for ChargedCluster

class PIDEfficiency {
public:
    PIDEfficiency(int nBins, double eMin, double eMax);

    void ProcessEvent(const std::vector<ChargedCluster>& clusters);

    void FinalizePlot(const std::string& outFileName, int pdgNumToPlot);

private:
    int nBins_;
    double eMin_, eMax_;

    TH1D* h_numPion_;
    TH1D* h_numPion_miss;
    TH1D* h_denPion_;

    TH1D* h_numProton_;
    TH1D* h_numProton_miss;
    TH1D* h_denProton_;

    // TH1D* h_numElectron_;
    // TH1D* h_numElectron_miss;
    // TH1D* h_denElectron_;
};

#endif