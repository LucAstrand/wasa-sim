#pragma once

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>

#include "Structures.hpp" 
#include "Clustering.hpp"
#include "TruePhotonCalc.hpp"
#include "PhotonMatch.hpp"

class Pi0Efficiency {
public:
    Pi0Efficiency(int nBins = 20, double xMin = 0, double xMax = 500);
    void ProcessEvent(const std::vector<TruePi0>& tpi0, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters);
    void FinalizePlot(const std::string& outFileName = "Pi0Efficiency.png");

private:
    double pi0Mass_;
    TH1D* h_num_;
    TH1D* h_den_;
    int nBins_;
    double xMin_, xMax_;

    bool IsReco(const TruePi0& tpi0, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters);
};
