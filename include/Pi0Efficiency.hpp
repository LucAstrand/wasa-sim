#pragma once

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include <vector>
#include <string>

#include "Structures.hpp"   // Cluster, TruePhoton
#include "Clustering.hpp"
#include "TruePhotonCalc.hpp"
#include "PhotonMatch.hpp"

class Pi0Efficiency {
public:
    Pi0Efficiency(double massLow = 120.0, double massHigh = 150.0, 
                  double pi0Mass = 134.977, int nBins = 20, double eMin = 0, double eMax = 500);

    // void ProcessEvent(const std::vector<Cluster>& clusters, const std::vector<TruePhoton>& truePhotons);
    void ProcessEvent(const std::vector<TruePi0>& tpi0, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters);
    void FinalizePlot(const std::string& outFileName = "Pi0Efficiency.png");

private:
    double massLow_, massHigh_;
    double pi0Mass_;
    TH1D* h_num_;
    TH1D* h_den_;
    int nBins_;
    double eMin_, eMax_;

    // double ComputeTruePi0Ekin(const std::vector<TruePhoton>& truePhotons);
    // int CountTruePi0s(const std::vector<TruePhoton>& truePhotons);
    // int CountRecoPi0s(const std::vector<Cluster>& clusters);
    bool IsReco(const TruePi0& tpi0, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters);
};
