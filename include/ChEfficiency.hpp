#pragma once

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include <vector>
#include <string>

#include "Structures.hpp"
#include "Clustering.hpp"
#include "TruePhotonCalc.hpp"
#include "PhotonMatch.hpp"
#include "PlotUtils.hpp"

enum class AccAxisTypeChPi { kEkin, kTheta, k2D, kTracks};

class ChEfficiency {
public:
    // Constructor for 1D plots (ekin or theta)
    ChEfficiency(const std::string& tag,
               AccAxisTypeChPi axisType,
               int nBins = 20, double xMin = 0.0, double xMax = 500.0);

    // Constructor for 2D plot (theta vs ekin)
    ChEfficiency(const std::string& tag,
               int nBinsEkin, double ekinMin, double ekinMax,
               int nBinsTheta, double thetaMin, double thetaMax);

    void ChPiProcessSignalEvent(const std::vector<TrueChPiInCal>& ChPiInCal,
                                const std::vector<primaryChPi>& primaryChPis,
                                const std::vector<ChargedCluster>& ChClusters);

    void FinalizePlot(const std::string& outFileName, PlotOptions opts);

private:
    std::string tag_;
    AccAxisTypeChPi axisType_;

    // 1D histograms
    TH1D* h_num_  = nullptr;
    TH1D* h_den_  = nullptr;

    // 2D histograms (theta vs ekin)
    TH2D* h2_num_ = nullptr;
    TH2D* h2_den_ = nullptr;

    // 1D binning
    int    nBins_ = 20;
    double xMin_  = 0.0;
    double xMax_  = 500.0;

    // 2D binning
    int    nBinsEkin_  = 20;
    double ekinMin_    = 0.0;
    double ekinMax_    = 500.0;
    int    nBinsTheta_ = 20;
    double thetaMin_   = 0.0;
    double thetaMax_   = 3.2;
};