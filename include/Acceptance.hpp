// #pragma once

// #include "TTree.h"
// #include "TFile.h"
// #include "TH1D.h"
// #include "TGraphAsymmErrors.h"
// #include "TLorentzVector.h"
// #include "TPaveText.h"
// #include <vector>
// #include <string>

// #include "Structures.hpp"   // Cluster, TruePhoton
// #include "Clustering.hpp"
// #include "TruePhotonCalc.hpp"
// #include "PhotonMatch.hpp"
// #include "PlotUtils.hpp"

// class Acceptance {
// public:
//     Acceptance(const std::string& tag, int nBins = 20, double xMin = 0, double xMax = 500);

//     void Pi0ProcessSignalEvent(const std::vector<TruePi0>& truePi0s, 
//                                const std::vector<primaryPi0>& primaryPi0s);

//     void ChPiProcessSignalEvent(const std::vector<TrueChPiInCal>& ChPiInCal, 
//                                 const std::vector<primaryChPi>& primaryChPis,
//                                 int type /* 0 = global, 1 = cond TPC, 2 = cond NO TPC*/);
    

//     void FinalizePlot(const std::string& outFileName, PlotOptions opts);

// private:
//     double massLow_, massHigh_;
//     double pi0Mass_;
//     TH1D* h_num_;
//     TH1D* h_den_;
//     TH1D* h_num_1;
//     TH1D* h_den_1;
//     TH1D* h_num_2;
//     TH1D* h_den_2;
//     int nBins_;
//     double xMin_, xMax_;
//     int countTot;

// };


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

enum class AccAxisType { kEkin, kTheta, k2D, kTracks};

class Acceptance {
public:
    // Constructor for 1D plots (ekin or theta)
    Acceptance(const std::string& tag,
               AccAxisType axisType,
               int nBins = 20, double xMin = 0.0, double xMax = 500.0);

    // Constructor for 2D plot (theta vs ekin)
    Acceptance(const std::string& tag,
               int nBinsEkin, double ekinMin, double ekinMax,
               int nBinsTheta, double thetaMin, double thetaMax);

    void Pi0ProcessSignalEvent(const std::vector<TruePi0>& truePi0s,
                               const std::vector<primaryPi0>& primaryPi0s);

    void ChPiProcessSignalEvent(const std::vector<TrueChPiInCal>& ChPiInCal,
                                const std::vector<primaryChPi>& primaryChPis,
                                int type /* 0 = global, 1 = cond TPC, 2 = cond NO TPC */);

    void FinalizePlot(const std::string& outFileName, PlotOptions opts);

private:
    std::string tag_;
    AccAxisType axisType_;

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