#pragma once

#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include <vector>
#include <string>

#include "Structures.hpp"   // Cluster, TruePhoton
#include "Clustering.hpp"
#include "TruePhotonCalc.hpp"
#include "PhotonMatch.hpp"

class Pi0Acceptance {
public:
    Pi0Acceptance(const std::string& tag, int nBins = 20, double xMin = 0, double xMax = 500);
    
    void ProcessGPSEvent(const std::vector<TruePi0>& truePi0s, 
                        const double xVar, // whatever vairable one chooses Eg Ekin, theta, eta...
                        int nPi0s = 1 // For the GPS generally we shoot one particle per event. 
                        );

    void ProcessSignalEvent(const std::vector<TruePi0>& truePi0s, 
                      const double xVar, // whatever vairable one chooses Eg Ekin, theta, eta...
                      int nPi0s);
    
    // void ProcessEventTwoHist(const std::vector<Cluster>& clusters,
    //                          const std::vector<TruePhoton>& TruePhotons,
    //                          const double primaryEkin,
    //                          const double primaryTheta);

    void FinalizePlot(const std::string& outFileName = "Pi0Efficiency.png");

private:
    double massLow_, massHigh_;
    double pi0Mass_;
    TH1D* h_num_;
    TH1D* h_den_;
    TH1D* h_num_1;
    TH1D* h_den_1;
    TH1D* h_num_2;
    TH1D* h_den_2;
    int nBins_;
    double xMin_, xMax_;
    int countTot;

};
