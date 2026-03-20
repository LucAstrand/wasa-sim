#ifndef EVENTLOOP_H
#define EVENTLOOP_H

#include "progressbar.hpp"

#include "TTree.h"

#include "BranchManager.hpp"
#include "DEDXTable.hpp"
#include "Calibration.hpp"
// analysis config deal with this
#include "AnalysisHistograms.hpp"
#include "SelectionHistograms.hpp"

struct AnalysisConfig {
    bool doPi0Analysis     = false;
    bool doChargedAnalysis = false;
    bool doTruthAnalysis   = false;
    bool doTruthAndMix     = false;
    bool doEventVariables  = false;
    bool isSignal          = true;   // false = background/cosmic
};

void RunSignalLoop(
    TTree* tree,
    BranchManagerVertex& brVtx,
    const DEDXTable& dedxTable,
    const ChargedKECalibration& calibration,
    const AnalysisConfig& cfg,
    Pi0Histograms* hPi0,
    TruthHistograms* hTruth,
    ChargedHistograms* hCharged,
    EventVarHistograms* hEvt,
    SelectionHistograms* hSel
);

void RunBackgroundLoop(
    TTree* tree,
    const DEDXTable& dedxTable,
    SelectionHistograms& hSel
);

#endif