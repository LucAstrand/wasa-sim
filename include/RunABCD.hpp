#ifndef RUNABCD_H
#define RUNABCD_H

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooFit.h"
#include "RooPlot.h"
#include "RooWorkspace.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "EventVariables.hpp"
#include <array>
#include <string>

struct ABCDResult {
    double N_A_sig   = 0;  // signal events in A
    double N_B_bkg   = 0;
    double N_C_bkg   = 0;
    double N_D_bkg   = 0;
    double N_A_bkg_true     = 0;  // true MC background in A
    double N_A_bkg_estimate = 0;  // ABCD estimate
    double closure          = 0;  // estimate / true
    double closureErr       = 0;
};

ABCDResult RunABCD(
    const std::vector<EventVariables>& sigEvents,
    const std::vector<EventVariables>& bkgEvents,
    // Variable accessors - pass lambdas so any variable pair works
    std::function<double(const EventVariables&)> getVar1,
    std::function<double(const EventVariables&)> getVar2,
    double threshold1,   // signal-like = above this
    double threshold2,   // signal-like = above this
    const std::string& var1Name,
    const std::string& var2Name,
    const std::string& outDir);

#endif