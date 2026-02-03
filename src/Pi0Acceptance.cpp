#include "Pi0Acceptance.hpp"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

#include "PlotUtils.hpp"

#include <iostream>

Pi0Acceptance::Pi0Acceptance(const std::string& tag, int nBins, double xMin, double xMax)
    : nBins_(nBins), xMin_(xMin), xMax_(xMax)
{
    h_num_ = new TH1D(Form("h_numPi0Acceptance_%s", tag.c_str()), ";True #pi^{0} E_{kin};Events", nBins_, xMin_, xMax_);
    h_den_ = new TH1D(Form("h_denPi0Acceptance_%s", tag.c_str()), ";GPS #pi^{0} E_{kin};Events", nBins_, xMin_, xMax_);
}

// Pi0Acceptance::Pi0Acceptance(double etaMin, double etaMax, int nBins)
//     : massLow_(0), massHigh_(0), pi0Mass_(0), // unused here
//       nBins_(nBins), eMin_(etaMin), eMax_(etaMax)
// {
//     h_num_ = new TH1D("h_num_eta", "Reconstructed #pi^{0};True #pi^{0} #eta;Events", nBins_, etaMin, etaMax);
//     h_den_ = new TH1D("h_den_eta", "All #pi^{0};True #pi^{0} #eta;Events", nBins_, etaMin, etaMax);
// }

// Pi0Acceptance::Pi0Acceptance(double thetaMin, double thetaMax, int nBins)
// {
//     h_num_1 = new TH1D("h_num_1", "Reco #pi^{0} (E_{kin}=50 MeV);#theta [rad];Events",
//                        nBins, thetaMin, thetaMax);
//     h_den_1 = new TH1D("h_den_1", "All #pi^{0} (E_{kin}=50 MeV);#theta [rad];Events",
//                        nBins, thetaMin, thetaMax);

//     h_num_2 = new TH1D("h_num_2", "Reco #pi^{0} (E_{kin}=500 MeV);#theta [rad];Events",
//                        nBins, thetaMin, thetaMax);
//     h_den_2 = new TH1D("h_den_2", "All #pi^{0} (E_{kin}=500 MeV);#theta [rad];Events",
//                        nBins, thetaMin, thetaMax);
// }


void Pi0Acceptance::ProcessGPSEvent(const std::vector<TruePi0>& truePi0s, const double xVar, int nPi0s) {

    if (nPi0s != 0) h_den_->Fill(xVar, nPi0s);
    if (truePi0s.size() != 0) h_num_->Fill(xVar, truePi0s.size());

}


void Pi0Acceptance::ProcessSignalEvent(const std::vector<TruePi0>& truePi0s, const double xVar, int nPi0s) {

    if (nPi0s != 0) h_den_->Fill(xVar, nPi0s);
    if (truePi0s.size() != 0) h_num_->Fill(xVar, truePi0s.size());

}

void Pi0Acceptance::FinalizePlot(const std::string& outFileName) {


    TGraphAsymmErrors* gEff = new TGraphAsymmErrors(h_num_, h_den_, "cl=0.683 b(1,1) mode"); // fraction in [0,1]

    // SCALE TO PERCENT 
    // Multiply y-values and Y-errors by 100 to convert fraction -> percent
    const int nPoints = gEff->GetN();
    for (int i = 0; i < nPoints; ++i) {
        double x, y;
        gEff->GetPoint(i, x, y);

        double exl = gEff->GetErrorXlow(i);
        double exh = gEff->GetErrorXhigh(i);
        double eyl = gEff->GetErrorYlow(i);
        double eyh = gEff->GetErrorYhigh(i);

        double y_new = y * 100.0;
        double eyl_new = eyl * 100.0;
        double eyh_new = eyh * 100.0;

        gEff->SetPoint(i, x, y_new);
        gEff->SetPointError(i, exl, exh, eyl_new, eyh_new);
    }

    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(1.0);
    gEff->SetMarkerColor(kBlack);
    gEff->SetLineColor(kBlack);

    gEff->GetYaxis()->SetTitle("Acceptance [%]");
    //vs Ekin
    gEff->GetXaxis()->SetTitle("GPS #pi^{0} E_{kin} [MeV]");
    //vs Eta 
    // gEff->GetXaxis()->SetTitle("GPS #pi^{0} pseudorapidity #eta");
    //vs Theta
    // gEff->GetXaxis()->SetTitle("GPS #pi^{0} #theta rad");


    gEff->GetYaxis()->SetRangeUser(0.0, 110.0);

    gStyle->SetOptStat(0);
    gEff->Draw("AP");

    PlotOptions opts;
    opts.topLatex = "#bf{Hibeam}  #it{Wasa full simulation}";
    opts.legendEntries = { "Acceptance" };
    opts.infoLines = {"GEANT4 #pi^{0} sample"};
    opts.addInfoPave = true;
    PlotGraph(gEff, outFileName.c_str(), opts);

    // cleanup
    delete gEff;
    delete h_num_;
    delete h_den_;
}
