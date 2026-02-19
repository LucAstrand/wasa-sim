#include "Pi0Efficiency.hpp"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

#include "PlotUtils.hpp"

#include <iostream>

Pi0Efficiency::Pi0Efficiency(int nBins, double xMin, double xMax)
    : nBins_(nBins), xMin_(xMin), xMax_(xMax)
{
    h_num_ = new TH1D("h_numPi0Efficiency", ";True #pi^{0} E_{kin};Events", nBins_, xMin_, xMax_);
    h_den_ = new TH1D("h_denPi0Efficiency", ";True #pi^{0} E_{kin};Events", nBins_, xMin_, xMax_);
}


bool Pi0Efficiency::IsReco(const TruePi0& tpi0, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters) {

    const auto& g1 = tpi0.photons[0];
    const auto& g2 = tpi0.photons[1];

    for (const auto& rpi0 : recoPi0s) {
        auto match = [&](const Cluster* c, const TruePhoton* g) {
            // std::cout << "Angle: " << c->p4.Vect().Angle(g->p4.Vect()) << std::endl;
            return c->p4.Vect().Angle(g->p4.Vect()) < 0.15;
        };

        if ( (match(rpi0.c1, g1) && match(rpi0.c2, g2)) ||
             (match(rpi0.c1, g2) && match(rpi0.c2, g1)) )
            return true;
    }
    return false;
}


void Pi0Efficiency::ProcessEvent(const std::vector<TruePi0>& truePi0s, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters) {
    for (const auto& tpi0 : truePi0s) {
        double E = tpi0.p4.E();
        h_den_->Fill(E);
        if (IsReco(tpi0, recoPi0s, clusters)) h_num_->Fill(E);
    }
}

void Pi0Efficiency::FinalizePlot(const std::string& outFileName) {

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

        // Set the scaled point and scaled errors back
        gEff->SetPoint(i, x, y_new);
        gEff->SetPointError(i, exl, exh, eyl_new, eyh_new);
    }

    gEff->GetXaxis()->SetTitle("True #pi^{0} E_{kin} [MeV]");
    gEff->GetYaxis()->SetTitle("Efficiency [%]");
    gEff->Draw("AP");

    PlotOptions opts;
    opts.topLatex = "#bf{Hibeam}  #it{Wasa full simulation}";
    // opts.legendEntries = { "Efficiency" };
    // opts.infoLines = {"GEANT4 #pi^{0} sample"}; // Add others
    // opts.addInfoPave = true;
    PlotGraph(gEff, outFileName.c_str(), opts);

    // cleanup
    delete gEff;
    delete h_num_;
    delete h_den_;
}
