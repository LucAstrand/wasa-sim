#include "PIDEfficiency.hpp"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"

#include "PlotUtils.hpp"

PIDEfficiency::PIDEfficiency(int nBins, double eMin, double eMax)
    : nBins_(nBins), eMin_(eMin), eMax_(eMax)
{
    h_numPion_ = new TH1D("h_numPion", ";True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_numPion_miss = new TH1D("h_numPion", ";True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_denPion_ = new TH1D("h_denPion", "All Pions;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);

    h_numProton_ = new TH1D("h_numProton", ";True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_numProton_miss = new TH1D("h_numProton", ";True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_denProton_ = new TH1D("h_denProton", "All Protons;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);

    h_numElectron_ = new TH1D("h_numElectron", ";True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_numElectron_miss = new TH1D("h_numElectron", ";True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_denElectron_ = new TH1D("h_denElectron", "All Electrons;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
}

void PIDEfficiency::ProcessEvent(const std::vector<ChargedCluster>& clusters)
{
    for (auto& c : clusters) {
        // Fill denominator first (all truth)
        switch (c.objectTruePDG) {
            case 211:  h_denPion_->Fill(c.objectTrueKE); break;
            case 2212: h_denProton_->Fill(c.objectTrueKE); break;
            case 11:   h_denElectron_->Fill(c.objectTrueKE); break;
            default: break;
        }

        // Fill numerator only if PID guess matches true PDG
        if ((c.objectTruePDG == 211 && c.pidGuess == PID::Pion) ||
            (c.objectTruePDG == 2212 && c.pidGuess == PID::Proton) ||
            (c.objectTruePDG == 11 && c.pidGuess == PID::Electron)) 
        {
            switch (c.objectTruePDG) {
                case 211:  h_numPion_->Fill(c.objectTrueKE); break;
                case 2212: h_numProton_->Fill(c.objectTrueKE); break;
                case 11:   h_numElectron_->Fill(c.objectTrueKE); break;
            }
        }
        else if ((c.objectTruePDG == 211 && c.pidGuess != PID::Pion) ||
                (c.objectTruePDG == 2212 && c.pidGuess != PID::Proton) ||
                (c.objectTruePDG == 11 && c.pidGuess != PID::Electron)) 
        {
            switch (c.objectTruePDG) {
                case 211:  h_numPion_miss->Fill(c.objectTrueKE); break;
                case 2212: h_numProton_miss->Fill(c.objectTrueKE); break;
                case 11:   h_numElectron_miss->Fill(c.objectTrueKE); break;
            }
        }
    }
}

void PIDEfficiency::FinalizePlot(const std::string& outFileName, int pdgNumToPlot)
{
    TGraphAsymmErrors* gEff = nullptr;
    TGraphAsymmErrors* gEffmiss = nullptr;

    switch (pdgNumToPlot) {

        case 211:
            gEff = new TGraphAsymmErrors(
                h_numPion_, h_denPion_, "cl=0.683 b(1,1) mode"
            );
            gEffmiss = new TGraphAsymmErrors(
                h_numPion_miss, h_denPion_, "cl=0.683 b(1,1) mode"
            );
            break;

        case 2212:
            gEff = new TGraphAsymmErrors(
                h_numProton_, h_denProton_, "cl=0.683 b(1,1) mode"
            );
            gEffmiss = new TGraphAsymmErrors(
                h_numProton_miss, h_denProton_, "cl=0.683 b(1,1) mode"
            );
            break;

        case 11:
            gEff = new TGraphAsymmErrors(
                h_numElectron_, h_denElectron_, "cl=0.683 b(1,1) mode"
            );
            gEffmiss = new TGraphAsymmErrors(
                h_numElectron_miss, h_denElectron_, "cl=0.683 b(1,1) mode"
            );
            break;

        default:
            std::cerr << "[PIDEfficiency] Unknown PDG " << pdgNumToPlot
                      << ", defaulting to pion\n";
            gEff = new TGraphAsymmErrors(
                h_numPion_, h_denPion_, "cl=0.683 b(1,1) mode"
            );
            gEffmiss = new TGraphAsymmErrors(
                h_numPion_miss, h_denPion_, "cl=0.683 b(1,1) mode"
            );
            break;
    }

    if (!gEff) return;


    gEff->GetXaxis()->SetTitle("True KE [MeV]");
    gEff->GetYaxis()->SetTitle("Efficiency [%]");
    gEff->GetYaxis()->SetRangeUser(0, 110);
    gEffmiss->GetXaxis()->SetTitle("True KE [MeV]");
    gEffmiss->GetYaxis()->SetTitle("miss-ID [%]");
    gEffmiss->GetYaxis()->SetRangeUser(0, 110);

    // scale to percent
    for (int i = 0; i < gEff->GetN(); ++i) {
        double x, y;
        gEff->GetPoint(i, x, y);
        gEff->SetPoint(i, x, 100.0 * y);
        gEff->SetPointError(
            i,
            gEff->GetErrorXlow(i),
            gEff->GetErrorXhigh(i),
            100.0 * gEff->GetErrorYlow(i),
            100.0 * gEff->GetErrorYhigh(i)
        );
    }
    // scale to percent
    for (int i = 0; i < gEff->GetN(); ++i) {
        double x, y;
        gEffmiss->GetPoint(i, x, y);
        gEffmiss->SetPoint(i, x, 100.0 * y);
        gEffmiss->SetPointError(
            i,
            gEffmiss->GetErrorXlow(i),
            gEffmiss->GetErrorXhigh(i),
            100.0 * gEffmiss->GetErrorYlow(i),
            100.0 * gEffmiss->GetErrorYhigh(i)
        );
    }

    PlotOptions opts;
    opts.legendEntries = { "#pi^{+} ID efficiency" };
    opts.topLatex = "#bf{Hibeam}  #it{Wasa full simulation}";
    opts.legendX1 = 0.65;
    opts.legendX2 = 0.98;

    PlotOptions optsmiss;
    optsmiss.legendEntries = { "#pi^{+} miss-ID percentage" };
    optsmiss.topLatex = "#bf{Hibeam}  #it{Wasa full simulation}";
    optsmiss.legendX1 = 0.65;
    optsmiss.legendX2 = 0.98;

    // PlotGraph(gEff, "pid_efficiency_pion.png", opts);
    // PlotGraph(gEffmiss, "miss_pid_percentage_pion.png", optsmiss);

    switch (pdgNumToPlot) {
        case 211:  PlotGraph(gEff, "pid_efficiency_pion.png", opts); PlotGraph(gEffmiss, "miss_pid_percentage_pion.png", optsmiss); break;
        case 2212: PlotGraph(gEff, "pid_efficiency_proton.png", opts); PlotGraph(gEffmiss, "miss_pid_percentage_proton.png", optsmiss); break;
        case 11:   PlotGraph(gEff, "pid_efficiency_electron.png", opts); PlotGraph(gEffmiss, "miss_pid_percentage_electron.png", optsmiss); break;
        default: break;
    }

    delete gEff;
    delete gEffmiss;

}
