#include "PIDEfficiency.hpp"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"

PIDEfficiency::PIDEfficiency(int nBins, double eMin, double eMax)
    : nBins_(nBins), eMin_(eMin), eMax_(eMax)
{
    h_numPion_ = new TH1D("h_numPion", "Pion PID Efficiency;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_denPion_ = new TH1D("h_denPion", "All Pions;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);

    h_numProton_ = new TH1D("h_numProton", "Proton PID Efficiency;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_denProton_ = new TH1D("h_denProton", "All Protons;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);

    h_numElectron_ = new TH1D("h_numElectron", "Electron PID Efficiency;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
    h_denElectron_ = new TH1D("h_denElectron", "All Electrons;True E_{kin} [MeV];Events", nBins_, eMin_, eMax_);
}

void PIDEfficiency::ProcessEvent(const std::vector<ChargedCluster>& clusters)
{
    for (auto& c : clusters) {
        std::cout << "cluster True E: " << c.objectTrueKE << std::endl;
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
    }
}

void PIDEfficiency::FinalizePlot(const std::string& outFileName)
{
    auto drawGraph = [&](TH1D* hNum, TH1D* hDen, const char* title, int color) {
        TCanvas* c = new TCanvas("c", title, 900, 650);
        TGraphAsymmErrors* gEff = new TGraphAsymmErrors(hNum, hDen, "cl=0.683 b(1,1) mode");
        gEff->SetMarkerStyle(20);
        gEff->SetMarkerSize(1.0);
        gEff->SetMarkerColor(color);
        gEff->SetLineColor(color);

        gEff->GetXaxis()->SetTitle("True E_{kin} [MeV]");
        gEff->GetYaxis()->SetTitle("Efficiency [%]");

        // scale to percent
        for (int i = 0; i < gEff->GetN(); ++i) {
            double x, y;
            gEff->GetPoint(i, x, y);
            gEff->SetPoint(i, x, y*100);
            gEff->SetPointError(i,
                                gEff->GetErrorXlow(i),
                                gEff->GetErrorXhigh(i),
                                gEff->GetErrorYlow(i)*100,
                                gEff->GetErrorYhigh(i)*100);
        }

        gEff->Draw("AP");

        TLatex l;
        l.SetNDC();
        l.SetTextFont(42);
        l.SetTextSize(0.045);
        l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

        c->SaveAs((std::string(outFileName) + "_" + title + ".png").c_str());

        delete gEff;
        delete c;
    };

    drawGraph(h_numPion_, h_denPion_, "PionEfficiency", kRed+1);
    drawGraph(h_numProton_, h_denProton_, "ProtonEfficiency", kBlue+1);
    drawGraph(h_numElectron_, h_denElectron_, "ElectronEfficiency", kGreen+2);
}
