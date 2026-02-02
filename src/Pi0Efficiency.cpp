#include "Pi0Efficiency.hpp"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

#include "PlotUtils.hpp"

#include <iostream>

Pi0Efficiency::Pi0Efficiency(double massLow, double massHigh, double pi0Mass, int nBins, double eMin, double eMax)
    : massLow_(massLow), massHigh_(massHigh), pi0Mass_(pi0Mass),
      nBins_(nBins), eMin_(eMin), eMax_(eMax)
{
    h_num_ = new TH1D("h_numPi0Efficiency", ";True #pi^{0} E_{kin};Events", nBins_, eMin_, eMax_);
    h_den_ = new TH1D("h_denPi0Efficiency", ";True #pi^{0} E_{kin};Events", nBins_, eMin_, eMax_);
}

// double Pi0Efficiency::ComputeTruePi0Ekin(const std::vector<TruePhoton>& truePhotons) {
//     if (truePhotons.size() < 2) return -1.0;
//     TLorentzVector pi0_true_p4 = truePhotons[0].p4 + truePhotons[1].p4;
//     double Ekin = (pi0_true_p4.E() - pi0Mass_);
//     return Ekin;
// }

// int Pi0Efficiency::CountTruePi0s(const std::vector<TruePhoton>& truePhotons) {
//     int count = 0;
//     for (size_t i = 0; i < truePhotons.size(); ++i) {
//         for (size_t j = i + 1; j < truePhotons.size(); ++j) {
//             count++;
//         }
//     }
//     return count;
// }

// int Pi0Efficiency::CountRecoPi0s(const std::vector<Cluster>& clusters) {
//     int count = 0;
//     for (size_t i = 0; i < clusters.size(); ++i) {
//         for (size_t j = i + 1; j < clusters.size(); ++j) {
//             TLorentzVector p = clusters[i].p4 + clusters[j].p4;
//             double m = p.M();
//             if (m > massLow_ && m < massHigh_) count++;
//         }
//     }
//     return count;
// }

bool Pi0Efficiency::IsReco(const TruePi0& tpi0, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters) {

    const auto& g1 = tpi0.photons[0];
    const auto& g2 = tpi0.photons[1];

    for (const auto& rpi0 : recoPi0s) {
        auto match = [&](const Cluster* c, const TruePhoton* g) {
            return c->p4.Vect()
                       .Angle(g->p4.Vect()) < 0.05;
        };

        if ( (match(rpi0.c1, g1) && match(rpi0.c2, g2)) ||
             (match(rpi0.c1, g2) && match(rpi0.c2, g1)) )
            return true;
    }
    return false;
}


// void Pi0Efficiency::ProcessEvent(const std::vector<Cluster>& clusters, const std::vector<TruePhoton>& truePhotons) {
void Pi0Efficiency::ProcessEvent(const std::vector<TruePi0>& truePi0s, const std::vector<Pi0Candidate>& recoPi0s, const std::vector<Cluster>& clusters) {

    // double Ekin_true = ComputeTruePi0Ekin(truePhotons);
    // if (Ekin_true < 0) return;

    // h_den_->Fill(Ekin_true);
    // int Nreco = CountRecoPi0s(clusters);
    // if (Nreco > 0) h_num_->Fill(Ekin_true);
    for (const auto& tpi0 : truePi0s) {
        double E = tpi0.p4.E();
        h_den_->Fill(E);
        if (IsReco(tpi0, recoPi0s, clusters)) h_num_->Fill(E);
    }
}

void Pi0Efficiency::FinalizePlot(const std::string& outFileName) {
    // TCanvas* c = new TCanvas("cEff", "Pi0 Efficiency", 900, 650);

    std::cout << "Total reco Pi0s: " << h_num_->GetEntries()
              << "Total true Pi0s: " << h_den_->GetEntries()
              << std::endl;

    // Create graph with binomial (asymmetric) errors
    TGraphAsymmErrors* gEff = new TGraphAsymmErrors(h_num_, h_den_, "cl=0.683 b(1,1) mode"); // fraction in [0,1]
    // gEff->SetTitle("");
    // gEff->SetName("Pi0EfficiencyGraph");

    // --- SCALE TO PERCENT (important) ---
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

    // gEff->SetMarkerStyle(20);
    // gEff->SetMarkerSize(1.0);
    // gEff->SetMarkerColor(kBlack);
    // gEff->SetLineColor(kBlack);
    // gEff->GetYaxis()->SetRangeUser(0.0, 110.0);
    
    // gStyle->SetOptStat(0);
    gEff->GetXaxis()->SetTitle("True #pi^{0} E_{kin} [MeV]");
    gEff->GetYaxis()->SetTitle("Efficiency [%]");
    gEff->Draw("AP");

    // auto legend = new TLegend(0.55, 0.75, 0.88, 0.88);
    // legend->AddEntry(gEff, "Efficiency", "p");
    // legend->SetFillStyle(0);   // no fill
    // legend->SetBorderSize(0);  // no border box
    // legend->Draw();

    // TLatex l;
    // l.SetNDC();
    // l.SetTextFont(42);
    // l.SetTextSize(0.045);
    // l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

    // c->SaveAs(outFileName.c_str());
    PlotOptions opts;
    opts.topLatex = "#bf{Hibeam}  #it{Wasa full simulation}";
    // opts.legendEntries = { "Efficiency" };
    // opts.infoLines = {"GEANT4 #pi^{0} sample"}; // Add others
    // opts.addInfoPave = true;
    PlotGraph(gEff, outFileName.c_str(), opts);
    // std::cout << "[Pi0Efficiency] Saved plot to " << outFileName << std::endl;

    // cleanup
    delete gEff;
    // delete c;
    delete h_num_;
    delete h_den_;
}
