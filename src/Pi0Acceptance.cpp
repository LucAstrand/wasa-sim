#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

#include "Pi0Acceptance.hpp"

#include <iostream>

// Pi0Acceptance::Pi0Acceptance(double massLow, double massHigh, double pi0Mass, int nBins, double eMin, double eMax)
//     : massLow_(massLow), massHigh_(massHigh), pi0Mass_(pi0Mass),
//       nBins_(nBins), eMin_(eMin), eMax_(eMax)
// {
//     h_num_ = new TH1D("h_num", "Reconstructed #pi^{0};True #pi^{0} E_{kin};Events", nBins_, eMin_, eMax_);
//     h_den_ = new TH1D("h_den", "All #pi^{0};GPS #pi^{0} E_{kin};Events", nBins_, eMin_, eMax_);
//     countTot = 0;
// }

Pi0Acceptance::Pi0Acceptance(double etaMin, double etaMax, int nBins)
    : massLow_(0), massHigh_(0), pi0Mass_(0), // unused here
      nBins_(nBins), eMin_(etaMin), eMax_(etaMax)
{
    h_num_ = new TH1D("h_num_eta", "Reconstructed #pi^{0};True #pi^{0} #eta;Events", nBins_, etaMin, etaMax);
    h_den_ = new TH1D("h_den_eta", "All #pi^{0};True #pi^{0} #eta;Events", nBins_, etaMin, etaMax);
}


double Pi0Acceptance::ComputeTruePi0Ekin(const std::vector<TruePhoton>& truePhotons) {
    if (truePhotons.size() < 2) return -1.0;
    TLorentzVector pi0_true_p4 = truePhotons[0].p4 + truePhotons[1].p4;
    double Ekin = (pi0_true_p4.E() - pi0Mass_);
    return Ekin;
}

int Pi0Acceptance::CountRecoPi0s(const std::vector<Cluster>& clusters) {
    int count = 0;
    for (size_t i = 0; i < clusters.size(); ++i) {
        for (size_t j = i + 1; j < clusters.size(); ++j) {
            
            //I don't think I should have this extra check... We should trust the clustering algorithm 
            //to give "proper" photons. If anything filtering is to be done before...
            
            // TLorentzVector p = clusters[i].p4 + clusters[j].p4;
            // double m = p.M();
            // if (m > massLow_ && m < massHigh_) {
            //     count++;
            //     countTot++;
            //     std::cout << "[CountRecoPi0s] RECO count: " << count << ", " << countTot << std::endl;
            // }

            //No extra filtering 
            count++;
            countTot++;
            std::cout << "[CountRecoPi0s] RECO count: " << count << ", " << countTot << std::endl;
        }
    }
    return count;
}

int Pi0Acceptance::CountTruePi0s(const std::vector<TruePhoton>& truePhotons) {
    int count = 0;
    for (size_t i = 0; i < truePhotons.size(); ++i) {
        for (size_t j = i + 1; j < truePhotons.size(); ++j) {
            
            //I don't think I should have this extra check... We should trust the clustering algorithm 
            //to give "proper" photons. If anything filtering is to be done before...
            
            // TLorentzVector p = clusters[i].p4 + clusters[j].p4;
            // double m = p.M();
            // if (m > massLow_ && m < massHigh_) {
            //     count++;
            //     countTot++;
            //     std::cout << "[CountRecoPi0s] RECO count: " << count << ", " << countTot << std::endl;
            // }

            //No extra filtering 
            count++;
            // countTot++;
            // std::cout << "[CountRecoPi0s] RECO count: " << count << ", " << countTot << std::endl;
        }
    }
    return count;
}

// void Pi0Acceptance::ProcessEvent(const std::vector<Cluster>& clusters,
//                                  const std::vector<TruePhoton>& truePhotons,
//                                  const double primaryEkin) {
                                    
//     // Use the generator-level kinetic energy directly
//     h_den_->Fill(primaryEkin);

//     // Reconstructed candidates
//     int Nreco = CountRecoPi0s(clusters);
//     // if (Nreco > 0)
//     //     h_num_->Fill(primaryEkin);
//     // True candidates                           
//     int Ntrue = CountTruePi0s(truePhotons);
//     if (Ntrue > 0) 
//         h_num_->Fill(primaryEkin);
// }

void Pi0Acceptance::ProcessEvent(const std::vector<Cluster>& clusters,
                                 const std::vector<TruePhoton>& truePhotons,
                                 double primaryEta)
{
    h_den_->Fill(primaryEta);

//     // Reconstructed candidates                           
    // int Nreco = CountRecoPi0s(clusters);
    // if (Nreco > 0) h_num_->Fill(primaryEta);
//     // True candidates                           
    int Ntrue = CountTruePi0s(truePhotons);
    if (Ntrue > 0) h_num_->Fill(primaryEta);
}

void Pi0Acceptance::FinalizePlot(const std::string& outFileName) {

    std::cout << "[FinilizePlot] Total RECO count: " << countTot << std::endl;

    TCanvas* c = new TCanvas("cEff", "Pi0 Acceptance", 900, 650);

    // Create graph with binomial (asymmetric) errors
    TGraphAsymmErrors* gEff = new TGraphAsymmErrors(h_num_, h_den_, "cl=0.683 b(1,1) mode"); // fraction in [0,1]
    gEff->SetTitle("");
    gEff->SetName("Pi0AccemptanceGraph");

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

    gEff->SetMarkerStyle(20);
    gEff->SetMarkerSize(1.0);
    gEff->SetMarkerColor(kBlack);
    gEff->SetLineColor(kBlack);

    //vs Ekin
    // gEff->GetXaxis()->SetTitle("GPS #pi^{0} E_{kin} [MeV]");
    // gEff->GetYaxis()->SetTitle("Acceptance [%]");
    //vs eta 
    gEff->GetXaxis()->SetTitle("GPS #pi^{0} pseudorapidity #eta");
    gEff->GetYaxis()->SetTitle("Acceptance [%]");

    gEff->GetYaxis()->SetRangeUser(0.0, 110.0);

    gStyle->SetOptStat(0);
    gEff->Draw("AP");

    auto legend = new TLegend(0.55, 0.75, 0.88, 0.88);
    legend->AddEntry(gEff, "Acceptance", "p");
    legend->Draw();

    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

    c->SaveAs(outFileName.c_str());
    // std::cout << "[Pi0Efficiency] Saved plot to " << outFileName << std::endl;

    // cleanup
    delete gEff;
    delete c;
}
