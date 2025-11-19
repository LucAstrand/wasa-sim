#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

#include "Pi0Acceptance.hpp"

#include <iostream>

Pi0Acceptance::Pi0Acceptance(double massLow, double massHigh, double pi0Mass, int nBins, double eMin, double eMax)
    : massLow_(massLow), massHigh_(massHigh), pi0Mass_(pi0Mass),
      nBins_(nBins), eMin_(eMin), eMax_(eMax)
{
    h_num_ = new TH1D("h_num", "Reconstructed #pi^{0};True #pi^{0} E_{kin};Events", nBins_, eMin_, eMax_);
    h_den_ = new TH1D("h_den", "All #pi^{0};GPS #pi^{0} E_{kin};Events", nBins_, eMin_, eMax_);
    countTot = 0;
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

void Pi0Acceptance::ProcessEventTwoHist(const std::vector<Cluster>& clusters,
                                 const std::vector<TruePhoton>& truePhotons,
                                 double primaryEkin,
                                 double primaryTheta)
{   
    if (primaryEkin == 50) {
    h_den_1->Fill(primaryTheta);

//     // Reconstructed candidates                           
    // int Nreco = CountRecoPi0s(clusters);
    // if (Nreco > 0) h_num_->Fill(primaryEta);
//     // True candidates                           
    int Ntrue = CountTruePi0s(truePhotons);
    if (Ntrue > 0) h_num_1->Fill(primaryTheta);
    }   

    if (primaryEkin == 500) {
    h_den_2->Fill(primaryTheta);

//     // Reconstructed candidates                           
    // int Nreco = CountRecoPi0s(clusters);
    // if (Nreco > 0) h_num_->Fill(primaryEta);
//     // True candidates                           
    int Ntrue = CountTruePi0s(truePhotons);
    if (Ntrue > 0) h_num_2->Fill(primaryTheta);
    }
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

    // auto legend = new TLegend(0.55, 0.75, 0.88, 0.88);
    // legend->AddEntry(gEff, "Acceptance", "p");
    // legend->SetFillStyle(0);   // no fill
    // legend->SetBorderSize(0);  // no border box
    // legend->Draw();
    
    TPaveText *info = new TPaveText(0.17, 0.65, 0.50, 0.85, "NDC");  // x1,y1,x2,y2 normalized coordinates
    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextFont(42);
    info->SetTextSize(0.04);
    // info->AddText("Hibeam Wasafull simulation");
    info->AddText("GEANT4 #pi^{0} sample");
    info->AddText("1000 simulated events");
    info->AddText("E_{kin} #in [1, 500] MeV");
    // info->AddText("E_{kin} = 100 MeV");
    info->Draw();

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

// void Pi0Acceptance::FinalizePlot(const std::string& outFileName)
// {
//     std::cout << "[FinalizePlot] Creating acceptance overlay..." << std::endl;

//     TCanvas* c = new TCanvas("cAcc", "Pi0 Acceptance Overlay", 900, 650);
//     gStyle->SetOptStat(0);

//     // --- Graph 1: E_kin = 50 MeV ---
//     TGraphAsymmErrors* gAcc1 = new TGraphAsymmErrors(h_num_1, h_den_1, "cl=0.683 b(1,1) mode");
//     gAcc1->SetName("Pi0Acceptance_50MeV");
//     gAcc1->SetTitle("");
//     gAcc1->SetMarkerStyle(20);
//     gAcc1->SetMarkerSize(1.0);
//     gAcc1->SetMarkerColor(kBlue+1);
//     gAcc1->SetLineColor(kBlue+1);

//     // --- Graph 2: E_kin = 500 MeV ---
//     TGraphAsymmErrors* gAcc2 = new TGraphAsymmErrors(h_num_2, h_den_2, "cl=0.683 b(1,1) mode");
//     gAcc2->SetName("Pi0Acceptance_500MeV");
//     gAcc2->SetMarkerStyle(21);
//     gAcc2->SetMarkerSize(1.0);
//     gAcc2->SetMarkerColor(kRed+1);
//     gAcc2->SetLineColor(kRed+1);

//     // --- Convert y values to percent for both graphs ---
//     auto scaleToPercent = [](TGraphAsymmErrors* g){
//         int n = g->GetN();
//         for (int i = 0; i < n; ++i) {
//             double x, y;
//             g->GetPoint(i, x, y);
//             double exl = g->GetErrorXlow(i);
//             double exh = g->GetErrorXhigh(i);
//             double eyl = g->GetErrorYlow(i);
//             double eyh = g->GetErrorYhigh(i);
//             g->SetPoint(i, x, y * 100.0);
//             g->SetPointError(i, exl, exh, eyl * 100.0, eyh * 100.0);
//         }
//     };
//     scaleToPercent(gAcc1);
//     scaleToPercent(gAcc2);

//     // --- Axes setup ---
//     gAcc1->GetXaxis()->SetTitle("GPS #pi^{0} #theta [rad]");
//     gAcc1->GetYaxis()->SetTitle("Acceptance [%]");
//     gAcc1->GetYaxis()->SetRangeUser(0.0, 110.0);

//     // --- Draw both ---
//     gAcc1->Draw("AP");
//     gAcc2->Draw("P SAME");

//     // --- Legend ---
//     // auto legend = new TLegend(0.55, 0.75, 0.88, 0.88);
//     auto legend = new TLegend(0.65, 0.75, 0.98, 0.88);
//     legend->AddEntry(gAcc1, "E_{kin} = 50 MeV", "p");
//     legend->AddEntry(gAcc2, "E_{kin} = 500 MeV", "p");
//     legend->SetFillStyle(0);
//     legend->SetBorderSize(0);
//     legend->SetTextSize(0.035);
//     legend->Draw();

//     // --- Info box ---
//     // TPaveText *info = new TPaveText(0.17, 0.65, 0.50, 0.85, "NDC");
//     // info->SetFillStyle(0);
//     // info->SetBorderSize(0);
//     // info->SetTextFont(42);
//     // info->SetTextSize(0.04);
//     // info->AddText("GEANT4 #pi^{0} sample");
//     // info->AddText("1000 simulated events");
//     // info->Draw();

//     // --- Label at top ---
//     TLatex l;
//     l.SetNDC();
//     l.SetTextFont(42);
//     l.SetTextSize(0.045);
//     l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

//     c->SaveAs(outFileName.c_str());
//     std::cout << "[FinalizePlot] Saved overlay to " << outFileName << std::endl;

//     delete gAcc1;
//     delete gAcc2;
//     delete c;
// }

