#include "PlotUtils.hpp"

#include "TStyle.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TF1.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TCanvas.h"

#include <filesystem>
// namespace fs = std::filesystem;

void SetPrettyStyle() {

    //Make sure plots folder exists, if not it creates it
    std::filesystem::create_directories("plots");

    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetLabelSize(0.045, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.4, "Y");
    // gStyle->SetTextSize(0.045);
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
}

void PrettyPi0MassPlot(TH1F* hPi0Mass, TString plotname, double fitMin, double fitMax) {
    SetPrettyStyle();

    TCanvas *c = new TCanvas("cPi0", "Pi0 Mass", 800, 600);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);

    // TF1 *fGaus = new TF1("fGaus", "gaus", 100, 170);
    TF1 *fGaus = new TF1("fGaus", "gaus", fitMin, fitMax);
    fGaus->SetLineColor(kRed+1);
    fGaus->SetLineWidth(2);
    // hPi0Mass->Fit(fGaus, "RQ");

    double maxBin = hPi0Mass->GetMaximum();
    double meanGuess = hPi0Mass->GetBinCenter(hPi0Mass->GetMaximumBin());
    double sigmaGuess = hPi0Mass->GetRMS();  // or ~5 MeV if you want
    fGaus->SetParameters(maxBin, meanGuess, sigmaGuess);
    hPi0Mass->Fit(fGaus, "RQ");

    double mean  = fGaus->GetParameter(1);
    double sigma = fGaus->GetParameter(2);
    double errMu = fGaus->GetParError(1);
    double errSi = fGaus->GetParError(2);

    hPi0Mass->SetLineColor(kBlack);
    hPi0Mass->SetLineWidth(2);
    hPi0Mass->GetXaxis()->SetTitle("M_{#gamma#gamma} [MeV]");
    hPi0Mass->GetYaxis()->SetTitle("Events");

    hPi0Mass->Draw("HIST");
    fGaus->Draw("SAME");

    TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(hPi0Mass, "Reconstructed #pi^{0} Invariant mass", "l");
    leg->AddEntry(fGaus, "Gaussian Fit:", "l");
    leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f MeV", mean, errMu), "");
    leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f MeV", sigma, errSi), "");
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);   // no fill
    leg->SetBorderSize(0);  // no border box
    leg->Draw();

    TPaveText *info = new TPaveText(0.17, 0.70, 0.50, 0.90, "NDC");  // x1,y1,x2,y2 normalized coordinates
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

    c->SaveAs("plots/" + plotname);
    // c->SaveAs("Pi0Mass_truth.png");

    // clean up
    delete leg;
    delete fGaus;
    delete c;
}

void TruthPi0MassPlot(TH1F* hPi0Mass, TString plotname) {
    SetPrettyStyle();
    gStyle->SetOptStat(1);

    TCanvas *c = new TCanvas("cPi0", "Pi0 Mass", 800, 600);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);

    TF1 *fGaus = new TF1("fGaus", "gaus", 100, 170);
    fGaus->SetLineColor(kRed+1);
    fGaus->SetLineWidth(2);
    // hPi0Mass->Fit(fGaus, "RQ");

    // double mean  = fGaus->GetParameter(1);
    // double sigma = fGaus->GetParameter(2);
    // double errMu = fGaus->GetParError(1);
    // double errSi = fGaus->GetParError(2);

    hPi0Mass->SetLineColor(kBlack);
    hPi0Mass->SetLineWidth(2);
    hPi0Mass->GetXaxis()->SetTitle("M_{#gamma#gamma} [MeV]");
    hPi0Mass->GetYaxis()->SetTitle("Events");

    hPi0Mass->Draw("HIST");
    // fGaus->Draw("SAME");
    Double_t mean = hPi0Mass->GetMean();

    TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(hPi0Mass, "Reconstructed #pi^{0} Invariant mass", "l");
    // leg->AddEntry(fGaus, "Gaussian Fit:", "l");
    leg->AddEntry((TObject*)0, Form("mean = %.3f MeV", mean), "");
    // leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f MeV", mean, errMu), "");
    // leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f MeV", sigma, errSi), "");
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);   // no fill
    leg->SetBorderSize(0);  // no border box
    leg->Draw();

    TPaveText *info = new TPaveText(0.17, 0.70, 0.50, 0.90, "NDC");  // x1,y1,x2,y2 normalized coordinates
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

    c->SaveAs("plots/" + plotname);
    // c->SaveAs("Pi0Mass_truth.png");

    // clean up
    // delete leg;
    delete fGaus;
    delete c;
}

void PrettyPi0NumClusterPlot(TH1F* hNCluster) {
    SetPrettyStyle();

    TCanvas *c = new TCanvas("cPi0", "Pi0 Mass", 800, 600);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);

    // TF1 *fGaus = new TF1("fGaus", "gaus", 100, 170);
    // fGaus->SetLineColor(kRed+1);
    // fGaus->SetLineWidth(2);
    // hPi0Mass->Fit(fGaus, "RQ");

    // double mean  = fGaus->GetParameter(1);
    // double sigma = fGaus->GetParameter(2);
    // double errMu = fGaus->GetParError(1);
    // double errSi = fGaus->GetParError(2);

    hNCluster->SetLineColor(kBlack);
    hNCluster->SetLineWidth(2);
    hNCluster->GetXaxis()->SetTitle("M_{#gamma#gamma} [MeV]");
    hNCluster->GetYaxis()->SetTitle("Events");

    hNCluster->Draw("HIST");
    // fGaus->Draw("SAME");

    int numEvts2Photon = hNCluster->GetBinContent(hNCluster->FindBin(2));

    TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(hNCluster, "Number of clusters per event", "l");
    // leg->AddEntry(fGaus, "Gaussian Fit:", "l");
    // leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f MeV", mean, errMu), "");
    // leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f MeV", sigma, errSi), "");

    leg->AddEntry((TObject*)0, Form("2 photon events = %d", numEvts2Photon), "");

    leg->SetTextSize(0.03);
    leg->Draw();

    TPaveText *info = new TPaveText(0.17, 0.70, 0.50, 0.90, "NDC");  // x1,y1,x2,y2 normalized coordinates
    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextFont(42);
    info->SetTextSize(0.04);
    // info->AddText("Hibeam Wasafull simulation");
    info->AddText("GEANT4 #pi^{0} sample");
    info->AddText("1000 simulated events");
    info->AddText("E_{kin} #in [1, 500] MeV");
    leg->SetFillStyle(0);   // no fill
    leg->SetBorderSize(0);  // no border box
    info->Draw();

    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

    // c->SaveAs("Pi0Mass_pretty.png");
    // c->SaveAs("Pi0Mass_truth.png");
    c->SaveAs("plots/ClusterNum.png");

    // clean up
    delete leg;
    // delete fGaus;
    delete c;
}

// void Pi0ClusterNumPlotEkin(TH1F* hNClusters_lowE, TH1F* hNClusters_midE, TH1F* hNClusters_highE) {
//     SetPrettyStyle();

//     TCanvas *c = new TCanvas("c", "Cluster multiplicity vs E_kin", 800, 600);
//     hNClusters_lowE->SetLineColor(kRed);
//     hNClusters_midE->SetLineColor(kBlue);
//     hNClusters_highE->SetLineColor(kGreen+2);

//     double maxY = std::max({
//     hNClusters_lowE->GetMaximum(),
//     hNClusters_midE->GetMaximum(),
//     hNClusters_highE->GetMaximum()
//     });

//     hNClusters_lowE->SetMaximum(1.2 * maxY);  

//     hNClusters_lowE->Draw();
//     hNClusters_midE->Draw("SAME");
//     hNClusters_highE->Draw("SAME");

//     auto legend = new TLegend(0.6,0.7,0.88,0.88);
//     legend->AddEntry(hNClusters_lowE, "E_{kin} < 200 MeV", "l");
//     legend->AddEntry(hNClusters_midE, "200-400 MeV", "l");
//     legend->AddEntry(hNClusters_highE, "E_{kin} > 400 MeV", "l");
//     legend->SetFillStyle(0);   // no fill
//     legend->SetBorderSize(0);  // no border box
//     legend->Draw();

//     // int numEvts2Photon = hNCluster->GetBinContent(hNCluster->FindBin(2));

//     // TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
//     // leg->AddEntry(hNCluster, "Number of clusters per event", "l");
//     // // leg->AddEntry(fGaus, "Gaussian Fit:", "l");
//     // // leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f MeV", mean, errMu), "");
//     // // leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f MeV", sigma, errSi), "");

//     // leg->AddEntry((TObject*)0, Form("2 photon events = %d", numEvts2Photon), "");

//     // leg->SetTextSize(0.03);
//     // leg->Draw();

//     TPaveText *info = new TPaveText(0.6,0.2,0.88,0.48, "NDC");  // x1,y1,x2,y2 normalized coordinates
//     info->SetFillStyle(0);
//     info->SetBorderSize(0);
//     info->SetTextFont(42);
//     info->SetTextSize(0.04);
//     // info->AddText("Hibeam Wasafull simulation");
//     info->AddText("GEANT4 #pi^{0} sample");
//     info->AddText("1000 simulated events");
//     info->AddText("E_{kin} #in [1, 500] MeV");
//     info->Draw();

//     TLatex l;
//     l.SetNDC();
//     l.SetTextFont(42);
//     l.SetTextSize(0.045);
//     l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

//     c->SaveAs("plots/ClusterNum_Ekin.png");

//     // clean up
//     delete legend;
//     delete c;
// }

void Pi0ClusterNumPlotEkin(TH1F* hNClusters_lowE, TH1F* hNClusters_highE) {
    SetPrettyStyle();

    TCanvas *c = new TCanvas("c", "Cluster multiplicity vs E_kin", 800, 600);
    hNClusters_lowE->SetLineColor(kRed);
    // hNClusters_midE->SetLineColor(kBlue);
    hNClusters_highE->SetLineColor(kGreen+2);

    double maxY = std::max({
    hNClusters_lowE->GetMaximum(),
    // hNClusters_midE->GetMaximum(),
    hNClusters_highE->GetMaximum()
    });

    hNClusters_lowE->SetMaximum(1.2 * maxY);  

    hNClusters_lowE->Draw();
    // hNClusters_midE->Draw("SAME");
    hNClusters_highE->Draw("SAME");

    auto legend = new TLegend(0.6,0.7,0.88,0.88);
    legend->AddEntry(hNClusters_lowE, "E_{kin} = 50 MeV", "l");
    // legend->AddEntry(hNClusters_midE, "200-400 MeV", "l");
    legend->AddEntry(hNClusters_highE, "E_{kin} = 500 MeV", "l");
    legend->SetFillStyle(0);   // no fill
    legend->SetBorderSize(0);  // no border box
    legend->Draw();

    // int numEvts2Photon = hNCluster->GetBinContent(hNCluster->FindBin(2));

    // TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    // leg->AddEntry(hNCluster, "Number of clusters per event", "l");
    // // leg->AddEntry(fGaus, "Gaussian Fit:", "l");
    // // leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f MeV", mean, errMu), "");
    // // leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f MeV", sigma, errSi), "");

    // leg->AddEntry((TObject*)0, Form("2 photon events = %d", numEvts2Photon), "");

    // leg->SetTextSize(0.03);
    // leg->Draw();

    TPaveText *info = new TPaveText(0.6,0.2,0.88,0.48, "NDC");  // x1,y1,x2,y2 normalized coordinates
    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextFont(42);
    info->SetTextSize(0.04);
    // info->AddText("Hibeam Wasafull simulation");
    info->AddText("GEANT4 #pi^{0} sample");
    info->AddText("1000 simulated events");
    info->AddText("E_{kin} #in [1, 500] MeV");
    info->Draw();

    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

    c->SaveAs("plots/ClusterNum_Ekin_50_150_MeV.png");

    // clean up
    delete legend;
    delete c;
}


void BasicHistPlot(TH1F* histogram) {
    SetPrettyStyle();

    TCanvas *c = new TCanvas("cPi0", "Pi0 Mass", 800, 600);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);

    // TF1 *fGaus = new TF1("fGaus", "gaus", 100, 170);
    // fGaus->SetLineColor(kRed+1);
    // fGaus->SetLineWidth(2);
    // hPi0Mass->Fit(fGaus, "RQ");

    // double mean  = fGaus->GetParameter(1);
    // double sigma = fGaus->GetParameter(2);
    // double errMu = fGaus->GetParError(1);
    // double errSi = fGaus->GetParError(2);

    histogram->SetLineColor(kBlack);
    histogram->SetLineWidth(2);
    histogram->GetXaxis()->SetTitle("M_{#gamma#gamma} [MeV]");
    histogram->GetYaxis()->SetTitle("Events");

    histogram->Draw("HIST");

    TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(histogram, "Counts", "l");
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);   // no fill
    leg->SetBorderSize(0);  // no border box
    leg->Draw();

    TPaveText *info = new TPaveText(0.17, 0.70, 0.50, 0.90, "NDC");  // x1,y1,x2,y2 normalized coordinates
    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextFont(42);
    info->SetTextSize(0.04);
    // info->AddText("Hibeam Wasafull simulation");
    info->AddText("GEANT4 #pi^{0} sample");
    info->AddText("1000 simulated events");
    info->AddText("E_{kin} #in [1, 500] MeV");
    info->Draw();

    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

    // c->SaveAs("Pi0Mass_pretty.png");
    // c->SaveAs("Pi0Mass_truth.png");
    c->SaveAs("plots/BasicPlot.png");

    // clean up
    delete leg;
    // delete fGaus;
    delete c;
}

void nSigmaPlot(TH1F* hNSigma, TString plotname, double fitMin, double fitMax) {
    SetPrettyStyle();

    TCanvas *c = new TCanvas("cPi0", "Pi0 Mass", 800, 600);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);

    // TF1 *fGaus = new TF1("fGaus", "gaus", 100, 170);
    TF1 *fGaus = new TF1("fGaus", "gaus", fitMin, fitMax);
    fGaus->SetLineColor(kRed+1);
    fGaus->SetLineWidth(2);
    // hPi0Mass->Fit(fGaus, "RQ");

    double maxBin = hNSigma->GetMaximum();
    double meanGuess = hNSigma->GetBinCenter(hNSigma->GetMaximumBin());
    double sigmaGuess = hNSigma->GetRMS();  // or ~5 MeV if you want
    fGaus->SetParameters(maxBin, meanGuess, sigmaGuess);
    hNSigma->Fit(fGaus, "RQ");

    double mean  = fGaus->GetParameter(1);
    double sigma = fGaus->GetParameter(2);
    double errMu = fGaus->GetParError(1);
    double errSi = fGaus->GetParError(2);

    hNSigma->SetLineColor(kBlack);
    hNSigma->SetLineWidth(2);
    hNSigma->GetXaxis()->SetTitle("n#sigma");
    hNSigma->GetYaxis()->SetTitle("Counts");

    hNSigma->Draw("HIST");
    fGaus->Draw("SAME");

    TLegend *leg = new TLegend(0.55, 0.7, 0.88, 0.88);
    leg->AddEntry(hNSigma, "n#sigma", "l");
    leg->AddEntry(fGaus, "Gaussian Fit:", "l");
    leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f MeV", mean, errMu), "");
    leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f MeV", sigma, errSi), "");
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);   // no fill
    leg->SetBorderSize(0);  // no border box
    leg->Draw();

    // TPaveText *info = new TPaveText(0.17, 0.70, 0.50, 0.90, "NDC");  // x1,y1,x2,y2 normalized coordinates
    // info->SetFillStyle(0);
    // info->SetBorderSize(0);
    // info->SetTextFont(42);
    // info->SetTextSize(0.04);
    // // info->AddText("Hibeam Wasafull simulation");
    // info->AddText("GEANT4 #pi^{0} sample");
    // info->AddText("1000 simulated events");
    // info->AddText("E_{kin} #in [1, 500] MeV");
    // // info->AddText("E_{kin} = 100 MeV");
    // info->Draw();

    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.16, 0.93, "#bf{Hibeam}  #it{Wasa full simulation}");

    c->SaveAs("plots/" + plotname);
    // c->SaveAs("Pi0Mass_truth.png");

    // clean up
    delete leg;
    delete fGaus;
    delete c;
}
