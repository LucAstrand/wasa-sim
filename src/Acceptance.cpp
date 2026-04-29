// #include "Acceptance.hpp"

// #include "TCanvas.h"
// #include "TLegend.h"
// #include "TStyle.h"
// #include "TGraphAsymmErrors.h"
// #include "TLatex.h"

// #include "PlotUtils.hpp"

// #include <iostream>

// Acceptance::Acceptance(const std::string& tag, int nBins, double xMin, double xMax)
//     : nBins_(nBins), xMin_(xMin), xMax_(xMax)
// {
//     h_num_ = new TH1D(Form("h_numPiAcceptance_%s", tag.c_str()), ";True #pi^{0} E_{kin};Events", nBins_, xMin_, xMax_);
//     h_den_ = new TH1D(Form("h_denPiAcceptance_%s", tag.c_str()), ";GPS #pi^{0} E_{kin};Events", nBins_, xMin_, xMax_);
// }

// void Acceptance::Pi0ProcessSignalEvent(const std::vector<TruePi0>& truePi0s, const std::vector<primaryPi0>& primaryPi0s) {

//     for (const auto& p : primaryPi0s) h_den_->Fill(p.p4.E() - p.p4.M());

//     std::unordered_map<int, double> genEkin;
//     genEkin.reserve(primaryPi0s.size());
//     for (const auto& p : primaryPi0s) genEkin.emplace(p.trackID, p.p4.E() - p.p4.M());

//     for (const auto& d : truePi0s) {
//         auto it = genEkin.find(d.trackID);
//         if (it != genEkin.end()) h_num_->Fill(it->second);
//     }
// }

// void Acceptance::ChPiProcessSignalEvent(const std::vector<TrueChPiInCal>& ChPiInCal, const std::vector<primaryChPi>& primaryChPis, int type) {
 
//     // w.r.t ekin 
//     // for (const auto& p : primaryChPis) h_den_->Fill(p.p4.E());
//     // std::unordered_map<int, double> genEkin;
//     // genEkin.reserve(primaryChPis.size());
//     // for (const auto& p : primaryChPis) genEkin.emplace(p.trackID, p.p4.E());

//     // for (const auto& d : ChPiInCal) {
//     //     if (type == 0) {
//     //         auto it = genEkin.find(d.trackID);
//     //         if (it != genEkin.end()) h_num_->Fill(it->second);
//     //     }
//     //     if (type == 1) {
//     //         if (!d.throughTPC) continue; // tricky logic here check not --> continue meaning we only do the ones that have d.throughTPC = true!!!
//     //         auto it = genEkin.find(d.trackID);
//     //         if (it != genEkin.end()) h_num_->Fill(it->second);
//     //     }
//     //     if (type == 2) {
//     //         if (d.throughTPC) continue;
//     //         auto it = genEkin.find(d.trackID);
//     //         if (it != genEkin.end()) h_num_->Fill(it->second);            
//     //     }
//     // }

//     // w.r.t theta
//     for (const auto& p : primaryChPis) h_den_->Fill(p.p4.Theta());
//     std::unordered_map<int, double> genTheta;
//     genTheta.reserve(primaryChPis.size());
//     for (const auto& p : primaryChPis) genTheta.emplace(p.trackID, p.p4.Theta());

//     for (const auto& d : ChPiInCal) {
//         if (type == 0) {
//             auto it = genTheta.find(d.trackID);
//             if (it != genTheta.end()) h_num_->Fill(it->second);
//         }
//         if (type == 1) {
//             if (!d.throughTPC) continue; // tricky logic here check not --> continue meaning we only do the ones that have d.throughTPC = true!!!
//             auto it = genTheta.find(d.trackID);
//             if (it != genTheta.end()) h_num_->Fill(it->second);
//         }
//         if (type == 2) {
//             if (d.throughTPC) continue;
//             auto it = genTheta.find(d.trackID);
//             if (it != genTheta.end()) h_num_->Fill(it->second);            
//         }
//     }
// }


// void Acceptance::FinalizePlot(const std::string& outFileName, PlotOptions opts) {


//     TGraphAsymmErrors* gEff = new TGraphAsymmErrors(h_num_, h_den_, "cl=0.683 b(1,1) mode"); // fraction in [0,1]

//     // SCALE TO PERCENT 
//     // Multiply y-values and Y-errors by 100 to convert fraction -> percent
//     const int nPoints = gEff->GetN();
//     for (int i = 0; i < nPoints; ++i) {
//         double x, y;
//         gEff->GetPoint(i, x, y);

//         double exl = gEff->GetErrorXlow(i);
//         double exh = gEff->GetErrorXhigh(i);
//         double eyl = gEff->GetErrorYlow(i);
//         double eyh = gEff->GetErrorYhigh(i);

//         double y_new = y * 100.0;
//         double eyl_new = eyl * 100.0;
//         double eyh_new = eyh * 100.0;

//         gEff->SetPoint(i, x, y_new);
//         gEff->SetPointError(i, exl, exh, eyl_new, eyh_new);
//     }

//     gEff->SetMarkerStyle(20);
//     gEff->SetMarkerSize(1.0);
//     gEff->SetMarkerColor(kBlack);
//     gEff->SetLineColor(kBlack);
//     gEff->GetYaxis()->SetTitle(opts.yAxisTitle.c_str());
//     gEff->GetXaxis()->SetTitle(opts.xAxisTitle.c_str());
//     gEff->GetYaxis()->SetRangeUser(0.0, 110.0);
//     gStyle->SetOptStat(0);
//     gEff->Draw("AP");

//     PlotGraph(gEff, outFileName.c_str(), opts);

//     // cleanup
//     delete gEff;
//     delete h_num_;
//     delete h_den_;
// }


#include "Acceptance.hpp"

#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TH2D.h"
#include "TColor.h"

#include "PlotUtils.hpp"
#include <iostream>
#include <unordered_map>

// ── 1D constructor ────────────────────────────────────────────────────────────
Acceptance::Acceptance(const std::string& tag,
                       AccAxisType axisType,
                       int nBins, double xMin, double xMax)
    : tag_(tag), axisType_(axisType),
      nBins_(nBins), xMin_(xMin), xMax_(xMax)
{
    h_num_ = new TH1D(Form("h_num_%s", tag.c_str()), "", nBins_, xMin_, xMax_);
    h_den_ = new TH1D(Form("h_den_%s", tag.c_str()), "", nBins_, xMin_, xMax_);
}

// ── 2D constructor ────────────────────────────────────────────────────────────
Acceptance::Acceptance(const std::string& tag,
                       int nBinsEkin, double ekinMin, double ekinMax,
                       int nBinsTheta, double thetaMin, double thetaMax)
    : tag_(tag), axisType_(AccAxisType::k2D),
      nBinsEkin_(nBinsEkin), ekinMin_(ekinMin), ekinMax_(ekinMax),
      nBinsTheta_(nBinsTheta), thetaMin_(thetaMin), thetaMax_(thetaMax)
{
    h2_num_ = new TH2D(Form("h2_num_%s", tag.c_str()), "",
                       nBinsEkin_, ekinMin_, ekinMax_,
                       nBinsTheta_, thetaMin_, thetaMax_);
    h2_den_ = new TH2D(Form("h2_den_%s", tag.c_str()), "",
                       nBinsEkin_, ekinMin_, ekinMax_,
                       nBinsTheta_, thetaMin_, thetaMax_);
}

// ── helpers to extract the variable we care about ────────────────────────────
namespace {
    double getVar(const TLorentzVector& p4, AccAxisType t) {
        return (t == AccAxisType::kEkin) ? (p4.E() - p4.M()) : p4.Theta();
    }
}

// ── Pi0 ───────────────────────────────────────────────────────────────────────
void Acceptance::Pi0ProcessSignalEvent(const std::vector<TruePi0>& truePi0s,
                                       const std::vector<primaryPi0>& primaryPi0s)
{
    if (axisType_ == AccAxisType::k2D) {
        // denominator
        for (const auto& p : primaryPi0s)
            h2_den_->Fill(p.p4.E() - p.p4.M(), p.p4.Theta());

        // build lookup: trackID -> primary
        std::unordered_map<int, const primaryPi0*> genMap;
        genMap.reserve(primaryPi0s.size());
        for (const auto& p : primaryPi0s) genMap.emplace(p.trackID, &p);

        for (const auto& d : truePi0s) {
            auto it = genMap.find(d.trackID);
            if (it != genMap.end())
                h2_num_->Fill(it->second->p4.E() - it->second->p4.M(),
                              it->second->p4.Theta());
        }
    } else {
        for (const auto& p : primaryPi0s)
            h_den_->Fill(getVar(p.p4, axisType_));

        std::unordered_map<int, double> genVar;
        genVar.reserve(primaryPi0s.size());
        for (const auto& p : primaryPi0s)
            genVar.emplace(p.trackID, getVar(p.p4, axisType_));

        for (const auto& d : truePi0s) {
            auto it = genVar.find(d.trackID);
            if (it != genVar.end()) h_num_->Fill(it->second);
        }
    }
}

// ── ChPi ──────────────────────────────────────────────────────────────────────
void Acceptance::ChPiProcessSignalEvent(const std::vector<TrueChPiInCal>& ChPiInCal,
                                        const std::vector<primaryChPi>& primaryChPis,
                                        int type)
{
    // Lambda to decide whether a hit passes the TPC condition for this type
    auto passesType = [&](const TrueChPiInCal& d) -> bool {
        if (type == 0) return true;
        if (type == 1) return d.throughTPC;
        if (type == 2) return !d.throughTPC;
        return false;
    };

    if (axisType_ == AccAxisType::k2D) {
        for (const auto& p : primaryChPis)
            h2_den_->Fill(p.p4.E() - p.p4.M(), p.p4.Theta());

        std::unordered_map<int, const primaryChPi*> genMap;
        genMap.reserve(primaryChPis.size());
        for (const auto& p : primaryChPis) genMap.emplace(p.trackID, &p);

        for (const auto& d : ChPiInCal) {
            if (!passesType(d)) continue;
            auto it = genMap.find(d.trackID);
            if (it != genMap.end())
                h2_num_->Fill(it->second->p4.E() - it->second->p4.M(),
                              it->second->p4.Theta());
        }
    } else if (axisType_ == AccAxisType::kTracks) {
        for (const auto& p : primaryChPis)
            h_den_->Fill(primaryChPis.size());
        std::unordered_map<int,int> genVar;
        genVar.reserve(primaryChPis.size());
        for (const auto& p : primaryChPis)
            genVar.emplace(p.trackID, primaryChPis.size());
        for (const auto& d: ChPiInCal) {
            if (!passesType(d)) continue;
            auto it = genVar.find(d.trackID);
            if (it != genVar.end()) h_num_->Fill(it->second);
        }

    } else {
        for (const auto& p : primaryChPis)
            h_den_->Fill(getVar(p.p4, axisType_));

        std::unordered_map<int, double> genVar;
        genVar.reserve(primaryChPis.size());
        for (const auto& p : primaryChPis)
            genVar.emplace(p.trackID, getVar(p.p4, axisType_));

        for (const auto& d : ChPiInCal) {
            if (!passesType(d)) continue;
            auto it = genVar.find(d.trackID);
            if (it != genVar.end()) h_num_->Fill(it->second);
        }
    }
}

// ── finalize ──────────────────────────────────────────────────────────────────
// void Acceptance::FinalizePlot(const std::string& outFileName, PlotOptions opts)
// {
//     gStyle->SetOptStat(0);

//     if (axisType_ == AccAxisType::k2D) {
//         SetPrettyStyle();
//         // Build acceptance TH2D by dividing bin-by-bin
//         TH2D* h2_acc = static_cast<TH2D*>(h2_num_->Clone(
//                             Form("h2_acc_%s", tag_.c_str())));
//         h2_acc->Divide(h2_den_);          // gives fraction [0,1]
//         h2_acc->Scale(100.0);             // convert to percent

//         h2_acc->GetXaxis()->SetTitle(opts.xAxisTitle.c_str());  // ekin
//         h2_acc->GetYaxis()->SetTitle(opts.yAxisTitle.c_str());  // theta
//         h2_acc->GetZaxis()->SetTitle("Acceptance (%)");
//         h2_acc->GetZaxis()->SetRangeUser(0.0, 100.0);

//         // Nice colour palette: kBird (default) or try kRainBow, kViridis
//         gStyle->SetPalette(kBird);
//         gStyle->SetNumberContours(99);

//         TCanvas* c = new TCanvas(Form("c2d_%s", tag_.c_str()), "", 800, 600);
//         c->SetRightMargin(0.15);   // room for the colour bar
//         h2_acc->Draw("COLZ");
//         std::cout << outFileName.c_str() << std::endl;
//         c->SaveAs(outFileName.c_str());

//         delete c;
//         delete h2_acc;
//         delete h2_num_;
//         delete h2_den_;

//     } else {
//         // ── original 1D path ──────────────────────────────────────────────
//         TGraphAsymmErrors* gEff = new TGraphAsymmErrors(
//                                       h_num_, h_den_, "cl=0.683 b(1,1) mode");

//         const int nPoints = gEff->GetN();
//         for (int i = 0; i < nPoints; ++i) {
//             double x, y;
//             gEff->GetPoint(i, x, y);
//             gEff->SetPoint(i, x, y * 100.0);
//             gEff->SetPointError(i,
//                 gEff->GetErrorXlow(i),  gEff->GetErrorXhigh(i),
//                 gEff->GetErrorYlow(i)  * 100.0,
//                 gEff->GetErrorYhigh(i) * 100.0);
//         }

//         gEff->SetMarkerStyle(20);
//         gEff->SetMarkerSize(1.0);
//         gEff->SetMarkerColor(kBlack);
//         gEff->SetLineColor(kBlack);
//         gEff->GetXaxis()->SetTitle(opts.xAxisTitle.c_str());
//         gEff->GetYaxis()->SetTitle(opts.yAxisTitle.c_str());
//         gEff->GetYaxis()->SetRangeUser(0.0, 110.0);
//         gEff->Draw("AP");

//         PlotGraph(gEff, outFileName.c_str(), opts);

//         delete gEff;
//         delete h_num_;
//         delete h_den_;
//     }
// }

void Acceptance::FinalizePlot(const std::string& outFileName, PlotOptions opts)
{
    gStyle->SetOptStat(0);

    if (axisType_ == AccAxisType::k2D) {
        TH2D* h2_acc = static_cast<TH2D*>(h2_num_->Clone(
                            Form("h2_acc_%s", tag_.c_str())));
        h2_acc->Divide(h2_den_);
        h2_acc->Scale(100.0);
        h2_acc->GetZaxis()->SetRangeUser(0.0, 100.0);

        // Set axis titles on the histogram itself so Plot2D picks them up
        h2_acc->GetXaxis()->SetTitle(opts.xAxisTitle.c_str());
        h2_acc->GetYaxis()->SetTitle(opts.yAxisTitle.c_str());
        h2_acc->GetZaxis()->SetTitle("Acceptance (%)");

        // Reuse Plot2D — it calls SavePlot which prepends "plots/" correctly
        opts.drawOption  = "COLZ";
        opts.isHeatmap   = true;
        opts.colorMap    = kBird;
        Plot2D(h2_acc, outFileName, opts);

        delete h2_acc;
        delete h2_num_;
        delete h2_den_;

    } else {
        // 1D path — unchanged, PlotGraph already calls SavePlot correctly
        TGraphAsymmErrors* gEff = new TGraphAsymmErrors(
                                      h_num_, h_den_, "cl=0.683 b(1,1) mode");
        const int nPoints = gEff->GetN();
        for (int i = 0; i < nPoints; ++i) {
            double x, y;
            gEff->GetPoint(i, x, y);
            gEff->SetPoint(i, x, y * 100.0);
            gEff->SetPointError(i,
                gEff->GetErrorXlow(i),  gEff->GetErrorXhigh(i),
                gEff->GetErrorYlow(i)  * 100.0,
                gEff->GetErrorYhigh(i) * 100.0);
        }

        gEff->SetMarkerStyle(20);
        gEff->SetMarkerSize(1.0);
        gEff->SetMarkerColor(kBlack);
        gEff->SetLineColor(kBlack);
        gEff->GetXaxis()->SetTitle(opts.xAxisTitle.c_str());
        gEff->GetYaxis()->SetTitle(opts.yAxisTitle.c_str());
        gEff->GetYaxis()->SetRangeUser(0.0, 110.0);

        PlotGraph(gEff, outFileName, opts);  // SavePlot called inside

        delete gEff;
        delete h_num_;
        delete h_den_;
    }
}