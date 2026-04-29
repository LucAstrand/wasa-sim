#include "ChEfficiency.hpp"

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
ChEfficiency::ChEfficiency(const std::string& tag,
                       AccAxisTypeChPi axisType,
                       int nBins, double xMin, double xMax)
    : tag_(tag), axisType_(axisType),
      nBins_(nBins), xMin_(xMin), xMax_(xMax)
{
    h_num_ = new TH1D(Form("h_num_%s", tag.c_str()), "", nBins_, xMin_, xMax_);
    h_den_ = new TH1D(Form("h_den_%s", tag.c_str()), "", nBins_, xMin_, xMax_);
}

// ── 2D constructor ────────────────────────────────────────────────────────────
ChEfficiency::ChEfficiency(const std::string& tag,
                       int nBinsEkin, double ekinMin, double ekinMax,
                       int nBinsTheta, double thetaMin, double thetaMax)
    : tag_(tag), axisType_(AccAxisTypeChPi::k2D),
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
    double getVar(const TLorentzVector& p4, AccAxisTypeChPi t) {
        return (t == AccAxisTypeChPi::kEkin) ? (p4.E() - p4.M()) : p4.Theta();
    }
}

// ── ChPi ──────────────────────────────────────────────────────────────────────
void ChEfficiency::ChPiProcessSignalEvent(
    const std::vector<TrueChPiInCal>& ChPiInCal,
    const std::vector<primaryChPi>& primaryChPis,
    const std::vector<ChargedCluster>& ChClusters)
{
    // --- Map: trackID -> generator p4 ---
    std::unordered_map<int, const primaryChPi*> p4Map;
    p4Map.reserve(primaryChPis.size());
    for (const auto& p : primaryChPis) {
        p4Map.emplace(p.trackID, &p);
    }

    // --- Set: reconstructed pion trackIDs (from reco objects) ---
    std::unordered_set<int> reconstructedIDs;
    reconstructedIDs.reserve(ChClusters.size());
    for (const auto& c : ChClusters) {
        if (std::abs(c.objectTruePDG) == 211) {
            reconstructedIDs.insert(c.trackID);
        }
    }

    // --- Precompute accepted truth pions (throughTPC) ---
    std::vector<const TrueChPiInCal*> accepted;
    accepted.reserve(ChPiInCal.size());

    for (const auto& d : ChPiInCal) {
        // if (d.throughTPC) {
            accepted.push_back(&d);
        // }
    }

    // ============================
    // === Fill histograms
    // ============================

    if (axisType_ == AccAxisTypeChPi::k2D) {

        for (const auto* d : accepted) {

            auto it = p4Map.find(d->trackID);
            if (it == p4Map.end()) continue;

            const auto& p4 = it->second->p4;

            double ekin  = p4.E() - p4.M();
            double theta = p4.Theta();

            h2_den_->Fill(ekin, theta);

            if (reconstructedIDs.count(d->trackID)) {
                h2_num_->Fill(ekin, theta);
            }
        }

    } else if (axisType_ == AccAxisTypeChPi::kTracks) {

        // multiplicity = number of accepted truth pions
        double mult = static_cast<double>(accepted.size());

        for (const auto* d : accepted) {

            h_den_->Fill(mult);

            if (reconstructedIDs.count(d->trackID)) {
                h_num_->Fill(mult);
            }
        }

    } else {

        // kEkin or kTheta
        for (const auto* d : accepted) {

            auto it = p4Map.find(d->trackID);
            if (it == p4Map.end()) continue;

            double val = getVar(it->second->p4, axisType_);

            h_den_->Fill(val);

            if (reconstructedIDs.count(d->trackID)) {
                h_num_->Fill(val);
            }
        }
    }
}


// void ChEfficiency::ChPiProcessSignalEvent(const std::vector<TrueChPiInCal>& ChPiInCal,
//                                           const std::vector<primaryChPi>& primaryChPis,
//                                           const std::vector<ChargedCluster>& ChClusters)
// {
//     // Build lookup: trackID -> primary (for p4 info) — used by k2D, kEkin, kTheta
//     std::unordered_map<int, const primaryChPi*> p4Map;
//     p4Map.reserve(primaryChPis.size());
//     for (const auto& p : primaryChPis) p4Map.emplace(p.trackID, &p);

//     // Build lookup: which trackIDs were reconstructed
//     std::unordered_set<int> reconstructedIDs;
//     reconstructedIDs.reserve(ChClusters.size());
//     for (const auto& c : ChClusters) {
//         if (c.objectTruePDG == 211 || c.objectTruePDG == -211)
//             reconstructedIDs.insert(c.trackID);
//     }
//     if (axisType_ == AccAxisTypeChPi::k2D) {
//         // den: all pions in acceptance, num: those also reconstructed
//         for (const auto& d : ChPiInCal) {
//             if (!d.throughTPC) continue;
//             auto it = p4Map.find(d.trackID);
//             if (it == p4Map.end()) continue;  // no p4 info, skip
//             double ekin  = it->second->p4.E() - it->second->p4.M();
//             double theta = it->second->p4.Theta();
//             h2_den_->Fill(ekin, theta);
//             if (reconstructedIDs.count(d.trackID))
//                 h2_num_->Fill(ekin, theta);
//         }

//     } else if (axisType_ == AccAxisTypeChPi::kTracks) {
//         // x-axis = multiplicity of pions in acceptance for this event
//         // fill den once per pion-in-acceptance, num once per reconstructed pion
//         // --> efficiency as function of event multiplicity
//         double mult = static_cast<double>(ChPiInCal.size());
//         for (const auto& d : ChPiInCal) {
//             if (!d.throughTPC) continue;
//             h_den_->Fill(mult);
//             if (reconstructedIDs.count(d.trackID))
//                 h_num_->Fill(mult);
//         }

//     } else {
//         // kEkin or kTheta: den = in acceptance, num = reconstructed
//         for (const auto& d : ChPiInCal) {
//             if (!d.throughTPC) continue;
//             auto it = p4Map.find(d.trackID);
//             if (it == p4Map.end()) continue;
//             double val = getVar(it->second->p4, axisType_);
//             h_den_->Fill(val);
//             if (reconstructedIDs.count(d.trackID))
//                 h_num_->Fill(val);
//         }
//     }
// }

void ChEfficiency::FinalizePlot(const std::string& outFileName, PlotOptions opts)
{
    gStyle->SetOptStat(0);

    if (axisType_ == AccAxisTypeChPi::k2D) {
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

// void ChEfficiency::FinalizePlot(const std::string& outFileName, PlotOptions opts)
// {
//     gStyle->SetOptStat(0);

//     std::cout << "=== FinalizePlot debug: " << tag_ << " ===" << std::endl;
//     for (int i = 1; i <= h_num_->GetNbinsX(); ++i) {
//         double num = h_num_->GetBinContent(i);
//         double den = h_den_->GetBinContent(i);
//         double center = h_den_->GetBinCenter(i);
//         std::cout << "  bin " << i
//                   << "  x=" << center
//                   << "  num=" << num
//                   << "  den=" << den
//                   << "  ratio=" << (den > 0 ? num/den : -1)
//                   << std::endl;
//     }

//     if (axisType_ == AccAxisTypeChPi::k2D) {
//         TH2D* h2_acc = static_cast<TH2D*>(h2_num_->Clone(
//                             Form("h2_acc_%s", tag_.c_str())));
//         h2_acc->Divide(h2_den_);
//         h2_acc->Scale(100.0);

//         // Suppress bins with no denominator statistics
//         for (int bx = 1; bx <= h2_acc->GetNbinsX(); ++bx)
//             for (int by = 1; by <= h2_acc->GetNbinsY(); ++by)
//                 if (h2_den_->GetBinContent(bx, by) == 0.0)
//                     h2_acc->SetBinContent(bx, by, 0.0);

//         h2_acc->GetZaxis()->SetRangeUser(0.0, 100.0);
//         h2_acc->GetXaxis()->SetTitle(opts.xAxisTitle.c_str());
//         h2_acc->GetYaxis()->SetTitle(opts.yAxisTitle.c_str());
//         h2_acc->GetZaxis()->SetTitle("Efficiency (%)");

//         opts.drawOption = "COLZ";
//         opts.isHeatmap  = true;
//         opts.colorMap   = kBird;
//         Plot2D(h2_acc, outFileName, opts);

//         delete h2_acc;
//         delete h2_num_;
//         delete h2_den_;

//     } else {
//         TGraphAsymmErrors* gEff = new TGraphAsymmErrors(
//                                       h_num_, h_den_, "cl=0.683 b(1,1) mode");
//         const int nPoints = gEff->GetN();
//         for (int i = 0; i < nPoints; ++i) {
//             double x, y;
//             gEff->GetPoint(i, x, y);

//             // Suppress points where denominator is empty —
//             // Bayesian prior would give ~50% for 0/0 bins otherwise
//             int bin = h_den_->FindBin(x);
//             if (h_den_->GetBinContent(bin) == 0.0) {
//                 gEff->SetPoint(i, x, -999.0);
//                 gEff->SetPointError(i, 0, 0, 0, 0);
//                 continue;
//             }

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
//         gEff->GetYaxis()->SetRangeUser(0.0, 110.0);  // -999 points invisible

//         PlotGraph(gEff, outFileName, opts);

//         delete gEff;
//         delete h_num_;
//         delete h_den_;
//     }
// }