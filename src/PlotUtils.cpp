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
#include <algorithm> 
#include <iostream>   


void SetPrettyStyle() {
    std::filesystem::create_directories("plots");
    std::filesystem::create_directories("plots/Neutral");
    std::filesystem::create_directories("plots/Charged");
    std::filesystem::create_directories("plots/Truth");
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
    gStyle->SetTextFont(42);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
}

std::unique_ptr<TCanvas> PlotCreateCanvas(const std::string& name, int width, int height) {
    auto c = std::make_unique<TCanvas>(name.c_str(), name.c_str(), width, height);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);
    return c;
}

void AddTopLatex(TCanvas* c, const std::string& text) {
    TLatex l;
    l.SetNDC();
    l.SetTextFont(42);
    l.SetTextSize(0.045);
    l.DrawLatex(0.16, 0.93, text.c_str());
}

// std::unique_ptr<TPaveText> PlotCreateInfoPave(const std::vector<std::string>& lines, double x1, double y1, double x2, double y2) {
TPaveText* PlotCreateInfoPave(const std::vector<std::string>& lines, double x1, double y1, double x2, double y2) {
    // auto info = std::make_unique<TPaveText>(x1, y1, x2, y2, "NDC");
    TPaveText *info = new TPaveText(x1, y1, x2, y2, "NDC");

    info->SetFillStyle(0);
    info->SetBorderSize(0);
    info->SetTextFont(42);
    info->SetTextSize(0.04);
    for (const auto& line : lines) {
        info->AddText(line.c_str());
    }
    return info;
}

// std::unique_ptr<TLegend> PlotCreateLegend(const std::vector<std::string>& entries, const std::vector<std::string>& extraLines, double x1, double y1, double x2, double y2, const std::vector<TObject*>& objects) {
//     auto leg = std::make_unique<TLegend>(x1, y1, x2, y2);
//     leg->SetTextSize(0.03);
//     leg->SetFillStyle(0);
//     leg->SetBorderSize(0);
//     for (size_t i = 0; i < entries.size(); ++i) {
//         leg->AddEntry(objects.empty() ? nullptr : objects[i], entries[i].c_str(), "l");
//     }
//     for (const auto& line : extraLines) {
//         leg->AddEntry((TObject*)0, line.c_str(), "");
//     }
//     return leg;
// }

std::unique_ptr<TLegend> PlotCreateLegend(const std::vector<std::string>& entries, const std::vector<std::string>& extraLines, double x1, double y1, double x2, double y2, const std::vector<TObject*>& objects, const std::string legendPlotOpts) {
    auto leg = std::make_unique<TLegend>(x1, y1, x2, y2);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    for (size_t i = 0; i < entries.size(); ++i) {
        leg->AddEntry(objects.empty() ? nullptr : objects[i], entries[i].c_str(), legendPlotOpts.c_str());
    }
    for (const auto& line : extraLines) {
        leg->AddEntry((TObject*)0, line.c_str(), "");
    }
    return leg;
}

// void PerformFitAndAddToLegend(TH1* hist, TLegend* leg, const PlotOptions& options) {
//     double fitMin = (options.fitMin == -999) ? hist->GetXaxis()->GetXmin() : options.fitMin;
//     double fitMax = (options.fitMax == -999) ? hist->GetXaxis()->GetXmax() : options.fitMax;
//     auto f = std::make_unique<TF1>("fit", options.fitFunction.c_str(), fitMin, fitMax);
//     f->SetLineColor(kRed + 1);
//     f->SetLineWidth(2);
//     double maxBin = hist->GetMaximum();
//     double meanGuess = hist->GetBinCenter(hist->GetMaximumBin());
//     double sigmaGuess = hist->GetRMS();
//     f->SetParameters(maxBin, meanGuess, sigmaGuess);
//     hist->Fit(f.get(), "RQ");  // Quiet, range
//     double mean = f->GetParameter(1);
//     double sigma = f->GetParameter(2);
//     double errMu = f->GetParError(1);
//     double errSi = f->GetParError(2);
//     f->Draw("SAME");  // Draw fit on current canvas
//     if (leg) {
//         leg->AddEntry(f.get(), "Gaussian Fit:", "l");
//         leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f", mean, errMu), "");
//         leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f", sigma, errSi), "");
//     }
// }

void PerformFitAndAddToLegend(TH1* hist, TLegend* leg, const PlotOptions& options)
{
    if (!hist || hist->GetEntries() < 5) return;

    double fitMin = (options.fitMin == -999)
                    ? hist->GetXaxis()->GetXmin()
                    : options.fitMin;
    double fitMax = (options.fitMax == -999)
                    ? hist->GetXaxis()->GetXmax()
                    : options.fitMax;

    // ROOT must own this
    TF1* f = new TF1(
        Form("fit_%s", hist->GetName()),
        options.fitFunction.c_str(),
        fitMin,
        fitMax
    );

    f->SetLineColor(kRed + 1);
    f->SetLineWidth(2);

    double maxBin = hist->GetMaximum();
    double meanGuess = hist->GetBinCenter(hist->GetMaximumBin());
    double sigmaGuess = hist->GetRMS();
    f->SetParameters(maxBin, meanGuess, sigmaGuess);

    // "S" stores fit result
    // "R" respects range
    // "Q" quiet
    hist->Fit(f, "SRQ");

    double mean = f->GetParameter(1);
    double sigma = f->GetParameter(2);
    double errMu = f->GetParError(1);
    double errSi = f->GetParError(2);

    f->Draw("SAME");

    if (leg) {
        leg->AddEntry(f, "Gaussian Fit:", "l");
        leg->AddEntry((TObject*)0, Form("#mu = %.1f #pm %.1f", mean, errMu), "");
        leg->AddEntry((TObject*)0, Form("#sigma = %.1f #pm %.1f", sigma, errSi), "");
    }
}


void SavePlot(TCanvas* c, const std::string& plotname) {
    c->SaveAs(("plots/" + plotname).c_str());
}

void Plot1D(const std::vector<TH1*>& hists, const std::vector<int>& colors, const std::string& plotname, const PlotOptions& options) {
    if (hists.empty() || hists[0] == nullptr) {
        std::cerr << "Error: No valid histograms provided." << std::endl;
        return;
    }
    SetPrettyStyle();
    gStyle->SetOptStat(options.showStats ? 1 : 0);
    auto c = PlotCreateCanvas("c1D_" + plotname);
    // Normalize if requested
    if (options.normalizeHists) {
        for (auto* h : hists) {
            if (h->Integral() > 0) h->Scale(1.0 / h->Integral());
        }
    }
    // Find max Y for scaling
    double maxY = 0;
    for (auto* h : hists) {
        maxY = std::max(maxY, h->GetMaximum());
    }
    // Draw first hist and set props
    hists[0]->SetLineColor(colors.empty() ? kBlack : colors[0]);
    hists[0]->SetLineWidth(2);
    hists[0]->SetMaximum(1.2 * maxY);  // Headroom
    hists[0]->Draw(options.drawOption.c_str());
    // Overlay others
    for (size_t i = 1; i < hists.size(); ++i) {
        hists[i]->SetLineColor(colors.size() > i ? colors[i] : kBlack);
        hists[i]->SetLineWidth(2);
        hists[i]->Draw(("SAME " + options.drawOption).c_str());
    }
    // Fit (on first hist only, for simplicity)
    std::unique_ptr<TLegend> leg;
    if (options.addLegend) {
        if (!options.legendEntries.empty() &&
            options.legendEntries.size() != hists.size()) {
            std::cerr << "Legend entries size mismatch!" << std::endl;
            return;
        }
        std::vector<TObject*> objs(hists.begin(), hists.end());
        leg = PlotCreateLegend(options.legendEntries, options.extraLegendLines, options.legendX1, options.legendY1, options.legendX2, options.legendY2, objs, options.legendDrawOpt);
    }
    if (options.doFit) {
        PerformFitAndAddToLegend(hists[0], leg.get(), options);
    }
    if (options.addLegend) leg->Draw();
    // Add info and top latex
    if (options.addInfoPave) {
        auto info = PlotCreateInfoPave(options.infoLines, options.infoX1, options.infoY1, options.infoX2, options.infoY2);
        info->Draw();
    }
    if (options.addTopLatex) AddTopLatex(c.get(), options.topLatex);
    SavePlot(c.get(), plotname);
}

void Plot2D(TH2F* hist, const std::string& plotname, const PlotOptions& options) {
    if (hist == nullptr) {
        std::cerr << "Error: No valid histogram provided." << std::endl;
        return;
    }
    SetPrettyStyle();
    gStyle->SetOptStat(options.showStats ? 1 : 0);
    auto c = PlotCreateCanvas("c2D_" + plotname);
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    hist->Draw(options.drawOption.c_str());  // e.g., "COLZ" for color map if you change default
    std::unique_ptr<TLegend> leg;
    if (options.addLegend) {
        leg = PlotCreateLegend(options.legendEntries, options.extraLegendLines, options.legendX1, options.legendY1, options.legendX2, options.legendY2, {} ,options.legendDrawOpt);
        leg->Draw();
    }
    if (options.addInfoPave) {
        auto info = PlotCreateInfoPave(options.infoLines, options.infoX1, options.infoY1, options.infoX2, options.infoY2);
        info->Draw();
    }
    if (options.addTopLatex) AddTopLatex(c.get(), options.topLatex);
    SavePlot(c.get(), plotname);
}

void Plot2DOverlay(
    const std::vector<TH2*>& hists,
    const std::vector<int>& colors,
    const std::string& plotname,
    const PlotOptions& options)
{
    if (hists.empty() || !hists[0]) return;

    SetPrettyStyle();
    auto c = PlotCreateCanvas("c2DOverlay_" + plotname);

    for (size_t i = 0; i < hists.size(); ++i) {
        auto* h = hists[i];
        h->SetMarkerColor(colors.size() > i ? colors[i] : kBlack);
        h->SetMarkerStyle(20 + i);
        h->SetMarkerSize(0.6);

        std::string drawOpt = (i == 0)
            ? options.drawOption
            : "SAME " + options.drawOption;

        h->Draw(drawOpt.c_str());
    }

    std::unique_ptr<TLegend> leg;
    if (options.addLegend) {
        std::vector<TObject*> objs(hists.begin(), hists.end());
        leg = PlotCreateLegend(options.legendEntries, options.extraLegendLines, options.legendX1, options.legendY1, options.legendX2, options.legendY2, objs, options.legendDrawOpt);
        leg->Draw();
    }

    if (options.addTopLatex)
        AddTopLatex(c.get(), options.topLatex);

    SavePlot(c.get(), plotname);
}


void PlotGraph(TGraph* graph, const std::string& plotname, const PlotOptions& options)
{
    if (!graph) {
        std::cerr << "Error: null graph passed to PlotGraph\n";
        return;
    }

    SetPrettyStyle();
    gStyle->SetOptStat(0);

    auto c = PlotCreateCanvas("cGraph_" + plotname);

    graph->SetLineWidth(2);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.0);

    graph->Draw("AP");

    std::unique_ptr<TLegend> leg;
    if (options.addLegend && !options.legendEntries.empty()) {
        leg = PlotCreateLegend(
            options.legendEntries,
            options.extraLegendLines,
            options.legendX1,
            options.legendY1,
            options.legendX2,
            options.legendY2,
            { graph },
            options.legendDrawOpt
        );
        leg->Draw();
    }

    if (options.addInfoPave) {
        auto info = PlotCreateInfoPave(
            options.infoLines,
            options.infoX1,
            options.infoY1,
            options.infoX2,
            options.infoY2
        );
        info->Draw();
    }

    if (options.addTopLatex)
        AddTopLatex(c.get(), options.topLatex);

    SavePlot(c.get(), plotname);
}
