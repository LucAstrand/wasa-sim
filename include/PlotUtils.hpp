#ifndef PLOTUTILS_H
#define PLOTUTILS_H

#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"  // For future if needed
#include "TEfficiency.h"  // For future if needed
#include <vector>
#include <string>
#include <memory>  // For unique_ptr
#include "TCanvas.h"   
#include "TPaveText.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TObject.h"

struct PlotOptions {
    bool doFit = false;
    std::string fitFunction = "gaus";  // e.g., "gaus", "pol2"
    double fitMin = -999;  // -999 means use hist range
    double fitMax = -999;
    bool showStats = false;  // For gStyle->SetOptStat
    std::vector<std::string> legendEntries;  // One per hist/object
    std::vector<std::string> extraLegendLines;  // e.g., "2 photon events = X"
    std::vector<std::string> infoLines;  // e.g., {"GEANT4 pi0 sample", "1000 events"}
    std::string topLatex = "#bf{Hibeam} #it{Wasa full simulation}";
    bool addTopLatex = true;
    bool addInfoPave = true;
    bool addLegend = true;
    std::string drawOption = "HIST";  // e.g., "HIST", "E", "COLZ" for 2D
    bool normalizeHists = false;  // Scale to unit area if true (for overlays)
    double legendX1 = 0.55, legendY1 = 0.7, legendX2 = 0.88, legendY2 = 0.88;  // Customizable position
    double infoX1 = 0.17, infoY1 = 0.70, infoX2 = 0.50, infoY2 = 0.90;
};

// Core functions
void Plot1D(const std::vector<TH1*>& hists, const std::vector<int>& colors, const std::string& plotname, const PlotOptions& options = PlotOptions());
void Plot2D(TH2F* hist, const std::string& plotname, const PlotOptions& options = PlotOptions());
void PlotGraph(TGraph* graph, const std::string& plotname, const PlotOptions& options);


// Helpers (exposed if you want to use standalone, else make static/private)
void SetPrettyStyle();
std::unique_ptr<TCanvas> PlotCreateCanvas(const std::string& name, int width = 800, int height = 600);
void AddTopLatex(TCanvas* c, const std::string& text);
std::unique_ptr<TPaveText> PlotCreateInfoPave(const std::vector<std::string>& lines, double x1, double y1, double x2, double y2);
std::unique_ptr<TLegend> PlotCreateLegend(const std::vector<std::string>& entries, const std::vector<std::string>& extraLines, double x1, double y1, double x2, double y2, const std::vector<TObject*>& objects = {});
void PerformFitAndAddToLegend(TH1* hist, TLegend* leg, const PlotOptions& options);
void SavePlot(TCanvas* c, const std::string& plotname);

#endif