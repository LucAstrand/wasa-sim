#ifndef SELECTIONHISTOGRAMS_H
#define SELECTIONHISTOGRAMS_H

#include "TH1F.h"
#include "TH2F.h"
#include "EventVariables.hpp"

#include <vector>


struct SelectionHistograms {
    TH1F* hNCharged       = nullptr;
    TH1F* hNNeutral       = nullptr;
    TH1F* hEtotal         = nullptr;
    TH1F* hECorrected     = nullptr;
    TH1F* hSphericity     = nullptr;
    TH1F* hMaxTrackAngle  = nullptr;
    TH1F* hVertexRadius   = nullptr;
    TH1F* hNPi0           = nullptr;

    // Container for cutflow
    std::vector<EventVariables> collectedEvents;

    void Book(const std::string& prefix) {
        hNCharged      = new TH1F((prefix+"_nCharged").c_str(),
                            ";N charged tracks;Events", 10, -0.5, 9.5);
        hNNeutral      = new TH1F((prefix+"_nNeutral").c_str(),
                            ";N neutral clusters;Events", 10, -0.5, 9.5);
        hEtotal        = new TH1F((prefix+"_Etotal").c_str(),
                            ";Total reco energy [MeV];Events", 100, 0, 3000);
        hECorrected    = new TH1F((prefix+"hECorrected").c_str(),
                            ";Total Corrected Energy [MeV];Events", 100, 0, 3000);
        hSphericity    = new TH1F((prefix+"_sphericity").c_str(),
                            ";Sphericity;Events", 50, 0, 1);
        hMaxTrackAngle = new TH1F((prefix+"_maxAngle").c_str(),
                            ";Max track opening angle [rad];Events", 60, 0, TMath::Pi());
        hVertexRadius  = new TH1F((prefix+"_vtxR").c_str(),
                            ";Vertex radius [cm];Events", 50, 0, 20);
        hNPi0          = new TH1F((prefix+"_nPi0").c_str(),
                            ";N #pi^{0} candidates;Events", 8, -0.5, 7.5);
    }

    void Fill(const EventVariables& ev) {
        // Fill Hists
        hNCharged      ->Fill(ev.nChargedTracks);
        hNNeutral      ->Fill(ev.nNeutralClusters);
        hEtotal        ->Fill(ev.totalRecoEnergy);
        hECorrected    ->Fill(ev.correctedEnergy);
        hSphericity    ->Fill(ev.sphericity);
        hMaxTrackAngle ->Fill(ev.maxTrackAngle);
        hVertexRadius  ->Fill(ev.vertexRadius);
        hNPi0          ->Fill(ev.nPi0Candidates);
        // Fill Container
        collectedEvents.push_back(ev);
    }   
    
    void PlotOverlay(const SelectionHistograms& bkg, 
                     const std::string& outDir) {
        PlotOptions opts;
        opts.addLegend = true;
        opts.legendEntries = {"Signal", "Background (cosmic)"};
    
        auto overlay = [&](TH1F* hS, TH1F* hB, const std::string& name) {
            if (!hS || !hB) return;
            TH1F* hSn = (TH1F*)hS->Clone(); hSn->Scale(1.0/hS->Integral());
            TH1F* hBn = (TH1F*)hB->Clone(); hBn->Scale(1.0/hB->Integral());
            Plot1D({hSn, hBn}, {kBlue+1, kRed+1}, outDir + name, opts);
            delete hSn; delete hBn;
        };
    
        overlay(hNCharged,      bkg.hNCharged,       "nCharged.png");
        overlay(hNNeutral,      bkg.hNNeutral,       "nNeutral.png");
        overlay(hEtotal,        bkg.hEtotal,         "Etotal.png");
        overlay(hECorrected,    bkg.hECorrected,     "ECorrected.png");
        overlay(hSphericity,    bkg.hSphericity,     "Sphericity.png");
        overlay(hMaxTrackAngle, bkg.hMaxTrackAngle,  "MaxTrackAngle.png");
        overlay(hVertexRadius,  bkg.hVertexRadius,   "VertexRadius.png");
        overlay(hNPi0,          bkg.hNPi0,           "NPi0.png");
    }
};

struct CorrelationMatrix {

    TH2F* hCorrelationMatrix = nullptr;
    std::vector<std::string> varNames = {
        "nCharged", "nNeutral", "Etotal", "Ecorrected",
        "sphericity", "maxTrackAngle", "vertexRadius", "nPi0"
    };
    std::vector<std::vector<double>> varValues;

    void BookCorrelation(const std::string& prefix) {
        int nVars = varNames.size();
        hCorrelationMatrix = new TH2F((prefix + "CorrelMatrix").c_str(), ";Variable;Variable", nVars,0,nVars, nVars,0,nVars);
        for (int i=0; i<nVars; ++i) {
            hCorrelationMatrix->GetXaxis()->SetBinLabel(i+1, varNames[i].c_str());
            hCorrelationMatrix->GetYaxis()->SetBinLabel(i+1, varNames[i].c_str());
        }
        varValues.resize(nVars);
        // hCorrelationMatrix->GetXaxis()->SetLabelOffset(-0.01);
        // hCorrelationMatrix->GetYaxis()->SetLabelOffset(0.001);
        // hCorrelationMatrix->GetXaxis()->LabelsOption("c");  // c = centered
        // hCorrelationMatrix->GetYaxis()->LabelsOption("c");
    }

    void FillCorrelation(const EventVariables& ev) {
        varValues[0].push_back(ev.nChargedTracks);
        varValues[1].push_back(ev.nNeutralClusters);
        varValues[2].push_back(ev.totalRecoEnergy);
        varValues[3].push_back(ev.correctedEnergy);
        varValues[4].push_back(ev.sphericity);
        varValues[5].push_back(ev.maxTrackAngle);
        varValues[6].push_back(ev.vertexRadius);
        varValues[7].push_back(ev.nPi0Candidates);
    }

    void ComputeAndPlotCorrelations(const std::string& outDir, const std::string& label) {

        int nVars = varNames.size();
        int nEvents = varValues[0].size();
        if (nEvents == 0) return;

        // Pearson correlation coefficient
        for (int i=0; i<nVars; ++i) {
            for(int j=0; j<nVars; ++j) {
                double meanI = 0, meanJ = 0;
                for (int e=0; e<nEvents; ++e) {
                    meanI += varValues[i][e];
                    meanJ += varValues[j][e];
                }
                meanI /= nEvents;
                meanJ /= nEvents;

                double cov = 0, varI = 0, varJ = 0;
                for (int e=0; e<nEvents; ++e) {
                    double di = varValues[i][e] - meanI;
                    double dj = varValues[j][e] - meanJ;
                    cov  += di * dj;
                    varI += di * di;
                    varJ += dj * dj;
                }

                double corr = 0.0;
                if (varI > 0 && varJ > 0) corr = cov / std::sqrt(varI*varJ);
                hCorrelationMatrix->SetBinContent(i+1, j+1, corr);
            }
        }
        // plot 
        PlotOptions opts;
        opts.drawOption = "COLZ TEXT"; // this draws the values in the cells, we shall see if it looks good :) 
        opts.isHeatmap = true;
        hCorrelationMatrix->SetMinimum(-1.0);
        hCorrelationMatrix->SetMaximum(1.0);
        Plot2D(hCorrelationMatrix, outDir + label + "_correlations.png", opts);
        delete hCorrelationMatrix;
    }

};

#endif