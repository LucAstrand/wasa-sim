#ifndef ANALYSISHISTOGRAMS_H
#define ANALYSISHISTOGRAMS_H


#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TProfile.h"
#include "TMath.h"
#include "Pi0Efficiency.hpp"
#include "ChEfficiency.hpp"
#include "Acceptance.hpp"
#include "PIDEfficiency.hpp"
#include "PlotUtils.hpp"
#include "Structures.hpp"

// ============================================================
struct Pi0Histograms {
    TH1F* hMass     = nullptr;
    TH1F* hMassPreSelection     = nullptr;
    TH2F* hppM_pre  = nullptr;
    TH2F* hppM_post = nullptr;
    Pi0Efficiency* effPlotter        = nullptr;
    // Acceptance*    accPlotter        = nullptr;
    // Acceptance*    accVsEta          = nullptr;
    // Acceptance*    accVsTheta        = nullptr;
    Acceptance* accTheta = nullptr;
    Acceptance* accEkin  = nullptr;
    Acceptance* acc2D    = nullptr;
    TH1F* hPhotonE       = nullptr;
    TH1F* hPhotonNum     = nullptr;

    void Book() {
        hMass      = new TH1F("hPi0Mass",
            ";M_{#gamma#gamma} [MeV];Events", 100, 1.5, 301.5);
        hMassPreSelection      = new TH1F("hPi0MassPreSelection",
            ";M_{#gamma#gamma} [MeV];Events", 100, 1.5, 301.5);
        hppM_pre   = new TH2F("hPi0ppM_pre",
            ";M_{#gamma#gamma} [MeV];#theta_{#gamma#gamma} [rad]",
            200, 0, 250, 200, 0, 4);
        hppM_post  = new TH2F("hPi0ppM_post",
            ";M_{#gamma#gamma} [MeV];#theta_{#gamma#gamma} [rad]",
            200, 0, 250, 200, 0, 4);
        effPlotter = new Pi0Efficiency(20, 1, 1000);
        // accPlotter = new Acceptance("nPiE",     100, 1, 1000);
        // accVsEta   = new Acceptance("nPiEta",   100, -10, 10);
        // accVsTheta = new Acceptance("nPiTheta", 60, 0, TMath::Pi());
        accTheta   = new Acceptance("npi_theta", AccAxisType::kTheta, 20, 0, 3.2);
        accEkin    = new Acceptance("npi_ekin", AccAxisType::kEkin, 20, 0, 1000);
        acc2D      = new Acceptance("chpi_2d", 25, 0, 1000, 20, 0, TMath::Pi());
        hPhotonE   = new TH1F("hPhotonE", ";E_{#gamma};Counts", 100, 0, 400);
        hPhotonNum = new TH1F("hPhotonNum", ";Multiplicity;Counts", 10, 0,10);
    }

    void Plot(int nentries, const std::string& outDir) {
        PlotOptions opts;
        opts.doFit = true; 
        opts.fitMin = 90; 
        opts.fitMax = 180; 
        opts.fitType = FitType::Gaussian;
        opts.drawBkgComponent = false;
        // opts.addInfoPave = true;
        // opts.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        Plot1D({hMass}, {kBlack}, outDir + "Pi0InvMass.pdf", opts);

        PlotOptions optsPreSelection;
        optsPreSelection.doFit = true; // opts.fitMin = 100; opts.fitMax = 170; use the full range?
        optsPreSelection.fitType = FitType::GaussPlusPhaseSpace;
        optsPreSelection.drawBkgComponent = true;
        // optsPreSelection.addInfoPave = true;
        // optsPreSelection.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        Plot1D({hMassPreSelection}, {kBlack}, outDir + "Pi0InvMassPreSelection.pdf", optsPreSelection);

        PlotOptions opts2D;
        opts2D.drawOption  = "SCAT";
        // opts2D.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        // Plot2D(hppM_pre,  "Neutral/Pi0ppM_pre.pdf",  opts2D);
        // Plot2D(hppM_post, "Neutral/Pi0ppM_post.pdf", opts2D);
        opts2D.legendEntries = {"All #gamma pairs", "Selected #gamma pairs"};
        opts2D.legendDrawOpt = "P";
        Plot2DOverlay({hppM_pre, hppM_post}, {kGray, kRed}, outDir + "Pi0ppM.pdf", opts2D);

        effPlotter->FinalizePlot(outDir + "Pi0_efficiency_vs_Ekin.pdf");

        PlotOptions optsAccEkin;
        optsAccEkin.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsAccEkin.addInfoPave = true;
        // optsAccEkin.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsAccEkin.xAxisTitle  = "Signal #pi^{0} E_{kin} [MeV]";
        optsAccEkin.yAxisTitle  = "Acceptance [%]";
        accEkin->FinalizePlot(outDir + "Pi0_acceptance_vs_Ekin.pdf", optsAccEkin);

        PlotOptions optsAccTheta;
        optsAccTheta.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsAccTheta.addInfoPave = true;
        // optsAccTheta.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsAccTheta.xAxisTitle  = "Signal #pi^{0} #theta [rad]";
        optsAccTheta.yAxisTitle  = "Acceptance [%]";
        accTheta->FinalizePlot(outDir + "Pi0_acceptance_vs_Theta.pdf", optsAccTheta);

        PlotOptions optsAcc2D;
        optsAcc2D.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsAcc2D.addInfoPave = true;
        // optsAcc2D.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsAcc2D.xAxisTitle  = "Signal #pi^{0} E_{kin} [MeV]";
        optsAcc2D.yAxisTitle  = "Signal #pi^{0} #theta [rad]";
        optsAcc2D.zAxisTitle  = "Acceptance (%)"; 
        acc2D->FinalizePlot(outDir + "Pi0_acceptance_2D.pdf", optsAcc2D);

        PlotOptions optsPhotonE;
        // optsPhotonE.addInfoPave = true;
        // optsPhotonE.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        Plot1D({hPhotonE}, {kBlack}, outDir + "PhotonEnergy.pdf", optsPhotonE);

        PlotOptions optsPhotonN;
        // optsPhotonN.addInfoPave = true;
        // optsPhotonN.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        Plot1D({hPhotonNum}, {kBlack}, outDir + "PhotonMultiplicity.pdf", optsPhotonN);
    }

    void Cleanup() {
        delete hMass; delete hppM_pre; delete hppM_post;
        delete effPlotter; delete accEkin;
        delete accTheta; delete acc2D;
        delete hPhotonE; delete hPhotonNum;
    }
};

// ============================================================
struct TruthHistograms {
    TH1F* hTrueMass          = nullptr;
    TH1F* h_tE_rA            = nullptr;
    TH1F* h_rE_tA            = nullptr;

    void Book() {
        hTrueMass = new TH1F("hPi0TrueMass",
            ";M_{#gamma#gamma} [MeV];Events", 100, 1.5, 301.5);
        h_tE_rA   = new TH1F("h_tE_rA",
            ";M_{#gamma#gamma} [MeV];Events", 100, 1.5, 301.5);
        h_rE_tA   = new TH1F("h_rE_tA",
            ";M_{#gamma#gamma} [MeV];Events", 100, 1.5, 301.5);
    }

    void Plot(const std::string& outDir) {
        PlotOptions opts;
        opts.addInfoPave = true;
        opts.legendEntries = {"Truth-level Invariant Mass"};
        if (hTrueMass) Plot1D({hTrueMass}, {kBlack}, 
                               outDir + "Pi0InvMassTT.pdf", opts);

        PlotOptions optsRA;
        optsRA.doFit = true; optsRA.fitMin = 80; optsRA.fitMax = 180;
        if (h_tE_rA) Plot1D({h_tE_rA}, {kBlack}, 
                             outDir + "Pi0InvMassRecoAngle.pdf", optsRA);

        PlotOptions optsTA;
        optsTA.legendEntries = {"Truth-level Invariant Mass"};
        if (h_rE_tA) Plot1D({h_rE_tA}, {kBlack},
                             outDir + "Pi0InvMassTruthAngle.pdf", optsTA);
    }

    void Cleanup() {
        delete hTrueMass; delete h_tE_rA; delete h_rE_tA;
    }
};

// ============================================================
struct ChargedHistograms {
    // PID
    TH1F* hNSigmaPion                     = nullptr;
    TH1F* hNSigmaProton                   = nullptr;
    TH1F* hNSigmaElectron                 = nullptr;
    TH1F* hNSigmaMuon                     = nullptr;
    TGraph* gdEdxPion                     = nullptr;
    TGraph* gdEdxProton                   = nullptr;
    // Pions
    TH2F* hdEdxVsE_cluster_Pion           = nullptr;
    TH2F* hdEdxVsE_true_Pion              = nullptr;
    TH1F* hdEdxTruePion                   = nullptr;
    TH1F* hdEdxSmearPion                  = nullptr;
    TH1F* hPionTheta                      = nullptr;
    TH1F* hPionThetaWASA                  = nullptr;
    TH1F* hPionCosTheta                   = nullptr;
    TH1F* hPionPhi                        = nullptr;
    // Protons
    TH2F* hdEdxVsE_cluster_Proton         = nullptr;
    TH2F* hdEdxVsE_true_Proton            = nullptr;
    TH1F* hdEdxTrueProton                 = nullptr;
    TH1F* hdEdxSmearProton                = nullptr;
    // Global
    TH2F* h2_Eres                         = nullptr;
    TH1F* hClusterE_Pion                  = nullptr;
    TH1F* hClusterE_Electron              = nullptr;
    // Cone diagnostics
    TH1F* hAngle_matched                  = nullptr;
    TH1F* hAngle_unmatched                = nullptr;
    TH2F* hAngleVsE                       = nullptr;
    // Acceptance
    Acceptance*    chAccGlobalEkin        = nullptr;
    Acceptance*    chAccCondTPCEkin       = nullptr;
    Acceptance*    chAccCondNoTPCEkin     = nullptr;
    Acceptance*    chAccGlobalTheta       = nullptr;
    Acceptance*    chAccCondTPCTheta      = nullptr;
    Acceptance*    chAccCondNoTPCTheta    = nullptr;
    // Acceptance*    accNChTracks      = nullptr;
    // Efficiency
    PIDEfficiency* pidEffPion             = nullptr;
    PIDEfficiency* pidEffProton           = nullptr;
    ChEfficiency*  chEffTheta             = nullptr;
    ChEfficiency*  chEffNChTracks         = nullptr;

    void Book() {
        hNSigmaPion   = new TH1F("hNSigmaPion",
            ";n#sigma;Counts", 100, -5, 5);
        hNSigmaProton = new TH1F("hNSigmaProton",
            ";n#sigma;Counts", 100, -5, 5);
        hNSigmaElectron = new TH1F("hNSigmaElectron",
            ";n#sigma;Counts", 100, -5, 5);
        hNSigmaMuon = new TH1F("hNSigmaMuon",
            ";n#sigma;Counts", 100, -5, 5);
        hPionTheta    = new TH1F("hPionTheta",
            ";#theta;Counts", 60, 0, TMath::Pi());
        hPionThetaWASA= new TH1F("hPionThetaWASA",
            ";#theta;Counts", 60, 0, TMath::Pi());
        hPionCosTheta = new TH1F("hPionCosTheta",
            ";cos #theta;Counts", 60, -1, 1);
        hPionPhi      = new TH1F("hPionPhi",
            ";#phi;Counts", 60, 0, TMath::Pi());
        hdEdxVsE_cluster_Pion = new TH2F("hdEdxVsE_cluster_Pion",
            ";E [MeV];dE/dx [MeV/cm]", 100, 0, 500, 100, 0, 0.1);
        hdEdxVsE_true_Pion    = new TH2F("hdEdxVsE_true_Pion",
            ";E [MeV];dE/dx [MeV/cm]", 200, 0, 500, 200, 0, 0.1);
        hdEdxTruePion  = new TH1F("hdEdxTruePion",
            ";dE/dx [MeV/cm];Counts", 100, 0, 0.025);
        hdEdxSmearPion = new TH1F("hdEdxSmearPion",
            ";dE/dx [MeV/cm];Counts", 100, 0, 0.025);
        hdEdxVsE_cluster_Proton = new TH2F("hdEdxVsE_cluster_Proton",
            ";E [MeV];dE/dx [MeV/cm]", 100, 0, 500, 100, 0, 0.1);
        hdEdxVsE_true_Proton    = new TH2F("hdEdxVsE_true_Proton",
            ";E [MeV];dE/dx [MeV/cm]", 200, 0, 500, 200, 0, 0.1);
        hdEdxTrueProton  = new TH1F("hdEdxTrueProton",
            ";dE/dx [MeV/cm];Counts", 100, 0, 0.025);
        hdEdxSmearProton = new TH1F("hdEdxSmearProton",
            ";dE/dx [MeV/cm];Counts", 100, 0, 0.025);
        h2_Eres   = new TH2F("h2_Eres",
            ";True KE [MeV];Energy Residual", 20, 0, 500, 100, -1.0, 1.0);
        hClusterE_Pion = new TH1F("hClusterE_Pion",
            ";Cluster E [MeV];Count", 100, 0, 600);
        hClusterE_Electron = new TH1F("hClusterE_Electron",
            ";Cluster E [MeV];Count", 100, 0, 600);
        hAngle_matched   = new TH1F("hHitAngle_matched",
            ";Angular distance [deg];Hits", 180, 0, 180);
        hAngle_unmatched = new TH1F("hHitAngle_unmatched",
            ";Angular distance [deg];Hits", 180, 0, 180);
        hAngleVsE = new TH2F("hHitAngleVsE",
            ";Hit Energy [MeV];Angular distance [deg]", 100, 0, 50, 180, 0, 180);
        pidEffPion     = new PIDEfficiency("PidPion", 20, 0, 500);
        pidEffProton   = new PIDEfficiency("PidProton", 20, 0, 500);
        // chAccGlobal    = new Acceptance("chPiEGlobal",    100, 1, 1000);
        // chAccCondTPC   = new Acceptance("chPiECondTPC",   100, 1, 1000);
        // chAccCondNoTPC = new Acceptance("chPiECondNoTPC", 100, 1, 1000);
        chAccGlobalEkin     = new Acceptance("chPiEGlobal", AccAxisType::kEkin,    100, 1, 1000);
        chAccCondTPCEkin    = new Acceptance("chPiECondTPC", AccAxisType::kEkin,   100, 1, 1000);
        chAccCondNoTPCEkin  = new Acceptance("chPiECondNoTPC", AccAxisType::kEkin, 100, 1, 1000);
        chAccGlobalTheta    = new Acceptance("chPiEGlobalTheta", AccAxisType::kTheta,    20, 0, 3.2);
        chAccCondTPCTheta   = new Acceptance("chPiECondTPCTheta", AccAxisType::kTheta,   20, 0, 3.2);
        chAccCondNoTPCTheta = new Acceptance("chPiECondNoTPCTheta", AccAxisType::kTheta, 20, 0, 3.2);
        // accNChTracks   = new Acceptance("chpi_chTracks", AccAxisType::kTracks, 9, -0.5, 8.5);
        chEffTheta          = new ChEfficiency("chEff_theta", AccAxisTypeChPi::kTheta, 20, 0, 3.2);
        chEffNChTracks      = new ChEfficiency("chEff_nChTracks", AccAxisTypeChPi::kTracks, 9, -0.5, 8.5);
    }

    void Plot(int nentries, const std::string& outDir) {
        PlotOptions opts_nSigma;
        opts_nSigma.doFit = true;
        opts_nSigma.fitMin = -3.0; 
        opts_nSigma.fitMax = 3.0; 
        opts_nSigma.fitType = FitType::Gaussian;
        opts_nSigma.drawBkgComponent = false;
        opts_nSigma.addLegend = true;
        // opts_nSigma.legendEntries = {"#pi^{#pm} hypothesis", "p hypothesis", "e hypothesis", "#mu hypothesis"};
        // Plot1D({hNSigmaPion, hNSigmaProton, hNSigmaElectron, hNSigmaMuon}, {kRed, kGray+1, kGray+2, kGray+3},
        //        outDir + "nSigmaPlots.pdf", opts_nSigma);
        opts_nSigma.legendEntries = {"#pi^{#pm} hypothesis", "p hypothesis", "e hypothesis"};
        Plot1D({hNSigmaPion, hNSigmaProton, hNSigmaElectron}, {kRed, kGray+2, kGray+3},
               outDir + "nSigmaPlots.pdf", opts_nSigma);
        // Plot1D({hNSigmaPion}, {kRed},
        //        outDir + "nSigmaPlots.pdf", opts_nSigma);
        

        Plot1D({hPionTheta},    {kBlack}, outDir + "primaryChPionTheta.pdf",    {});
        Plot1D({hPionThetaWASA},{kBlack}, outDir + "primaryChPionThetaWASA.pdf",{});
        Plot1D({hPionCosTheta}, {kBlack}, outDir + "primaryChPionCosTheta.pdf", {});
        Plot1D({hPionPhi},      {kBlack}, outDir + "primaryChPionPhi.pdf",      {});

        PlotOptions opts2D;
        opts2D.setLogX = true;
        opts2D.drawOption  = "SCAT";
        opts2D.legendEntries = {"#pi^{#pm}", "p"};
        opts2D.legendDrawOpt = "P";
        Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton},
                      {kRed, kBlack},
                      outDir + "dedx_vs_E_overlay.pdf", opts2D);

        // PlotOptions opts2D;
        // opts2D.setLogX = true;
        // opts2D.drawOption  = "SCAT";
        // opts2D.legendEntries = {"#pi^{#pm} smeared", "#pi^{#pm} truth", "p smeared", "p truth"};
        // opts2D.legendDrawOpt = "P";
        // Plot2DOverlayGraph({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton},
        //                    {gdEdxPion, gdEdxProton},
        //                    {kRed, kRed},
        //                    {kGray, kGray},
        //                    outDir + "dedx_vs_E_overlay.pdf",
        //                    opts2D);



        // Plot2DWithBands(
        // {hdEdxVsE_true_Pion, hdEdxVsE_cluster_Pion, hdEdxVsE_true_Proton, hdEdxVsE_cluster_Proton},
        // {kP6Violet, kP6Gray, kP6Yellow, kP6Red},
        // outDir + "dEdx_overlay.pdf",
        // opts2D,
        // {true, false, true, false},  // shade truth hists only
        // 0.15);                        // 15% resolution

        TProfile* pEres = h2_Eres->ProfileX(
            Form("pEres_%s", h2_Eres->GetName()));
        Plot1D({pEres}, {kBlack}, outDir + "energy_residual.pdf", {});

        PlotOptions opts_dEdx;
        opts_dEdx.addLegend = true;
        opts_dEdx.legendEntries = {"True #pi^{#pm}", "Smeared #pi^{#pm}",
                                   "True p",         "Smeared p"};
        Plot1D({hdEdxTruePion, hdEdxSmearPion, hdEdxTrueProton, hdEdxSmearProton},
               {kBlack, kGray, kRed, kRed+1},
               outDir + "dEdxPlots.pdf", opts_dEdx);

        pidEffPion->FinalizePlot(outDir + "plots/PIDEfficiencyPion", 211);
        pidEffPion->FinalizePlot(outDir + "plots/PIDEfficiencyProton", 2212);

        // Plot1D({hClusterE_Pion}, {kBlack}, outDir + "ClusterE_pion.pdf", {});
        // Plot1D({hClusterE_Electron}, {kBlack}, outDir + "ClusterE_Electron.pdf", {});
        PlotOptions optsClusterE;
        optsClusterE.setLogX = true;
        Plot1D({hClusterE_Pion, hClusterE_Electron}, {kBlack, kRed}, outDir + "ClusterE_PionElectron.pdf", optsClusterE);


        PlotOptions optsAccEkin;
        optsAccEkin.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsAccEkin.addInfoPave = true;
        // optsAccEkin.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsAccEkin.xAxisTitle  = "Signal #pi^{#pm} E_{kin} [MeV]";
        optsAccEkin.yAxisTitle  = "Acceptance [%]";
        chAccGlobalEkin->FinalizePlot(outDir + "chPi_acceptance_Global_Ekin.pdf",    optsAccEkin);
        chAccCondTPCEkin->FinalizePlot(outDir + "chPi_acceptance_CondTPC_Ekin.pdf",   optsAccEkin);
        chAccCondNoTPCEkin->FinalizePlot(outDir + "chPi_acceptance_CondNoTPC_Ekin.pdf", optsAccEkin);

        PlotOptions optsAccTheta;
        optsAccTheta.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsAccTheta.addInfoPave = true;
        // optsAccTheta.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsAccTheta.xAxisTitle  = "Signal #pi^{#pm} #theta [rad]";
        optsAccTheta.yAxisTitle  = "Acceptance [%]";
        chAccGlobalTheta->FinalizePlot(outDir + "chPi_acceptance_Global_Theta.pdf",    optsAccTheta);
        chAccCondTPCTheta->FinalizePlot(outDir + "chPi_acceptance_CondTPC_Theta.pdf",   optsAccTheta);
        chAccCondNoTPCTheta->FinalizePlot(outDir + "chPi_acceptance_CondNoTPC_Theta.pdf", optsAccTheta);

        // PlotOptions optsAccNChTracks;
        // optsAccNChTracks.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsAccNChTracks.infoLines   = {"Signal dataset"};
        // optsAccNChTracks.addInfoPave = true;
        // optsAccNChTracks.xAxisTitle  = "Number of charged tracks";
        // optsAccNChTracks.yAxisTitle  = "Acceptance [%]";
        // accNChTracks->FinalizePlot(outDir + "chPi_acceptance_vs_nTracks.pdf", optsAccNChTracks);

        PlotOptions optsEffTheta;
        optsEffTheta.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsEffTheta.addInfoPave = true;
        // optsEffTheta.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsEffTheta.xAxisTitle  = "Signal #pi^{#pm} #theta [rad]";
        optsEffTheta.yAxisTitle  = "Efficiency [%]";
        chEffTheta->FinalizePlot(outDir + "chPi_eff_vs_theta.pdf", optsEffTheta);

        PlotOptions optsEffNChTracks;
        optsEffNChTracks.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        // optsEffNChTracks.addInfoPave = true;
        // optsEffNChTracks.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        optsEffNChTracks.xAxisTitle  = "Number of charged tracks";
        optsEffNChTracks.yAxisTitle  = "Efficiency [%]";
        chEffNChTracks->FinalizePlot(outDir + "chPi_eff_vs_nTracks.pdf", optsEffNChTracks);

        PlotOptions opts_cone;
        opts_cone.addLegend = true;
        opts_cone.legendEntries = {"Matched hits", "Unmatched hits"};
        hAngle_matched  ->Scale(1.0 / hAngle_matched  ->Integral());
        hAngle_unmatched->Scale(1.0 / hAngle_unmatched->Integral());
        Plot1D({hAngle_matched, hAngle_unmatched}, {kBlack, kRed},
               outDir + "HitAngleDiagnostic.pdf", opts_cone);
        Plot2D(hAngleVsE, outDir + "HitAngleVsE.pdf", {});
    }

    void Cleanup() {
        delete hNSigmaPion; delete hNSigmaProton;
        delete hPionTheta; delete hPionCosTheta; delete hPionPhi;
        delete hdEdxVsE_cluster_Pion; delete hdEdxVsE_true_Pion;
        delete hdEdxTruePion; delete hdEdxSmearPion;
        delete hdEdxVsE_cluster_Proton; delete hdEdxVsE_true_Proton;
        delete hdEdxTrueProton; delete hdEdxSmearProton;
        delete h2_Eres; delete hClusterE_Pion; delete hClusterE_Electron;
        delete hAngle_matched; delete hAngle_unmatched; delete hAngleVsE;
        delete pidEffPion;
        delete pidEffProton;
        delete chAccGlobalEkin; delete chAccCondTPCEkin; delete chAccCondNoTPCEkin;
        delete chAccGlobalTheta; delete chAccCondTPCTheta; delete chAccCondNoTPCTheta;
        delete chEffNChTracks;
    }
};

// ============================================================
struct EventVarHistograms {
    TH1F* hEvis             = nullptr;
    TH1F* hEcorrected       = nullptr;
    TH2F* hErecoVsEtrue     = nullptr;
    TH2F* hEvisVsEtrue      = nullptr;
    TH1F* hDiffErecoEtrue   = nullptr;
    TH2I* hNPionMult        = nullptr;
    TH2I* hChPionMult       = nullptr;

    void Book() {
        hEvis           = new TH1F("hEventTotE",
            ";Total Energy [MeV];Events", 100, 0, 3000);
        hEcorrected     = new TH1F("hEventCorrectedTotE",
            ";Total Corrected Energy [MeV];Events", 100, 0, 3000);
        hErecoVsEtrue   = new TH2F("hErecoVsEtrue",
            ";E_{true} [MeV];E_{reco} [MeV]", 100, 0, 3000, 100, 0, 3000);
        hEvisVsEtrue    = new TH2F("hEvisVsEtrue",
            ";E_{true} [MeV];E_{vis} [MeV]", 100, 0, 3000, 100, 0, 3000);
        hDiffErecoEtrue = new TH1F("hERecoVsETrue",
            ";E_{reco} - E_{true} [MeV];Counts", 100, -1000, 1000);
        hNPionMult      = new TH2I("hNPionMultiplicity",
            ";True #pi^{0} Mult.;Reco #pi^{0} Mult.;Counts", 8, -0.5, 7.5, 8, -0.5, 7.5);
        hChPionMult     = new TH2I("hChPionMultiplicity",
            ";True #pi^{#pm} Mult.;Reco #pi^{#pm} Mult.;Counts", 9, -0.5, 8.5, 9, -0.5, 8.5);
    }

    void Fill(double eVis, double eReco, double eTrue,
              int truePi0Mult, int recoPi0Mult,
              int trueChPiMult, int recoChPiMult) {
        if (hEvis)           hEvis          ->Fill(eVis);
        if (hEcorrected)     hEcorrected    ->Fill(eReco);
        if (hErecoVsEtrue)   hErecoVsEtrue  ->Fill(eTrue, eReco);
        if (hEvisVsEtrue)    hEvisVsEtrue   ->Fill(eTrue, eVis);
        if (hDiffErecoEtrue) hDiffErecoEtrue->Fill(eReco - eTrue);
        if (hNPionMult)      hNPionMult     ->Fill(truePi0Mult,  recoPi0Mult);
        if (hChPionMult)     hChPionMult    ->Fill(trueChPiMult, recoChPiMult);
    }

    void Plot(const std::string& outDir) {
        PlotOptions opts;
        opts.addLegend = true;
        opts.legendEntries = {"E_{reco}", "E_{vis}"};
        Plot1D({hEcorrected, hEvis}, {kBlack, kRed},
               outDir + "EventTotE.pdf", opts);

        PlotOptions opts2D;
        opts2D.overlayProfileX = true;
        opts2D.profileColor    = kRed;
        Plot2D(hEvisVsEtrue,  outDir + "EvisVsEtrue.pdf",  opts2D);
        Plot2D(hErecoVsEtrue, outDir + "ErecoVsEtrue.pdf", opts2D);

        Plot1D({hDiffErecoEtrue}, {kBlack}, outDir + "DiffErecoEtrue.pdf", {});

        PlotOptions optsHeatmap;
        optsHeatmap.drawOption = "COLZ";
        optsHeatmap.isHeatmap  = true;
        optsHeatmap.zAxisTitle = "Counts";
        Plot2D(hNPionMult,  outDir + "NPionMultiplicity.pdf",  optsHeatmap);
        Plot2D(hChPionMult, outDir + "ChPionMultiplicity.pdf", optsHeatmap);
    }

    void Cleanup() {
        delete hEvis; delete hEcorrected;
        delete hErecoVsEtrue; delete hEvisVsEtrue; delete hDiffErecoEtrue;
        delete hNPionMult; delete hChPionMult;
    }
};

#endif