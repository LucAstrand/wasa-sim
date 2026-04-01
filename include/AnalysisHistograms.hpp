#ifndef ANALYSISHISTOGRAMS_H
#define ANALYSISHISTOGRAMS_H


#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TProfile.h"
#include "TMath.h"
#include "Pi0Efficiency.hpp"
#include "Acceptance.hpp"
#include "PIDEfficiency.hpp"
#include "PlotUtils.hpp"
#include "Structures.hpp"

// ============================================================
struct Pi0Histograms {
    TH1F* hMass     = nullptr;
    TH2F* hppM_pre  = nullptr;
    TH2F* hppM_post = nullptr;
    Pi0Efficiency* effPlotter        = nullptr;
    Acceptance*    accPlotter        = nullptr;
    Acceptance*    accVsEta          = nullptr;
    Acceptance*    accVsTheta        = nullptr;

    void Book() {
        hMass     = new TH1F("hPi0Mass",
            ";M_{#gamma#gamma} [MeV];Events", 100, 1.5, 301.5);
        hppM_pre  = new TH2F("hPi0ppM_pre",
            ";M_{#gamma#gamma} [MeV];#theta_{#gamma#gamma} [rad]",
            200, 0, 250, 200, 0, 4);
        hppM_post = new TH2F("hPi0ppM_post",
            ";M_{#gamma#gamma} [MeV];#theta_{#gamma#gamma} [rad]",
            200, 0, 250, 200, 0, 4);
        effPlotter = new Pi0Efficiency(4, 1, 500);
        accPlotter = new Acceptance("nPiE",     100, 1, 1000);
        accVsEta   = new Acceptance("nPiEta",   100, -10, 10);
        accVsTheta = new Acceptance("nPiTheta", 60, 0, TMath::Pi());
    }

    void Plot(int nentries, const std::string& outDir) {
        PlotOptions opts;
        opts.doFit = true; // opts.fitMin = 100; opts.fitMax = 170; use the full range?
        opts.addInfoPave = true;
        opts.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        Plot1D({hMass}, {kBlack}, outDir + "Pi0InvMass.png", opts);

        PlotOptions opts2D;
        opts2D.drawOption  = "SCAT";
        opts2D.infoLines = {"Signal Dataset", Form("%d Events", nentries)};
        // Plot2D(hppM_pre,  "Neutral/Pi0ppM_pre.png",  opts2D);
        // Plot2D(hppM_post, "Neutral/Pi0ppM_post.png", opts2D);
        Plot2DOverlay({hppM_pre, hppM_post}, {kGray, kRed}, outDir + "Pi0ppM.png", opts2D);

        effPlotter->FinalizePlot(outDir + "Pi0_efficiency_vs_Ekin.png");

        PlotOptions optsAcc;
        optsAcc.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        optsAcc.infoLines   = {"Signal dataset"};
        optsAcc.addInfoPave = true;
        optsAcc.xAxisTitle  = "Signal #pi^{0} E_{kin} [MeV]";
        optsAcc.yAxisTitle  = "Acceptance [%]";
        accPlotter->FinalizePlot(outDir + "Pi0_acceptance_vs_Ekin.png", optsAcc);
    }

    void Cleanup() {
        delete hMass; delete hppM_pre; delete hppM_post;
        delete effPlotter; delete accPlotter;
        delete accVsEta; delete accVsTheta;
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
                               outDir + "Pi0InvMassTT.png", opts);

        PlotOptions optsRA;
        optsRA.doFit = true; optsRA.fitMin = 80; optsRA.fitMax = 180;
        if (h_tE_rA) Plot1D({h_tE_rA}, {kBlack}, 
                             outDir + "Pi0InvMassRecoAngle.png", optsRA);

        PlotOptions optsTA;
        optsTA.legendEntries = {"Truth-level Invariant Mass"};
        if (h_rE_tA) Plot1D({h_rE_tA}, {kBlack},
                             outDir + "Pi0InvMassTruthAngle.png", optsTA);
    }

    void Cleanup() {
        delete hTrueMass; delete h_tE_rA; delete h_rE_tA;
    }
};

// ============================================================
struct ChargedHistograms {
    // PID
    TH1F* hNSigmaPion   = nullptr;
    TH1F* hNSigmaProton = nullptr;
    // Pions
    TH2F* hdEdxVsE_cluster_Pion = nullptr;
    TH2F* hdEdxVsE_true_Pion    = nullptr;
    TH1F* hdEdxTruePion         = nullptr;
    TH1F* hdEdxSmearPion        = nullptr;
    TH1F* hPionTheta            = nullptr;
    TH1F* hPionCosTheta         = nullptr;
    TH1F* hPionPhi              = nullptr;
    // Protons
    TH2F* hdEdxVsE_cluster_Proton = nullptr;
    TH2F* hdEdxVsE_true_Proton    = nullptr;
    TH1F* hdEdxTrueProton         = nullptr;
    TH1F* hdEdxSmearProton        = nullptr;
    // Global
    TH2F* h2_Eres       = nullptr;
    TH1F* hClusterE     = nullptr;
    // Cone diagnostics
    TH1F* hAngle_matched   = nullptr;
    TH1F* hAngle_unmatched = nullptr;
    TH2F* hAngleVsE        = nullptr;
    // Acceptance
    PIDEfficiency* pidEff            = nullptr;
    Acceptance*    chAccGlobal       = nullptr;
    Acceptance*    chAccCondTPC      = nullptr;
    Acceptance*    chAccCondNoTPC    = nullptr;

    void Book() {
        hNSigmaPion   = new TH1F("hNSigmaPion",
            ";n#sigma;Counts", 100, -5, 5);
        hNSigmaProton = new TH1F("hNSigmaProton",
            ";n#sigma;Counts", 100, -5, 5);
        hPionTheta    = new TH1F("hPionTheta",
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
        hClusterE = new TH1F("hClusterE",
            ";Cluster E [MeV];Count", 100, 0, 500);
        hAngle_matched   = new TH1F("hHitAngle_matched",
            ";Angular distance [deg];Hits", 180, 0, 180);
        hAngle_unmatched = new TH1F("hHitAngle_unmatched",
            ";Angular distance [deg];Hits", 180, 0, 180);
        hAngleVsE = new TH2F("hHitAngleVsE",
            ";Hit Energy [MeV];Angular distance [deg]", 100, 0, 50, 180, 0, 180);
        pidEff         = new PIDEfficiency(20, 0, 500);
        chAccGlobal    = new Acceptance("chPiEGlobal",    100, 1, 1000);
        chAccCondTPC   = new Acceptance("chPiECondTPC",   100, 1, 1000);
        chAccCondNoTPC = new Acceptance("chPiECondNoTPC", 100, 1, 1000);
    }

    void Plot(const std::string& outDir) {
        PlotOptions opts_nSigma;
        opts_nSigma.doFit = true;
        opts_nSigma.addLegend = true;
        opts_nSigma.legendEntries = {"#pi^{#pm} hypothesis", "p hypothesis"};
        Plot1D({hNSigmaPion, hNSigmaProton}, {kRed+1, kBlue+1},
               outDir + "nSigmaPlots.png", opts_nSigma);

        Plot1D({hPionTheta},    {kBlack}, outDir + "primaryChPionTheta.png",    {});
        Plot1D({hPionCosTheta}, {kBlack}, outDir + "primaryChPionCosTheta.png", {});
        Plot1D({hPionPhi},      {kBlack}, outDir + "primaryChPionPhi.png",      {});

        PlotOptions opts2D;
        opts2D.drawOption  = "SCAT";
        opts2D.legendEntries = {"#pi^{#pm}", "p"};
        opts2D.legendDrawOpt = "P";
        Plot2DOverlay({hdEdxVsE_true_Pion, hdEdxVsE_true_Proton},
                      {kBlack, kRed},
                      outDir + "dedx_vs_E_overlay_true.png", opts2D);
        Plot2DOverlay({hdEdxVsE_cluster_Pion, hdEdxVsE_cluster_Proton},
                      {kBlack, kRed},
                      outDir + "dedx_vs_E_overlay_cluster.png", opts2D);

        TProfile* pEres = h2_Eres->ProfileX(
            Form("pEres_%s", h2_Eres->GetName()));
        Plot1D({pEres}, {kBlack}, outDir + "energy_residual.png", {});

        PlotOptions opts_dEdx;
        opts_dEdx.addLegend = true;
        opts_dEdx.legendEntries = {"True #pi^{#pm}", "Smeared #pi^{#pm}",
                                   "True p",         "Smeared p"};
        Plot1D({hdEdxTruePion, hdEdxSmearPion, hdEdxTrueProton, hdEdxSmearProton},
               {kBlack, kGray, kRed, kRed+1},
               outDir + "dEdxPlots.png", opts_dEdx);

        pidEff->FinalizePlot(outDir + "plots/PIDEfficiency", 211);

        Plot1D({hClusterE}, {kBlack}, outDir + "ClusterE_pion.png", {});

        PlotOptions optsAcc;
        optsAcc.topLatex    = "#bf{Hibeam}  #it{Wasa full simulation}";
        optsAcc.infoLines   = {"Signal dataset"};
        optsAcc.addInfoPave = true;
        optsAcc.xAxisTitle  = "Signal #pi^{#pm} E_{kin} [MeV]";
        optsAcc.yAxisTitle  = "Acceptance [%]";
        chAccGlobal   ->FinalizePlot(outDir + "chPi_acceptance_Global.png",    optsAcc);
        chAccCondTPC  ->FinalizePlot(outDir + "chPi_acceptance_CondTPC.png",   optsAcc);
        chAccCondNoTPC->FinalizePlot(outDir + "chPi_acceptance_CondNoTPC.png", optsAcc);

        PlotOptions opts_cone;
        opts_cone.addLegend = true;
        opts_cone.legendEntries = {"Matched hits", "Unmatched hits"};
        hAngle_matched  ->Scale(1.0 / hAngle_matched  ->Integral());
        hAngle_unmatched->Scale(1.0 / hAngle_unmatched->Integral());
        Plot1D({hAngle_matched, hAngle_unmatched}, {kBlack, kRed},
               outDir + "HitAngleDiagnostic.png", opts_cone);
        Plot2D(hAngleVsE, outDir + "HitAngleVsE.png", {});
    }

    void Cleanup() {
        delete hNSigmaPion; delete hNSigmaProton;
        delete hPionTheta; delete hPionCosTheta; delete hPionPhi;
        delete hdEdxVsE_cluster_Pion; delete hdEdxVsE_true_Pion;
        delete hdEdxTruePion; delete hdEdxSmearPion;
        delete hdEdxVsE_cluster_Proton; delete hdEdxVsE_true_Proton;
        delete hdEdxTrueProton; delete hdEdxSmearProton;
        delete h2_Eres; delete hClusterE;
        delete hAngle_matched; delete hAngle_unmatched; delete hAngleVsE;
        delete pidEff;
        delete chAccGlobal; delete chAccCondTPC; delete chAccCondNoTPC;
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
            ";True #pi^{0} Mult.;Reco #pi^{0} Mult.", 8, -0.5, 7.5, 8, -0.5, 7.5);
        hChPionMult     = new TH2I("hChPionMultiplicity",
            ";True #pi^{#pm} Mult.;Reco #pi^{#pm} Mult.", 9, -0.5, 8.5, 9, -0.5, 8.5);
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
               outDir + "EventTotE.png", opts);

        PlotOptions opts2D;
        opts2D.overlayProfileX = true;
        opts2D.profileColor    = kRed;
        Plot2D(hEvisVsEtrue,  outDir + "EvisVsEtrue.png",  opts2D);
        Plot2D(hErecoVsEtrue, outDir + "ErecoVsEtrue.png", opts2D);

        Plot1D({hDiffErecoEtrue}, {kBlack}, outDir + "DiffErecoEtrue.png", {});

        PlotOptions optsHeatmap;
        optsHeatmap.drawOption = "COLZ";
        optsHeatmap.isHeatmap  = true;
        Plot2D(hNPionMult,  outDir + "NPionMultiplicity.png",  optsHeatmap);
        Plot2D(hChPionMult, outDir + "ChPionMultiplicity.png", optsHeatmap);
    }

    void Cleanup() {
        delete hEvis; delete hEcorrected;
        delete hErecoVsEtrue; delete hEvisVsEtrue; delete hDiffErecoEtrue;
        delete hNPionMult; delete hChPionMult;
    }
};

#endif