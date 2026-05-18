/**
 * sec_beamline_gaps.C
 *
 * Visualises the beamline opening regions of the WASA SEC barrel calorimeter.
 *
 * Gap edge definition: geometric polar angle to crystal centre position,
 *   theta_i = arccos( dz[i] / sqrt(r[i]^2 + dz[i]^2) )
 * This gives:
 *   Forward  gap edge (ring 16): theta_fwd_nom = 39.06 deg
 *   Backward gap edge (ring  0): theta_bwd_nom = 141.03 deg
 *
 * For a vertex displaced to radius rho_v on the target disk, the effective
 * gap edge angle is computed from the near-side and far-side approach to the
 * gap-edge ring in the (rho, z) plane:
 *   drho_near = |R - rho_v|,  drho_far = R + rho_v
 *   theta(drho, dz) = arccos( dz / sqrt(drho^2 + dz^2) )
 *   theta_fwd(rho_v) = min(theta_near, theta_far)   [forward ring]
 *   theta_bwd(rho_v) = max(theta_near, theta_far)   [backward ring]
 *
 * Coordinate convention (from Sec::BuildBarrel):
 *   theta = polar angle from +z (beam) axis  [0 = forward, 180 = backward]
 *   phi   = azimuthal angle around beam axis (48 crystals x 7.5 deg)
 *
 * Output:
 *   sec_gap_edges_vs_rho.pdf
 *
 * Usage:
 *   root -l -b -q AcceptanceStudy.C
 */

#include "TCanvas.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStyle.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

// ---------------------------------------------------------------------------
// SEC barrel geometry from Sec::BuildBarrel()
// r[i]  = radial distance of ring inner-face centre from beam axis [cm]
// dz[i] = z position of ring inner-face centre [cm]
// Ring 0  = most backward  (dz = -39.28 cm, geometric theta = 141.03 deg)
// Ring 16 = most forward   (dz = +41.81 cm, geometric theta =  39.06 deg)
// ---------------------------------------------------------------------------
static const int    kNRings = 17;
static const double kR[kNRings] = {
    31.77, 37.75, 42.15, 44.57, 46.38,
    47.75, 48.05, 48.15, 48.0,  48.15,
    48.05, 47.75, 47.27, 46.65, 42.86,
    38.50, 33.93
};
static const double kDZ[kNRings] = {
    -39.28, -36.16, -32.10, -27.12, -21.69,
    -15.91,  -9.74,  -3.60,   2.50,   8.60,
     14.74,  20.91,  27.15,  33.50,  37.47,
     40.07,  41.81
};

static const int kFwdRing = 16;  // most forward  ring -> forward  gap edge
static const int kBwdRing =  0;  // most backward ring -> backward gap edge

// Target: disk of diameter 40 cm, centred at origin
static const double kTargetRadius = 20.0; // cm
static const double kTargetZ      =  0.0; // cm

static const int kNVertices = 500000;

// ---------------------------------------------------------------------------
// Polar angle [deg] from transverse separation drho and longitudinal dz.
// drho >= 0.  Result in [0, 180].
// ---------------------------------------------------------------------------
inline double PolarAngle(double drho, double dz)
{
    double dist = TMath::Sqrt(drho*drho + dz*dz);
    if(dist < 1e-10) return 0.;
    return TMath::RadToDeg() * TMath::ACos(dz / dist);
}

// ---------------------------------------------------------------------------
// Geometric (nominal) polar angle of ring i as seen from the origin.
// ---------------------------------------------------------------------------
inline double NominalTheta(int i)
{
    return PolarAngle(kR[i], kDZ[i]);
}

// ---------------------------------------------------------------------------
// Effective gap edge angles for a vertex at (rho_v, vz).
// Near-side and far-side refer to the two azimuthal extremes in the rho-z plane.
// ---------------------------------------------------------------------------
void GetGapEdges(double rho_v, double vz,
                 double& theta_fwd, double& theta_bwd)
{
    // Forward gap edge: minimum theta to ring 16
    {
        double dz        = kDZ[kFwdRing] - vz;
        double drho_near = TMath::Abs(kR[kFwdRing] - rho_v);
        double drho_far  = kR[kFwdRing] + rho_v;
        theta_fwd = TMath::Min(PolarAngle(drho_near, dz),
                               PolarAngle(drho_far,  dz));
    }

    // Backward gap edge: maximum theta to ring 0
    {
        double dz        = kDZ[kBwdRing] - vz;
        double drho_near = TMath::Abs(kR[kBwdRing] - rho_v);
        double drho_far  = kR[kBwdRing] + rho_v;
        theta_bwd = TMath::Max(PolarAngle(drho_near, dz),
                               PolarAngle(drho_far,  dz));
    }
}

// ---------------------------------------------------------------------------
// Closed TGraph band for filling between two curves.
// ---------------------------------------------------------------------------
TGraph* MakeBand(const std::vector<double>& x,
                 const std::vector<double>& ylo,
                 const std::vector<double>& yhi)
{
    int n = (int)x.size();
    TGraph* g = new TGraph(2 * n);
    for(int i = 0;   i < n; i++) g->SetPoint(i,           x[i], ylo[i]);
    for(int i = n-1; i >= 0; i--) g->SetPoint(2*n-1-i,   x[i], yhi[i]);
    return g;
}

// ---------------------------------------------------------------------------
void AcceptanceStudy()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    const double kThetaFwdNom = NominalTheta(kFwdRing);  // 39.06 deg
    const double kThetaBwdNom = NominalTheta(kBwdRing);  // 141.03 deg

    std::cout << "\n=== Nominal gap edges (geometric theta from origin) ===\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Forward  gap edge (ring 16): " << kThetaFwdNom << " deg\n";
    std::cout << "  Backward gap edge (ring  0): " << kThetaBwdNom << " deg\n\n";

    // -----------------------------------------------------------------------
    // MC: sample vertices uniformly over target disk, record gap edges
    // -----------------------------------------------------------------------
    const int kRhoBins = 50;
    std::vector<std::vector<double>> fwd_per_rho(kRhoBins), bwd_per_rho(kRhoBins);

    TRandom3 rng(42);
    for(int iv = 0; iv < kNVertices; iv++){
        double rho_v = kTargetRadius * TMath::Sqrt(rng.Uniform()); // uniform in area
        double tfwd, tbwd;
        GetGapEdges(rho_v, kTargetZ, tfwd, tbwd);
        int ib = TMath::Min((int)(rho_v / kTargetRadius * kRhoBins), kRhoBins - 1);
        fwd_per_rho[ib].push_back(tfwd);
        bwd_per_rho[ib].push_back(tbwd);
    }

    // Bin statistics: min, median, max
    std::vector<double> rho_c(kRhoBins);
    std::vector<double> fwd_min(kRhoBins), fwd_med(kRhoBins), fwd_max(kRhoBins);
    std::vector<double> bwd_min(kRhoBins), bwd_med(kRhoBins), bwd_max(kRhoBins);

    for(int ib = 0; ib < kRhoBins; ib++){
        rho_c[ib] = (ib + 0.5) * kTargetRadius / kRhoBins;
        auto sf = fwd_per_rho[ib];
        auto sb = bwd_per_rho[ib];
        if(sf.empty()) continue;
        std::sort(sf.begin(), sf.end());
        std::sort(sb.begin(), sb.end());
        int n = (int)sf.size();
        fwd_min[ib] = sf.front();  fwd_med[ib] = sf[n/2];  fwd_max[ib] = sf.back();
        bwd_min[ib] = sb.front();  bwd_med[ib] = sb[n/2];  bwd_max[ib] = sb.back();
    }

    double global_fwd_min = *std::min_element(fwd_min.begin(), fwd_min.end());
    double global_fwd_max = *std::max_element(fwd_max.begin(), fwd_max.end());
    double global_bwd_min = *std::min_element(bwd_min.begin(), bwd_min.end());
    double global_bwd_max = *std::max_element(bwd_max.begin(), bwd_max.end());

    std::cout << "=== Gap edge envelope over full target ===\n";
    std::cout << "  Forward  gap edge: [" << global_fwd_min << ", " << global_fwd_max << "] deg\n";
    std::cout << "  Backward gap edge: [" << global_bwd_min << ", " << global_bwd_max << "] deg\n";
    std::cout << "  Conservative calorimeter coverage: ["
              << global_fwd_max << ", " << global_bwd_min << "] deg\n\n";

    // -----------------------------------------------------------------------
    // Sanity check: at rho_v = 0 the median should equal the nominal theta
    // -----------------------------------------------------------------------
    std::cout << "=== Sanity check: first rho bin median vs nominal ===\n";
    std::cout << "  fwd median (rho~0): " << fwd_med[0]
              << " deg  (nominal: " << kThetaFwdNom << " deg)\n";
    std::cout << "  bwd median (rho~0): " << bwd_med[0]
              << " deg  (nominal: " << kThetaBwdNom << " deg)\n\n";

    // -----------------------------------------------------------------------
    // Plot
    // -----------------------------------------------------------------------
    TCanvas* c = new TCanvas("c", "SEC beamline gap edges", 900, 620);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.05);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.05);

    TH1D* hframe = new TH1D("hframe", "", 1, 0., kTargetRadius);
    hframe->SetMinimum(0.);
    hframe->SetMaximum(180.);
    hframe->GetXaxis()->SetTitle("#rho_{v} [cm]");
    hframe->GetYaxis()->SetTitle("#theta [deg]");
    hframe->GetXaxis()->SetTitleSize(0.050);
    hframe->GetYaxis()->SetTitleSize(0.050);
    hframe->GetXaxis()->SetLabelSize(0.042);
    hframe->GetYaxis()->SetLabelSize(0.042);
    hframe->Draw();

    // // Light background shading for the gap regions
    // TBox* fwdGapBox = new TBox(0., 0., kTargetRadius, global_fwd_max);
    // fwdGapBox->SetFillColorAlpha(kRed, 0.10);
    // fwdGapBox->SetLineColor(0);
    // fwdGapBox->Draw("SAME");

    // TBox* bwdGapBox = new TBox(0., global_bwd_min, kTargetRadius, 180.);
    // bwdGapBox->SetFillColorAlpha(kBlue, 0.10);
    // bwdGapBox->SetLineColor(0);
    // bwdGapBox->Draw("SAME");

        // Light background shading for the gap regions
    TBox* fwdGapBox = new TBox(0., 18.43, kTargetRadius, global_fwd_max);
    fwdGapBox->SetFillColorAlpha(kRed, 0.10);
    fwdGapBox->SetLineColor(0);
    fwdGapBox->Draw("SAME");

    TBox* bwdGapBox = new TBox(0., global_bwd_min, kTargetRadius, 163.32);
    bwdGapBox->SetFillColorAlpha(kBlack, 0.10);
    bwdGapBox->SetLineColor(0);
    bwdGapBox->Draw("SAME");

    // Min-max bands
    TGraph* gFwdBand = MakeBand(rho_c, fwd_min, fwd_max);
    gFwdBand->SetFillColorAlpha(kRed, 0.30);
    gFwdBand->SetLineColor(0);
    gFwdBand->Draw("F SAME");

    TGraph* gBwdBand = MakeBand(rho_c, bwd_min, bwd_max);
    gBwdBand->SetFillColorAlpha(kBlack, 0.30);
    gBwdBand->SetLineColor(0);
    gBwdBand->Draw("F SAME");

    // Median lines
    TGraph* gFwdMed = new TGraph(kRhoBins, rho_c.data(), fwd_med.data());
    gFwdMed->SetLineColor(kRed);
    gFwdMed->SetLineWidth(2);
    gFwdMed->Draw("L SAME");

    TGraph* gBwdMed = new TGraph(kRhoBins, rho_c.data(), bwd_med.data());
    gBwdMed->SetLineColor(kBlack);
    gBwdMed->SetLineWidth(2);
    gBwdMed->Draw("L SAME");

    // Nominal reference lines (constant, rho_v = 0 value)
    TLine* lFwdNom = new TLine(0., kThetaFwdNom, kTargetRadius, kThetaFwdNom);
    lFwdNom->SetLineColor(kBlack);
    lFwdNom->SetLineStyle(2);
    lFwdNom->SetLineWidth(2);
    lFwdNom->Draw();

    TLine* lBwdNom = new TLine(0., kThetaBwdNom, kTargetRadius, kThetaBwdNom);
    lBwdNom->SetLineColor(kBlack);
    lBwdNom->SetLineStyle(2);
    lBwdNom->SetLineWidth(2);
    lBwdNom->Draw();

    // // Labels
    // TLatex lat;
    // lat.SetNDC(false);
    // lat.SetTextSize(0.037);
    // lat.SetTextColor(kRed + 1);
    // lat.DrawLatex(0.5, kThetaFwdNom + 1.8,
    //               Form("#theta_{fwd}^{nom} = %.1f#circ", kThetaFwdNom));
    // lat.SetTextColor(kBlue + 1);
    // lat.DrawLatex(0.5, kThetaBwdNom - 4.0,
    //               Form("#theta_{bwd}^{nom} = %.1f#circ", kThetaBwdNom));
    // lat.SetTextSize(0.033);
    // lat.SetTextColor(kRed);
    // lat.DrawLatex(10.5, 5., "forward beamline opening");
    // lat.SetTextColor(kBlue + 2);
    // lat.DrawLatex(10.5, 175., "backward beamline opening");

        // Labels
    TLatex lat;
    lat.SetNDC(false);
    lat.SetTextSize(0.037);
    lat.SetTextColor(kRed + 1);
    lat.DrawLatex(0.5, kThetaFwdNom + 3.8,
                  Form("#theta_{fwd}^{nom} #approx %.0f#circ", kThetaFwdNom));
    lat.SetTextColor(kBlack);
    lat.DrawLatex(0.5, kThetaBwdNom - 8.0,
                  Form("#theta_{bwd}^{nom} #approx %.0f#circ", kThetaBwdNom));
    lat.SetTextSize(0.033);
    lat.SetTextColor(kRed);
    lat.DrawLatex(10.5, kThetaFwdNom + 3.8, "forward beamline opening");
    lat.SetTextColor(kBlack);
    lat.DrawLatex(10.5, kThetaBwdNom - 8.0, "backward beamline opening");

    // Legend
    TLegend* leg = new TLegend(0.35, 0.38, 0.80, 0.60);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(gFwdMed,  "fwd gap edge #theta_{min}(#rho_{v}), median", "L");
    leg->AddEntry(gFwdBand, "fwd gap edge, min-max range",                  "F");
    leg->AddEntry(gBwdMed,  "bwd gap edge #theta_{max}(#rho_{v}), median", "L");
    leg->AddEntry(gBwdBand, "bwd gap edge, min-max range",                  "F");
    leg->AddEntry(lFwdNom,  "nominal (#rho_{v} = 0)",                       "L");
    leg->Draw();

    c->SaveAs("sec_gap_edges_vs_rho.pdf");
    std::cout << "Saved: sec_gap_edges_vs_rho.pdf\n\n";
}