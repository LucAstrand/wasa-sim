#include "RunABCD.hpp"

ABCDResult RunABCD(
    const std::vector<EventVariables>& sigEvents,
    const std::vector<EventVariables>& bkgEvents,
    // Variable accessors - pass lambdas so any variable pair works
    std::function<double(const EventVariables&)> getVar1,
    std::function<double(const EventVariables&)> getVar2,
    double threshold1,   // signal-like = above this
    double threshold2,   // signal-like = above this
    const std::string& var1Name,
    const std::string& var2Name,
    const std::string& outDir)
{
    ABCDResult result;

    // Count events in each region for background
    int N_A_bkg = 0, N_B_bkg = 0, N_C_bkg = 0, N_D_bkg = 0;
    for (const auto& ev : bkgEvents) {
        bool highV1 = getVar1(ev) > threshold1;
        bool highV2 = getVar2(ev) > threshold2;
        if      ( highV1 &&  highV2) N_A_bkg++;  // signal region
        else if ( highV1 && !highV2) N_B_bkg++;
        else if (!highV1 &&  highV2) N_D_bkg++;
        else                         N_C_bkg++;
    }

    // Count signal events in signal region
    int N_A_sig = 0;
    for (const auto& ev : sigEvents) {
        bool highV1 = getVar1(ev) > threshold1;
        bool highV2 = getVar2(ev) > threshold2;
        if (highV1 && highV2) N_A_sig++;
    }

    result.N_A_sig          = N_A_sig;
    result.N_B_bkg          = N_B_bkg;
    result.N_C_bkg          = N_C_bkg;
    result.N_D_bkg          = N_D_bkg;
    result.N_A_bkg_true     = N_A_bkg;
    result.N_A_bkg_estimate = (N_C_bkg > 0) ? 
        (double)N_B_bkg * N_D_bkg / N_C_bkg : 0.0;
    result.closure          = (N_A_bkg > 0) ?
        result.N_A_bkg_estimate / N_A_bkg : 0.0;

    // Poisson uncertainty on ABCD estimate:
    // sigma^2(N_A_est) = N_A_est^2 * (1/N_B + 1/N_C + 1/N_D)
    if (N_B_bkg > 0 && N_C_bkg > 0 && N_D_bkg > 0) {
        result.closureErr = result.N_A_bkg_estimate *
            std::sqrt(1.0/N_B_bkg + 1.0/N_C_bkg + 1.0/N_D_bkg);
    }

    // Print results
    std::cout << "\n=== ABCD Results: " 
              << var1Name << " vs " << var2Name << " ===\n";
    std::cout << "Thresholds: " << var1Name << " > " << threshold1
              << ", " << var2Name << " > " << threshold2 << "\n";
    std::cout << "Region counts (background):\n";
    std::cout << "  C=" << N_C_bkg << "  D=" << N_D_bkg << "\n";
    std::cout << "  B=" << N_B_bkg << "  A=" << N_A_bkg 
              << " (signal region)\n";
    std::cout << "ABCD estimate in A: " << result.N_A_bkg_estimate
              << " +/- " << result.closureErr << "\n";
    std::cout << "True MC bkg in A:   " << result.N_A_bkg_true << "\n";
    std::cout << "Signal in A:        " << result.N_A_sig << "\n";
    std::cout << "Closure (est/true): " << result.closure << "\n";

    // ---- RooFit visualisation ----
    // Create RooFit variables
    double v1min = 0, v1max = 0, v2min = 0, v2max = 0;
    for (const auto& ev : bkgEvents) {
        v1min = std::min(v1min, getVar1(ev));
        v1max = std::max(v1max, getVar1(ev));
        v2min = std::min(v2min, getVar2(ev));
        v2max = std::max(v2max, getVar2(ev));
    }
    for (const auto& ev : sigEvents) {
        v1min = std::min(v1min, getVar1(ev));
        v1max = std::max(v1max, getVar1(ev));
        v2min = std::min(v2min, getVar2(ev));
        v2max = std::max(v2max, getVar2(ev));
    }
    // Add 10% margin
    double m1 = 0.1*(v1max-v1min), m2 = 0.1*(v2max-v2min);
    v1min -= m1; v1max += m1;
    v2min -= m2; v2max += m2;

    RooRealVar var1(var1Name.c_str(), var1Name.c_str(), v1min, v1max);
    RooRealVar var2(var2Name.c_str(), var2Name.c_str(), v2min, v2max);

    // Build RooDataSets
    RooDataSet dsSig("dsSig", "Signal", RooArgSet(var1, var2));
    RooDataSet dsBkg("dsBkg", "Background", RooArgSet(var1, var2));

    for (const auto& ev : sigEvents) {
        var1.setVal(getVar1(ev));
        var2.setVal(getVar2(ev));
        dsSig.add(RooArgSet(var1, var2));
    }
    for (const auto& ev : bkgEvents) {
        var1.setVal(getVar1(ev));
        var2.setVal(getVar2(ev));
        dsBkg.add(RooArgSet(var1, var2));
    }

    // Plot 2D scatter with ABCD region boundaries
    TCanvas* c = new TCanvas("cABCD", "ABCD", 800, 700);
    c->SetRightMargin(0.12);

    // Create 2D histogram from datasets for visualisation
    TH2F* h2Sig = new TH2F("h2SigABCD", 
        (";"+var1Name+";"+var2Name).c_str(),
        50, v1min, v1max, 50, v2min, v2max);
    TH2F* h2Bkg = new TH2F("h2BkgABCD",
        (";"+var1Name+";"+var2Name).c_str(),
        50, v1min, v1max, 50, v2min, v2max);

    for (const auto& ev : sigEvents) 
        h2Sig->Fill(getVar1(ev), getVar2(ev));
    for (const auto& ev : bkgEvents) 
        h2Bkg->Fill(getVar1(ev), getVar2(ev));

    // Draw background as base
    h2Bkg->SetMarkerColor(kRed);
    h2Bkg->SetMarkerStyle(20);
    h2Bkg->SetMarkerSize(0.5);
    h2Bkg->Draw("SCAT");

    // Overlay signal
    h2Sig->SetMarkerColor(kBlue);
    h2Sig->SetMarkerStyle(20);
    h2Sig->SetMarkerSize(0.5);
    h2Sig->Draw("SCAT SAME");

    // Draw ABCD boundary lines
    TLine* lV1 = new TLine(threshold1, v2min, threshold1, v2max);
    TLine* lV2 = new TLine(v1min, threshold2, v1max, threshold2);
    lV1->SetLineColor(kBlack); lV1->SetLineWidth(2); 
    lV1->SetLineStyle(2);  // dashed
    lV2->SetLineColor(kBlack); lV2->SetLineWidth(2); 
    lV2->SetLineStyle(2);
    lV1->Draw(); lV2->Draw();

    // Label the regions
    TLatex latex;
    latex.SetTextSize(0.05);
    latex.SetTextAlign(22);
    double midLow1  = 0.5*(v1min + threshold1);
    double midHigh1 = 0.5*(threshold1 + v1max);
    double midLow2  = 0.5*(v2min + threshold2);
    double midHigh2 = 0.5*(threshold2 + v2max);

    latex.DrawLatex(midHigh1, midHigh2, 
        Form("A (SR)\nN_{sig}=%d\nN_{bkg}=%d", N_A_sig, N_A_bkg));
    latex.DrawLatex(midHigh1, midLow2,  
        Form("B\nN=%d", N_B_bkg));
    latex.DrawLatex(midLow1,  midLow2,  
        Form("C\nN=%d", N_C_bkg));
    latex.DrawLatex(midLow1,  midHigh2, 
        Form("D\nN=%d", N_D_bkg));

    // ABCD estimate box
    TLegend* leg = new TLegend(0.15, 0.75, 0.55, 0.92);
    leg->SetBorderSize(1);
    leg->AddEntry(h2Sig, "Signal", "P");
    leg->AddEntry(h2Bkg, "Background", "P");
    leg->AddEntry((TObject*)nullptr, 
        Form("ABCD est: %.1f #pm %.1f", 
             result.N_A_bkg_estimate, result.closureErr), "");
    leg->AddEntry((TObject*)nullptr,
        Form("True bkg: %.0f", result.N_A_bkg_true), "");
    leg->AddEntry((TObject*)nullptr,
        Form("Closure: %.2f", result.closure), "");
    leg->Draw();

    c->SaveAs((outDir + "ABCD_" + var1Name + "_vs_" + 
               var2Name + ".png").c_str());

    delete c; delete h2Sig; delete h2Bkg;
    delete lV1; delete lV2;

    return result;
}