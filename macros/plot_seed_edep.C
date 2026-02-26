#include <string.h>

void plot_seed_edep(const char* fname="TEST.root",
                    const char* tname="hibeam",
                    int nbins=200,
                    double xmin=0.0,
                    double xmax=600.0)
{
  // Open file & get tree
  TFile* f = TFile::Open(fname);
  if (!f || f->IsZombie()) {
    Error("plot_seed_edep","Cannot open file %s", fname);
    return;
  }
  TTree* t = (TTree*)f->Get(tname);
  if (!t) {
    Error("plot_seed_edep","Cannot find tree %s", tname);
    return;
  }

  // Vector branches
  std::vector<int>*    seedType = nullptr;
  std::vector<double>* seedEdep = nullptr;

  t->SetBranchAddress("seedType",      &seedType);
  t->SetBranchAddress("seedEdepEmcal", &seedEdep);

  // Histograms
  TH1D* hCompt = new TH1D("hCompt",
      "EMCAL Edep per Seed;Edep (MeV);Seeds",
      nbins, xmin, xmax);
  TH1D* hPion  = new TH1D("hPion",
      "EMCAL Edep per Seed;Edep (MeV);Seeds",
      nbins, xmin, xmax);

  hCompt->SetLineColor(kBlue+1);
  hCompt->SetLineWidth(2);
  hPion->SetLineColor(kRed+1);
  hPion->SetLineWidth(2);

  // Loop over events
  const Long64_t n = t->GetEntries();
  for (Long64_t ie=0; ie<n; ++ie) {
    t->GetEntry(ie);

    if (ie < 5) {
      printf("evt %lld: nSeeds=%zu\n", ie, seedEdep->size());
    }
    if (!seedType || !seedEdep) continue;
    if (seedType->size() != seedEdep->size()) continue;

    // Loop over seeds in this event
    for (size_t i=0; i<seedType->size(); ++i) {
      const int    type = seedType->at(i);     // 1=compt e, 2=pion
      const double edep = seedEdep->at(i);     

      // if (edep > 100) printf("Edep > 100 MeV: %f \n", edep);

      if (type == 1)      hCompt->Fill(edep);
      else if (type == 2) hPion->Fill(edep);
    }
  }

  // printf("Entries in Pion hist %f \n", hPion->GetEntries());

  // Draw
  TCanvas* c = new TCanvas("c","Per-seed EMCAL Edep",900,700);

  //beauty options :) 
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
  std::string text = {"#bf{Hibeam} #it{Wasa full simulation}"};
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize(0.045);
  l.DrawLatex(0.16, 0.93, text.c_str());


  hCompt->Draw("HIST");
  hPion->Draw("HIST SAME");

  auto leg = new TLegend(0.62,0.70,0.88,0.88);
  leg->AddEntry(hCompt, "Compton e seeds", "l");
  leg->AddEntry(hPion,  "Charged pion seeds", "l");
  leg->Draw();

  c->SetLogy();  // often helpful for tails; comment out if you prefer linear

  c->SaveAs("seed_edep_emcal.png");
}
