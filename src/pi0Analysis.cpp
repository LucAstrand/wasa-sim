#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

#include "Structures.hpp"
#include "Clustering.hpp"
#include "PlotUtils.hpp"
#include "EventDisplay.hpp"
#include "Utils.hpp"
#include "TruePhotonCalc.hpp"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

    SetPrettyStyle();

    TFile *f = TFile::Open(argv[1]);
    if (!f || f->IsZombie()) return 1;

    TTree *t = (TTree *)f->Get("digitizedHits");
    if (!t) { std::cerr << "Tree digitizedHits not found" << std::endl; return 1; }

    // Single consolidated vectors (matches preprocessed output—no per-ring arrays!)
    std::vector<double> *centerXs = nullptr, *centerYs = nullptr, *centerZs = nullptr, *energies = nullptr;
    t->SetBranchAddress("centerX", &centerXs);
    t->SetBranchAddress("centerY", &centerYs);
    t->SetBranchAddress("centerZ", &centerZs);
    t->SetBranchAddress("energy", &energies);

    // Primary vertex (vector, but size=1 per event)
    std::vector<double> *primaryX=nullptr,*primaryY=nullptr,*primaryZ=nullptr;
    t->SetBranchAddress("PrimaryPosX",&primaryX);
    t->SetBranchAddress("PrimaryPosY",&primaryY);
    t->SetBranchAddress("PrimaryPosZ",&primaryZ);

    // Truth level info 
    std::vector<double> *truthPosX=nullptr, *truthPosY=nullptr, *truthPosZ=nullptr, *truthE=nullptr;
    t->SetBranchAddress("truthPosX", &truthPosX);
    t->SetBranchAddress("truthPosY", &truthPosY);
    t->SetBranchAddress("truthPosZ", &truthPosZ);
    t->SetBranchAddress("truthE", &truthE);


    Long64_t nentries = t->GetEntries();
    TH1F *hPi0Mass = new TH1F("hPi0Mass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);
    TH1F *hClusterE = new TH1F("hClusterE",";Cluster E [MeV];Count",100,0,500);
    TH1F *hNClusters = new TH1F("hNClusters",";N_{clusters};Events",20,0,20);

    // TH1F *hSingleClusterE = new TH1F("hSingleClusterE",";Cluster E [MeV];Count",100,0,500);

    TH1F *hPi0TrueMass = new TH1F("hPi0TrueMass",";M_{#gamma#gamma} [MeV];Events",100,1.5,301.5);

    for (Long64_t ievt=0; ievt<nentries; ++ievt) {
        t->GetEntry(ievt);

        // Get the event's vertex
        if (!primaryX||primaryX->empty()) continue;
        TVector3 vertex((*primaryX)[0],(*primaryY)[0],(*primaryZ)[0]);

        std::vector<Hit> hits;
        if (!energies || energies->empty()) {
            std::cout << "Event " << ievt << ": No hits (skipping)" << std::endl;
            continue;
        }
        size_t nHits = energies->size();
        // std::cout << "Event " << ievt << ": " << nHits << " total hits across all rings" << std::endl;
        for (size_t k=0; k<nHits; ++k) {
            hits.push_back({(*centerXs)[k], (*centerYs)[k], (*centerZs)[k], (*energies)[k]});
        }

        std::vector<Hit> trueHits;
        size_t nTrueHits = truthE->size();
        std::cout << "[DEBUG] Event " << ievt << ": " << nTrueHits << " total true hits" << std::endl;
        for (size_t k=0; k<nTrueHits; ++k) {
            trueHits.push_back({(*truthPosX)[k], (*truthPosY)[k], (*truthPosZ)[k], (*truthE)[k]});
        }

        auto truePhotons = TruePhotonBuilder(trueHits, vertex);

        if (truePhotons.size() == 2) {
            TLorentzVector diphoton = truePhotons[0].p4 + truePhotons[1].p4;
            // std::cout << "[DIPHOTON] " << diphoton.M() << std::endl;
            hPi0TrueMass->Fill(diphoton.M());
        }

        // DEBUG print bin contents 
        // hPi0TrueMass->Print("all");


        std::vector<Cluster> clusters;
        double dEta = 0.30/ 2; 
        double dPhi = 0.30/2;
        double E_seed = 20.00;
        double E_neighbor = 0.03;
        int winSize = 3; 
        
        double totalE_Evt = 0;

        clusters = SlidingWindowClusterHits(hits, vertex, dEta, dPhi, E_seed, E_neighbor, winSize);

        // Apply cluster energy threshold
        clusters.erase(std::remove_if(clusters.begin(), clusters.end(),
                                    [](const Cluster &c){ return c.p4.E() < 50.0; }),
                    clusters.end());

        for (size_t ci=0; ci<clusters.size(); ++ci) {
            hClusterE->Fill(clusters[ci].p4.E());
        }
        hNClusters->Fill(clusters.size());

        for (size_t a=0;a<clusters.size();++a) {
            for (size_t b=a+1;b<clusters.size();++b) {
                TLorentzVector pi0 = clusters[a].p4 + clusters[b].p4;
                hPi0Mass->Fill(pi0.M());
            }
        }
    }


    // RECO-RECO
    // PrettyPi0MassPlot(hPi0Mass, "Pi0Mass_Clustered.png");
    
    // TRUTH-TRUTH
    // TruthPi0MassPlot(hPi0TrueMass, "Pi0Mass_Truth.png");
    
    // CLUSTER NUM &/or DEBUG PLOTS
    // PrettyPi0NumClusterPlot(hNClusters);
    // BasicHistPlot(hSingleClusterE);

    //
    // PrettyPi0MassPlot(hPi0Mass, "Pi0Mass_Clustered.png");


    // delete c1; 
    delete hPi0Mass; delete hClusterE; delete hNClusters; delete hPi0TrueMass; 
    f->Close(); delete f;
    return 0;
}