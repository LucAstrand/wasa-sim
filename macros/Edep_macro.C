
void EDep_macro(const char* filename)
{ 
    //Open the file
    auto *f = TFile::Open(filename);
    if (!f || f->IsZombie()){
        std::cerr << "Error: Cannot open file" << filename << std::endl;
        return;
    }
    //Extract the tree
    auto t = dynamic_cast<TTree*>(f->Get("hibeam"));
    if (!t) {
        std::cerr << "Error: TTRee hibeam not found in" << filename << std::endl;
        return;
    }
    //set the branches
    const int n_seces{16};
    std::vector<double>* SECE_EDep[n_seces];
    std::vector<int>* SECE_Pdg[n_seces];

    for (int i=0; i<n_seces; i++) {
        SECE_EDep[i] = nullptr;
        SECE_Pdg[i] = nullptr;
        std::string EDepName = "SECE" + std::to_string(i) + "_EDep";
        std::string PdgName = "SECE" + std::to_string(i) + "_PDG";

        t->SetBranchAddress(EDepName.c_str(), &SECE_EDep[i]);
        t->SetBranchAddress(PdgName.c_str(), &SECE_Pdg[i]);
    }

    //Make the histograms
    TH1F *hem = new TH1F("EMCAL", "Photon Energy", 50,0,500);
    hem->SetDirectory(0);
    TH1F *hpdg = new TH1F("PDG", "PDS", 50,10,400);

    //Fill the histograms

    double sum = 0;
    double energy_dep = 0;

    for (int j=0;j < t-> GetEntries(); j++) {

        t->GetEntry(j);
        for (int k=0; k < n_seces; k++) {
            if (size(*SECE_EDep[k]) == 0) continue;
            sum = std::reduce((*SECE_EDep[k]).begin(), (*SECE_EDep[k]).end());
            energy_dep += sum;
        }
        if (energy_dep != 0) {
            hem->Fill(energy_dep);
        }
        energy_dep = 0;
    }

    hem->Draw();
    for (int i = 0; i < n_seces; ++i) delete SECE_EDep[i];
    delete f;
}

