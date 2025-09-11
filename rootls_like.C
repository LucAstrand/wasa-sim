void rootls_like(const char *filename)
{
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    TList *keys = f->GetListOfKeys();
    if (!keys) {
        std::cerr << "No keys found in file" << std::endl;
        return;
    }

    keys->Sort();

    for (int i = 0; i < keys->GetEntries(); ++i) {
        TKey *key = (TKey *) keys->At(i);
        TDatime dt = key->GetDatime();

        printf("%-20s %-15s %04d-%02d-%02d %02d:%02d\n",
               key->GetName(),
               key->GetClassName(),
               dt.GetYear(), dt.GetMonth(), dt.GetDay(),
               dt.GetHour(), dt.GetMinute());

        // If it's a TTree, inspect its branches
        if (strcmp(key->GetClassName(), "TTree") == 0) {
            TTree *t = (TTree *)key->ReadObj();
            TObjArray *branches = t->GetListOfBranches();

            for (int j = 0; j < branches->GetEntries(); ++j) {
                TBranch *br = (TBranch *)branches->At(j);
                printf("  └── %-30s (%s)\n",
                       br->GetName(),
                       br->GetClassName());
            }
        }
    }

    delete f;
}
