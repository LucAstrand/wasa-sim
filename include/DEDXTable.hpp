#ifndef DEDXTABLE_H
#define DEDXTABLE_H

#include "TFile.h"
#include "TGraph.h"
#include <map>
#include <string>
#include <stdexcept>
#include <iostream>


class DEDXTable {

public:
    // DEDXTable(const DEDXTable&) = delete;
    // DEDXTable& operator=(const DEDXTable&) = delete;

    DEDXTable(const std::string& rootFile) {
        fFile = TFile::Open(rootFile.c_str(), "READ");
        if (!fFile || fFile->IsZombie())
            throw std::runtime_error("Cannot open DEDX table: " + rootFile);
        
        // Map PDG to graph name in file
        std::map<int, std::string> nameMap = {
            { 211,  "dedx_pion_plus"},
            {-211,  "dedx_pion_minus"},
            {2212,  "dedx_proton"},
            {  13,  "dedx_muon_minus"},
            { -13,  "dedx_muon_plus"},
            {  11,  "dedx_electron"},
            { -11,  "dedx_positron"}
        };

        for (auto& [pdg, name] : nameMap) {
            TGraph* g = (TGraph*)fFile->Get(name.c_str());
            if (!g) {
                std::cerr << "[DEDXTable] WARNING: missing " << name << std::endl;
                continue;
            }
            // g->SetDirectory(nullptr);
            fGraphs[pdg] = g;
        }
        fFile->Close();
    }

    ~DEDXTable() {
        for (auto& [pdg, g] : fGraphs) delete g;
    }

    // Get method --> returns dEdx in MeV/cm for given PDG and KE in MeV
    // returns -1 if PDG is not found
    double GetDEDX(int pdg, double KE_MeV) const {
        auto it = fGraphs.find(pdg);
        if (it == fGraphs.end() || !it->second) {
            std::cerr << "[DEDXTable] No table for PDG " << pdg << std::endl;
            return -1;
        }
        return it->second->Eval(KE_MeV);
    }

    // For convenience
    double Pion (double KE) const {return GetDEDX(211, KE); }
    double Proton (double KE) const {return GetDEDX(2212, KE); }
    double Muon (double KE) const {return GetDEDX(13, KE); }
    double Electron (double KE) const {return GetDEDX(11, KE); }


private:
    TFile* fFile = nullptr;
    std::map<int, TGraph*> fGraphs;
};


#endif