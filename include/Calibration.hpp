#ifndef CALIBRATION_H
#define CALIBRATION_H

#include "TH3D.h"
#include "TTree.h"
#include "TFile.h"

#include "Structures.hpp"
#include "RecoEvent.hpp"

class ChargedKECalibration {
public:
    ChargedKECalibration();
    ChargedKECalibration(const std::string& filename);

    void FillCalibration(int Nch, double Eem, double KEtrue);
    void Finalize();
    double GetMeanKE(int Nch, double Eem) const;
    double GetRMSKE(int Nch, double Eem) const;

    void Save(const std::string& filename) const;

private:
    TH3D* h_ke = nullptr; //(Nch, Eem, KEtrue)

};

void DoCalibration(
    TTree* t,
    Long64_t nentries,
    const std::vector<int>& pi0_per_event,
    const std::vector<double>* centerXs,
    const std::vector<double>* centerYs,
    const std::vector<double>* centerZs,
    const std::vector<double>* energies,
    const std::vector<double>* primaryX,
    const std::vector<double>* primaryY,
    const std::vector<double>* primaryZ,
    const std::vector<double>* primaryEkin,
    const std::vector<double>* TPC_firstPosX,
    const std::vector<double>* TPC_firstPosY,
    const std::vector<double>* TPC_firstPosZ,
    const std::vector<double>* TPC_lastPosX,
    const std::vector<double>* TPC_lastPosY,
    const std::vector<double>* TPC_lastPosZ,
    const std::vector<double>* TPC_TrueKE,
    const std::vector<double>*    TPC_pdg,
    const std::vector<double>* TPC_dEdx,
    const std::vector<double>* TPC_smearedEdep,
    const std::vector<double>* TPC_PathLength,
    const std::string& outFile
);

#endif