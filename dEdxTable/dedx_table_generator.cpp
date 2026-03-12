// dedx_table_generator.cc
// Compile with:
// g++ dedx_table_generator.cc -o dedx_table_generator \
//     $(geant4-config --cflags --libs) \
//     $(root-config --cflags --libs)

#include <iostream>
#include <vector>
#include <string>
#include <cmath>

// GEANT4
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4VModularPhysicsList.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4EmCalculator.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

// ROOT
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"

// ----------------------------------------------------------------
// Minimal detector construction - just a box of TPC gas
// G4EmCalculator needs *something* to be constructed
// ----------------------------------------------------------------
class MinimalDetector : public G4VUserDetectorConstruction {
public:
    MinimalDetector(G4Material* tpcMat) : fTPCMaterial(tpcMat) {}

    G4VPhysicalVolume* Construct() override {
        G4NistManager* nist = G4NistManager::Instance();
        G4Material* vacuum  = nist->FindOrBuildMaterial("G4_Galactic");

        // World volume - must be vacuum/low density
        G4Box* worldBox = new G4Box("world", 2*m, 2*m, 2*m);
        G4LogicalVolume* worldLV = new G4LogicalVolume(
            worldBox, vacuum, "world");
        G4VPhysicalVolume* worldPV = new G4PVPlacement(
            nullptr, G4ThreeVector(), worldLV, "world", nullptr, false, 0);

        // TPC gas volume placed inside world
        // This ensures EM tables are built for the TPC material
        G4Box* tpcBox = new G4Box("tpc", 1*m, 1*m, 1*m);
        G4LogicalVolume* tpcLV = new G4LogicalVolume(
            tpcBox, fTPCMaterial, "tpc");
        new G4PVPlacement(
            nullptr, G4ThreeVector(), tpcLV, "tpc", worldLV, false, 0);

        return worldPV;
    }

private:
    G4Material* fTPCMaterial;
};

class MinimalActionInit : public G4VUserActionInitialization {
public:
    void Build() const override {}  // nothing needed
};

// Then add before Initialize():

// ----------------------------------------------------------------
// Minimal physics list - only need EM physics for dEdx tables
// ----------------------------------------------------------------
class MinimalPhysics : public G4VModularPhysicsList {
public:
    MinimalPhysics() {
        RegisterPhysics(new G4EmStandardPhysics());
    }
};

// ----------------------------------------------------------------
// Main
// ----------------------------------------------------------------
int main(int argc, char** argv)
{
    std::string outputFile  = "dedx_tables.root";
    std::string materialName = "Ar80CO2";

    if (argc > 1) outputFile    = argv[1];
    if (argc > 2) materialName  = argv[2];

    std::cout << "Output file:  " << outputFile   << std::endl;
    std::cout << "Material:     " << materialName << std::endl;

    // ----------------------------------------------------------------
    // MUST create run manager FIRST before touching any G4 classes
    // ----------------------------------------------------------------
    auto* runManager = G4RunManagerFactory::CreateRunManager(
        G4RunManagerType::Serial);

    if (!runManager) {
        std::cerr << "ERROR: Failed to create run manager" << std::endl;
        return 1;
    }
    std::cout << "Run manager created OK" << std::endl;

    // std::cout << "Getting NistManager..." << std::endl;
    // G4NistManager* nist = G4NistManager::Instance();
    // G4Material* tpcMaterial = nullptr;

    // std::cout << "NistManager OK" << std::endl;

    // std::cout << "Building G4_Ar..." << std::endl;
    // G4Material* ar = nist->FindOrBuildMaterial("G4_Ar");
    // std::cout << "G4_Ar OK: " << (ar ? "valid" : "null") << std::endl;

    // std::cout << "Building G4_CARBON_DIOXIDE..." << std::endl;
    // G4Material* co2 = nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
    // std::cout << "CO2 OK: " << (co2 ? "valid" : "null") << std::endl;

    // std::cout << "Building Ar80CO2..." << std::endl;
    // tpcMaterial = new G4Material("Ar80CO2", 1.823e-3 * g/cm3, 2);
    // std::cout << "G4Material created OK" << std::endl;
    // tpcMaterial->AddMaterial(ar,  0.80);
    // std::cout << "Added Ar OK" << std::endl;
    // tpcMaterial->AddMaterial(co2, 0.20);
    // std::cout << "Added CO2 OK" << std::endl;

    // ----------------------------------------------------------------
    // NOW build material
    // ----------------------------------------------------------------
    G4NistManager* nist = G4NistManager::Instance();

    G4Material* tpcMaterial = nullptr;

    if (materialName == "Ar80CO2") {
        tpcMaterial = new G4Material("Ar80CO2",
                                     1.823e-3 * g/cm3,
                                     2);
        tpcMaterial->AddMaterial(nist->FindOrBuildMaterial("G4_Ar"),  0.80);
        tpcMaterial->AddMaterial(nist->FindOrBuildMaterial("G4_CARBON_DIOXIDE"), 0.20);
    }
    else {
        tpcMaterial = nist->FindOrBuildMaterial(materialName);
    }

    if (!tpcMaterial) {
        std::cerr << "ERROR: Could not build material: "
                  << materialName << std::endl;
        delete runManager;
        return 1;
    }

    std::cout << "Material density: "
              << tpcMaterial->GetDensity()/(g/cm3)
              << " g/cm3" << std::endl;

    // ----------------------------------------------------------------
    // Initialise with minimal detector and physics
    // ----------------------------------------------------------------
    runManager->SetUserInitialization(new MinimalDetector(tpcMaterial));
    runManager->SetUserInitialization(new MinimalPhysics());
    runManager->SetUserInitialization(new MinimalActionInit());
    std::cout << "Initializing (building EM tables)..." << std::endl;
    runManager->Initialize();
    std::cout << "Initialized OK - starting table generation" << std::endl;
    runManager->BeamOn(0);
    std::cout << "BeamOn(0) completed - processes fully instantiated" << std::endl;

    // ----------------------------------------------------------------
    // Particles to tabulate
    // ----------------------------------------------------------------
    struct ParticleEntry {
        std::string g4name;
        std::string label;
        int         pdg;
    };

    std::vector<ParticleEntry> particles = {
        {"pi+",    "pion_plus",   211},
        {"pi-",    "pion_minus", -211},
        {"proton", "proton",     2212},
        {"mu-",    "muon_minus",   13},
        {"mu+",    "muon_plus",   -13},
        {"e-",     "electron",    11},
        {"e+",     "positron",   -11},
    };

    // ----------------------------------------------------------------
    // Energy range: log-spaced from 0.1 MeV to 10 GeV
    // ----------------------------------------------------------------
    const int    nPoints = 500;
    const double Emin    = 0.01;      // MeV
    const double Emax    = 1000.0;  // MeV

    // ----------------------------------------------------------------
    // Compute and write tables
    // ----------------------------------------------------------------
    G4EmCalculator emCalc;
    G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();

    TFile* outFile = TFile::Open(outputFile.c_str(), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "ERROR: Cannot open output file: " 
                  << outputFile << std::endl;
        return 1;
    }

    for (auto& entry : particles) {
        G4ParticleDefinition* particle = partTable->FindParticle(entry.g4name);
        if (!particle) {
            std::cerr << "WARNING: Particle not found: " 
                      << entry.g4name << std::endl;
            continue;
        }

        TGraph* gTotal    = new TGraph(nPoints);
        TGraph* gElec     = new TGraph(nPoints);
        TGraph* gNuclear  = new TGraph(nPoints);

        gTotal  ->SetName(("dedx_"         + entry.label).c_str());
        gElec   ->SetName(("dedx_elec_"    + entry.label).c_str());
        gNuclear->SetName(("dedx_nuclear_" + entry.label).c_str());

        gTotal  ->SetTitle((entry.label + " total;KE [MeV];dEdx [MeV/cm]").c_str());
        gElec   ->SetTitle((entry.label + " electronic;KE [MeV];dEdx [MeV/cm]").c_str());
        gNuclear->SetTitle((entry.label + " nuclear;KE [MeV];dEdx [MeV/cm]").c_str());

        for (int i = 0; i < nPoints; i++) {
            double KE = Emin * std::pow(Emax/Emin, (double)i/(nPoints-1));
            
            // if (i == 0) {
            //     std::cout << "First ComputeTotalDEDX call: "
            //             << entry.g4name << " KE=" << KE << " MeV" << std::endl;
            // }

            double dedx_total = emCalc.ComputeTotalDEDX(
                KE * MeV, particle, tpcMaterial) / (MeV/cm);

            // if (i == 0) {
            //     std::cout << "First dedx_total = " << dedx_total << std::endl;
            // }

            double dedx_elec = emCalc.ComputeElectronicDEDX(
                KE * MeV, particle, tpcMaterial) / (MeV/cm);

            double dedx_nuclear = emCalc.ComputeNuclearDEDX(
                KE * MeV, particle, tpcMaterial) / (MeV/cm);

            gTotal  ->SetPoint(i, KE, dedx_total);
            gElec   ->SetPoint(i, KE, dedx_elec);
            gNuclear->SetPoint(i, KE, dedx_nuclear);
        }

        outFile->cd();
        gTotal  ->Write();
        gElec   ->Write();
        gNuclear->Write();

        std::cout << "Written tables for " << entry.label 
                  << " (pdg=" << entry.pdg << ")" << std::endl;
    }

    // // Also write a metadata tree so the analysis knows what material was used
    // TTree* meta = new TTree("metadata", "DEDX table metadata");
    // double density = tpcMaterial->GetDensity() / (g/cm3);
    // std::string matName = tpcMaterial->GetName();
    // meta->Branch("material",  &matName);
    // meta->Branch("density",   &density);
    // meta->Branch("Emin_MeV",  &Emin);
    // meta->Branch("Emax_MeV",  &Emax);
    // meta->Branch("nPoints",   const_cast<int*>(&nPoints));
    // meta->Fill();
    // meta->Write();

    outFile->Close();

    std::cout << "Done. Tables written to: " << outputFile << std::endl;

    delete runManager;
    return 0;
}