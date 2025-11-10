#define DEBUG  // Comment this out to silence debug prints

#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <vector>
#include <iostream>
#include <TString.h>  // For TString and Form

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.root output.root" << std::endl;
        return 1;
    }

    // Hardcode barrel geometry arrays from Geobuilder code
    std::vector<double> r = {31.77, 37.75, 42.15, 44.57, 46.38, 47.75, 48.05, 48.15, 48.0, 48.15, 48.05, 47.75, 47.27, 46.65, 42.86, 38.50, 33.93};
    std::vector<double> dz = {-39.28, -36.16, -32.10, -27.12, -21.69, -15.91, -9.74, -3.60, 2.50, 8.60, 14.74, 20.91, 27.15, 33.50, 37.47, 40.07, 41.81};
    const size_t nRings = 17;
    const size_t nCrystalsPerRing = 48;

    // Open input file and tree
    TFile* inFile = new TFile(argv[1], "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Failed to open input file: " << argv[1] << std::endl;
        return 1;
    }
    TTree* inTree = (TTree*)inFile->Get("hibeam");
    if (!inTree) {
        std::cerr << "Failed to find 'hibeam' tree in input file." << std::endl;
        return 1;
    }

    // --- Precompute all hardcoded detector cell centers ---
    std::vector<double> detCenterX, detCenterY, detCenterZ;
    for (size_t ring = 0; ring < nRings; ++ring) {
        for (size_t copyNo = 0; copyNo < nCrystalsPerRing; ++copyNo) {
            double phi = 3.75 + copyNo * 7.5;
            if (phi > 180.0) phi -= 360.0;
            double phiRad = phi * TMath::DegToRad();

            double cx = r[ring] * TMath::Sin(phiRad);
            double cy = -r[ring] * TMath::Cos(phiRad);
            double cz = dz[ring];

            detCenterX.push_back(cx);
            detCenterY.push_back(cy);
            detCenterZ.push_back(cz);
        }
    }
    std::cout << "Precomputed " << detCenterX.size() << " detector cell centers." << std::endl;

    std::vector<int>* LayerID[nRings];
    std::vector<double>* edep[nRings];
    // std::vector<double>* truthPosX[nRings];
    // std::vector<double>* truthPosY[nRings];
    // std::vector<double>* truthPosZ[nRings];      
    for (int i = 0; i < nRings; ++i) {
        LayerID[i] = nullptr;
        edep[i] = nullptr;
        // truthPosX[i] = nullptr;
        // truthPosY[i] = nullptr;
        // truthPosZ[i] = nullptr;
    }
    bool allBranchesSet = true;
    for (int i = 0; i < nRings; ++i) {
        TString sec = Form("SECE%d", i);
        TString layerBranch = sec + "_LayerID";
        TString edepBranch = sec + "_EDep";
        if (inTree->GetBranch(layerBranch.Data()) == nullptr || inTree->GetBranch(edepBranch.Data()) == nullptr) {
            std::cerr << "Missing branch: " << layerBranch.Data() 
                      << " or " << edepBranch.Data() 
                      << std::endl;
            allBranchesSet = false;
        } else {
            inTree->SetBranchAddress(layerBranch.Data(), &LayerID[i]);
            inTree->SetBranchAddress(edepBranch.Data(), &edep[i]);
#ifdef DEBUG
            std::cout << "Set branch for ring " << i << ": " << layerBranch.Data() << " and " << edepBranch.Data() << std::endl;
#endif
        }
    }
    if (!allBranchesSet) {
        std::cerr << "Not all branches found—check tree structure with inTree->Print()." << std::endl;
        return 1;
    }

    // Set addresses for primary vertex branches (vectors per event)
    std::vector<double>* primaryPosX = nullptr;
    std::vector<double>* primaryPosY = nullptr;
    std::vector<double>* primaryPosZ = nullptr;
    std::vector<double>* primaryEkin = nullptr;

    std::vector<double>* primaryMomX = nullptr;
    std::vector<double>* primaryMomY = nullptr;
    std::vector<double>* primaryMomZ = nullptr;

    std::vector<double>* TruePhotonX = nullptr;
    std::vector<double>* TruePhotonY = nullptr;
    std::vector<double>* TruePhotonZ = nullptr;
    std::vector<double>* TruePhotonE = nullptr;

    if (!inTree->GetBranch("PrimaryPosX") || !inTree->GetBranch("PrimaryPosY") || !inTree->GetBranch("PrimaryPosZ")) {
        std::cerr << "Missing one or more PrimaryPos* branches—check tree structure." << std::endl;
        return 1;
    }
    inTree->SetBranchAddress("PrimaryPosX", &primaryPosX);
    inTree->SetBranchAddress("PrimaryPosY", &primaryPosY);
    inTree->SetBranchAddress("PrimaryPosZ", &primaryPosZ);
    inTree->SetBranchAddress("PrimaryEkin", &primaryEkin);
    inTree->SetBranchAddress("TruePhotonX", &TruePhotonX);
    inTree->SetBranchAddress("PrimaryMomX", &primaryMomX);
    inTree->SetBranchAddress("PrimaryMomY", &primaryMomY);
    inTree->SetBranchAddress("PrimaryMomZ", &primaryMomZ);
    
    inTree->SetBranchAddress("TruePhotonY", &TruePhotonY);
    inTree->SetBranchAddress("TruePhotonZ", &TruePhotonZ);
    inTree->SetBranchAddress("TruePhotonE", &TruePhotonE);
    

    // Quick test read of first entry to catch issues early
    if (inTree->GetEntries() > 0) {
        Long64_t testEntry = 0;
        std::cout << "Testing first entry..." << std::endl;
        inTree->GetEntry(testEntry);
        if (primaryPosX && primaryPosY && primaryPosZ && !primaryPosX->empty()) {
            std::cout << "Test: Primary vertex for event 0: (" << (*primaryPosX)[0] << ", " << (*primaryPosY)[0] << ", " << (*primaryPosZ)[0] << ")" << std::endl;
        } else {
            std::cout << "Test: Primary vertex for event 0: Empty or not set" << std::endl;
        }
        int testHits = 0;
        for (int ring = 0; ring < nRings; ++ring) {
            if (LayerID[ring] && !LayerID[ring]->empty()) {
                std::cout << "Test: Ring " << ring << " LayerID size=" << LayerID[ring]->size() 
                          << ", first value=" << (*LayerID[ring])[0] << std::endl;
                testHits += LayerID[ring]->size();
            }
            if (edep[ring] && !edep[ring]->empty()) {
                std::cout << "Test: Ring " << ring << " EDep size=" << edep[ring]->size() 
                          << ", first value=" << (*edep[ring])[0] << std::endl;
            }
        }
        std::cout << "Test: First event has " << testHits << " total hits across rings." << std::endl;
        if (testHits == 0) {
            std::cerr << "WARNING: First event has no hits—check if data is in expected branches!" << std::endl;
        }
    } else {
        std::cerr << "Input tree has 0 entries—nothing to process." << std::endl;
        return 1;
    }

    // Create output file and tree (one entry per event, with vectors for hits + vertex vectors)
    TFile* outFile = new TFile(argv[2], "RECREATE");
    TTree* outTree = new TTree("digitizedHits", "Digitized calorimeter");
    std::vector<double> centerXs, centerYs, centerZs, energies;
    std::vector<int> ringNos, copyNos;
    std::vector<double> outPrimaryPosX, outPrimaryPosY, outPrimaryPosZ, outPrimaryEkin, outPrimaryMomX, outPrimaryMomY, outPrimaryMomZ;
    std::vector<double> outTruePhotonX, outTruePhotonY, outTruePhotonZ, outTruePhotonE;
    // std::vector<double> residualX, residualY, residualZ; // --- Create output tree branches (add residuals) ---
    outTree->Branch("centerX", &centerXs);
    outTree->Branch("centerY", &centerYs);
    outTree->Branch("centerZ", &centerZs);
    outTree->Branch("energy", &energies);
    outTree->Branch("ringNo", &ringNos);  // Optional: preserve ringNo
    outTree->Branch("copyNo", &copyNos);
    outTree->Branch("PrimaryPosX", &outPrimaryPosX);
    outTree->Branch("PrimaryPosY", &outPrimaryPosY);
    outTree->Branch("PrimaryPosZ", &outPrimaryPosZ);
    outTree->Branch("PrimaryEkin", &outPrimaryEkin);
    outTree->Branch("PrimaryMomX", &outPrimaryMomX);
    outTree->Branch("PrimaryMomY", &outPrimaryMomY);
    outTree->Branch("PrimaryMomZ", &outPrimaryMomZ);

    outTree->Branch("truthPosX", &outTruePhotonX);
    outTree->Branch("truthPosY", &outTruePhotonY);
    outTree->Branch("truthPosZ", &outTruePhotonZ);
    outTree->Branch("truthE", &outTruePhotonE);
    // outTree->Branch("residualX", &residualX);
    // outTree->Branch("residualY", &residualY);
    // outTree->Branch("residualZ", &residualZ);

    // Process each entry (event)
    Long64_t nEntries = inTree->GetEntries();
    std::cout << "Processing " << nEntries << " events..." << std::endl;
    Long64_t totalHits = 0;
    int emptyEvents = 0;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        if (entry % 100 == 0) std::cout << "Event " << entry << std::endl;
        inTree->GetEntry(entry);

        // // Copy primary vertex vectors for this event
        // if (primaryPosX) outPrimaryPosX = *primaryPosX;
        // else outPrimaryPosX.clear();
        // if (primaryPosY) outPrimaryPosY = *primaryPosY;
        // else outPrimaryPosY.clear();
        // if (primaryPosZ) outPrimaryPosZ = *primaryPosZ;
        // else outPrimaryPosZ.clear();

        // // Copy truth-level hit positions 
        // if (TruePhotonX) outTruePhotonX = *TruePhotonX;
        // else outTruePhotonX.clear();
        // if (TruePhotonY) outTruePhotonY = *TruePhotonY;
        // else outTruePhotonY.clear();
        // if (TruePhotonZ) outTruePhotonZ = *TruePhotonZ;
        // else outTruePhotonZ.clear();

        // for (int ring = 0; ring < nRings; ++ring) {
        //     if (!truthPosX[ring] || truthPosX[ring]->empty()) continue;

        //     // Append all hits in this ring
        //     outTruthPosX.insert(outTruthPosX.end(), truthPosX[ring]->begin(), truthPosX[ring]->end());
        //     outTruthPosY.insert(outTruthPosY.end(), truthPosY[ring]->begin(), truthPosY[ring]->end());
        //     outTruthPosZ.insert(outTruthPosZ.end(), truthPosZ[ring]->begin(), truthPosZ[ring]->end());
        // }

        outPrimaryPosX = primaryPosX ? *primaryPosX : std::vector<double>();
        outPrimaryPosY = primaryPosY ? *primaryPosY : std::vector<double>();
        outPrimaryPosZ = primaryPosZ ? *primaryPosZ : std::vector<double>();
        outPrimaryEkin = primaryEkin ? *primaryEkin : std::vector<double>();
        outPrimaryMomX = primaryMomX ? *primaryMomX : std::vector<double>();
        outPrimaryMomY = primaryMomY ? *primaryMomY : std::vector<double>();
        outPrimaryMomZ = primaryMomZ ? *primaryMomZ : std::vector<double>();

        outTruePhotonX = TruePhotonX ? *TruePhotonX : std::vector<double>();
        outTruePhotonY = TruePhotonY ? *TruePhotonY : std::vector<double>();
        outTruePhotonZ = TruePhotonZ ? *TruePhotonZ : std::vector<double>();
        outTruePhotonE = TruePhotonE ? *TruePhotonE : std::vector<double>();


        // // Clear per-event residuals
        // residualX.clear();
        // residualY.clear();
        // residualZ.clear();

        // // --- Loop over all TruePhotons in this event ---
        // if (TruePhotonX && !TruePhotonX->empty()) {
        //     for (size_t i = 0; i < TruePhotonX->size(); ++i) {
        //         double tx = (*TruePhotonX)[i];
        //         double ty = (*TruePhotonY)[i];
        //         double tz = (*TruePhotonZ)[i];

        //         double minDist2 = 1e12;
        //         double closestCX = 0, closestCY = 0, closestCZ = 0;

        //         // Find nearest detector cell center
        //         for (size_t j = 0; j < detCenterX.size(); ++j) {
        //             double dx = tx - detCenterX[j];
        //             double dy = ty - detCenterY[j];
        //             double dz_ = tz - detCenterZ[j];
        //             double dist2 = dx*dx + dy*dy + dz_*dz_;
        //             if (dist2 < minDist2) {
        //                 minDist2 = dist2;
        //                 closestCX = detCenterX[j];
        //                 closestCY = detCenterY[j];
        //                 closestCZ = detCenterZ[j];
        //             }
        //         }

        //         // Compute residuals with respect to the nearest cell
        //         residualX.push_back(tx - closestCX);
        //         residualY.push_back(ty - closestCY);
        //         residualZ.push_back(tz - closestCZ);
        //     }
        // }

        // outTree->Fill();
    // }

        // Clear vectors for this event
        centerXs.clear();
        centerYs.clear();
        centerZs.clear();
        energies.clear();
        ringNos.clear();
        copyNos.clear();

        int eventHits = 0;
        for (int ring = 0; ring < nRings; ++ring) {
            if (!LayerID[ring] || !edep[ring] || LayerID[ring]->empty()) {
#ifdef DEBUG
                std::cout << "Debug: Skipping empty ring " << ring << " in event " << entry << std::endl;
#endif
                continue;
            }
            if (LayerID[ring]->size() != edep[ring]->size()) {
                std::cerr << "Mismatch in vector sizes for ring " << ring << " in event " << entry 
                          << " (LayerID: " << LayerID[ring]->size() << ", EDep: " << edep[ring]->size() << ")" << std::endl;
                continue;
            }
#ifdef DEBUG
            std::cout << "Debug: Ring " << ring << " in event " << entry << " has " << LayerID[ring]->size() << " hits" << std::endl;
#endif
            for (size_t j = 0; j < LayerID[ring]->size(); ++j) {
                int copyNo = (*LayerID[ring])[j];  // Direct assign: vector<int>
                double energy = (*edep[ring])[j];

                // Compute phi (matches builder code)
                double phi = 3.75 + copyNo * 7.5;  // Assuming copyNo is 0-based; adjust if 1-based
                if (phi > 180.0) phi -= 360.0;

                // Convert to radians
                double phiRad = phi * TMath::DegToRad();

                // Compute center
                double centerX = r[ring] * TMath::Sin(phiRad);
                double centerY = -r[ring] * TMath::Cos(phiRad);  // Negative as per original code
                double centerZ = dz[ring];

                // Add to vectors
                centerXs.push_back(centerX);
                centerYs.push_back(centerY);
                centerZs.push_back(centerZ);
                energies.push_back(energy);
                ringNos.push_back(ring);
                copyNos.push_back(copyNo);

                eventHits++;
                totalHits++;
            }
        }

        // Fill the output tree once per event (with all hits in vectors + vertex vectors)
        outTree->Fill();
#ifdef DEBUG
        if (eventHits == 0) {
            std::cout << "Debug: Event " << entry << " had 0 hits (empty vectors filled)" << std::endl;
            emptyEvents++;
        } else {
            std::cout << "Debug: Event " << entry << " filled with " << eventHits << " hits (vectors size=" << energies.size() << ")" << std::endl;
        }
#endif
    }
#ifdef DEBUG
    std::cout << "Debug: " << emptyEvents << "/" << nEntries << " events were empty." << std::endl;
#endif

    // Write and clean up
    outFile->cd();
    outTree->Write();
    std::cout << "Output tree has " << outTree->GetEntries() << " entries (" << totalHits << " hits processed) from " << nEntries << " events." << std::endl;
    if (totalHits == 0) {
        std::cerr << "ERROR: No hits processed—double-check branch types/names or simulation output." << std::endl;
    }
    outFile->Close();
    inFile->Close();
    delete outFile;
    delete inFile;
    return 0;
}