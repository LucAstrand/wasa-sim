#include "EventVariables.hpp"

// EventVariables ComputeEventVariables(
//     const int& eventNumber,
//     const RecoEvent& reco,
//     const TVector3& vertex,
//     const ChargedKECalibration& calibration
// ) {
EventVariables ComputeEventVariables(
    const int& eventNumber,
    const RecoEvent& reco,
    const Vtx& vertex,
    const ChargedKECalibration& calibration
) {
    EventVariables ev;

    if (vertex.has) ev.hasVertex = true;
    
    //count related variables
    for (const auto& cl : reco.chargedClusters) {
        if (cl.isOrphanElectron || cl.isUsedInConversion) continue;
        ev.nChargedTracks++;
        ev.chargedEnergy += cl.totalEnergy; // cl.EdepSmeared; IDK do we want this?
        if (cl.pidGuess == PID::Pion) ev.nPionCandidates++;
        if (cl.pidGuess == PID::Proton) ev.nProtonCandidates++;
    }
    ev.nNeutralClusters = reco.clusters.size();
    ev.nTotalObjects = ev.nChargedTracks + ev.nNeutralClusters;
    ev.EM_energy = reco.EM_energy;
    double calibratedKE = calibration.GetMeanKE(reco.chargedClusters.size(), ev.EM_energy);
    ev.correctedEnergy = ev.EM_energy + calibratedKE;
    ev.totalRecoEnergy = reco.EM_energy; // + ev.chargedEnergy; // this is just not correct... charged energy is included in EM energy
    ev.nPi0Candidates = reco.nPionMultiplicity;

    //proxy sphericity using direction * totalE instead of momentum vector --> this should still give us spatial information about the event

    TMatrixDSym S(3);
    double norm = 0.0;
    auto addToSphericity = [&](const TVector3& p) {
        double w = p.Mag();
        if (w < 1e-9) return;
        S(0,0) += w * p.X() * p.X(); S(0,1) += w * p.X() * p.Y();
        S(0,2) += w * p.X() * p.Z(); S(1,1) += w * p.Y() * p.Y();
        S(1,2) += w * p.Y() * p.Z(); S(2,2) += w * p.Z() * p.Z();
        norm += w * p.Mag2(); 
    };
    //Charged --> Here is where we use the proxy
    for (const auto& cl : reco.chargedClusters) {
        if (cl.isOrphanElectron || cl.isUsedInConversion) continue;
        addToSphericity(cl.direction * cl.totalEnergy);
    }
    //Neutral --> No proxy since we have full p4 vector
    for (const auto& cl : reco.clusters) {
        addToSphericity(cl.p4.Vect());
    }
    S(1,0)=S(0,1); S(2,0)=S(0,2); S(2,1)=S(1,2);
    if (norm > 0) {
        S *= (1.0/norm);
        TVectorD eigenVals(3);
        S.EigenVectors(eigenVals);
        std::vector<double> lam = {eigenVals[0], eigenVals[1], eigenVals[2]};
        std::sort(lam.begin(), lam.end(), std::greater<double>());
        ev.sphericity = 1.5 * (lam[1] + lam[2]); 
    }

    // // max opening angle between tracks --> IDK if this is a really good metric
    // const auto& ccs = reco.chargedClusters;
    // for (size_t i = 0; i < ccs.size(); ++i) {
    //     if (ccs[i].isOrphanElectron || ccs[i].isUsedInConversion) continue;
    //     for (size_t j = i+1; j < ccs.size(); ++j) {
    //         if (ccs[j].isOrphanElectron || ccs[j].isUsedInConversion) continue;
    //         double angle = ccs[i].direction.Angle(ccs[j].direction);
    //         ev.maxTrackAngle = std::max(ev.maxTrackAngle, angle);
    //     }
    // }

    // Vertex radius
    ev.vertexRadius = vertex.vertexVec.Perp(); // the perpendicular component to the Z-axis --> radius on the disk.



    return ev;

}