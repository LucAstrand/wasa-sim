#ifndef SELECTIONCUTS_H
#define SELECTIONCUTS_H

// #include "EventVariables.hpp"

// struct SelectionCuts {
//     bool hasVertex            = true; 
//     int minChargedTracks      = 2;
//     double minTotalEnergy     = 300.0; //--> MeV
//     int minPi0Candidates      = 1;     // could also have a min pi+- candidates etc...
//     // double minSphericity      = 0.3;   // we need to establish if this, since we use proxy momentum is okay! 

// };

// bool PassesSelection(const EventVariables& ev,
//                      const SelectionCuts& cuts)
// {
//     if (ev.hasVertex        != cuts.hasVertex) return false;
//     if (ev.nChargedTracks   < cuts.minChargedTracks) return false;
//     if (ev.totalRecoEnergy  < cuts.minTotalEnergy) return false;
//     if (ev.nPi0Candidates   < cuts.minPi0Candidates) return false;
//     // if (ev.sphericity       < cuts.minSphericity) return false;
//     return true;
// }

#include "EventVariables.hpp"

struct SelectionCuts {
    bool   requireVertex     = true;
    int    minChargedTracks  = 2;
    double minTotalEnergy    = 300.0; // MeV
    int    minPi0Candidates  = 1;
    // double minSphericity  = 0.3;
};

// ── Named cut predicates ──────────────────────────────────────────────────────
// Each lambda tests ONE condition only, independently.
// The cutflow engine applies them cumulatively — no manual chaining here.

using CutPredicate = std::function<bool(const EventVariables&)>;

std::vector<std::pair<std::string, CutPredicate>>
BuildCutList(const SelectionCuts& c)
{
    return {
        { "hasVertex",
          [&](const EventVariables& ev) {
              return !c.requireVertex || ev.hasVertex; }},

        { "nCharged >= " + std::to_string(c.minChargedTracks),
          [&](const EventVariables& ev) {
              return ev.nChargedTracks >= c.minChargedTracks; }},

        { "Etotal > " + std::to_string((int)c.minTotalEnergy) + " MeV",
          [&](const EventVariables& ev) {
              return ev.totalRecoEnergy > c.minTotalEnergy; }},

        { "nPi0 >= " + std::to_string(c.minPi0Candidates),
          [&](const EventVariables& ev) {
              return ev.nPi0Candidates >= c.minPi0Candidates; }},

        // { "sphericity > " + std::to_string(c.minSphericity),
        //   [&](const EventVariables& ev) {
        //       return ev.sphericity > c.minSphericity; }},
    };
}

// PassesSelection just runs all cuts in sequence — single source of truth
bool PassesSelection(const EventVariables& ev, const SelectionCuts& cuts)
{
    for (auto& [name, fn] : BuildCutList(cuts))
        if (!fn(ev)) return false;
    return true;
}

// struct CutflowResult {
//     std::string cutName;
//     int nPassSig   = 0;
//     int nPassBkg   = 0;
//     double effSig  = 0.0;  // fraction of total signal passing
//     double effBkg  = 0.0;  // fraction of total background passing
// };

// struct Cutflow {
//     std::vector<CutflowResult> cuts;
//     int nTotalSig = 0;
//     int nTotalBkg = 0;

//     void Init(int totalSig, int totalBkg) {
//         nTotalSig = totalSig;
//         nTotalBkg = totalBkg;
//     }

//     void Evaluate(const std::vector<EventVariables>& sigEvents,
//                   const std::vector<EventVariables>& bkgEvents,
//                   const SelectionCuts& c) {

//         std::vector<std::pair<std::string, std::function<bool(const EventVariables&)>>> cutList = {

//             {"hasVertex != " + std::to_string(c.hasVertex),
//                 [&](const EventVariables& ev){ 
//                     return ev.hasVertex != c.hasVertex; }},  

//             {"nCharged >= " + std::to_string(c.minChargedTracks),
//                 [&](const EventVariables& ev){ 
//                     return ev.hasVertex != c.hasVertex &&
//                            ev.nChargedTracks >= c.minChargedTracks; }},  

//             {"Etotal > " + std::to_string((int)c.minTotalEnergy) + " MeV",
//                 [&](const EventVariables& ev){ 
//                     return ev.hasVertex != c.hasVertex &&
//                            ev.nChargedTracks >= c.minChargedTracks &&
//                            ev.totalRecoEnergy > c.minTotalEnergy; }},


//             {"nPi0 >= " + std::to_string(c.minPi0Candidates),
//                 [&](const EventVariables& ev){ 
//                     return ev.hasVertex != c.hasVertex &&
//                            ev.nChargedTracks >= c.minChargedTracks &&
//                            ev.totalRecoEnergy > c.minTotalEnergy &&
//                            ev.nPi0Candidates > c.minPi0Candidates; }},
                        
//             // {"sphericity > " + std::to_string(c.minSphericity),
//             //     [&](const EventVariables& ev){ 
//             //         return PassesSelection(ev, c); }};  // final = all cuts
//             };

//         cuts.clear();
//         for (auto& [name, cutFn] : cutList) {
//             CutflowResult res;
//             res.cutName = name;
//             for (auto& ev : sigEvents) if (cutFn(ev)) res.nPassSig++;
//             for (auto& ev : bkgEvents) if (cutFn(ev)) res.nPassBkg++;
//             res.effSig = (nTotalSig > 0) ? 100.0 * res.nPassSig / nTotalSig : 0.0;
//             res.effBkg = (nTotalBkg > 0) ? 100.0 * res.nPassBkg / nTotalBkg : 0.0;
//             cuts.push_back(res);
//         }
//     }

    struct CutflowResult {
        std::string cutName;
        int    nPassSig = 0, nPassBkg = 0;
        double effSig   = 0.0, effBkg  = 0.0;
    };

    struct Cutflow {
        std::vector<CutflowResult> cuts;
        int nTotalSig = 0, nTotalBkg = 0;

        void Init(int totalSig, int totalBkg) {
            nTotalSig = totalSig;
            nTotalBkg = totalBkg;
        }

        void Evaluate(const std::vector<EventVariables>& sigEvents,
                    const std::vector<EventVariables>& bkgEvents,
                    const SelectionCuts& c)
        {
            auto cutList = BuildCutList(c);
            cuts.clear();

            for (size_t i = 0; i < cutList.size(); ++i) {
                CutflowResult res;
                res.cutName = cutList[i].first;

                // Cumulative: event passes if it passes cuts 0..i
                auto passesCumulative = [&](const EventVariables& ev) {
                    for (size_t j = 0; j <= i; ++j)
                        if (!cutList[j].second(ev)) return false;
                    return true;
                };

                for (auto& ev : sigEvents) if (passesCumulative(ev)) res.nPassSig++;
                for (auto& ev : bkgEvents) if (passesCumulative(ev)) res.nPassBkg++;
                res.effSig = (nTotalSig > 0) ? 100.0 * res.nPassSig / nTotalSig : 0.0;
                res.effBkg = (nTotalBkg > 0) ? 100.0 * res.nPassBkg / nTotalBkg : 0.0;
                cuts.push_back(res);
            }
        }

    void Print() const {
        std::cout << "\n";
        std::cout << std::string(80, '=') << "\n";
        std::cout << std::left
                  << std::setw(35) << "Cut"
                  << std::setw(12) << "N_sig"
                  << std::setw(12) << "eff_sig(%)"
                  << std::setw(12) << "N_bkg"
                  << std::setw(12) << "eff_bkg(%)"
                  << "\n";
        std::cout << std::string(80, '-') << "\n";

        // Print total before cuts
        std::cout << std::left
                  << std::setw(35) << "No cuts"
                  << std::setw(12) << nTotalSig
                  << std::setw(12) << "100.0"
                  << std::setw(12) << nTotalBkg
                  << std::setw(12) << "100.0"
                  << "\n";

        for (const auto& c : cuts) {
            std::cout << std::left
                      << std::setw(35) << c.cutName
                      << std::setw(12) << c.nPassSig
                      << std::setw(12) << std::fixed << std::setprecision(1) 
                      << c.effSig
                      << std::setw(12) << c.nPassBkg
                      << std::setw(12) << std::fixed << std::setprecision(1) 
                      << c.effBkg
                      << "\n";
        }
        std::cout << std::string(80, '=') << "\n\n";
    }

    // Also save as a ROOT TH1 for visual cutflow plot
    void PlotCutflow(const std::string& outFile) const {
        int nCuts = cuts.size() + 1;  // +1 for "no cuts"

        TH1F* hSig = new TH1F("hCutflowSig", 
            ";Cut;Efficiency [%]", nCuts, 0, nCuts);
        TH1F* hBkg = new TH1F("hCutflowBkg", 
            ";Cut;Efficiency [%]", nCuts, 0, nCuts);

        hSig->SetBinContent(1, 100.0);
        hBkg->SetBinContent(1, 100.0);
        hSig->GetXaxis()->SetBinLabel(1, "No cuts");
        hBkg->GetXaxis()->SetBinLabel(1, "No cuts");

        for (size_t i = 0; i < cuts.size(); ++i) {
            hSig->SetBinContent(i+2, cuts[i].effSig);
            hBkg->SetBinContent(i+2, cuts[i].effBkg);
            hSig->GetXaxis()->SetBinLabel(i+2, cuts[i].cutName.c_str());
            hBkg->GetXaxis()->SetBinLabel(i+2, cuts[i].cutName.c_str());
        }

        PlotOptions opts;
        opts.addLegend = true;
        opts.legendEntries = {"Signal", "Background"};
        hSig->GetXaxis()->LabelsOption("v");  // vertical labels
        hSig->GetXaxis()->SetLabelSize(0.03);
        Plot1D({hSig, hBkg}, {kBlack, kRed}, outFile, opts);

        delete hSig; delete hBkg;
    }
};


#endif