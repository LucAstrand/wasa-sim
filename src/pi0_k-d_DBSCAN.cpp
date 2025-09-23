#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <queue>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TStyle.h"


// ------------------ Hit & Cluster ------------------
struct Hit {
    double x, y, z, e;
};

struct Cluster {
    std::vector<int> hitIndices;
    TVector3 centroid;
    TLorentzVector p4;
};

// ------------------ KD-Tree ------------------
struct KDNode {
    int index;
    int axis;
    std::unique_ptr<KDNode> left;
    std::unique_ptr<KDNode> right;
};

struct KDTree {
    const std::vector<Hit> *hits;
    std::unique_ptr<KDNode> root;
    double zWeight;

    KDTree(const std::vector<Hit> &h, double zW=1.0) : hits(&h), zWeight(zW) {
        std::vector<int> indices(h.size());
        std::iota(indices.begin(), indices.end(), 0);
        root = build(indices, 0);
    }

    double dist2(const Hit &a, const Hit &b) const {
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        double dz = a.z - b.z;
        return dx*dx + dy*dy + zWeight*dz*dz;
    }

    std::unique_ptr<KDNode> build(std::vector<int> &indices, int depth) {
        if (indices.empty()) return nullptr;
        int axis = depth % 3;

        auto comp = [&](int lhs, int rhs) {
            if (axis == 0) return (*hits)[lhs].x < (*hits)[rhs].x;
            if (axis == 1) return (*hits)[lhs].y < (*hits)[rhs].y;
            return (*hits)[lhs].z < (*hits)[rhs].z;
        };
        std::sort(indices.begin(), indices.end(), comp);

        int mid = indices.size()/2;
        auto node = std::make_unique<KDNode>();
        node->index = indices[mid];
        node->axis = axis;

        std::vector<int> left(indices.begin(), indices.begin()+mid);
        std::vector<int> right(indices.begin()+mid+1, indices.end());

        node->left = build(left, depth+1);
        node->right = build(right, depth+1);
        return node;
    }

    void radiusSearch(const KDNode *node, const Hit &query, double eps, std::vector<int> &results) const {
        if (!node) return;
        const Hit &hit = (*hits)[node->index];
        double eps2 = eps*eps;
        if (dist2(hit, query) <= eps2) results.push_back(node->index);

        double qval = (node->axis == 0 ? query.x : node->axis == 1 ? query.y : query.z);
        double hval = (node->axis == 0 ? hit.x   : node->axis == 1 ? hit.y   : hit.z);
        double diff = qval - hval;

        if (diff <= eps) radiusSearch(node->left.get(), query, eps, results);
        if (diff >= -eps) radiusSearch(node->right.get(), query, eps, results);
    }

    std::vector<int> neighbors(int idx, double eps) const {
        std::vector<int> results;
        radiusSearch(root.get(), (*hits)[idx], eps, results);
        results.erase(std::remove(results.begin(), results.end(), idx), results.end());
        return results;
    }
};

// ------------------ DBSCAN with KD-Tree ------------------
std::vector<Cluster> DBSCANClusterHitsWithKDTree(
    const std::vector<Hit>& hits,
    const TVector3& vertex,
    double eps=12.0,               // neighborhood radius [cm]
    int minPts=3,                  // min neighbors for core point
    double zWeight=1.0,            // scale factor for z-axis
    double minClusterEnergy=30.0,  // discard low-energy clusters
    double minHitEnergyConsidered=0.5) // ignore very small hits
{
    std::vector<Cluster> clusters;
    int N = hits.size();
    if (N == 0) return clusters;

    KDTree tree(hits, zWeight);

    enum class State { Unvisited, Visited, Clustered };
    std::vector<State> state(N, State::Unvisited);

    for (int i = 0; i < N; ++i) {
        if (hits[i].e < minHitEnergyConsidered) continue;
        if (state[i] != State::Unvisited) continue;

        state[i] = State::Visited;
        auto neighbors = tree.neighbors(i, eps);

        if ((int)neighbors.size() + 1 < minPts) {
            continue; // noise
        }

        // start new cluster
        Cluster c;
        std::queue<int> q;
        q.push(i);

        while (!q.empty()) {
            int cur = q.front();
            q.pop();

            if (state[cur] == State::Clustered) continue;
            state[cur] = State::Clustered;
            c.hitIndices.push_back(cur);

            auto neigh = tree.neighbors(cur, eps);
            if ((int)neigh.size() + 1 >= minPts) {
                for (int nb : neigh) {
                    if (state[nb] == State::Unvisited) {
                        state[nb] = State::Visited;
                        q.push(nb);
                    }
                }
            }
        }

        // --- Energy-weighted centroid and momentum ---
        double sumE = 0, meanX=0, meanY=0, meanZ=0;
        TVector3 momentum(0,0,0);
        for (int idx : c.hitIndices) {
            const auto &h = hits[idx];
            double e = h.e;
            sumE += e;
            meanX += h.x * e;
            meanY += h.y * e;
            meanZ += h.z * e;

            TVector3 v(h.x - vertex.X(), h.y - vertex.Y(), h.z - vertex.Z());
            if (v.Mag2() > 1e-12) {
                v = v.Unit();
                momentum += e * v;
            }
        }

        if (sumE < minClusterEnergy) continue; // discard noise cluster

        meanX /= sumE;
        meanY /= sumE;
        meanZ /= sumE;

        c.centroid = TVector3(meanX, meanY, meanZ);
        c.p4.SetPxPyPzE(momentum.X(), momentum.Y(), momentum.Z(), momentum.Mag());

        clusters.push_back(c);
    }

    return clusters;
}

int main(int argc, char **argv) {  // Fixed main signature
    if (argc < 2) {
        printf("Usage: %s <filename>\n", argv[0]);
        return 1;
    }
    
    const char *filename = argv[1];
    
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) return 1;

    TTree *t = (TTree*)f->Get("hibeam");
    if (!t) {
        printf("Tree 'hibeam' not found\n");
        return 1;
    }

    // For each section, we will have separate vectors
    const int NSEC = 16;
    std::vector<double> *posX[NSEC], *posY[NSEC], *posZ[NSEC], *edep[NSEC];

    // Initialize pointers to nullptr
    for (int i = 0; i < NSEC; ++i) {
        posX[i] = nullptr;
        posY[i] = nullptr;
        posZ[i] = nullptr;
        edep[i] = nullptr;
    }

    // Dynamically set branch addresses
    for (int i = 0; i < NSEC; ++i) {
        TString sec = Form("SECE%d", i);
        t->SetBranchAddress(sec + "_GlobalPosition_X", &posX[i]);
        t->SetBranchAddress(sec + "_GlobalPosition_Y", &posY[i]);
        t->SetBranchAddress(sec + "_GlobalPosition_Z", &posZ[i]);
        t->SetBranchAddress(sec + "_EDep", &edep[i]);
    }

    std::vector<double> *primaryX = nullptr, *primaryY = nullptr, *primaryZ = nullptr;
    t->SetBranchAddress("PrimaryPosX", &primaryX);
    t->SetBranchAddress("PrimaryPosY", &primaryY);
    t->SetBranchAddress("PrimaryPosZ", &primaryZ);

    Long64_t nentries = t->GetEntries();
    // TH1F *hPi0Mass = new TH1F("hPi0Mass", "Reconstructed #pi^{0} mass;M_{#gamma#gamma} [GeV];Events", 100, 0, 0.3);
    TH1F *hPi0Mass = new TH1F("hPi0Mass", "Reconstructed #pi^{0} mass;M_{#gamma#gamma} [MeV];Events", 100, 0, 300);


    for (Long64_t i = 0; i < nentries; ++i) {
        // std::cout << "Event: " << i << " / " << nentries << std::endl;
        t->GetEntry(i);

        // Check if primary vectors are valid and non-empty
        if (!primaryX || !primaryY || !primaryZ || 
            primaryX->empty() || primaryY->empty() || primaryZ->empty()) {
            continue;
        }

        std::vector<Hit> hits;
        for (int sec = 0; sec < NSEC; ++sec) {
            if (!edep[sec] || !posX[sec] || !posY[sec] || !posZ[sec]) continue;
            
            for (size_t k = 0; k < edep[sec]->size(); ++k) {
                double E = (*edep[sec])[k];
                // I dont think this is correct, like the hits are very low energy, what we probably really wanted to do is filter on clusters not here
                // if (E < 25) continue; 
                hits.push_back({(*posX[sec])[k], (*posY[sec])[k], (*posZ[sec])[k], E});
            }
        }

        TVector3 vertex((*primaryX)[0], (*primaryY)[0], (*primaryZ)[0]);

        double eps = 15.0;                 // neighborhood radius (cm)
        int    minPts = 3;                 // minimum hits for a core point (neighbours) --> We avoid isolated hits which are probs noise.
        double zWeight = 1.0;              // scale for z-axis
        double minClusterEnergy = 50.0;    // MeV
        double minHitEnergyConsidered = 0.0;

        auto clusters = DBSCANClusterHitsWithKDTree(
                            hits, vertex, eps, minPts, zWeight,
                            minClusterEnergy, minHitEnergyConsidered);

        // Some info to print about the clustering process 
        std::cout << "Found " << clusters.size() << " clusters in event " << i << std::endl;
        for (size_t ci=0; ci<clusters.size(); ++ci) {
            const auto &c = clusters[ci];
            std::cout << " cluster " << ci 
                    << " E = " << c.p4.E()
                    << " nHits = " << c.hitIndices.size()
                    << " pos = (" << c.centroid.X() << "," 
                                << c.centroid.Y() << "," 
                                << c.centroid.Z() << ")\n";
        }

        // Combine photon pairs into pi0 candidates
        for (size_t a = 0; a < clusters.size(); ++a) {
            for (size_t b = a+1; b < clusters.size(); ++b) {
                double angle = clusters[a].p4.Vect().Angle(clusters[b].p4.Vect()) * 180/TMath::Pi();
                std::cout << "Opening angle: " << angle << " deg\n";
                TLorentzVector pi0 = clusters[a].p4 + clusters[b].p4;
                hPi0Mass->Fill(pi0.M());
            }
        }

    }

    TCanvas *c1 = new TCanvas("c1", "Pi0 Mass Distribution", 800, 600);

    // --- Optional: automatically zoom in to region around the peak ---
    int peakBin = hPi0Mass->GetMaximumBin();
    double peakCenter = hPi0Mass->GetBinCenter(peakBin);

    // Zoom to window around peak (here I just use +-50 MeV as a first guess)
    hPi0Mass->GetXaxis()->SetRangeUser(peakCenter - 50, peakCenter + 50);

    // Draw & fit
    hPi0Mass->Draw();
    hPi0Mass->Fit("gaus", "Q"); 

    gStyle->SetOptStat(0); // turn off default stats box
    gStyle->SetOptFit(0);  // turn off auto-fit box

    // Get the fit object and parameters
    TF1 *fit = hPi0Mass->GetFunction("gaus");
    if (fit) {
        double mu = fit->GetParameter(1);
        double sigma = fit->GetParameter(2);
        TLegend *legend = new TLegend(0.55, 0.65, 0.85, 0.85);
        legend->AddEntry(hPi0Mass, "pi^{0} invariant mass", "l");
        legend->AddEntry((TObject*)0, Form("#mu = %.2f MeV", mu), "");
        legend->AddEntry((TObject*)0, Form("#sigma = %.2f MeV", sigma), "");
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
        legend->Draw();
    }

    c1->SaveAs("pi0_mass.png");



    
    // Clean up
    delete c1;
    delete hPi0Mass;
    f->Close();
    delete f;
    
    return 0;
}