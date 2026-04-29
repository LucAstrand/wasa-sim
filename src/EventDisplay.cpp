#include "EventDisplay.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <unordered_set>

// ---------- SafeSetBranch ----------
template <typename T>
bool SafeSetBranch(TTree* tree, const char* branchName, T*& ptr) {
    if (tree->GetListOfBranches()->FindObject(branchName)) {
        tree->SetBranchAddress(branchName, &ptr);
        return true;
    } else {
        std::cout << "[Warning] Branch " << branchName << " not found. Skipping." << std::endl;
        return false;
    }
}

// ---------- Distinct per-cluster colors ----------
static Color_t MakeDistinctColor(int i, float sat=1.0f, float val=1.0f)
{
    // Golden ratio hue stepping for visually distinct colors
    double hue = fmod(i * 0.61803398875, 1.0);

    float r=0, g=0, b=0;

    float h = (float)(hue * 6.0);
    float c = val * sat;
    float x = c * (1.0f - fabsf(fmodf(h, 2.0f) - 1.0f));
    float m = val - c;

    float rp=0, gp=0, bp=0;
    if      (0 <= h && h < 1) { rp = c; gp = x; bp = 0; }
    else if (1 <= h && h < 2) { rp = x; gp = c; bp = 0; }
    else if (2 <= h && h < 3) { rp = 0; gp = c; bp = x; }
    else if (3 <= h && h < 4) { rp = 0; gp = x; bp = c; }
    else if (4 <= h && h < 5) { rp = x; gp = 0; bp = c; }
    else                      { rp = c; gp = 0; bp = x; }

    r = rp + m; g = gp + m; b = bp + m;

    Color_t idx = TColor::GetFreeColorIndex();
    new TColor(idx, r, g, b); // register in ROOT
    return idx;
}

static inline void MakePerpBasis(const TVector3& dir, TVector3& u, TVector3& v)
{
    // Pick a vector not parallel to dir
    TVector3 a = (std::fabs(dir.Z()) < 0.9) ? TVector3(0,0,1) : TVector3(0,1,0);
    u = dir.Cross(a).Unit();
    v = dir.Cross(u).Unit();
}

// Draw the same cone used by your clustering: angle between (hit - vertex) and direction < thetaMax
TEveElement* DrawMatchCone(const TVector3& apex,
                           const TVector3& dirIn,
                           double thetaMax,      // radians
                           double L,             // length in your detector units
                           Color_t col = kCyan,
                           int nCircle = 48,
                           int nRays   = 12,
                           double alpha = 0.7)   // (0..1) line transparency feel; TEveLine doesn't do true alpha
{
    auto group = new TEveElementList("MatchCone");
    group->SetMainColor(col);

    TVector3 dir = dirIn;
    if (dir.Mag() <= 0) return group;
    dir = dir.Unit();

    TVector3 u, v;
    MakePerpBasis(dir, u, v);

    // Base circle at distance L along axis
    const TVector3 baseC = apex + L * dir;
    const double rBase   = L * std::tan(thetaMax);

    auto base = new TEveLine("ConeBase");
    base->SetMainColor(col);
    base->SetLineWidth(2);

    for (int i = 0; i <= nCircle; ++i) {
        double phi = 2.0 * TMath::Pi() * (double)i / (double)nCircle;
        TVector3 p = baseC + rBase * (std::cos(phi) * u + std::sin(phi) * v);
        base->SetNextPoint(p.X(), p.Y(), p.Z());
    }
    group->AddElement(base);

    // Rays from apex to a few points on the base circle
    for (int k = 0; k < nRays; ++k) {
        double phi = 2.0 * TMath::Pi() * (double)k / (double)nRays;
        TVector3 p = baseC + rBase * (std::cos(phi) * u + std::sin(phi) * v);

        auto ray = new TEveLine(Form("ConeRay_%d", k));
        ray->SetMainColor(col);
        ray->SetLineStyle(2);
        ray->SetLineWidth(1);
        ray->SetNextPoint(apex.X(), apex.Y(), apex.Z());
        ray->SetNextPoint(p.X(), p.Y(), p.Z());
        group->AddElement(ray);
    }

    // Optional: draw axis line for clarity
    auto axis = new TEveLine("ConeAxis");
    axis->SetMainColor(col);
    axis->SetLineWidth(2);
    axis->SetNextPoint(apex.X(), apex.Y(), apex.Z());
    axis->SetNextPoint(baseC.X(), baseC.Y(), baseC.Z());
    group->AddElement(axis);

    return group;
}


// ---------- Draw a reconstructed cluster consistently ----------
TEveElementList* DrawRecoNCluster(const Cluster& c, int cid, Color_t col,
                                bool drawBBox=true, bool drawHitBoxes=false)
{
    auto group = new TEveElementList(Form("NCluster_%d", cid));
    group->SetMainColor(col);

    // ---- Hits: one TEvePointSet per cluster (fast + consistent) ----
    auto hitsPS = new TEvePointSet(Form("NClusterHits_%d", cid));
    hitsPS->SetMarkerStyle(4);
    hitsPS->SetMarkerSize(1.0);
    hitsPS->SetMarkerColor(col);

    for (const Hit* h : c.hits) {
        if (!h) continue;
        hitsPS->SetNextPoint(h->x, h->y, h->z);
    }
    group->AddElement(hitsPS);

    // ---- Centroid ----
    auto cent = new TEvePointSet(Form("Centroid_%d", cid));
    cent->SetMarkerStyle(20);
    cent->SetMarkerSize(1.6);
    cent->SetMarkerColor(col);
    cent->SetNextPoint(c.centroid.X(), c.centroid.Y(), c.centroid.Z());
    group->AddElement(cent);

    // ---- Momentum arrow ----
    TVector3 cc = c.centroid;
    TVector3 dir = c.p4.Vect();
    if (dir.Mag() > 0) dir = dir.Unit();

    double E = c.p4.E();
    double arrowLen = 4.0 + std::log(std::max(1e-6, E));
    TVector3 vec = arrowLen * dir;

    auto arrow = new TEveArrow(vec.X(), vec.Y(), vec.Z(), cc.X(), cc.Y(), cc.Z());
    arrow->SetMainColor(col);
    // arrow->SetLineWidth(2);
    group->AddElement(arrow);

    // ---- Optional: bounding box from cluster hits ----
    if (drawBBox && !c.hits.empty()) {
        double xmin=1e9, xmax=-1e9, ymin=1e9, ymax=-1e9, zmin=1e9, zmax=-1e9;
        bool any=false;

        for (const Hit* h : c.hits) {
            if (!h) continue;
            any = true;
            xmin = std::min(xmin, h->x); xmax = std::max(xmax, h->x);
            ymin = std::min(ymin, h->y); ymax = std::max(ymax, h->y);
            zmin = std::min(zmin, h->z); zmax = std::max(zmax, h->z);
        }

        if (any) {
            auto box = new TEveBox(Form("BBox_%d", cid));
            box->SetMainColor(col);
            box->SetMainTransparency(75);

            box->SetVertex(0, xmin, ymin, zmin);
            box->SetVertex(1, xmax, ymin, zmin);
            box->SetVertex(2, xmax, ymax, zmin);
            box->SetVertex(3, xmin, ymax, zmin);
            box->SetVertex(4, xmin, ymin, zmax);
            box->SetVertex(5, xmax, ymin, zmax);
            box->SetVertex(6, xmax, ymax, zmax);
            box->SetVertex(7, xmin, ymax, zmax);

            group->AddElement(box);
        }
    }

    // ---- Optional: little boxes per hit (cluster "shape") ----
    // Looks nice if hits correspond to cells. Can be heavy if many hits.
    if (drawHitBoxes && !c.hits.empty()) {
        auto bs = new TEveBoxSet(Form("HitBoxes_%d", cid));
        bs->SetMainColor(col);
        bs->SetMainTransparency(40);
        bs->Reset(TEveBoxSet::kBT_AABox, true, 64);

        const float half = 0.5f; // tune this to your cell size
        for (const Hit* h : c.hits) {
            if (!h) continue;
            bs->AddBox(h->x - half, h->y - half, h->z - half,
                         2*half, 2*half, 2*half);
        }
        bs->RefitPlex();
        group->AddElement(bs);
    }

    return group;
}

TEveElementList* DrawRecoCCluster(const ChargedCluster& c, int cid, Color_t col,
                                bool drawBBox=true, bool drawHitBoxes=false)
{
    std::string particle_label = PIDToString(c.pidGuess);
    auto group = new TEveElementList(Form("CCluster_%s_%d", particle_label.c_str(), cid));
    group->SetMainColor(col);

    // ---- Hits: one TEvePointSet per cluster (fast + consistent) ----
    auto hitsPS = new TEvePointSet(Form("CClusterHits_%d", cid));
    hitsPS->SetMarkerStyle(4);
    hitsPS->SetMarkerSize(1.0);
    hitsPS->SetMarkerColor(col);

    for (const Hit* h : c.hits) {
        if (!h) continue;
        if (h->owner == HitOwner::Charged) {
            hitsPS->SetNextPoint(h->x, h->y, h->z);
        }
    }
    group->AddElement(hitsPS);

    // // ---- Centroid ----
    // auto cent = new TEvePointSet(Form("Direction_%d", cid));
    // cent->SetMarkerStyle(20);
    // cent->SetMarkerSize(1.6);
    // cent->SetMarkerColor(col);
    // cent->SetNextPoint(c.direction.X(), c.direction.Y(), c.direction.Z());
    // group->AddElement(cent);

    // ---- Momentum arrow ----
    // TVector3 cc = c.centroid;
    TVector3 dir = c.direction;
    TVector3 exitPoint = c.TPCExitPoint;
    if (dir.Mag() > 0) dir = dir.Unit();

    double E = c.totalEnergy;
    double arrowLen = 4.0 + std::log(std::max(1e-6, E));
    TVector3 vec = arrowLen * dir;

    auto arrow = new TEveArrow(vec.X(), vec.Y(), vec.Z(), exitPoint.X(), exitPoint.Y(), exitPoint.Z());
    arrow->SetMainColor(col);
    // arrow->SetLineWidth(2);
    group->AddElement(arrow);

    const TVector3 axis = dir.Unit();
    double theta = 25*TMath::DegToRad();      // MUST be radians
    double L = 100.0;             // choose: track length or detector size

    auto coneEl = DrawMatchCone(exitPoint, axis, theta, L, kMagenta);
    group->AddElement(coneEl);


    return group;
}

void DrawPionsAsArrows(const std::vector<primaryChPi>& pions,
                       double x0 = 0, double y0 = 0, double z0 = 0,
                       double scale = 1.0)
{
  // Make sure TEve is running
  if (!gEve) TEveManager::Create();

  auto* list = new TEveElementList("Primary charged pions");
  gEve->AddElement(list);

  for (const auto& pi : pions) {
    // Use momentum as arrow direction
    const double px = pi.p4.Px();
    const double py = pi.p4.Py();
    const double pz = pi.p4.Pz();

    // Normalize direction; use magnitude separately if you like
    const double p = TMath::Sqrt(px*px + py*py + pz*pz);
    if (p <= 0) continue;

    const double dx = (px / p) * (p * scale);
    const double dy = (py / p) * (p * scale);
    const double dz = (pz / p) * (p * scale);

    auto* a = new TEveArrow(dx, dy, dz, x0, y0, z0);
    a->SetName(Form("pi trackID=%d", pi.trackID));

    // Optional styling
    a->SetMainColor(kRed);   // or color by charge/pT/etc.
    a->SetTubeR(0.001);       // shaft radius (world units)
    a->SetConeR(0.0001);
    a->SetConeL(0.0001);

    list->AddElement(a);
  }

  gEve->Redraw3D(kTRUE);
}

TEveElementList* MakePrimaryTracks(
    const std::vector<double>* primaryX,
    const std::vector<double>* primaryY,
    const std::vector<double>* primaryZ,
    const std::vector<double>* primaryPx,
    const std::vector<double>* primaryPy,
    const std::vector<double>* primaryPz,
    const std::vector<double>* primaryEkin,
    const std::vector<int>*    primaryPDG,
    const std::vector<int>*    primaryTrackID,
    double trackLength = 200.0  // cm, how far to draw the arrow
)
{
    TEveElementList* group = new TEveElementList("Primary Tracks");

    if (!primaryX || primaryX->empty()) return group;

    int nPrimaries = primaryX->size();

    for (int i = 0; i < nPrimaries; i++) {

        double x0 = primaryX->at(i);
        double y0 = primaryY->at(i);
        double z0 = primaryZ->at(i);

        double px = primaryPx->at(i);
        double py = primaryPy->at(i);
        double pz = primaryPz->at(i);
        double pmag = std::sqrt(px*px + py*py + pz*pz);

        int pdg = primaryPDG->at(i);
        int tid = primaryTrackID->at(i);
        double ekin = primaryEkin->at(i);

        // Normalised direction
        double dx = 0, dy = 0, dz = 0;
        if (pmag > 0) {
            dx = px / pmag;
            dy = py / pmag;
            dz = pz / pmag;
        }

        // End point of the drawn track segment
        double x1 = x0 + dx * trackLength;
        double y1 = y0 + dy * trackLength;
        double z1 = z0 + dz * trackLength;

        // Use TEveLine for the track
        TEveLine* line = new TEveLine();
        line->SetNextPoint(x0, y0, z0);
        line->SetNextPoint(x1, y1, z1);
        line->SetLineWidth(2);

        // Color by PDG type
        Color_t color = kGray;
        std::string pname = "unknown";
        if      (pdg ==  211) { color = kBlue;    pname = "pi+";     }
        else if (pdg == -211) { color = kCyan;    pname = "pi-";     }
        else if (pdg == 2212) { color = kRed;     pname = "proton";  }
        else if (pdg ==   13) { color = kYellow;  pname = "mu-";     }
        else if (pdg ==  -13) { color = kOrange;  pname = "mu+";     }
        else if (pdg ==   11) { color = kMagenta; pname = "e-";      }
        else if (pdg ==  -11) { color = kPink;    pname = "e+";      }
        else if (pdg ==  111) { color = kGreen;   pname = "pi0";     }
        else if (pdg == 2112) { color = kWhite;   pname = "neutron"; }

        line->SetLineColor(color);

        // Descriptive name shown in the Eve browser panel
        char name[256];
        snprintf(name, sizeof(name), 
                 "%s [TID=%d] Ekin=%.1f MeV", 
                 pname.c_str(), tid, ekin);
        line->SetName(name);

        // Optional: add a marker at the vertex position
        TEvePointSet* vtxMarker = new TEvePointSet(1);
        vtxMarker->SetNextPoint(x0, y0, z0);
        vtxMarker->SetMarkerStyle(4);  // open circle
        vtxMarker->SetMarkerSize(1.5);
        vtxMarker->SetMarkerColor(color);
        vtxMarker->SetName("vertex");

        // Group the line and its vertex marker together
        TEveElementList* trackGroup = new TEveElementList(name);
        trackGroup->AddElement(line);
        trackGroup->AddElement(vtxMarker);

        group->AddElement(trackGroup);
    }

    return group;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <data.root> <geometry.root>\n";
        return 1;
    }
    const char* dataFile = argv[1];
    const char* geoFile  = argv[2];

    TGeoManager::Import(geoFile, "geometry");
    if (!gGeoManager) {
        std::cerr << "ERROR: Could not load geometry!\n";
        return 1;
    }

    new TApplication("app", &argc, argv);
    TEveManager::Create();

    // Make detector geometry semi-transparent
    TGeoVolume* topVolume = gGeoManager->GetTopVolume();
    int nd = topVolume->GetNdaughters();
    for (int i = 0; i < nd; i++) {
        TGeoNode* node = topVolume->GetNode(i);
        TGeoVolume* vol = node->GetVolume();
        vol->SetTransparency(80);
        vol->SetLineColor(kGray + 1);
    }

    auto top = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
    top->SetVisLevel(4);
    gEve->AddGlobalElement(top);

    TFile* f = TFile::Open(dataFile);
    if (!f) { std::cerr << "ERROR: cannot open datafile\n"; return 1; }

    TTree* t = (TTree*)f->Get("digitizedHits");
    if (!t) { std::cerr << "ERROR: TTree not found\n"; return 1; }

    // Hit cell centers + energy
    std::vector<double> *centerXs = nullptr, *centerYs = nullptr, *centerZs = nullptr, *energies = nullptr;
    SafeSetBranch(t, "centerX", centerXs);
    SafeSetBranch(t, "centerY", centerYs);
    SafeSetBranch(t, "centerZ", centerZs);
    SafeSetBranch(t, "energy", energies);

    // Primary vertex
    std::vector<double> *primaryX = nullptr, *primaryY = nullptr, *primaryZ = nullptr, *primaryEkin = nullptr;
    std::vector<double> *primaryPx = nullptr, *primaryPy = nullptr, *primaryPz = nullptr;
    std::vector<int> *primaryPDG = nullptr, *primaryTrackID = nullptr;
    SafeSetBranch(t, "PrimaryPosX", primaryX);
    SafeSetBranch(t, "PrimaryPosY", primaryY);
    SafeSetBranch(t, "PrimaryPosZ", primaryZ);
    SafeSetBranch(t, "PrimaryEkin", primaryEkin);
    SafeSetBranch(t, "PrimaryMomX", primaryPx);
    SafeSetBranch(t, "PrimaryMomY", primaryPy);
    SafeSetBranch(t, "PrimaryMomZ", primaryPz);
    SafeSetBranch(t, "PrimaryPDG", primaryPDG);
    SafeSetBranch(t, "PrimaryTrackID", primaryTrackID);

    // Truth info (photons)
    std::vector<double> *truePhotonPosX = nullptr, *truePhotonPosY = nullptr, *truePhotonPosZ = nullptr, *truePhotonE = nullptr;
    std::vector<int> *truePhotonTrackID = nullptr, *truePhotonParentID = nullptr;
    SafeSetBranch(t, "truthPosX", truePhotonPosX);
    SafeSetBranch(t, "truthPosY", truePhotonPosY);
    SafeSetBranch(t, "truthPosZ", truePhotonPosZ);
    SafeSetBranch(t, "truthE", truePhotonE);
    SafeSetBranch(t, "TruePhotonTrackID", truePhotonTrackID);
    SafeSetBranch(t, "TruePhotonParentID", truePhotonParentID);

    std::vector<double>* TruePhotonCreationX = nullptr;
    std::vector<double>* TruePhotonCreationY = nullptr;
    std::vector<double>* TruePhotonCreationZ = nullptr;
    std::vector<double>* TruePhotonEndX = nullptr;
    std::vector<double>* TruePhotonEndY = nullptr;
    std::vector<double>* TruePhotonEndZ = nullptr;

    SafeSetBranch(t, "TruePhotonCreationX", TruePhotonCreationX);
    SafeSetBranch(t, "TruePhotonCreationY", TruePhotonCreationY);
    SafeSetBranch(t, "TruePhotonCreationZ", TruePhotonCreationZ);
    SafeSetBranch(t, "TruePhotonEndX", TruePhotonEndX);
    SafeSetBranch(t, "TruePhotonEndY", TruePhotonEndY);
    SafeSetBranch(t, "TruePhotonEndZ", TruePhotonEndZ);

    // TPC info
    // std::vector<double> *TPC_trackID = nullptr, *TPC_Edep = nullptr, *TPC_smearedEdep = nullptr;
    // std::vector<double> *TPC_firstPosX = nullptr, *TPC_firstPosY = nullptr, *TPC_firstPosZ = nullptr;
    // std::vector<double> *TPC_lastPosX = nullptr, *TPC_lastPosY = nullptr, *TPC_lastPosZ = nullptr;
    // std::vector<double> *TPC_PathLength = nullptr, *TPC_dEdx = nullptr, *TPC_TrueKE = nullptr, *TPC_pdg = nullptr;
    // SafeSetBranch(t, "TPC_trackID", TPC_trackID);
    // SafeSetBranch(t, "TPC_Edep", TPC_Edep);
    // SafeSetBranch(t, "TPC_smearedEdep", TPC_smearedEdep);
    // SafeSetBranch(t, "TPC_firstPosX", TPC_firstPosX);
    // SafeSetBranch(t, "TPC_firstPosY", TPC_firstPosY);
    // SafeSetBranch(t, "TPC_firstPosZ", TPC_firstPosZ);
    // SafeSetBranch(t, "TPC_lastPosX", TPC_lastPosX);
    // SafeSetBranch(t, "TPC_lastPosY", TPC_lastPosY);
    // SafeSetBranch(t, "TPC_lastPosZ", TPC_lastPosZ);
    // SafeSetBranch(t, "TPC_PathLength", TPC_PathLength);
    // SafeSetBranch(t, "TPC_dEdx", TPC_dEdx);
    // SafeSetBranch(t, "TPC_TrueKE", TPC_TrueKE);
    // SafeSetBranch(t, "TPC_pdg", TPC_pdg);
    std::vector<int>    *TPC_trackID     = nullptr;
    std::vector<int>    *TPC_pdg         = nullptr;
    std::vector<double> *TPC_TrueKE      = nullptr;
    std::vector<double> *TPC_firstPosX   = nullptr;
    std::vector<double> *TPC_firstPosY   = nullptr;
    std::vector<double> *TPC_firstPosZ   = nullptr;
    std::vector<double> *TPC_lastPosX    = nullptr;
    std::vector<double> *TPC_lastPosY    = nullptr;
    std::vector<double> *TPC_lastPosZ    = nullptr;
    std::vector<double> *TPC_smearedDedx = nullptr;
    std::vector<double> *TPC_theoryDedx  = nullptr;

    SafeSetBranch(t, "TPC_trackID",    TPC_trackID);
    SafeSetBranch(t, "TPC_pdg",        TPC_pdg);
    SafeSetBranch(t, "TPC_TrueKE",     TPC_TrueKE);
    SafeSetBranch(t, "TPC_firstPosX",  TPC_firstPosX);
    SafeSetBranch(t, "TPC_firstPosY",  TPC_firstPosY);
    SafeSetBranch(t, "TPC_firstPosZ",  TPC_firstPosZ);
    SafeSetBranch(t, "TPC_lastPosX",   TPC_lastPosX);
    SafeSetBranch(t, "TPC_lastPosY",   TPC_lastPosY);
    SafeSetBranch(t, "TPC_lastPosZ",   TPC_lastPosZ);
    SafeSetBranch(t, "TPC_smearedDedx", TPC_smearedDedx);
    SafeSetBranch(t, "TPC_theoryDedx", TPC_theoryDedx);

    Long64_t nEvents = 1; // or: t->GetEntries();

    std::vector<primaryChPi> primaryChPis;
    DEDXTable dedxTABLE("dedx_tables_Ar80CO2.root");

    for (Long64_t ev = 0; ev < nEvents; ev++) {
        t->GetEntry(ev);

        // Build event list
        auto eventList = new TEveElementList(Form("Event_%lld", ev));
        gEve->AddElement(eventList);

        // ---------- Event color (only used for truth/charged tracks here) ----------
        // Keep eventColor distinct, but clusters will have their own colors.
        Color_t eventColor = MakeDistinctColor((int)ev, 0.6f, 1.0f);

        // ---------- Draw photons (truth) ----------
        if (TruePhotonCreationX && TruePhotonEndX) {
            for (size_t i = 0; i < TruePhotonCreationX->size(); i++) {
                auto photon = new TEveLine(Form("Photon_%lld_%zu", ev, i));
                photon->SetLineWidth(3);
                photon->SetMainColor(eventColor);

                photon->SetPoint(0,
                    (*TruePhotonCreationX)[i],
                    (*TruePhotonCreationY)[i],
                    (*TruePhotonCreationZ)[i]);

                photon->SetPoint(1,
                    (*TruePhotonEndX)[i],
                    (*TruePhotonEndY)[i],
                    (*TruePhotonEndZ)[i]);

                eventList->AddElement(photon);
            }
        }

        // ---------- Draw charged tracks ----------
        if (TPC_firstPosX && TPC_lastPosX) {
            for (size_t i = 0; i < TPC_firstPosX->size(); i++) {
                auto charged = new TEveLine(Form("charged_%lld_%zu", ev, i));
                charged->SetLineWidth(3);
                charged->SetMainColor(eventColor);

                charged->SetPoint(0,
                    (*TPC_firstPosX)[i],
                    (*TPC_firstPosY)[i],
                    (*TPC_firstPosZ)[i]);

                charged->SetPoint(1,
                    (*TPC_lastPosX)[i],
                    (*TPC_lastPosY)[i],
                    (*TPC_lastPosZ)[i]);

                eventList->AddElement(charged);
            }
        }

        // ---------- Build hits ----------
        std::vector<Hit> hits;
        if (!centerXs || !centerYs || !centerZs || !energies) continue;

        const size_t nHits = energies->size();
        hits.reserve(nHits);
        for (size_t k = 0; k < nHits; ++k) {
            hits.push_back({(*centerXs)[k], (*centerYs)[k], (*centerZs)[k], (*energies)[k]});
        }

        // Primary vertex
        TVector3 vertex(0,0,0);
        if (primaryX && !primaryX->empty()) {
            vertex = TVector3((*primaryX)[0], (*primaryY)[0], (*primaryZ)[0]);
        }

        size_t nPrimaries = primaryPDG->size();
        // reset the vetors as they are filled on an event level
        primaryChPis.clear();
        for (size_t i=0; i<nPrimaries; ++i) {
            
            //init event value
            // chPi_per_event.push_back(0);

            if (std::abs((*primaryPDG)[i]) == 211) { // pi+-
                // chPi_per_event[ievt]++;
                primaryChPis.push_back({(*primaryTrackID)[i], TLorentzVector(TVector3((*primaryPx)[i], (*primaryPy)[i], (*primaryPz)[i]), (*primaryEkin)[i])});
            }
        }

        // ---------- Build chargedTracks ----------
        std::vector<ChargedTrack> chargedTracks;
        if (TPC_firstPosX && TPC_lastPosX && TPC_TrueKE && TPC_pdg) {
            chargedTracks.reserve(TPC_pdg->size());
            for (size_t k = 0; k < TPC_pdg->size(); ++k) {
                TVector3 first((*TPC_firstPosX)[k], (*TPC_firstPosY)[k], (*TPC_firstPosZ)[k]);
                TVector3 last((*TPC_lastPosX)[k],  (*TPC_lastPosY)[k],  (*TPC_lastPosZ)[k]);
                TVector3 dir = last - first;

                // chargedTracks.push_back({
                //     (*TPC_trackID)[k],
                //     vertex,
                //     last,
                //     dir,
                //     (*TPC_TrueKE)[k],
                //     (*TPC_pdg)[k],
                //     (*TPC_dEdx)[k],
                //     (*TPC_smearedEdep)[k],
                //     (*TPC_PathLength)[k],
                //     0.0,   // dEdxTheory placeholder
                //     0.15   // resolution
                // });
                ChargedTrack trk;
                trk.id            = k;
                trk.vertex        = vertex;
                trk.exitPoint     = TVector3((*TPC_lastPosX)[k], (*TPC_lastPosY)[k], (*TPC_lastPosZ)[k]);
                trk.direction     = trk.exitPoint
                                - TVector3((*TPC_firstPosX)[k], (*TPC_firstPosY)[k], (*TPC_firstPosZ)[k]);
                trk.direction     = trk.direction.Unit();
                trk.TrueKE        = (*TPC_TrueKE)[k];
                trk.TruePDG       = (*TPC_pdg)[k];
                trk.smearedDedx   = (*TPC_smearedDedx)[k];
                trk.theoryDedx    = (*TPC_theoryDedx)[k];
                trk.resolution    = 0.15;
                chargedTracks.push_back(trk);
            }
        }

        // ---------- Reconstruct ----------
        RecoEvent reco = ReconstructEvent(hits, chargedTracks, vertex, dedxTABLE);

        // ---------- Draw ALL hits as gray background ----------
        {
            auto allHits = new TEvePointSet(Form("AllHits_%lld", ev));
            allHits->SetMarkerStyle(4);
            allHits->SetMarkerSize(0.6);
            allHits->SetMarkerColor(kGray + 1);
            for (size_t i = 0; i < centerXs->size(); ++i) {
                allHits->SetNextPoint(centerXs->at(i), centerYs->at(i), centerZs->at(i));
            }
            eventList->AddElement(allHits);
        }

        // ---------- candidate selection ----------
        double param_h = 135;
        double param_k = 1.6;
        double param_a = 45;
        double param_b = 1.55;

        std::vector<Pi0Candidate> candidates;
        std::vector<Pi0Candidate> selected;

        for (size_t a = 0; a < reco.clusters.size(); ++a) {
            for (size_t b = a + 1; b < reco.clusters.size(); ++b) {
                const auto& g1 = reco.clusters[a].p4;
                const auto& g2 = reco.clusters[b].p4;

                double mgg = (g1 + g2).M();
                double theta = openingAngle(g1, g2);

                double relipse =
                    std::pow((mgg - param_h), 2) / std::pow(param_a, 2) +
                    std::pow((theta - param_k), 2) / std::pow(param_b, 2);

                if (relipse > 1) continue;

                Pi0Candidate cand;
                cand.c1 = &reco.clusters[a];
                cand.c2 = &reco.clusters[b];
                cand.mgg = mgg;
                cand.theta = theta;
                cand.p4 = g1 + g2;
                cand.score = relipse;

                candidates.push_back(cand);
            }
        }

        std::sort(candidates.begin(), candidates.end(),
            [](const Pi0Candidate& a, const Pi0Candidate& b) { return a.score < b.score; });

        std::unordered_set<const Cluster*> used;
        for (const auto& cand : candidates) {
            if (used.count(cand.c1) || used.count(cand.c2)) continue;
            used.insert(cand.c1);
            used.insert(cand.c2);
            selected.push_back(cand);
        }

        // ---------- Helper: cluster pointer -> index ----------
        auto clusterIndex = [&](const Cluster* cptr)->int {
            for (size_t i = 0; i < reco.clusters.size(); ++i) {
                if (&reco.clusters[i] == cptr) return (int)i;
            }
            return -1;
        };

        // ---------- DRAWING OPTIONS ----------
        const bool drawAllNClusters = false;     // set false if you only want selected pi0 clusters
        const bool drawBBox         = true;
        const bool drawHitBoxes     = false;    // heavy if many hits; try true for a few clusters
        const bool drawAllCClusters = true;
        const bool drawPrimaryPions = true;

        if (drawAllNClusters) {
            // Draw every cluster with its own unique color
            for (size_t ci = 0; ci < reco.clusters.size(); ++ci) {
                Color_t ccol = MakeDistinctColor((int)ci);
                auto clEve = DrawRecoNCluster(reco.clusters[ci], (int)ci, ccol, drawBBox, drawHitBoxes);
                eventList->AddElement(clEve);
            }
        } else {
            // Draw only selected pi0 clusters, each with its own unique color
            for (const Pi0Candidate& pi0 : selected) {
                int i1 = clusterIndex(pi0.c1);
                int i2 = clusterIndex(pi0.c2);

                Color_t col1 = MakeDistinctColor(i1 >= 0 ? i1 : 1000);
                Color_t col2 = MakeDistinctColor(i2 >= 0 ? i2 : 1001);

                // For selected, you might want hit boxes ON to see the shape:
                eventList->AddElement(DrawRecoNCluster(*pi0.c1, i1, col1, /*drawBBox=*/true, /*drawHitBoxes=*/true));
                eventList->AddElement(DrawRecoNCluster(*pi0.c2, i2, col2, /*drawBBox=*/true, /*drawHitBoxes=*/true));
            }
        }

        if (drawAllCClusters) {
            for (size_t ch = 0; ch < reco.chargedClusters.size(); ++ch) {
                Color_t chcol = MakeDistinctColor((int)ch);
                auto chclEve = DrawRecoCCluster(reco.chargedClusters[ch], (int)ch, chcol, drawBBox, drawHitBoxes);
                eventList->AddElement(chclEve);
            }
        }

        if (drawPrimaryPions) {

            DrawPionsAsArrows(primaryChPis, vertex.X(), vertex.Y(), vertex.Z(), 1);
        }        

        TEveElementList* primaries = MakePrimaryTracks(
            primaryX, primaryY, primaryZ,
            primaryPx, primaryPy, primaryPz,
            primaryEkin, primaryPDG, primaryTrackID
        );
        gEve->AddElement(primaries);
        // gEve->Redraw3D(kTRUE);

        // Done with this event
    }

    gEve->Redraw3D(kTRUE);
    gApplication->Run(kTRUE);
    return 0;
}