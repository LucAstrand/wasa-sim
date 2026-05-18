#include "EventDisplay.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "TEveProjectionManager.h"
#include "TEveProjectionBases.h"

// ============================================================
//  Utility: safe TTree branch binding
// ============================================================
template <typename T>
bool SafeSetBranch(TTree* tree, const char* name, T*& ptr)
{
    if (tree->GetListOfBranches()->FindObject(name)) {
        tree->SetBranchAddress(name, &ptr);
        return true;
    }
    std::cout << "[Warning] Branch '" << name << "' not found.\n";
    return false;
}

// ============================================================
//  Colour palette
// ============================================================
namespace Palette
{
    static Color_t Reg(float r, float g, float b)
    {
        Color_t idx = TColor::GetFreeColorIndex();
        new TColor(idx, r, g, b);
        return idx;
    }

    static std::vector<Color_t> MakeClusterPalette()
    {
        return {
            Reg(0.22f, 0.62f, 0.98f),   // azure
            Reg(0.98f, 0.52f, 0.10f),   // amber
            Reg(0.22f, 0.88f, 0.44f),   // green
            Reg(0.95f, 0.18f, 0.28f),   // coral
            Reg(0.78f, 0.42f, 0.95f),   // violet
            Reg(0.08f, 0.82f, 0.82f),   // teal
            Reg(0.97f, 0.88f, 0.12f),   // gold
            Reg(0.92f, 0.50f, 0.72f),   // rose
            Reg(0.38f, 0.92f, 0.72f),   // mint
            Reg(0.96f, 0.64f, 0.38f),   // peach
        };
    }

    static Color_t ForPDG(int pdg)
    {
        switch (pdg) {
            case  211: return Reg(0.20f, 0.60f, 1.00f);
            case -211: return Reg(0.00f, 0.85f, 0.85f);
            case 2212: return Reg(1.00f, 0.25f, 0.25f);
            case   13: return Reg(1.00f, 0.85f, 0.10f);
            case  -13: return Reg(1.00f, 0.65f, 0.10f);
            case   11: return Reg(0.90f, 0.20f, 0.80f);
            case  -11: return Reg(0.95f, 0.55f, 0.85f);
            case  111: return Reg(0.30f, 1.00f, 0.30f);
            case 2112: return Reg(0.80f, 0.80f, 0.80f);
            case   22: return Reg(1.00f, 1.00f, 0.50f);
            default:   return Reg(0.55f, 0.55f, 0.55f);
        }
    }

    static std::string NameForPDG(int pdg)
    {
        switch (pdg) {
            case  211: return "pi+";
            case -211: return "pi-";
            case 2212: return "p";
            case   13: return "mu-";
            case  -13: return "mu+";
            case   11: return "e-";
            case  -11: return "e+";
            case  111: return "pi0";
            case 2112: return "n";
            case   22: return "gamma";
            default:   return Form("PDG%d", pdg);
        }
    }
}

// ============================================================
//  Geometry helpers
// ============================================================

// Walk the current gGeoManager mother stack upward looking for
// the first node whose volume name starts with 'prefix'.
// Must be called immediately after FindNode() so the stack is valid.
static TGeoNode* FindAncestorByPrefix(const char* prefix)
{
    int level = gGeoManager->GetLevel();
    for (int lv = 0; lv <= level; ++lv) {
        TGeoNode* node = gGeoManager->GetMother(lv);
        if (!node) continue;
        const char* vname = node->GetVolume()->GetName();
        if (std::string(vname).rfind(prefix, 0) == 0)
            return node;
    }
    return nullptr;
}

// Retrieve the world-space centre of the current node
// (call right after FindNode, before moving the navigator).
static TVector3 CurrentNodeWorldCentre()
{
    double local[3] = {0, 0, 0};
    double world[3] = {0, 0, 0};
    gGeoManager->LocalToMaster(local, world);
    return TVector3(world[0], world[1], world[2]);
}

// Reset every volume in the tree to a dim, semi-transparent state.
static void ResetGeometryColors(Color_t dimCol, int transp)
{
    if (!gGeoManager) return;
    TGeoVolume* top = gGeoManager->GetTopVolume();
    if (!top) return;
    std::function<void(TGeoVolume*, int)> walk = [&](TGeoVolume* vol, int d) {
        if (!vol || d > 8) return;
        vol->SetLineColor(dimCol);
        vol->SetFillColor(dimCol);
        vol->SetTransparency(transp);
        for (int i = 0; i < vol->GetNdaughters(); ++i)
            walk(vol->GetNode(i)->GetVolume(), d + 1);
    };
    walk(top, 0);
}

// ============================================================
//  Energy-spike parameters  (tune to your detector units / MeV)
// ============================================================
static constexpr double SPIKE_SCALE  = 0.20;   // length per MeV  [cm/MeV]
static constexpr double SPIKE_MIN    = 2.0;    // minimum spike length [cm]
static constexpr double SPIKE_MAX    = 55.0;   // maximum spike length [cm]

// ============================================================
//  Core cluster drawing
//  For each hit:
//    1. FindNode → walk stack to SECE cell node
//    2. Recolour that geometry volume
//    3. Accumulate energy per cell
//    4. Draw one radial spike per cell, length ~ sum(E_hits)
//    5. Place a centroid star at the energy-weighted centre
// ============================================================
static TEveElementList* DrawClusterCells(
    const std::vector<Hit*>& hits,
    int                       cid,
    const char*               label,
    Color_t                   col,
    bool                      chargedOnly = false)
{
    auto group = new TEveElementList(label);
    group->SetMainColor(col);

    if (!gGeoManager) return group;
    std::unordered_map<TGeoNode*, double>   cellE;
    std::unordered_map<TGeoNode*, TVector3> cellCtr;
    // hits that couldn't be matched to a SECE node
    std::vector<const Hit*> unmatched;

    for (Hit* h : hits) {
        if (!h) continue;
        if (chargedOnly && h->owner != HitOwner::Charged) continue;
        TGeoNode* cellNode = gGeoManager->FindNode(h->x, h->y, h->z);
        if (!cellNode) {
            unmatched.push_back(h);
            continue;
        }
        if (std::string(cellNode->GetName()).rfind("SECE", 0) != 0) {
            unmatched.push_back(h);
            continue;
        }
        cellE[cellNode] += h->e;

        if (cellCtr.find(cellNode) == cellCtr.end()) {
            gGeoManager->FindNode(h->x, h->y, h->z);
            double loc[3] = {0,0,0}, wld[3] = {0,0,0};
            gGeoManager->LocalToMaster(loc, wld);
            cellCtr[cellNode] = TVector3(wld[0], wld[1], wld[2]);

            // Restore navigator
            gGeoManager->FindNode(h->x, h->y, h->z);
        }
    }

    // ── Recolour cells + draw spikes ───────────────────────
    auto spikeGroup = new TEveElementList(Form("Spikes_%d", cid));
    spikeGroup->SetMainColor(col);

    TVector3 wCentroid(0, 0, 0);
    double   wSum = 0;

    for (auto& kv : cellE) {
        TGeoNode*       node   = kv.first;
        double          eSum   = kv.second;
        const TVector3& ctr    = cellCtr[node];

        // Radial direction outward from origin
        TVector3 rhat = ctr;
        if (rhat.Mag() > 0) rhat = rhat.Unit();
        else rhat.SetXYZ(0, 0, 1);
        double len =
            std::clamp(
                6.0 + 18.0 * std::log10(1.0 + eSum),
                6.0,
                90.0);
        int width =
            std::clamp(
                (int)(2 + 1.8 * std::log10(1.0 + eSum)),
                2,
                12);
        TVector3 base = ctr + 1.5 * rhat;
        TVector3 tip = ctr + len * rhat;

        auto outer = new TEveLine();
        outer->SetLineWidth(width + 4);
        outer->SetLineColor(kBlack);
        outer->SetMainTransparency(70);
        outer->SetNextPoint(base.X(), base.Y(), base.Z());
        outer->SetNextPoint(tip.X(),  tip.Y(),  tip.Z());
        auto spike = new TEveLine(Form("spike_%.0fMeV", eSum));
        spike->SetLineColor(col);
        spike->SetLineWidth(width);
        spike->SetLineStyle(1);
        // spike->SetMainTransparency(15);
        spike->SetNextPoint(base.X(), base.Y(), base.Z());
        spike->SetNextPoint(tip.X(),  tip.Y(),  tip.Z());
        spikeGroup->AddElement(outer);
        spikeGroup->AddElement(spike);

        // Accumulate energy-weighted centroid
        wCentroid += eSum * ctr;
        wSum += eSum;
    }
    group->AddElement(spikeGroup);

    // ── Fallback points for unmatched hits ─────────────────
    if (!unmatched.empty()) {
        auto ps = new TEvePointSet(Form("UnmatchedHits_%d", cid));
        ps->SetMarkerStyle(kFullDotSmall);
        ps->SetMarkerSize(0.9);
        ps->SetMarkerColor(col);
        for (const Hit* h : unmatched) {
            ps->SetNextPoint(h->x, h->y, h->z);
            wCentroid += TVector3(h->x, h->y, h->z);
            wSum += 1.0;
        }
        group->AddElement(ps);
    }

    // ── Centroid star ──────────────────────────────────────
    if (wSum > 0) {
        wCentroid *= 1.0 / wSum;
        auto centPS = new TEvePointSet(Form("Centroid_%d", cid));
        centPS->SetMarkerStyle(kFullStar);
        centPS->SetMarkerSize(2.8);
        centPS->SetMarkerColor(col);
        centPS->SetNextPoint(wCentroid.X(), wCentroid.Y(), wCentroid.Z());
        group->AddElement(centPS);
    }

    return group;
}

// ── Public wrappers ────────────────────────────────────────

static TEveElementList* DrawNeutralCluster(const Cluster& c, int cid, Color_t col)
{
    std::string lbl = Form("NCluster_%d  E=%.0f MeV  nhits=%zu",
                           cid, c.p4.E(), c.hits.size());
    return DrawClusterCells(c.hits, cid, lbl.c_str(), col, false);
}

static TEveElementList* DrawChargedCluster(const ChargedCluster& c, int cid, Color_t col)
{
    std::string pid = PIDToString(c.pidGuess);
    std::string lbl = Form("CCluster_%s_%d  E=%.0f MeV",
                           pid.c_str(), cid, c.totalEnergy);
    return DrawClusterCells(c.hits, cid, lbl.c_str(), col, true);
}

// ============================================================
//  TPC reconstructed tracks  (straight lines, no B-field)
// ============================================================
static TEveElementList* DrawTPCTracks(
    const std::vector<double>* fX, const std::vector<double>* fY, const std::vector<double>* fZ,
    const std::vector<double>* lX, const std::vector<double>* lY, const std::vector<double>* lZ,
    const std::vector<int>*    pdg,
    const std::vector<double>* trueKE,
    const std::vector<double>* smearedDedx)
{
    auto group = new TEveElementList("TPC Tracks (reco)");
    if (!fX || fX->empty()) return group;

    for (int i = 0; i < (int)fX->size(); ++i) {
        TVector3 p0((*fX)[i], (*fY)[i], (*fZ)[i]);
        TVector3 p1((*lX)[i], (*lY)[i], (*lZ)[i]);

        int    thisPDG = pdg       ? (*pdg)[i]        : 0;
        double ke      = trueKE    ? (*trueKE)[i]     : 0.0;
        double dedx    = smearedDedx ? (*smearedDedx)[i] : -1.0;

        Color_t     col   = Palette::ForPDG(thisPDG);
        std::string pname = Palette::NameForPDG(thisPDG);

        std::ostringstream lbl;
        lbl << "TPC " << pname
            << "  Ekin=" << std::fixed << std::setprecision(0) << ke << " MeV";
        if (dedx >= 0)
            lbl << "  dE/dx=" << std::setprecision(2) << dedx << " MeV/cm";

        auto trkGrp = new TEveElementList(lbl.str().c_str());
        trkGrp->SetMainColor(col);

        auto line = new TEveLine("segment");
        line->SetLineColor(col);
        line->SetLineWidth(3);
        line->SetLineStyle(1);
        line->SetNextPoint(p0.X(), p0.Y(), p0.Z());
        line->SetNextPoint(p1.X(), p1.Y(), p1.Z());
        trkGrp->AddElement(line);

        auto psA = new TEvePointSet("entry");
        psA->SetMarkerStyle(kOpenCircle);
        psA->SetMarkerSize(2.0);
        psA->SetMarkerColor(col);
        psA->SetNextPoint(p0.X(), p0.Y(), p0.Z());
        trkGrp->AddElement(psA);

        auto psB = new TEvePointSet("exit");
        psB->SetMarkerStyle(kFullCircle);
        psB->SetMarkerSize(1.6);
        psB->SetMarkerColor(col);
        psB->SetNextPoint(p1.X(), p1.Y(), p1.Z());
        trkGrp->AddElement(psB);

        group->AddElement(trkGrp);
    }
    return group;
}

// ============================================================
//  Primary truth tracks  (long-dashed, PDG-coloured)
// ============================================================
static TEveElementList* DrawPrimaryTracks(
    const std::vector<double>* pX,  const std::vector<double>* pY,  const std::vector<double>* pZ,
    const std::vector<double>* mX,  const std::vector<double>* mY,  const std::vector<double>* mZ,
    const std::vector<double>* ekin,
    const std::vector<int>*    pdg,
    const std::vector<int>*    trackID,
    double trackLen = 180.0)
{
    auto group = new TEveElementList("Primary Tracks (truth)");
    if (!pX || pX->empty()) return group;

    for (int i = 0; i < (int)pX->size(); ++i) {
        double px = (*mX)[i], py = (*mY)[i], pz = (*mZ)[i];
        double pmag = std::sqrt(px*px + py*py + pz*pz);
        if (pmag <= 0) continue;

        TVector3 start((*pX)[i], (*pY)[i], (*pZ)[i]);
        TVector3 dir(px/pmag, py/pmag, pz/pmag);
        TVector3 end = start + trackLen * dir;

        int    thisPDG = (*pdg)[i];
        int    tid     = (*trackID)[i];
        double ke      = (*ekin)[i];

        Color_t     col   = Palette::ForPDG(thisPDG);
        std::string pname = Palette::NameForPDG(thisPDG);

        std::ostringstream lbl;
        lbl << "Truth " << pname << " [TID=" << tid << "]"
            << "  Ekin=" << std::fixed << std::setprecision(0) << ke << " MeV";

        auto trkGrp = new TEveElementList(lbl.str().c_str());
        trkGrp->SetMainColor(col);

        auto line = new TEveLine("truth_seg");
        line->SetLineColor(col);
        line->SetLineWidth(2);
        line->SetLineStyle(7);   // long dash  → clearly different from reco
        line->SetNextPoint(start.X(), start.Y(), start.Z());
        line->SetNextPoint(end.X(),   end.Y(),   end.Z());
        trkGrp->AddElement(line);

        auto psVtx = new TEvePointSet("vtx");
        psVtx->SetMarkerStyle(kOpenDiamond);
        psVtx->SetMarkerSize(2.2);
        psVtx->SetMarkerColor(col);
        psVtx->SetNextPoint(start.X(), start.Y(), start.Z());
        trkGrp->AddElement(psVtx);

        group->AddElement(trkGrp);
    }
    return group;
}

// ============================================================
//  True photon trajectories (dotted)
// ============================================================
static TEveElementList* DrawTruePhotons(
    const std::vector<double>* cx, const std::vector<double>* cy, const std::vector<double>* cz,
    const std::vector<double>* ex, const std::vector<double>* ey, const std::vector<double>* ez,
    const std::vector<double>* phPid)
{
    auto group = new TEveElementList("True Photons (truth)");
    if (!cx || cx->empty()) return group;

    Color_t col = TColor::GetColor(1.0f, 1.0f, 0.55f);
    for (size_t i = 0; i < cx->size(); ++i) {
        double photonParentID = (*phPid)[i];
        auto line = new TEveLine(Form("gamma_%zu_%f", i, photonParentID));
        line->SetLineColor(col);
        line->SetLineWidth(2);
        line->SetLineStyle(3);
        line->SetNextPoint((*cx)[i], (*cy)[i], (*cz)[i]);
        line->SetNextPoint((*ex)[i], (*ey)[i], (*ez)[i]);
        group->AddElement(line);
    }
    return group;
}

// ============================================================
//  Background hit cloud (very dim – just shows detector coverage)
// ============================================================
static TEvePointSet* DrawAllHits(
    const std::vector<double>* xs,
    const std::vector<double>* ys,
    const std::vector<double>* zs)
{
    auto ps = new TEvePointSet("All Calo Hits");
    ps->SetMarkerStyle(kFullDotSmall);
    ps->SetMarkerSize(0.28);
    ps->SetMarkerColor(TColor::GetColor(0.24f, 0.24f, 0.26f));
    if (xs)
        for (size_t i = 0; i < xs->size(); ++i)
            ps->SetNextPoint((*xs)[i], (*ys)[i], (*zs)[i]);
    return ps;
}

// ============================================================
//  Primary vertex
// ============================================================
static TEvePointSet* DrawVertex(const TVector3& v)
{
    auto ps = new TEvePointSet("Primary Vertex");
    ps->SetMarkerStyle(kFullStar);
    ps->SetMarkerSize(1);
    ps->SetMarkerColor(kWhite);
    ps->SetNextPoint(v.X(), v.Y(), v.Z());
    return ps;
}

// ============================================================
//  pi0 candidate connectors
// ============================================================
static TEveElementList* DrawPi0Candidates(
    const std::vector<Pi0Candidate>& selected,
    const std::vector<Cluster>&      allClusters)
{
    auto group = new TEveElementList("pi0 Candidates");
    Color_t col = TColor::GetColor(1.0f, 1.0f, 0.55f);

    for (size_t k = 0; k < selected.size(); ++k) {
        const auto& cand = selected[k];
        std::ostringstream lbl;
        lbl << "pi0 #" << k
            << "  m_gg=" << std::fixed << std::setprecision(1) << cand.mgg << " MeV"
            << "  theta=" << std::setprecision(1)
            << cand.theta * TMath::RadToDeg() << " deg";

        auto candGrp = new TEveElementList(lbl.str().c_str());

        auto conn = new TEveLine("connector");
        conn->SetLineColor(col);
        conn->SetLineWidth(2);
        conn->SetLineStyle(5);
        conn->SetNextPoint(cand.c1->centroid.X(), cand.c1->centroid.Y(), cand.c1->centroid.Z());
        conn->SetNextPoint(cand.c2->centroid.X(), cand.c2->centroid.Y(), cand.c2->centroid.Z());
        candGrp->AddElement(conn);

        group->AddElement(candGrp);
    }
    return group;
}

// ============================================================
//  MAIN
// ============================================================
int main(int argc, char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <data.root> <geometry.root>\n";
        return 1;
    }

    TGeoManager::Import(argv[2], "geometry");
    if (!gGeoManager) {
        std::cerr << "ERROR: cannot load geometry\n";
        return 1;
    }

    new TApplication("app", &argc, argv);
    TEveManager::Create();

    // auto projMgrRPhi = new TEveProjectionManager();
    // projMgrRPhi->SetProjection(TEveProjection::kPT_RPhi);
    // gEve->AddToListTree(projMgrRPhi, kTRUE);

    // auto viewRPhi = gEve->SpawnNewViewer("RPhi View");
    // auto sceneRPhi = gEve->SpawnNewScene("RPhi Scene");

    // viewRPhi->AddScene(sceneRPhi);
    // projMgrRPhi->SetCurrentCamera(viewRPhi->GetGLViewer()->GetCurrentCamera());
    auto projMgrRPhi = new TEveProjectionManager();
    projMgrRPhi->SetProjection(TEveProjection::kPT_RPhi);
    gEve->AddToListTree(projMgrRPhi, kFALSE);

    auto projScene = gEve->SpawnNewScene("RPhi Projected Scene");
    // projMgrRPhi->SetProjModel(projScene);

    auto viewRPhi = gEve->SpawnNewViewer("RPhi View");
    viewRPhi->AddScene(projScene);

    // Dark background
    if (auto gl = gEve->GetDefaultViewer()->GetGLViewer())
        gl->SetClearColor(TColor::GetColor(0.04f, 0.04f, 0.06f));

    // All geometry starts dim
    Color_t dimCol = TColor::GetColor(0.28f, 0.30f, 0.33f);
    ResetGeometryColors(dimCol, 90);

    // Add geometry to scene; vis level 5 reaches the SECE cell level
    auto geoTop = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
    geoTop->SetVisLevel(5);
    gEve->AddGlobalElement(geoTop);
    // auto gse  = (TEveGeoShapeExtract*) gGeoManager->Get("Gentle");
    projMgrRPhi->AddElement(geoTop);

    // AddCoordinateAxes();

    // Open data
    std::cout << "datafile name: "<< argv[1] << "geometry filename: "<< argv[2] << std::endl;
    // TFile* f = TFile::Open(argv[1]);
    TFile* f = TFile::Open("../build/processedData/processed_fullSignal_2026-04-29.root");
    if (!f || f->IsZombie()) { std::cerr << "ERROR: cannot open data file\n"; return 1; }
    TTree* t = dynamic_cast<TTree*>(f->Get("digitizedHits"));
    if (!t) { std::cerr << "ERROR: TTree 'digitizedHits' not found\n"; return 1; }

    // Branches
    std::vector<double> *cX=nullptr,*cY=nullptr,*cZ=nullptr,*cE=nullptr;
    SafeSetBranch(t,"centerX",cX); SafeSetBranch(t,"centerY",cY);
    SafeSetBranch(t,"centerZ",cZ); SafeSetBranch(t,"energy", cE);

    std::vector<double> *primX=nullptr,*primY=nullptr,*primZ=nullptr,*primEkin=nullptr;
    std::vector<double> *primPx=nullptr,*primPy=nullptr,*primPz=nullptr;
    std::vector<int>    *primPDG=nullptr,*primTID=nullptr;
    SafeSetBranch(t,"PrimaryPosX",primX);   SafeSetBranch(t,"PrimaryPosY",primY);
    SafeSetBranch(t,"PrimaryPosZ",primZ);   SafeSetBranch(t,"PrimaryEkin",primEkin);
    SafeSetBranch(t,"PrimaryMomX",primPx);  SafeSetBranch(t,"PrimaryMomY",primPy);
    SafeSetBranch(t,"PrimaryMomZ",primPz);  SafeSetBranch(t,"PrimaryPDG", primPDG);
    SafeSetBranch(t,"PrimaryTrackID",primTID);

    std::vector<double> *phCX=nullptr,*phCY=nullptr,*phCZ=nullptr;
    std::vector<double> *phEX=nullptr,*phEY=nullptr,*phEZ=nullptr;
    std::vector<double> *phPid=nullptr;
    SafeSetBranch(t,"TruePhotonCreationX",phCX); SafeSetBranch(t,"TruePhotonCreationY",phCY);
    SafeSetBranch(t,"TruePhotonCreationZ",phCZ); SafeSetBranch(t,"TruePhotonEndX",phEX);
    SafeSetBranch(t,"TruePhotonEndY",phEY);      SafeSetBranch(t,"TruePhotonEndZ",phEZ);
    SafeSetBranch(t, "TruePhotonParentID",phPid);

    std::vector<int>    *tpcPDG=nullptr,*tpcTID=nullptr;
    std::vector<double> *tpcKE=nullptr;
    std::vector<double> *tpcFX=nullptr,*tpcFY=nullptr,*tpcFZ=nullptr;
    std::vector<double> *tpcLX=nullptr,*tpcLY=nullptr,*tpcLZ=nullptr;
    std::vector<double> *tpcDS=nullptr,*tpcDT=nullptr;
    SafeSetBranch(t,"TPC_pdg",tpcPDG);         SafeSetBranch(t,"TPC_trackID",tpcTID);
    SafeSetBranch(t,"TPC_TrueKE",tpcKE);
    SafeSetBranch(t,"TPC_firstPosX",tpcFX);    SafeSetBranch(t,"TPC_firstPosY",tpcFY);
    SafeSetBranch(t,"TPC_firstPosZ",tpcFZ);    SafeSetBranch(t,"TPC_lastPosX",tpcLX);
    SafeSetBranch(t,"TPC_lastPosY",tpcLY);     SafeSetBranch(t,"TPC_lastPosZ",tpcLZ);
    SafeSetBranch(t,"TPC_smearedDedx",tpcDS);  SafeSetBranch(t,"TPC_theoryDedx",tpcDT);

    auto palette = Palette::MakeClusterPalette();
    DEDXTable dedxTable("dedx_tables_Ar80CO2.root");

    const double pi0_h=135,pi0_a=45,pi0_k=1.6,pi0_b=1.55;

    Long64_t nEvents = std::min((Long64_t)2, t->GetEntries());

    for (Long64_t ev = 0; ev < nEvents; ++ev) {

        // Reset geometry highlights before each new event
        ResetGeometryColors(dimCol, 90);

        t->GetEntry(ev);

        auto eventList = new TEveElementList(Form("Event %lld", ev));
        gEve->AddElement(eventList);

        TVector3 vertex(0,0,0);
        if (primX && !primX->empty())
            vertex.SetXYZ((*primX)[0], (*primY)[0], (*primZ)[0]);

        // Draw order: back to front
        eventList->AddElement(DrawAllHits(cX, cY, cZ));
        eventList->AddElement(DrawTruePhotons(phCX,phCY,phCZ,phEX,phEY,phEZ,phPid));
        eventList->AddElement(DrawPrimaryTracks(primX,primY,primZ,
                                                primPx,primPy,primPz,
                                                primEkin,primPDG,primTID));
        eventList->AddElement(DrawTPCTracks(tpcFX,tpcFY,tpcFZ,
                                            tpcLX,tpcLY,tpcLZ,
                                            tpcPDG,tpcKE,tpcDS));
        eventList->AddElement(DrawVertex(vertex));

        // Build hits for reconstruction
        if (!cX || cX->empty()) continue;
        std::vector<Hit> hits;
        hits.reserve(cX->size());
        for (size_t k = 0; k < cX->size(); ++k)
            hits.push_back({(*cX)[k],(*cY)[k],(*cZ)[k],(*cE)[k]});

        std::vector<ChargedTrack> chargedTracks;
        if (tpcFX && tpcPDG) {
            for (size_t k = 0; k < tpcPDG->size(); ++k) {
                ChargedTrack trk;
                trk.id        = (int)k;
                trk.vertex    = vertex;
                trk.exitPoint = TVector3((*tpcLX)[k],(*tpcLY)[k],(*tpcLZ)[k]);
                trk.direction = (trk.exitPoint -
                    TVector3((*tpcFX)[k],(*tpcFY)[k],(*tpcFZ)[k])).Unit();
                trk.TrueKE      = (*tpcKE)[k];
                trk.TruePDG     = (*tpcPDG)[k];
                trk.smearedDedx = tpcDS ? (*tpcDS)[k] : 0.0;
                trk.theoryDedx  = tpcDT ? (*tpcDT)[k] : 0.0;
                trk.resolution  = 0.15;
                chargedTracks.push_back(trk);
            }
        }

        RecoEvent reco = ReconstructEvent(hits, chargedTracks, vertex, dedxTable);

        // pi0 selection
        std::vector<Pi0Candidate> candidates;
        for (size_t a = 0; a < reco.clusters.size(); ++a)
            for (size_t b = a+1; b < reco.clusters.size(); ++b) {
                const auto& g1=reco.clusters[a].p4, &g2=reco.clusters[b].p4;
                double mgg=( g1+g2).M(), theta=openingAngle(g1,g2);
                double r2=std::pow((mgg-pi0_h)/pi0_a,2)+std::pow((theta-pi0_k)/pi0_b,2);
                if (r2>1.0) continue;
                Pi0Candidate c; c.c1=&reco.clusters[a]; c.c2=&reco.clusters[b];
                c.mgg=mgg; c.theta=theta; c.p4=g1+g2; c.score=r2;
                candidates.push_back(c);
            }
        std::sort(candidates.begin(),candidates.end(),
            [](const Pi0Candidate& a,const Pi0Candidate& b){return a.score<b.score;});
        std::vector<Pi0Candidate> selected;
        std::unordered_set<const Cluster*> used;
        for (const auto& c:candidates) {
            if (used.count(c.c1)||used.count(c.c2)) continue;
            used.insert(c.c1); used.insert(c.c2);
            selected.push_back(c);
        }

        // Neutral clusters
        auto neutGrp = new TEveElementList("Neutral Clusters (ECAL)");
        for (size_t ci=0; ci<reco.clusters.size(); ++ci)
            neutGrp->AddElement(
                DrawNeutralCluster(reco.clusters[ci],(int)ci,
                                   palette[ci%palette.size()]));
        eventList->AddElement(neutGrp);

        // Charged clusters
        auto chGrp = new TEveElementList("Charged Clusters");
        for (size_t ci=0; ci<reco.chargedClusters.size(); ++ci)
            chGrp->AddElement(
                DrawChargedCluster(reco.chargedClusters[ci],(int)ci,
                                   palette[(ci+5)%palette.size()]));
        eventList->AddElement(chGrp);

        // pi0 connectors
        eventList->AddElement(DrawPi0Candidates(selected, reco.clusters));
        // projMgrRPhi->ImportElements(eventList);
        projMgrRPhi->ImportElements(eventList, projScene);
    }
    gEve->Redraw3D(kTRUE);
    gApplication->Run(kTRUE);
    return 0;
}

