#include "EventDisplay.hpp"

TEveElementList* DrawCluster(const Cluster& c, int id, Color_t color)
{
    auto group = new TEveElementList(Form("Cluster_%d", id));
    group->SetMainColor(color);

    // ---------------- Centroid point ----------------
    TEvePointSet* ps = new TEvePointSet();
    ps->SetNextPoint(c.centroid.X(), c.centroid.Y(), c.centroid.Z());
    ps->SetMarkerStyle(20);
    ps->SetMarkerSize(1.5);
    ps->SetMarkerColor(color);
    group->AddElement(ps);

    // ---------------- Cluster box ----------------
    double L = 3.0;  // half-size in cm
    TVector3 cc = c.centroid;

    TEveBox* box = new TEveBox();
    box->SetMainColor(color);
    box->SetMainTransparency(60);  // semi-transparent

    box->SetVertex(0, cc.X()-L, cc.Y()-L, cc.Z()-L);
    box->SetVertex(1, cc.X()+L, cc.Y()-L, cc.Z()-L);
    box->SetVertex(2, cc.X()+L, cc.Y()+L, cc.Z()-L);
    box->SetVertex(3, cc.X()-L, cc.Y()+L, cc.Z()-L);
    box->SetVertex(4, cc.X()-L, cc.Y()-L, cc.Z()+L);
    box->SetVertex(5, cc.X()+L, cc.Y()-L, cc.Z()+L);
    box->SetVertex(6, cc.X()+L, cc.Y()+L, cc.Z()+L);
    box->SetVertex(7, cc.X()-L, cc.Y()+L, cc.Z()+L);

    group->AddElement(box);

    // ---------------- Momentum arrow ----------------
    TVector3 dir = c.p4.Vect().Unit();
    double arrowLen = 4.0 + log(c.p4.E());

    // Convert (direction × length) → vector
    TVector3 vec = arrowLen * dir;

    // TEveArrow(vec, origin)
    TEveArrow* arrow = new TEveArrow(
        vec.X(), vec.Y(), vec.Z(),
        cc.X(), cc.Y(), cc.Z()
    );

    // Color as needed
    arrow->SetMainColor(color);
    // arrow->SetTubeRadius(0.3);

    return group;
}

TEveElementList* BuildClusterBoundingBox(const Cluster& c,
                                 const std::vector<double>& centerX,
                                 const std::vector<double>& centerY,
                                 const std::vector<double>& centerZ,
                                 Color_t col)
{
    auto group = new TEveElementList("Cluster");
    group->SetMainColor(col);
    TVector3 cc = c.centroid;


    // ---------------- Centroid point ----------------
    TEvePointSet* ps = new TEvePointSet();
    ps->SetNextPoint(c.centroid.X(), c.centroid.Y(), c.centroid.Z());
    ps->SetMarkerStyle(20);
    ps->SetMarkerSize(1.5);
    ps->SetMarkerColor(col);
    group->AddElement(ps);

    // // ---------------- Cluster box ----------------
    // double xmin=1e9, xmax=-1e9;
    // double ymin=1e9, ymax=-1e9;
    // double zmin=1e9, zmax=-1e9;

    // for (int idx : c.hitIndices) {
    //     double x = centerX[idx];
    //     double y = centerY[idx];
    //     double z = centerZ[idx];

    //     xmin = std::min(xmin, x);
    //     xmax = std::max(xmax, x);
    //     ymin = std::min(ymin, y);
    //     ymax = std::max(ymax, y);
    //     zmin = std::min(zmin, z);
    //     zmax = std::max(zmax, z);
    // }

    // TEveBox* box = new TEveBox();
    // box->SetMainColor(col);
    // box->SetMainTransparency(60);

    // box->SetVertex(0, xmin, ymin, zmin);
    // box->SetVertex(1, xmax, ymin, zmin);
    // box->SetVertex(2, xmax, ymax, zmin);
    // box->SetVertex(3, xmin, ymax, zmin);
    // box->SetVertex(4, xmin, ymin, zmax);
    // box->SetVertex(5, xmax, ymin, zmax);
    // box->SetVertex(6, xmax, ymax, zmax);
    // box->SetVertex(7, xmin, ymax, zmax);

    // group->AddElement(box);
    
    // ---------------- Momentum arrow ----------------
    TVector3 dir = c.p4.Vect().Unit();
    double arrowLen = 4.0 + log(c.p4.E());
    
    // Convert (direction × length) → vector
    TVector3 vec = arrowLen * dir;
    
    // TEveArrow(vec, origin)
    TEveArrow* arrow = new TEveArrow(
        vec.X(), vec.Y(), vec.Z(),
        cc.X(), cc.Y(), cc.Z()
    );
    
    // Color as needed
    arrow->SetMainColor(col);
    group->AddElement(arrow);
    
    // return box;
    return group;
    
}

int main(int argc, char **argv) {
    if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " <data.root> <geometry.root> \n";
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

//    // camera
//    auto s = gEve->SpawnNewScene("Projected Event");
//    gEve->GetDefaultViewer()->AddScene(s);
//    auto v = gEve->GetDefaultGLViewer();
//    v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
//    TGLOrthoCamera &cam = (TGLOrthoCamera &)v->CurrentCamera();
//    cam.SetZoomMinMax(0.2, 20);
 
//    // projections
//    auto mng = new TEveProjectionManager(TEveProjection::kPT_RPhi);
//    s->AddElement(mng);
//    auto axes = new TEveProjectionAxes(mng);
//    axes->SetTitle("R-PhiProjection");
//    s->AddElement(axes);
//    gEve->AddToListTree(axes, kTRUE);
//    gEve->AddToListTree(mng, kTRUE);

    // // =====================================================
    // // 1. Create two scenes: one for 3D, one for projections
    // // =====================================================
    // TEveScene* scene3D  = gEve->SpawnNewScene("3D Scene");
    // TEveScene* scene2D  = gEve->SpawnNewScene("Projection Scene");

    // // =====================================================
    // // 2. Create two viewers (two separate windows/tabs)
    // // =====================================================
    // TEveViewer* viewer3D = gEve->SpawnNewViewer("3D View");
    // TEveViewer* viewer2D = gEve->SpawnNewViewer("Projected View");

    // // Attach scenes to viewers
    // viewer3D->AddScene(scene3D);
    // viewer2D->AddScene(scene2D);

    // // =====================================================
    // // 3. Configure cameras
    // // =====================================================

    // // --- 3D viewer uses standard perspective camera ---
    // TGLViewer* glv3D = viewer3D->GetGLViewer();
    // glv3D->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);

    // // --- 2D viewer uses R–Phi orthographic camera ---
    // TGLViewer* glv2D = viewer2D->GetGLViewer();
    // glv2D->SetCurrentCamera(TGLViewer::kCameraOrthoXOY);
    // TGLOrthoCamera& cam2D = (TGLOrthoCamera&) glv2D->CurrentCamera();
    // cam2D.SetZoomMinMax(0.2, 20);

    // // =====================================================
    // // 4. Create projection manager & attach only to 2D scene
    // // =====================================================
    // TEveProjectionManager* projMgr =
    //     new TEveProjectionManager(TEveProjection::kPT_RPhi);
    // scene2D->AddElement(projMgr);

    // TEveProjectionAxes* projAxes = new TEveProjectionAxes(projMgr);
    // projAxes->SetTitle("R-φ Projection");
    // scene2D->AddElement(projAxes);

    // gEve->AddToListTree(projMgr, kTRUE);
    // gEve->AddToListTree(projAxes, kTRUE);

    

    TGeoVolume* topVolume = gGeoManager->GetTopVolume();
    int nd = topVolume->GetNdaughters();

    for (int i = 0; i < nd; i++) {
        TGeoNode* node = topVolume->GetNode(i);
        TGeoVolume* vol = node->GetVolume();

        vol->SetTransparency(80);     // 0-100
        vol->SetLineColor(kGray + 1); // optional: lighter color
    }


    auto top = new TEveGeoTopNode(gGeoManager, gGeoManager->GetTopNode());
    top->SetVisLevel(4);
    gEve->AddGlobalElement(top);
    // scene3D->AddElement(top);

    TFile* f = TFile::Open(dataFile);
    if (!f) { std::cerr << "ERROR: cannot open datafile\n"; return 1; }

    // TTree* t = (TTree*)f->Get("hibeam");
    TTree* t = (TTree*)f->Get("digitizedHits");
    if (!t) { std::cerr << "ERROR: TTree not found\n"; return 1; }

    std::vector<double>* PrimaryEkin = nullptr;
    std::vector<double>* PrimaryTime = nullptr;
    std::vector<double>* PrimaryPosX = nullptr;
    std::vector<double>* PrimaryPosY = nullptr;
    std::vector<double>* PrimaryPosZ = nullptr;
    std::vector<double>* PrimaryMomX = nullptr;
    std::vector<double>* PrimaryMomY = nullptr;
    std::vector<double>* PrimaryMomZ = nullptr;
    std::vector<double>* centerX = nullptr;
    std::vector<double>* centerZ = nullptr;
    std::vector<double>* centerY = nullptr;
    std::vector<double>* energies = nullptr;
    std::vector<double>* TruePhotonCreationX = nullptr;
    std::vector<double>* TruePhotonCreationY = nullptr;
    std::vector<double>* TruePhotonCreationZ = nullptr;
    std::vector<double>* TruePhotonEndX = nullptr;
    std::vector<double>* TruePhotonEndY = nullptr;
    std::vector<double>* TruePhotonEndZ = nullptr;

    t->SetBranchAddress("PrimaryEkin",    &PrimaryEkin);
    t->SetBranchAddress("PrimaryTime",    &PrimaryTime);
    t->SetBranchAddress("PrimaryPosX",    &PrimaryPosX);
    t->SetBranchAddress("PrimaryPosY",    &PrimaryPosY);
    t->SetBranchAddress("PrimaryPosZ",    &PrimaryPosZ);
    t->SetBranchAddress("PrimaryMomX",    &PrimaryMomX);
    t->SetBranchAddress("PrimaryMomY",    &PrimaryMomY);
    t->SetBranchAddress("PrimaryMomZ",    &PrimaryMomZ);

    t->SetBranchAddress("centerX",    &centerX);
    t->SetBranchAddress("centerY",    &centerY);
    t->SetBranchAddress("centerZ",    &centerZ);
    t->SetBranchAddress("energy", &energies);


    t->SetBranchAddress("TruePhotonCreationX",    &TruePhotonCreationX);
    t->SetBranchAddress("TruePhotonCreationY",    &TruePhotonCreationY);
    t->SetBranchAddress("TruePhotonCreationZ",    &TruePhotonCreationZ);
    t->SetBranchAddress("TruePhotonEndX",    &TruePhotonEndX);
    t->SetBranchAddress("TruePhotonEndY",    &TruePhotonEndY);
    t->SetBranchAddress("TruePhotonEndZ",    &TruePhotonEndZ);

    // --------------------------------------------------------
    // 4. Detector subsystems 
    // --------------------------------------------------------
    std::vector<std::string> dets = {
        "SECE0",
        "SECE1",
        "SECE2",
        "SECE3",
        "SECE4",
        "SECE5",
        "SECE6",
        "SECE7",
        "SECE8",
        "SECE9",
        "SECE10",
        "SECE11",
        "SECE12",
        "SECE13",
        "SECE14",
        "SECE15",
        "SECE16",
    };

    // Hit branches
    std::map<std::string, std::vector<double>*> HitX, HitY, HitZ, HitE;

    for (auto& D : dets) {
        HitX[D] = HitY[D] = HitZ[D] = HitE[D] = nullptr;

        t->SetBranchAddress((D+"_Position_X").c_str(), &HitX[D]);
        t->SetBranchAddress((D+"_Position_Y").c_str(), &HitY[D]);
        t->SetBranchAddress((D+"_Position_Z").c_str(), &HitZ[D]);
        t->SetBranchAddress((D+"_EDep").c_str(),        &HitE[D]);
    }

    Long64_t nEvents = t->GetEntries();

    for (Long64_t ev = 0; ev < nEvents; ev++) {
        t->GetEntry(ev);

        // ============ 1. Unique color per event ============
        // gStyle->SetPalette(kAurora);
        // int idx = ev % gStyle->GetNumberOfColors();
        // Color_t eventColor = gStyle->GetColorPalette(idx);
        // std::cout << eventColor << std::endl;

        // Pick a hue we can vary per event:
        double hue = fmod(ev * 0.61803398875, 1.0);

        // Convert HSV → RGB manually (this version works)
        float r, g, b;
        {
            float h = hue * 6;
            float c = 1.0f;
            float x = (1 - fabs(fmod(h, 2) - 1)) * c;

            if      (0 <= h && h < 1) { r = c; g = x; b = 0; }
            else if (1 <= h && h < 2) { r = x; g = c; b = 0; }
            else if (2 <= h && h < 3) { r = 0; g = c; b = x; }
            else if (3 <= h && h < 4) { r = 0; g = x; b = c; }
            else if (4 <= h && h < 5) { r = x; g = 0; b = c; }
            else                      { r = c; g = 0; b = x; }
        }

        // Normalize
        float R = r;
        float G = g;
        float B = b;

        // Allocate a real color index in ROOT + OpenGL
        Color_t colIndex = TColor::GetFreeColorIndex();
        new TColor(colIndex, R, G, B);   // <-- THIS MAKES IT REAL

        // Use this color
        Color_t eventColor = colIndex;

        // ============ 2. Group for this event ===============
        TEveElementList* eventList =
            new TEveElementList(Form("Event_%lld", ev));
        eventList->SetMainColor(eventColor);
        gEve->AddElement(eventList);
        // scene3D->AddElement(eventList);


        // ============ 3. Draw photons =======================
        for (size_t i = 0; i < TruePhotonCreationX->size(); i++) {
            auto photon = new TEveLine(Form("Photon_%zu", i));
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

        // ============ 4. Build clusters ======================
        TVector3 vertex((*PrimaryPosX)[0],(*PrimaryPosY)[0],(*PrimaryPosZ)[0]);
        std::vector<Hit> hits;
        size_t nHits = energies->size();
        // std::cout << "Event " << ievt << ": " << nHits << " total hits across all rings" << std::endl;
        for (size_t k=0; k<nHits; ++k) {
            hits.push_back({(*centerX)[k], (*centerY)[k], (*centerZ)[k], (*energies)[k]});
        }
        // std::vector<Cluster> clusters = runClustering(centerX, centerY, centerZ);
        std::vector<Cluster> clusters;
        double dEta = 0.10/2;
        double dPhi = 0.10/2;
        double E_seed = 15.00;
        double E_neighbor = 0.03;
        int winSize = 3;

        double totalE_Evt = 0;
        clusters = SlidingWindowClusterHits(hits, vertex, dEta, dPhi, E_seed, E_neighbor, winSize);

        // ============ 5. Draw clusters =======================
        // for (size_t ci = 0; ci < clusters.size(); ci++) {
        //     TEveElementList* cl = DrawCluster(clusters[ci], ci, eventColor);
        //     eventList->AddElement(cl);
        // }

        for (Cluster c : clusters) {
            TEveElementList* clusterElements = BuildClusterBoundingBox(
                c, *centerX, *centerY, *centerZ, eventColor
            );
            eventList->AddElement(clusterElements);
        }   

        for (size_t i = 0; i < centerX->size(); i++) {
        if (!centerX) continue;

            auto ps = new TEvePointSet(Form("Hit_%lld_%zu", ev, i));
            ps->SetMarkerStyle(4);
            ps->SetMarkerSize(1.0);
            ps->SetMarkerColor(eventColor);

            ps->SetNextPoint(centerX->at(i), centerY->at(i), centerZ->at(i));

            eventList->AddElement(ps);
        }

        // ============ 6. Update TEve ========================
        // projMgr->ImportElements(eventList, scene2D);
    }
    gEve->Redraw3D(kTRUE);

    gApplication->Run(kTRUE);
    return 0;
}
