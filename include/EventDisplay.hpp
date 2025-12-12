#ifndef EVENTDISPLAY_H
#define EVENTDISPLAY_H

#include <iostream>
#include <vector>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGeoManager.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEveLine.h"
#include "TEvePointSet.h"
#include "TEveArrow.h"
#include "TEveBox.h"
#include "TApplication.h"
#include "TColor.h"
#include "TStyle.h"
// #include "TEveProjectionAxes.h"
#include "TEveProjectionManager.h"
#include "TEveProjections.h"

#include <TApplication.h>
#include <TFile.h>
#include <TTree.h>
#include <TGeoManager.h>
#include <TColor.h>
#include <TStyle.h>

#include <TEveManager.h>
#include <TEveBrowser.h>
// #include <TEveGeoTopNode.h>
#include <TEveArrow.h>
// #include <TEveElementList.h>
#include <TEveLine.h>
#include <TEvePointSet.h>

#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>
#include <TEveViewer.h>
#include <TEveScene.h>
#include <TGLViewer.h>

#include <TGClient.h>
#include <TGFrame.h>
#include <TGButton.h>
#include <TGLabel.h>
// #include <TGPictureButton.h>
#include <TGPicture.h>
#include <TGText.h>
// #include <TGTextButton.h>
#include <TSystem.h>
#include <TRootBrowser.h>
#include <TROOT.h>

#include <TGTab.h>


#include "Clustering.hpp"

TEveElementList* DrawCluster(const Cluster& c, int id, Color_t color);

TEveElementList* BuildClusterBoundingBox(const Cluster& c,
                                 const std::vector<double>& centerX,
                                 const std::vector<double>& centerY,
                                 const std::vector<double>& centerZ,
                                 Color_t col);


#endif
