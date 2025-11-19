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

#include "Clustering.hpp"

TEveElementList* DrawCluster(const Cluster& c, int id, Color_t color);

TEveElementList* BuildClusterBoundingBox(const Cluster& c,
                                 const std::vector<double>& centerX,
                                 const std::vector<double>& centerY,
                                 const std::vector<double>& centerZ,
                                 Color_t col);


#endif
