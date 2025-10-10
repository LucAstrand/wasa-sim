# Program Title
**File:** `slidingWindowClustering.cpp`  
**Author:** Lucas Åstrand 
**Date:** 2025-10-XX  

---

## Overview

This program is dedicated to the reconstruction of $\pi^0$ particles simulated by the HibeamG4 simulation framework. 

The program accepts root file input, including energy depositions and calorimeter cell centers. A Sliding window - style clustering algorithm is employed, clustering hits in $(\eta - \phi)$ space.

The final invariant mass plot of the reconstructed $\pi^0$ is returned. 

---

## High-Level Description

The following list breaks down the overall functionality of the program. 

1. Parse input file, extracting hits and their energy deposition
2. A Clustering algorithm is called, event-wise, on the collected hits.
3. Cluster characteristics are built and later stored for invariant mass computation using built in ROOT funcitonality.
4. A histogram is populated, with each event's reconstructed invariant mass. 

---

## Program Structure

**Files involved:**
- `slidingWindowClstering.cpp`

**Main Components:**
- `main()` - coordinates execution.
- `Hit`, `Custer` - structs defining needed datastructures.
- `SlidingWindowClusterHits()` - performs core clustering algorithm.
- `SplitClusters()` - performs cluster splitting if formed cluster surpass an energy threshold.

---

## Detailed Function Descriptions

### `int main()`

**Purpose:**  
Entry point of the program; sets up data, calls core routines, and handles output.

**Inputs:**  
1. `root-data-file`

**Outputs:**  
Exit code `0` on success. Saves invariant mass plot. 

### `std::vector<Cluster> SlidingWindowClusterHits()`

**Purpose:**
Main clustering algorithm in $(\eta, \phi)$ space, inspired by calorimeter clustering in high-energy physics.

Detailed Steps:

1. Build Towers: Convert each hit to $(\eta, \phi)$ and group energies into discrete towers of width $(\Delta\eta, \Delta\phi)$.

2. Find Seeds: Identify towers with energy > E_seed and verify that they are local maxima.

3. Sort Seeds: Order by decreasing tower energy.

4. Form Clusters: For each seed, define a sliding window of size winSize in both $\eta$ and $\phi$. Sum all unassigned towers with energy > E_neighbor inside the window.

5. Compute Cluster Properties: Energy-weighted centroid in $(x, y, z)$. Momentum direction from vertex to cluster. Four-vector $(p_x, p_y, p_z, E)$ for physics analysis.

**Outputs:**
A list of reconstructed clusters with energies and hit indices.

### `std::vector<Cluster> SplitClusters()`

**Purpose:** 
Splits clusters with energy above a configurable threshold (split_E_thresh) into sub-clusters, mimicking electromagnetic shower deconvolution.

Algorithm Overview:

1. Energy Check: If cluster energy < threshold, keep as is.

2. Tower Map Rebuild: Re-project cluster hits into $(\eta, \phi)$ tower grid.

3. Find Sub-Seeds: Identify local energy maxima above E_subseed.

4. Initialize Sub-Clusters: Create new clusters centered on sub-seeds.

5. Iterative Reassignment: Reassign hits fractionally to sub-clusters using an exponential distance weighting $f_i = e^{-r_i/r_0}$ until centroids stabilize or maximum iterations reached.

6. Finalize: Keep sub-clusters with energy > 50 MeV; if no valid split, keep original cluster.