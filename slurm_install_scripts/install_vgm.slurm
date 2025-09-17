#!/bin/bash
#SBATCH -J install_vgm
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -A lu2025-2-51
#SBATCH --output=l_%x-%j.out
#SBATCH --error=l_%x-%j.err

# Load required modules
module load GCCcore/13.2.0
module load Xerces-C++/3.2.5
module load expat/2.5.0
module load Qt5/5.15.13

# Source ROOT and Geant4 environments
source /home/lucast/git/wasa-sim/root_install/bin/thisroot.sh
source /home/lucast/git/wasa-sim/geant4-11.3.2-install/bin/geant4.sh

# Clone VGM
git clone https://github.com/vmc-project/vgm.git

# Create build directory
mkdir vgm_build
cd vgm_build

# Configure with cmake
cmake -DCMAKE_INSTALL_PREFIX=/home/lucast/git/wasa-sim/vgm_install \
      -DCMAKE_PREFIX_PATH="$Qt5_DIR;$CMAKE_PREFIX_PATH" \
      ../vgm

# Build and install
make -j$(nproc)
make install

echo "[INFO] VGM installation completed!"