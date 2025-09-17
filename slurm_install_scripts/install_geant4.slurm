#!/bin/bash
#SBATCH -J install_geant
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -A lu2025-2-51
#SBATCH --output=l_%x-%j.out
#SBATCH --error=l_%x-%j.err

module load GCCcore/13.2.0
module load X11/20231019
module load expat/2.5.0
module load Xerces-C++/3.2.5
module load Qt5/5.15.13

wget https://github.com/Geant4/geant4/archive/refs/tags/v11.3.2.tar.gz
tar -xvf v11.3.2.tar.gz

mkdir geant4-v11.3.2-build
cd geant4-v11.3.2-build

cmake -DCMAKE_INSTALL_PREFIX=../geant4-11.3.2-install \
      -DGEANT4_INSTALL_DATA=ON \
      -DGEANT4_USE_GDML=ON \
      -DGEANT4_USE_OPENGL_X11=ON \
      -DGEANT4_USE_QT=ON \
      ../geant4-11.3.2

make -j$(nproc)
make install
