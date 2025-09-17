#!/bin/bash
#SBATCH -J install_root
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -A lu2025-2-51
#SBATCH --output=l_%x-%j.out
#SBATCH --error=l_%x-%j.err

module purge

module load GCCcore/13.2.0
module load X11/20231019
module load Python/3.11.5
module load OpenSSL/1.1

git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git root_src
mkdir root_build root_install
cd root_build
cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src
cmake -D geombuilder=on
cmake --build . -- install -j4 
