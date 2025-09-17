#!/bin/bash
#SBATCH -J install_hibeam_g4
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH -A lu2025-2-51
#SBATCH --output=l_%x-%j.out
#SBATCH --error=l_%x-%j.err

# Load modules
module load GCCcore/13.2.0
module load expat/2.5.0
module load Xerces-C++/3.2.5
module load Qt5/5.15.13

# Source ROOT and Geant4 environments
source /home/lucast/git/wasa-sim/root_install/bin/thisroot.sh
source /home/lucast/git/wasa-sim/geant4-11.3.2-install/bin/geant4.sh

# Set environment variables
export WASA_INSTALL_DIR="$(pwd)/wasa_install"
export VGM_DIR="/home/lucast/git/wasa-sim/vgm_install"

# Explicitly set Geant4_DIR if not set by geant4.sh
if [ -z "$Geant4_DIR" ]; then
    export Geant4_DIR="/home/lucast/git/wasa-sim/geant4-11.3.2-install/lib64/cmake/Geant4"
fi

# Create installation directory
mkdir -p "$WASA_INSTALL_DIR"

# -------------------- hibeam_g4_geobuilder --------------------
cd "$WASA_INSTALL_DIR"
if [ ! -d hibeam_g4_geobuilder ]; then
    git clone https://git.esss.dk/nnbar-sim/hibeam_g4_geobuilder.git
fi
cd hibeam_g4_geobuilder
# Fix CMake version
if [ -f CMakeLists.txt ]; then
    sed -i 's/^cmake_minimum_required.*$/cmake_minimum_required(VERSION 3.10)/' CMakeLists.txt
fi
mkdir -p build && cd build
cmake -DVGM_DIR=$VGM_DIR/lib/cmake/VGM \
      -DGeant4_DIR=$Geant4_DIR \
      -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH \
      ..
make -j$(nproc)
# Run geobuilder
if [ -f "./geobuilder" ]; then
    ./geobuilder --barrel --save_as=wasa_geometry.root || echo "[WARN] Geobuilder completed with warnings"
else
    echo "[ERROR] Geobuilder executable not found!"
fi
cd ../..

# -------------------- hibeam_g4 --------------------
cd "$WASA_INSTALL_DIR"
if [ ! -d hibeam_g4 ]; then  
    git clone https://git.esss.dk/nnbar-sim/hibeam_g4.git
fi
cd hibeam_g4
# Fix CMake version
if [ -f CMakeLists.txt ]; then
    sed -i 's/cmake_minimum_required(VERSION 3\.0)/cmake_minimum_required(VERSION 3.10)/' CMakeLists.txt
    sed -i 's/cmake_minimum_required(VERSION 2\.[0-9])/cmake_minimum_required(VERSION 3.10)/' CMakeLists.txt
fi
mkdir -p build && cd build
cmake -DVGM_DIR=$VGM_DIR/lib/cmake/VGM \
      -DGeant4_DIR=$Geant4_DIR \
      -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH \
      ..
make -j$(nproc)

# -------------------- COPY GEOMETRY FILE --------------------
echo "[INFO] Copying geometry file..."
GEOMETRY_FILES=(
    "$WASA_INSTALL_DIR/hibeam_g4_geobuilder/build/wasa_geometry.root"
    "$WASA_INSTALL_DIR/hibeam_g4_geobuilder/build/geometry.root"
    "$WASA_INSTALL_DIR/hibeam_g4_geobuilder/wasa_geometry.root"
)
for geometry_file in "${GEOMETRY_FILES[@]}"; do
    if [ -f "$geometry_file" ]; then
        cp "$geometry_file" "$WASA_INSTALL_DIR/hibeam_g4/build/"
        echo "[INFO] Copied geometry file: $(basename $geometry_file)"
        break
    fi
done
if [ ! -f "$WASA_INSTALL_DIR/hibeam_g4/build/wasa_geometry.root" ]; then
    echo "[WARN] Geometry file not found in expected locations."
    find "$WASA_INSTALL_DIR" -name "*geometry*.root" -type f
fi

# -------------------- PATCH CONFIG FILES --------------------
echo "[INFO] Patching config files..."
for f in $(find "$WASA_INSTALL_DIR/hibeam_g4" -name "*.config" 2>/dev/null); do
    temp_file="${f}.tmp"
    while IFS= read -r line; do
        if [[ "$line" == "Detectors TPC,TARGET,SECE0,SECE1,SECE2,SECE3,SECE4,SECE5,SECE6,SECE7,SECE8,SECE9,SECE10,SECE11,SECE12,SECE13,SECE14,SECE15,SECE16" ]]; then
            echo "##$line" >> "$temp_file"
            echo "Detectors SECE0,SECE1,SECE2,SECE3,SECE4,SECE5,SECE6,SECE7,SECE8,SECE9,SECE10,SECE11,SECE12,SECE13,SECE14,SECE15,SECE16" >> "$temp_file"
        elif [[ "$line" == "# MCPL_Inputfile mcpl_file.mcpl" ]]; then
            echo "# Source MCPL" >> "$temp_file"
            echo "$line" >> "$temp_file"
        elif [[ "$line" == "Geometry_Namefile geometry_file.root" ]]; then
            echo "Geometry_Namefile wasa_geometry.root" >> "$temp_file"
        else
            echo "$line" >> "$temp_file"
        fi
    done < "$f"
    if ! grep -q "WriteTree" "$temp_file"; then
        echo "" >> "$temp_file"
        echo "WriteTree 1" >> "$temp_file"
    fi
    if ! grep -q "WriteHistograms" "$temp_file"; then
        echo "WriteHistograms 0" >> "$temp_file"
    fi
    mv "$temp_file" "$f"
    sed -i '/GeometryPath/d' "$f"
done

echo "[INFO] HIBEAM_G4 installation and configuration completed!"