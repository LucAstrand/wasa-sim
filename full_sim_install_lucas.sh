#!/bin/bash
set -e

# -------------------- USER CONFIG --------------------
WASA_INSTALL_DIR="$(pwd)/wasa_full"
export WASA_INSTALL_DIR

# -------------------- LOAD REQUIRED MODULES --------------------
module purge
# module load GCC/11.3.0
# module load GCC/12.2.0

# module load GCC/12.3.0
# module load OpenMPI/4.1.5
# module load ROOT/6.30.06

# module load GCCcore/11.3.0
# module load GCCcore/12.2.0
# module load GCCcore/12.3.0
# module load X11/20230603

module load GCCcore/13.2.0
module load X11/20231019
module load Python/3.11.5

module load OpenSSL/1.1
module load expat/2.5.0
module load Xerces-C++/3.2.5
module load Qt5/5.15.13

export CC=$(which gcc)
export CXX=$(which g++)

# -------------------- CLEAN OLD INSTALL AND MAKE NEW DIR --------------------
echo "[INFO] Cleaning previous installation..."
rm -rf "$WASA_INSTALL_DIR"
mkdir -p "$WASA_INSTALL_DIR"
cd "$WASA_INSTALL_DIR"

# -------------------- ROOT BUILD --------------------

git clone --branch latest-stable --depth=1 https://github.com/root-project/root.git root_src
mkdir root_build root_install
cd root_build
cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src
cmake -D geombuilder=on
cmake --build . -- install -j4 

# -------------------- GEANT4 BUILD --------------------
if [ ! -d "$WASA_INSTALL_DIR/geant4-install" ]; then 
    echo "[INFO] Building GEANT4..."
    wget https://github.com/Geant4/geant4/archive/refs/tags/v11.3.2.tar.gz
    tar -xzf v11.3.2.tar.gz
    mkdir geant4-build && cd geant4-build

    cmake -DGEANT4_INSTALL_DATA=ON \
          -DGEANT4_USE_QT=ON \
          -DGEANT4_USE_OPENGL_X11=ON \
          -DGEANT4_USE_GDML=ON \
          -DCMAKE_INSTALL_PREFIX=$WASA_INSTALL_DIR/geant4-install \
          ../geant4-11.3.2

    make -j$(nproc)
    make install
    cd $WASA_INSTALL_DIR
else
    echo "[INFO] GEANT4 already installed, skipping build."
fi

# Source Geant4 environment
source $WASA_INSTALL_DIR/geant4-install/bin/geant4.sh

# -------------------- VGM BUILD --------------------
echo "[INFO] Cloning and building VGM..."
if [ ! -d vgm ]; then 
    git clone https://github.com/vmc-project/vgm.git
fi
cd vgm

# Fix CMake version
if [ -f CMakeLists.txt ]; then
    sed -i 's/cmake_minimum_required(VERSION 2.8)/cmake_minimum_required(VERSION 3.10)/' CMakeLists.txt
    sed -i 's/cmake_minimum_required(VERSION 2.8.12)/cmake_minimum_required(VERSION 3.10)/' cmake/VGMRequiredPackages.cmake
fi

mkdir -p build install
cd build

cmake -DCMAKE_INSTALL_PREFIX=../install \
      -DGeant4_DIR=$Geant4_DIR \
      -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH \
      -DCMAKE_CXX_COMPILER=$CXX \
      -DCMAKE_C_COMPILER=$CC \
      ..

make -j$(nproc)
make install
cd ../..

# Set VGM env
export VGM_DIR=$WASA_INSTALL_DIR/vgm/install
export LD_LIBRARY_PATH=$VGM_DIR/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=$VGM_DIR:$CMAKE_PREFIX_PATH

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

cmake -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
      -DVGM_DIR=$VGM_DIR/lib/cmake/VGM \
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

# -------------------- FINISH --------------------
echo ""
echo "================================================================================"
echo " INSTALLATION COMPLETE"
echo "================================================================================"
echo "Environment set up under: $WASA_INSTALL_DIR"
echo "You are already using the ROOT module (no need to source thisroot.sh)."
echo "Geant4 environment has been sourced for this session."
echo "To run simulations in a new shell, reload modules and re-source geant4.sh:"
echo "   module load ROOT/6.30.06"
echo "   source $WASA_INSTALL_DIR/geant4-install/bin/geant4.sh"
echo ""
echo "Then run:"
echo "   cd $WASA_INSTALL_DIR/hibeam_g4/build"
echo "   ./hibeam_g4 -m run_gps.mac -c example.config output.root"
echo "================================================================================"
