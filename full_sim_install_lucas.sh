#!/bin/bash

# -------------------- USER CONFIG --------------------
WASA_INSTALL_DIR="$(pwd)/wasa_full"
export WASA_INSTALL_DIR
# echo "export WASA_INSTALL_DIR=$WASA_INSTALL_DIR" >> ~/.bashrc

ENV_NAME="hibeam_env"

module load Anaconda3/2024.02-1

#Commented to avoid deleting stuff when it is not needed! For a fresh install un-comment.
# # -------------------- CLEAN OLD ENV --------------------
# echo "Cleaning previous installation..."
# conda clean --all -y
# conda env remove -n $ENV_NAME -y 2>/dev/null
# rm -rf ~/.conda/envs/$ENV_NAME 2>/dev/null
# rm -rf ~/.conda/pkgs/* 2>/dev/null
# rm -rf "$WASA_INSTALL_DIR" 2>/dev/null

# # -------------------- CREATE ENV --------------------
# echo "Creating conda environment '$ENV_NAME'..."
# conda create -n $ENV_NAME -c conda-forge root cmake make gcc gxx git wget curl python=3.10 xerces-c expat -y

# Activate conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

export CC=$CONDA_PREFIX/bin/gcc
export CXX=$CONDA_PREFIX/bin/g++

# -------------------- INSTALL PATH SETUP --------------------
mkdir -p "$WASA_INSTALL_DIR"
cd "$WASA_INSTALL_DIR"

# -------------------- LOAD REQUIRED MODULES --------------------
module purge
module load GCC/11.3.0
module load GCC/12.2.0
module load intel-compilers/2022.2.1
module load CMake/3.24.3
module load Boost/1.81.0
module load Xerces-C++/3.2.4
module load OpenMPI/4.1.4
module load Qt5/5.15.7

# -------------------- GEANT4 BUILD --------------------
if [ ! -d "$WASA_INSTALL_DIR/geant4-install" ]; then 
    echo "Building GEANT4"
    wget https://gitlab.cern.ch/geant4/geant4/-/archive/v11.2.2/geant4-v11.2.2.tar.gz
    tar -xzf geant4-v11.2.2.tar.gz
    mkdir geant4-build && cd geant4-build

    cmake -DGEANT4_INSTALL_DATA=ON \
        -DGEANT4_USE_QT=ON \
        -DGEANT4_USE_OPENGL_X11=ON \
        -DGEANT4_USE_XM=OFF \
        -DGEANT4_USE_GDML=ON \
        -DCMAKE_INSTALL_PREFIX=$WASA_INSTALL_DIR/geant4-install \
        ../geant4-v11.2.2

    make -j$(nproc)
    make install
    cd $WASA_INSTALL_DIR

else
        echo "GEANT4 already installed, skipping build"
fi

# -------------------- SET GEANT4 ENV --------------------
GEANT4_CONFIG_DIR=$(find $WASA_INSTALL_DIR/geant4-install -name Geant4Config.cmake -printf "%h")
if [ -z "$GEANT4_CONFIG_DIR" ]; then
    echo "Could not find Geant4Config.cmake!"
    exit 1
fi

export Geant4_DIR=$GEANT4_CONFIG_DIR
export LD_LIBRARY_PATH=$WASA_INSTALL_DIR/geant4-install/lib:$WASA_INSTALL_DIR/geant4-install/lib64:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=$WASA_INSTALL_DIR/geant4-install:$CMAKE_PREFIX_PATH

# Save to bashrc
# {
#     echo "export Geant4_DIR=$Geant4_DIR"
#     echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
#     echo "export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH"
# } >> ~/.bashrc

# -------------------- VGM --------------------
echo "Cloning and building VGM..."
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

# Set VGM env
export VGM_DIR=$WASA_INSTALL_DIR/vgm/install
export LD_LIBRARY_PATH=$VGM_DIR/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=$VGM_DIR:$CMAKE_PREFIX_PATH

# Save VGM to bashrc
# {
#     echo "export VGM_DIR=$VGM_DIR"
#     echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
#     echo "export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH"
# } >> ~/.bashrc

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
    ./geobuilder --barrel --save_as=wasa_geometry.root || echo "Geobuilder completed with warnings"
else
    echo "Geobuilder executable not found!"
fi

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
echo "Copying geometry file..."
# Try multiple possible locations for geometry file
GEOMETRY_FILES=(
    "$WASA_INSTALL_DIR/hibeam_g4_geobuilder/build/wasa_geometry.root"
    "$WASA_INSTALL_DIR/hibeam_g4_geobuilder/build/geometry.root"
    "$WASA_INSTALL_DIR/hibeam_g4_geobuilder/wasa_geometry.root"
)

for geometry_file in "${GEOMETRY_FILES[@]}"; do
    if [ -f "$geometry_file" ]; then
        cp "$geometry_file" "$WASA_INSTALL_DIR/hibeam_g4/build/"
        echo "Copied geometry file: $(basename $geometry_file)"
        break
    fi
done

if [ ! -f "$WASA_INSTALL_DIR/hibeam_g4/build/wasa_geometry.root" ]; then
    echo "Geometry file not found in expected locations"
    echo "Searching for any geometry files:"
    find "$WASA_INSTALL_DIR" -name "*geometry*.root" -type f
fi

# -------------------- PATCH CONFIG FILES --------------------
echo "Patching config files..."
for f in $(find "$WASA_INSTALL_DIR/hibeam_g4" -name "*.config" 2>/dev/null); do
    echo "Patching config file: $f"
    
    # Create a temporary file for clean modifications
    temp_file="${f}.tmp"
    
    # Process the file line by line with proper modifications
    while IFS= read -r line; do
        # Comment out the old Detectors line and add new one
        if [[ "$line" == "Detectors TPC,TARGET,SECE0,SECE1,SECE2,SECE3,SECE4,SECE5,SECE6,SECE7,SECE8,SECE9,SECE10,SECE11,SECE12,SECE13,SECE14,SECE15,SECE16" ]]; then
            echo "##$line" >> "$temp_file"
            echo "Detectors SECE0,SECE1,SECE2,SECE3,SECE4,SECE5,SECE6,SECE7,SECE8,SECE9,SECE10,SECE11,SECE12,SECE13,SECE14,SECE15,SECE16" >> "$temp_file"
        # Fix the MCPL comments order - Source MCPL should come first
        elif [[ "$line" == "# MCPL_Inputfile mcpl_file.mcpl" ]]; then
            echo "# Source MCPL" >> "$temp_file"
            echo "$line" >> "$temp_file"
        # Fix Geometry_Namefile to use just the filename without path
        elif [[ "$line" == "Geometry_Namefile geometry_file.root" ]]; then
            echo "Geometry_Namefile wasa_geometry.root" >> "$temp_file"
        # Handle any other lines
        else
            echo "$line" >> "$temp_file"
        fi
    done < "$f"
    
    # Add WriteTree and WriteHistograms if not present (with proper newlines)
    if ! grep -q "WriteTree" "$temp_file"; then
        echo "" >> "$temp_file"
        echo "WriteTree 1" >> "$temp_file"
    fi
    
    if ! grep -q "WriteHistograms" "$temp_file"; then
        echo "WriteHistograms 0" >> "$temp_file"
    fi
    
    # Replace the original file with the patched version
    mv "$temp_file" "$f"
    
    # Remove GeometryPath line if it was incorrectly added previously
    sed -i '/GeometryPath/d' "$f"
    
    echo "Successfully patched: $f"
done

# -------------------- FINISH --------------------
echo ""
echo "   Installation completed!"
echo ""
echo "   Environment variables set:"
echo "   WASA_INSTALL_DIR: $WASA_INSTALL_DIR"
echo "   Geant4_DIR: $Geant4_DIR"
echo "   VGM_DIR: $VGM_DIR"
echo "   LD_LIBRARY_PATH includes Geant4 and VGM libs"
echo ""
echo "   To use this environment, run:"
echo "   module load Anaconda3/2024.02-1"
echo "   conda activate $ENV_NAME"
echo "   source ~/.bashrc"
echo ""

echo ""
echo "  If any step failed, check the specific component:"
echo "   - VGM: $WASA_INSTALL_DIR/vgm/build/CMakeCache.txt"
echo "   - Geobuilder: $WASA_INSTALL_DIR/hibeam_g4_geobuilder/build/CMakeCache.txt"
echo "   - Main sim: $WASA_INSTALL_DIR/hibeam_g4/build/CMakeCache.txt"
echo ""
echo ""
echo "================================================================================"
echo "  TO RUN A SIMULATION:"
echo "================================================================================"
#echo "   cd $WASA_INSTALL_DIR/hibeam_g4/build"
echo "   ./hibeam_g4 -m run_gps.mac -c example.config output.root"
echo "================================================================================"
echo ""








