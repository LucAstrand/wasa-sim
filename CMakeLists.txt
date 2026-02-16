cmake_minimum_required(VERSION 3.10)
project(pi0Reco LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fno-omit-frame-pointer")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic -no-pie")



# ROOT must be sourced before running cmake
find_package(ROOT REQUIRED COMPONENTS
    Core
    RIO
    Hist
    Tree
    MathCore
    Graf
    Gpad
    Gui
    Physics
)

message(STATUS "ROOTSYS = $ENV{ROOTSYS}")

# ------------------------------------------------------------------
# Main analysis executable
# ------------------------------------------------------------------
add_executable(Analysis
    src/Analysis.cpp
    src/ParticleID.cpp
    src/Clustering.cpp
    src/PhotonMatch.cpp
    src/Pi0Acceptance.cpp
    src/Pi0Efficiency.cpp
    src/PlotUtils.cpp
    src/TruePhotonCalc.cpp
    src/PIDEfficiency.cpp
    src/RecoEvent.cpp 
    src/Calibration.cpp 
)

target_include_directories(Analysis PRIVATE include)

add_subdirectory(mcpl)

target_link_libraries(Analysis
    mcpl
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Tree
    ROOT::MathCore
    ROOT::Graf
    ROOT::Gpad
    ROOT::Gui
    ROOT::Physics
)

# ------------------------------------------------------------------
# EventDisplay executable
# ------------------------------------------------------------------
add_executable(EventDisplay
    src/EventDisplay.cpp
    src/Clustering.cpp
    src/ParticleID.cpp
    src/RecoEvent.cpp
)

target_include_directories(EventDisplay PRIVATE include)

target_link_libraries(EventDisplay
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Tree
    ROOT::MathCore
    ROOT::Physics

    ROOT::Geom
    ROOT::Graf
    ROOT::Graf3d
    ROOT::Gpad
    ROOT::Gui
    ROOT::Eve
)
