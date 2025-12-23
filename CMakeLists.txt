# cmake_minimum_required(VERSION 3.10)
# project(pi0Reco)

# # --- ROOT setup ---
# find_package(ROOT REQUIRED COMPONENTS RIO Hist Tree MathCore Gpad Graf)

# # Include ROOT macros
# include(${ROOT_USE_FILE})

# # --- Project include path ---
# include_directories(${ROOT_INCLUDE_DIRS} include)

# # --- Source files ---
# file(GLOB SOURCES src/*.cpp)

# # --- Create executable ---
# add_executable(pi0Analysis ${SOURCES})

# # --- Link to ROOT libraries ---
# target_link_libraries(pi0Analysis ${ROOT_LIBRARIES})


# cmake_minimum_required(VERSION 3.10)
# project(pi0Reco LANGUAGES CXX)

# # --- ROOT setup with EVE ---
# find_package(ROOT REQUIRED COMPONENTS 
#     RIO Hist Tree MathCore 
#     Gpad Graf Geom
#     Eve EG RGL          # ← These three are essential for TEve*
# )

# include(${ROOT_USE_FILE})

# # --- Include directories ---
# include_directories(${ROOT_INCLUDE_DIRS} include)

# # ------------------------------------------------------------------
# # Main analysis executable
# # ------------------------------------------------------------------
# file(GLOB MAIN_SOURCES src/*.cpp)
# list(REMOVE_ITEM MAIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/EventDisplay.cpp")

# add_executable(pi0Analysis ${MAIN_SOURCES})
# target_link_libraries(pi0Analysis ${ROOT_LIBRARIES})

# # ------------------------------------------------------------------
# # Standalone EventDisplay executable
# # ------------------------------------------------------------------
# add_executable(EventDisplay src/EventDisplay.cpp src/Clustering.cpp)
# target_include_directories(EventDisplay PRIVATE ${ROOT_INCLUDE_DIRS} include)
# target_link_libraries(EventDisplay ${ROOT_LIBRARIES})   # ← Now includes Eve, EG, RGL

# # Optional: same C++ standard
# target_compile_features(EventDisplay PRIVATE cxx_std_17)


# cmake_minimum_required(VERSION 3.10)
# project(pi0Reco LANGUAGES CXX)

# # --- ROOT setup with EVE and GUI ---
# find_package(ROOT REQUIRED COMPONENTS 
#     RIO Hist Tree MathCore 
#     Gpad Graf Geom Gui      # <-- Gui added for TGMainFrame, TGButton, etc.
#     Eve EG RGL              # TEve components
# )

# include(${ROOT_USE_FILE})

# # --- Include directories ---
# include_directories(${ROOT_INCLUDE_DIRS} include)

# # ------------------------------------------------------------------
# # Main analysis executable
# # ------------------------------------------------------------------
# file(GLOB MAIN_SOURCES src/*.cpp)
# list(REMOVE_ITEM MAIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/EventDisplay.cpp")

# add_executable(pi0Analysis ${MAIN_SOURCES})
# target_link_libraries(pi0Analysis ${ROOT_LIBRARIES})

# # ------------------------------------------------------------------
# # Standalone EventDisplay executable
# # ------------------------------------------------------------------
# add_executable(EventDisplay src/EventDisplay.cpp src/Clustering.cpp)
# target_include_directories(EventDisplay PRIVATE ${ROOT_INCLUDE_DIRS} include)

# # Link ROOT libraries, explicitly include Gui
# target_link_libraries(EventDisplay 
#     ${ROOT_LIBRARIES} 
#     ROOT::Gui               # required for TGMainFrame, TGButton, TGPictureButton
# )

# # Optional: set C++ standard
# target_compile_features(EventDisplay PRIVATE cxx_std_17)


cmake_minimum_required(VERSION 3.10)
project(pi0Reco LANGUAGES CXX)

# --- ROOT setup ---
find_package(ROOT REQUIRED COMPONENTS 
    RIO Hist Tree MathCore 
    Gpad Graf Geom Gui
    Eve EG RGL
)

include(${ROOT_USE_FILE})

# --- Include directories ---
include_directories(${ROOT_INCLUDE_DIRS} include)

# ------------------------------------------------------------------
# Main analysis executable - explicitly list sources
# ------------------------------------------------------------------
add_executable(pi0Analysis
    src/ParticleID.cpp
    src/pi0Analysis.cpp          # adjust names as needed
    src/Clustering.cpp
    src/PhotonMatch.cpp
    src/Pi0Acceptance.cpp
    src/Pi0Efficiency.cpp
    src/PlotUtils.cpp
    src/TruePhotonCalc.cpp
    # add any other .cpp files that belong to the main analysis here
)

target_link_libraries(pi0Analysis ${ROOT_LIBRARIES})

# Optional: set C++ standard for main exe too
target_compile_features(pi0Analysis PRIVATE cxx_std_17)

# ------------------------------------------------------------------
# Standalone EventDisplay executable
# ------------------------------------------------------------------
add_executable(EventDisplay 
    src/EventDisplay.cpp 
    src/Clustering.cpp
    src/ParticleID.cpp
)

target_include_directories(EventDisplay PRIVATE ${ROOT_INCLUDE_DIRS} include)

target_link_libraries(EventDisplay 
    ${ROOT_LIBRARIES} 
    ROOT::Gui
)

target_compile_features(EventDisplay PRIVATE cxx_std_17)