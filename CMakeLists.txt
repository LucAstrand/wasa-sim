cmake_minimum_required(VERSION 3.10)
project(pi0Reco)

# --- ROOT setup ---
find_package(ROOT REQUIRED COMPONENTS RIO Hist Tree MathCore Gpad Graf)

# Include ROOT macros
include(${ROOT_USE_FILE})

# --- Project include path ---
include_directories(${ROOT_INCLUDE_DIRS} include)

# --- Source files ---
file(GLOB SOURCES src/*.cpp)

# --- Create executable ---
add_executable(pi0Analysis ${SOURCES})

# --- Link to ROOT libraries ---
target_link_libraries(pi0Analysis ${ROOT_LIBRARIES})
