cmake_minimum_required(VERSION 3.10)
project(HoppedProjectileSim LANGUAGES CXX)

# Add source files
set(SOURCES
    HoppedProjectileSim/dllmain.cpp
    HoppedProjectileSim/Simulation.cpp
)

# Create static library
add_library(HoppedProjectileSim STATIC ${SOURCES})

# Include header files
target_include_directories(HoppedProjectileSim PUBLIC ${PROJECT_SOURCE_DIR}/HoppedProjectileSim)
