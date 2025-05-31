#!/bin/bash

# Define directories
LLVM_DIR="/opt/homebrew/opt/llvm"
OMP_DIR="/opt/homebrew/opt/libomp"
IM_DIR="/opt/homebrew/Cellar/imagemagick/7.1.1-47"
FFTW_DIR="/opt/homebrew/Cellar/fftw/3.3.10_2"

g++ -std=c++17 mainparticlemesh.cpp \
  src/particlemesh.cpp src/core.cpp src/simplesimulation.cpp src/particlemesh_thread.cpp \
  -I./src \
  -I"$IM_DIR/include/ImageMagick-7" \
  -I"$FFTW_DIR/include" \
  -L"$LLVM_DIR/lib" \
  -L"$OMP_DIR/lib" \
  -L"$IM_DIR/lib" \
  -L"$FFTW_DIR/lib" \
  -lMagick++-7.Q16HDRI -lMagickCore-7.Q16HDRI \
  -lfftw3 -lfftw3_threads -lpthread \
  -DMAGICKCORE_HDRI_ENABLE=1 -DMAGICKCORE_QUANTUM_DEPTH=16 \
  -o nbody_sim_particle_mesh

