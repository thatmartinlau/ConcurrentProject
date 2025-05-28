#!/bin/bash

g++ -std=c++17 mainparticlemesh.cpp \
  src/particlemesh.cpp src/core.cpp src/simplesimulation.cpp src/particlemesh_thread.cpp  \
  -I./src \
  -I/opt/homebrew/Cellar/imagemagick/7.1.1-47/include/ImageMagick-7 \
  -I/opt/homebrew/Cellar/fftw/3.3.10_2/include \
  -L/opt/homebrew/Cellar/imagemagick/7.1.1-47/lib \
  -L/opt/homebrew/Cellar/fftw/3.3.10_2/lib \
  -lMagick++-7.Q16HDRI -lMagickCore-7.Q16HDRI -lfftw3 \
  -DMAGICKCORE_HDRI_ENABLE=1 -DMAGICKCORE_QUANTUM_DEPTH=16 \
  -o nbody_sim_particle_mesh
