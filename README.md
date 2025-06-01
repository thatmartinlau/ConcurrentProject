# ConcurrentProject

N-body simulation using various methods.

**Project members**: Martin LAU, Ziyue QIU (William), Oscar PEYRON.

## Project Prerequisites

- **C++17** compiler (e.g., g++ 7+)
- **ImageMagick++** (for image creation)
- **ffmpeg** (for video creation)
- **(Optional)** NVIDIA CUDA Toolkit (for GPU mode)
- **fftw3** Fast Fourier transform
- **Make** (optional, for convenience)


## Building 

### CPU-only executable (naive & Barnes–Hut)

```bash
# 1. Compile sources
g++ -std=c++17 -O3 -I./src \
    -c src/core.cpp src/simplesimulation.cpp src/barneshutt.cpp main.cpp mainparticlemesh.cpp

# 2. Link into executable
g++ -O3 \
    core.o simplesimulation.o barneshutt.o main.o \
    -o nbody_cpu $(Magick++-config --cppflags --cxxflags --ldflags --libs)
```

### GPU-enabled executable (brute-force PRAM style)

```bash
nvcc -std=c++14 -O3 -Xcompiler "-std=c++17 -I./src" -DUSE_CUDA \
    src/core.cpp src/simplesimulation.cpp src/barneshutt.cpp main.cpp mainparticlemesh.cpp \
    -o nbody_gpu $(Magick++-config --cppflags --cxxflags --ldflags --libs) -lcudart
```

## Running

All executables accept a `-method=` flag to select the simulation algorithm:

- `naive`    — simple \(O(N^2)\) CPU simulation
- `barneshut`— Barnes–Hut \(O(N \log N)\) CPU simulation
- `barneshut_gpu`      — Barnes–Hut brute‑force \(O(N^2)\) GPU simulation (if built with `USE_CUDA`)

```bash
# Naive CPU
./nbody_cpu -method=naive

# Barnes–Hut CPU
./nbody_cpu -method=barneshut

# Barnes–Hut GPU brute‑force
./nbody_gpu -method=barneshut_gpu
```

Each run will produce a GIF (e.g., `naive_simulation.gif`, `barneshut_simulation.gif`, or `gpu_simulation.gif`) in the working directory.
