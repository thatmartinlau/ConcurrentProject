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


## Express build and run

```bash
# Basic simulation with default parameters
make nbody

# Customize simulation parameters
make nbody NTHREADS=13 VISUALIZE=true PRINT_TELEMETRY=true EXPORT_CSV=true

```

Available Make Parameters

    NTHREADS: Number of threads to use (default: 5)
    VISUALIZE: Generate visualization (default: false)
    PRINT_TELEMETRY: Print simulation data to console (default: false)
    EXPORT_CSV: Export simulation data to CSV (default: false)

Make Targets

    nbody: Run basic N-body simulation
    nbody_particle: Run particle mesh simulation
    bs: Run Barnes-Hut serial implementation
    bsp: Run Barnes-Hut parallel implementation (13 threads)


## Manual Build 

### CPU-only executable (naive & Barnes–Hut)

```bash
# 1. Compile sources
g++ -std=c++17 -O3 -I./src \
    -c src/core.cpp src/simplesimulation.cpp src/barneshutt.cpp main.cpp 

# 2. Link into executable
g++ -O3 \
    core.o simplesimulation.o barneshutt.o main.o \
    -o nbody_cpu $(Magick++-config --cppflags --cxxflags --ldflags --libs)
```

### CPU-only executable (Particle mesh)
```bash
# 1. Compile sources
g++ -std=c++17 -O3 -I./src \
  -c src/core.cpp src/simplesimulation.cpp src/particlemesh_thread.cpp src/particlemesh.cpp mainparticlemesh.cpp

# 2. Link into executable
g++ -O3 \
    core.o simplesimulation.o particlemesh_thread.o particlemesh.o mainparticlemesh.o \
    -o nbodyparticlemesh_cpu $(Magick++-config --cppflags --cxxflags --ldflags --libs)
```


### GPU-enabled executable (brute-force PRAM style)

```bash
nvcc -std=c++14 -O3 -Xcompiler "-std=c++17 -I./src" -DUSE_CUDA \
    src/core.cpp src/simplesimulation.cpp src/barneshutt.cpp main.cpp  \
    -o nbody_gpu $(Magick++-config --cppflags --cxxflags --ldflags --libs) -lcudart
```

## Running

All executables accept a `-method=` flag to select the simulation algorithm:

- `naive`    — simple \(O(N^2)\) CPU simulation
- `barneshut`— Barnes–Hut \(O(N \log N)\) CPU simulation

```bash
# Naive CPU
./nbody_cpu -method=naive

# Barnes–Hut sequential
./nbody_cpu -method=barneshut

# Barnes–Hut parallel
./nbody_cpu -method=barneshut --parallel

#Particle-mesh sequential
./nbodyparticlemesh_cpu -method=particlemesh

#Particle-mesh parallel
./nbodyparticlemesh_cpu -method=particlemesh_thread
```

Each run will produce a GIF (e.g., `naive_simulation.gif`, `barneshut_simulation.gif`, or `gpu_simulation.gif` or `particle_mesh.gif` or `particle_mesh_thread.gif`) in the working directory.
