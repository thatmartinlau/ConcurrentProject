#include "core.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>


__global__ void 


void particle_mesh_simulation_CUDA(System &universe, double dt, int grid_size){

    const size_t BLOCKS_NUM = 128; 
    const size_t THREADS_PER_BLOCK = 1024;
    const size_t TOTAL_THREADS = BLOCKS_NUM * THREADS_PER_BLOCK; 

    Body* d_bodies;
    cudaMalloc(&d_bodies, num_bodies * sizeof(Body));
    cudaMemcpy(d_bodies, h_bodies, num_bodies * sizeof(Body), cudaMemcpyHostToDevice);

    double* d_grid;
    cudaMalloc(&d_grid, grid_size * grid_size * sizeof(double));
    cudaMemset(d_grid, 0, grid_size * grid_size * sizeof(double));


}