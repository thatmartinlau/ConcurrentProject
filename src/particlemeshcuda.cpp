#include "core.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <cufft.h>

__global__ void assign_mass_to_grid(Body* bodies, double* grid, int num_bodies, int grid_size, double min_x, double min_y, double cell_size) {
int i = blockIdx.x * blockDim.x + threadIdx.x;
if (i >= num_bodies) return;

double x = bodies[i].pos.x;
double y = bodies[i].pos.y;
int gx = (int)((x - min_x) / cell_size);
int gy = (int)((y - min_y) / cell_size);

if (gx >= 0 && gx < grid_size && gy >= 0 && gy < grid_size) {
atomicAdd(&grid[gy * grid_size + gx], bodies[i].m);
}
}


__global__ void solve_poisson_in_freq_space(cufftDoubleComplex* freq_data, 
    int grid_size, double domain_size_x, double domain_size_y) {
int i = blockIdx.x * blockDim.x + threadIdx.x;
if (i >= grid_size * grid_size) return;

int kx = i % grid_size;
int ky = i / grid_size;

double fx = 2.0 * M_PI * (kx < grid_size / 2 ? kx : kx - grid_size) / domain_size_x;
double fy = 2.0 * M_PI * (ky < grid_size / 2 ? ky : ky - grid_size) / domain_size_y;
double k_sq = fx * fx + fy * fy;

if (k_sq > 0) {
freq_data[i].x /= -k_sq;
freq_data[i].y /= -k_sq;
}
}

__global__ void compute_forces(Body* bodies, double* potential,
    int num_bodies, int grid_size, double min_x, double min_y, double cell_size) {
int i = blockIdx.x * blockDim.x + threadIdx.x;
if (i >= num_bodies) return;

double x = bodies[i].pos.x;
double y = bodies[i].pos.y;
int gx = (int)((x - min_x) / cell_size);
int gy = (int)((y - min_y) / cell_size);

double fx = 0.0, fy = 0.0;
if (gx > 0 && gx < grid_size - 1 && gy > 0 && gy < grid_size - 1) {
fx = (potential[gy * grid_size + (gx + 1)] - potential[gy * grid_size + (gx - 1)]) / (2.0 * cell_size);
fy = (potential[(gy + 1) * grid_size + gx] - potential[(gy - 1) * grid_size + gx]) / (2.0 * cell_size);
}

bodies[i].acc.x = -fx;
bodies[i].acc.y = -fy;
}


__global__ void update_bodies(Body* bodies, int num_bodies, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_bodies) return;

    bodies[i].vel.x += dt * bodies[i].acc.x;
    bodies[i].vel.y += dt * bodies[i].acc.y;
    bodies[i].pos.x += dt * bodies[i].vel.x;
    bodies[i].pos.y += dt * bodies[i].vel.y;
}


__global__ void normalize_potential(double* potential, int length, double norm_factor) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= length) return;
    potential[i] /= norm_factor;
}

void particle_mesh_simulation_CUDA(System &universe, double dt, int grid_size, int num_bodies){

    const size_t THREADS_PER_BLOCK = 256;
    const size_t BLOCKS_NUM_BODIES = (num_bodies + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    const size_t BLOCKS_NUM_GRID = ((grid_size * grid_size) + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;

    Body* d_bodies;
    cudaMalloc(&d_bodies, num_bodies * sizeof(Body));
    cudaMemcpy(d_bodies, h_bodies, num_bodies * sizeof(Body), cudaMemcpyHostToDevice);

    double* d_grid;
    cudaMalloc(&d_grid, grid_size * grid_size * sizeof(double));
    cudaMemset(d_grid, 0, grid_size * grid_size * sizeof(double));

    auto [min_pos, max_pos] = universe.exposeBounds();
    double domain_size_x = max_pos.data[0] - min_pos.data[0];
    double domain_size_y = max_pos.data[1] - min_pos.data[1];
    double cell_size = domain_size_x / grid_size;

    // Assign masses to grid
    assign_mass_to_grid<<<BLOCKS_NUM_BODIES, THREADS_PER_BLOCK>>>(d_bodies, d_grid, num_bodies, grid_size, min_pos.data[0], min_pos.data[1], cell_size);
    cudaDeviceSynchronize();

    // Allocate frequency space data
    cufftDoubleComplex* d_freq_data;
    cudaMalloc(&d_freq_data, grid_size * grid_size * sizeof(cufftDoubleComplex));

    // Create FFT plan for real-to-complex transform
    cufftHandle plan;
    cufftPlan2d(&plan, grid_size, grid_size, CUFFT_D2Z);

    // Execute forward FFT (real to complex)
    cufftExecD2Z(plan, d_grid, d_freq_data);

    // Solve Poisson in frequency space
    solve_poisson_in_freq_space<<<BLOCKS_NUM_GRID, THREADS_PER_BLOCK>>>(d_freq_data, grid_size, domain_size_x, domain_size_y);
    cudaDeviceSynchronize();

    // Execute inverse FFT (complex to real)
    cufftExecZ2D(plan, d_freq_data, d_grid);

    // Normalize potential (divide by grid_size^2)
    // You can do this with a CUDA kernel or cudaMemcpy + CPU loop. Let's do kernel here:
    int total_grid_points = grid_size * grid_size;
    normalize_potential<<<BLOCKS_NUM_GRID, THREADS_PER_BLOCK>>>(d_grid, total_grid_points, grid_size * grid_size);
    cudaDeviceSynchronize();

    // Compute forces on bodies
    compute_forces<<<BLOCKS_NUM_BODIES, THREADS_PER_BLOCK>>>(d_bodies, d_grid, num_bodies, grid_size, min_pos.data[0], min_pos.data[1], cell_size);
    cudaDeviceSynchronize();

    // Update bodies positions and velocities
    update_bodies<<<BLOCKS_NUM_BODIES, THREADS_PER_BLOCK>>>(d_bodies, num_bodies, dt);
    cudaDeviceSynchronize();

    // Copy updated bodies back to host
    cudaMemcpy(h_bodies, d_bodies, num_bodies * sizeof(Body), cudaMemcpyDeviceToHost);

    // Cleanup
    cufftDestroy(plan);
    cudaFree(d_grid);
    cudaFree(d_freq_data);
    cudaFree(d_bodies);
}   