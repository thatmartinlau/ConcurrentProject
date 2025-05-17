#include "particlemesh.hpp"
#include <cmath>
#include <algorithm>

ParticleMeshSystem::ParticleMeshSystem() {
    density.resize(GRID_SIZE, std::vector<double>(GRID_SIZE, 0.0));
    potential_freq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * GRID_SIZE * (GRID_SIZE / 2 + 1));
    density_freq   = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * GRID_SIZE * (GRID_SIZE / 2 + 1));
    potential      = (double*) fftw_malloc(sizeof(double) * GRID_SIZE * GRID_SIZE);

    forward_plan  = fftw_plan_dft_r2c_2d(GRID_SIZE, GRID_SIZE, &density[0][0], density_freq, FFTW_MEASURE);
    backward_plan = fftw_plan_dft_c2r_2d(GRID_SIZE, GRID_SIZE, potential_freq, potential, FFTW_MEASURE);
}

ParticleMeshSystem::~ParticleMeshSystem() {
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(potential_freq);
    fftw_free(density_freq);
    fftw_free(potential);
}

void ParticleMeshSystem::simulate(double dt) {
    computeDensity();
    solvePotential();
    computeForcesAndUpdate(dt);
    recordPositions();
}

void ParticleMeshSystem::computeDensity() {
    for (auto& row : density)
        std::fill(row.begin(), row.end(), 0.0);

    for (const auto& body : bodies) {
        int x_idx = std::min(std::max(0, int((body.coordinates.data[0] + 1) * 0.5 * GRID_SIZE)), GRID_SIZE - 1);
        int y_idx = std::min(std::max(0, int((body.coordinates.data[1] + 1) * 0.5 * GRID_SIZE)), GRID_SIZE - 1);
        density[x_idx][y_idx] += body.m;
    }

    fftw_execute(forward_plan);
}

void ParticleMeshSystem::solvePotential() {
    for (int i = 0; i < GRID_SIZE; ++i) {
        for (int j = 0; j <= GRID_SIZE / 2; ++j) {
            int k = i * (GRID_SIZE / 2 + 1) + j;
            double kx = (i <= GRID_SIZE / 2) ? i : i - GRID_SIZE;
            double ky = j;
            double denom = kx * kx + ky * ky + SOFTENING;

            if (denom != 0) {
                potential_freq[k][0] = -density_freq[k][0] / denom * (4 * M_PI * G);
                potential_freq[k][1] = -density_freq[k][1] / denom * (4 * M_PI * G);
            } else {
                potential_freq[k][0] = 0;
                potential_freq[k][1] = 0;
            }
        }
    }

    fftw_execute(backward_plan);

    for (int i = 0; i < GRID_SIZE * GRID_SIZE; ++i)
        potential[i] /= (GRID_SIZE * GRID_SIZE);
}

void ParticleMeshSystem::computeForcesAndUpdate(double dt) {
    for (auto& body : bodies) {
        int x_idx = std::min(std::max(1, int((body.coordinates.data[0] + 1) * 0.5 * GRID_SIZE)), GRID_SIZE - 2);
        int y_idx = std::min(std::max(1, int((body.coordinates.data[1] + 1) * 0.5 * GRID_SIZE)), GRID_SIZE - 2);

        int idx = x_idx + y_idx * GRID_SIZE;

        double dx = (potential[idx + 1] - potential[idx - 1]) * 0.5;
        double dy = (potential[idx + GRID_SIZE] - potential[idx - GRID_SIZE]) * 0.5;

        body.acceleration = Vector(-dx, -dy);
        body.update(dt);
    }
}
