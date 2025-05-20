#include "src/core.hpp"
#include "src/simplesimulation.hpp"
#include "src/barneshutt.hpp"
#include "src/particlemesh.hpp"

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <thread>
#include <mutex>

#define M_PI 3.14159265358979323846
#define BIG_G 6.67e-11



void particle_mesh_simulation_parallel(System &universe, double dt, int grid_size, size_t num_threads) {
    universe.telemetry.clear();

     // Define the grid
    auto [min_pos, max_pos] = universe.exposeBounds();
    double cell_size = (max_pos.data[0] - min_pos.data[0]) / grid_size;

    // Allocate shared data
    std::vector<std::vector<double>> grid(grid_size, std::vector<double>(grid_size, 0.0));
    std::vector<std::vector<double>> potential(grid_size, std::vector<double>(grid_size, 0.0));

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_plan forward_plan = fftw_plan_dft_2d(grid_size, grid_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_2d(grid_size, grid_size, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    std::mutex grid_mutex;

    for (int step = 0; step < STEP_COUNT; ++step) {
        // Zero the grid
        for (auto &row : grid) {
            std::fill(row.begin(), row.end(), 0.0);
        }

        // === 1. Parallel Mass Assignment ===
        std::vector<std::thread> threads;
        size_t chunk_size = universe.bodies.size() / num_threads;

        for (size_t t = 0; t < num_threads; ++t) {
            size_t start = t * chunk_size;
            size_t end = (t == num_threads - 1) ? universe.bodies.size() : start + chunk_size;

            threads.emplace_back([&, start, end]() {
                std::vector<std::vector<double>> local_grid(grid_size, std::vector<double>(grid_size, 0.0));
                for (size_t i = start; i < end; ++i) {
                    const auto &body = universe.bodies[i];
                    int gx = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                    int gy = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);
                    if (gx >= 0 && gx < grid_size && gy >= 0 && gy < grid_size) {
                        local_grid[gx][gy] += body.m;
                    }
                }

                std::lock_guard<std::mutex> lock(grid_mutex);
                for (int i = 0; i < grid_size; ++i)
                    for (int j = 0; j < grid_size; ++j)
                        grid[i][j] += local_grid[i][j];
            });
        }
        for (auto &t : threads) t.join();
        threads.clear();

        // === 2. FFT and Potential Computation (single-threaded, FFTW is not thread-safe) ===
        for (int i = 0; i < grid_size; ++i)
            for (int j = 0; j < grid_size; ++j) {
                in[i * grid_size + j][0] = grid[i][j];
                in[i * grid_size + j][1] = 0.0;
            }

        fftw_execute(forward_plan);

        for (int i = 0; i < grid_size; ++i) {
            int kx_index = (i <= grid_size / 2) ? i : i - grid_size;
            double kx = 2.0 * M_PI * kx_index / (grid_size * cell_size);
            for (int j = 0; j < grid_size; ++j) {
                int ky_index = (j <= grid_size / 2) ? j : j - grid_size;
                double ky = 2.0 * M_PI * ky_index / (grid_size * cell_size);
                double k2 = kx * kx + ky * ky;
                int idx = i * grid_size + j;
                if (k2 > 0) {
                    out[idx][0] *= -BIG_G / k2;
                    out[idx][1] *= -BIG_G / k2;
                } else {
                    out[idx][0] = out[idx][1] = 0.0;
                }
            }
        }

        fftw_execute(backward_plan);

        for (int i = 0; i < grid_size; ++i)
            for (int j = 0; j < grid_size; ++j)
                potential[i][j] = in[i * grid_size + j][0] / (grid_size * grid_size);

        // === 3. Parallel Force Calculation and Update ===
        threads.clear();
        for (size_t t = 0; t < num_threads; ++t) {
            size_t start = t * chunk_size;
            size_t end = (t == num_threads - 1) ? universe.bodies.size() : start + chunk_size;

            threads.emplace_back([&, start, end]() {
                for (size_t i = start; i < end; ++i) {
                    auto &body = universe.bodies[i];
                    int gx = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                    int gy = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);

                    if (gx >= 1 && gx < grid_size - 1 && gy >= 1 && gy < grid_size - 1) {
                        double fx = -(potential[gx + 1][gy] - potential[gx - 1][gy]) / (2.0 * cell_size);
                        double fy = -(potential[gx][gy + 1] - potential[gx][gy - 1]) / (2.0 * cell_size);
                        body.acceleration.data[0] = fx;
                        body.acceleration.data[1] = fy;
                    } else {
                        body.acceleration.data[0] = 0.0;
                        body.acceleration.data[1] = 0.0;
                    }
                    body.update(dt);
                }
            });
        }
        for (auto &t : threads) t.join();

        // === Save telemetry ===
        std::vector<Vector> positions;
        for (auto &body : universe.bodies)
            positions.push_back(body.coordinates);
        universe.telemetry.push_back(std::move(positions));
    }

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(in);
    fftw_free(out);
}
