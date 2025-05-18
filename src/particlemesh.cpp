
#include "core.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>

void particle_mesh_simulation(System &universe, double dt, int grid_size) {
            universe.telemetry.clear();
        
            // Define the grid
            auto [min_pos, max_pos] = universe.exposeBounds();
            double cell_size = (max_pos.data[0] - min_pos.data[0]) / grid_size;
        
            // Initialize the grid
            std::vector<std::vector<double>> grid(grid_size, std::vector<double>(grid_size, 0.0));
            std::vector<std::vector<double>> potential(grid_size, std::vector<double>(grid_size, 0.0));
        
            // FFTW setup
            fftw_complex *in, *out;
            fftw_plan plan;
            in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
            plan = fftw_plan_dft_2d(grid_size, grid_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        
            for (int step = 0; step < STEP_COUNT; ++step) {
                // Reset grid
                for (auto &row : grid) {
                    std::fill(row.begin(), row.end(), 0.0);
                }
        
                // Assign masses to the grid
                for (const auto &body : universe.bodies) {
                    int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                    int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);
        
                    if (grid_x >= 0 && grid_x < grid_size && grid_y >= 0 && grid_y < grid_size) {
                        grid[grid_x][grid_y] += body.m;
                    }
                }
        
                // Solve Poisson's equation using FFT
                for (int i = 0; i < grid_size; ++i) {
                    for (int j = 0; j < grid_size; ++j) {
                        in[i * grid_size + j][0] = grid[i][j];
                        in[i * grid_size + j][1] = 0.0;
                    }
                }
        
                fftw_execute(plan);
        
                // Calculate potential (simplified)
                for (int i = 0; i < grid_size; ++i) {
                    for (int j = 0; j < grid_size; ++j) {
                        double kx = 2.0 * M_PI * (i < grid_size / 2 ? i : i - grid_size) / (max_pos.data[0] - min_pos.data[0]);
                        double ky = 2.0 * M_PI * (j < grid_size / 2 ? j : j - grid_size) / (max_pos.data[1] - min_pos.data[1]);
                        double k_squared = kx * kx + ky * ky;
                        if (k_squared > 0) {
                            out[i * grid_size + j][0] /= -k_squared;
                            out[i * grid_size + j][1] /= -k_squared;
                        }
                    }
                }
        
                fftw_execute(fftw_plan_dft_2d(grid_size, grid_size, out, in, FFTW_BACKWARD, FFTW_ESTIMATE));
        
                for (int i = 0; i < grid_size; ++i) {
                    for (int j = 0; j < grid_size; ++j) {
                        potential[i][j] = in[i * grid_size + j][0] / (grid_size * grid_size);
                    }
                }
        
                // Calculate forces from potential
                for (auto &body : universe.bodies) {
                    int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                    int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);
        
                    if (grid_x >= 0 && grid_x < grid_size && grid_y >= 0 && grid_y < grid_size) {
                        // Calculate gradient of potential (force)
                        double fx = 0.0, fy = 0.0;
                        if (grid_x > 0 && grid_x < grid_size - 1) {
                            fx = (potential[grid_x + 1][grid_y] - potential[grid_x - 1][grid_y]) / (2.0 * cell_size);
                        }
                        if (grid_y > 0 && grid_y < grid_size - 1) {
                            fy = (potential[grid_x][grid_y + 1] - potential[grid_x][grid_y - 1]) / (2.0 * cell_size);
                        }
                        body.acceleration.data[0] = -fx;
                        body.acceleration.data[1] = -fy;
                    }
                }
        
                // Update positions
                std::vector<Vector> positions;
                for (auto &body : universe.bodies) {
                    body.update(dt);
                    positions.push_back(body.coordinates);
                }
                universe.telemetry.push_back(positions);
            }
        
            fftw_destroy_plan(plan);
            fftw_free(in);
            fftw_free(out);
        }
        
