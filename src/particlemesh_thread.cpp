#include "core.hpp"


#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fftw3.h>
#include <thread>
#include <mutex>

#define M_PI 3.14159265358979323846
#define BIG_G 6.67e-11
#define DEBUG true

void particle_mesh_simulation_parallel(System &universe, double dt, int grid_size, size_t num_threads, double R) {
    universe.telemetry.clear();

    std::vector<Vector> initial_positions;
    for (const auto& body : universe.bodies)
        initial_positions.push_back(body.coordinates);
    universe.telemetry.push_back(initial_positions);

    Vector min_pos(-R, -R);
    double cell_size = 2 * R / grid_size;
    double chunk_width = 2 * R / num_threads;

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_plan forward_plan = fftw_plan_dft_2d(grid_size, grid_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_2d(grid_size, grid_size, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    std::vector<std::vector<double>> potential(grid_size, std::vector<double>(grid_size, 0.0));
    std::vector<std::vector<double>> grid(grid_size, std::vector<double>(grid_size, 0.0));

    for (int step = 0; step < STEP_COUNT; ++step) {
        for (auto &row : grid)
            std::fill(row.begin(), row.end(), 0.0);

        // === 1. Parallel Mass Assignment (chunk-based) ===
        std::vector<std::thread> threads;

        for (size_t chunk_x = 0; chunk_x < num_threads; ++chunk_x) {
            for (size_t chunk_y = 0; chunk_y < num_threads; ++chunk_y) {
                threads.emplace_back([&, chunk_x, chunk_y]() {
                    double x0 = -R + chunk_x * chunk_width;
                    double x1 = x0 + chunk_width;
                    double y0 = -R + chunk_y * chunk_width;
                    double y1 = y0 + chunk_width;

                    std::vector<std::vector<double>> local_grid(grid_size, std::vector<double>(grid_size, 0.0));

                    for (const auto& body : universe.bodies) {
                        double x = body.coordinates.data[0];
                        double y = body.coordinates.data[1];

                        if (x >= x0 && x < x1 && y >= y0 && y < y1) {
                            int gx = static_cast<int>((x - min_pos.data[0]) / cell_size);
                            int gy = static_cast<int>((y - min_pos.data[1]) / cell_size);

                            if (gx >= 0 && gx < grid_size && gy >= 0 && gy < grid_size)
                                local_grid[gx][gy] += body.m;
                        }
                    }

                    // Merge local grid into global grid
                    for (int i = 0; i < grid_size; ++i) {
                        for (int j = 0; j < grid_size; ++j) {
                            if (local_grid[i][j] != 0.0) {
                                // Lock-free because each thread writes to disjoint regions
                                grid[i][j] += local_grid[i][j];
                            }
                        }
                    }
                });
            }
        }

        for (auto &t : threads) t.join();
        threads.clear();

        std::cout<<"managed first parallelization!"<< "\n";
        // === 2. FFT and Potential Computation (single-threaded, FFTW is not thread-safe) ===
        for (int i = 0; i < grid_size; ++i){
            for (int j = 0; j < grid_size; ++j) {
                in[i * grid_size + j][0] = grid[i][j];
                in[i * grid_size + j][1] = 0.0;
            }
        }

        fftw_execute(forward_plan);

        for (int i = 0; i < grid_size; ++i) {
            int kx_index = (i <= grid_size / 2) ? i : i - grid_size;
            double kx = 2.0 * M_PI * kx_index / (grid_size * cell_size);

            for (int j = 0; j < grid_size; ++j) {
                int ky_index = (j <= grid_size / 2) ? j : j - grid_size;
                double ky = 2.0 * M_PI * ky_index / (grid_size * cell_size);

                double k_squared = kx * kx + ky * ky;
                int idx = i * grid_size + j;
                if (k_squared > 0) {
                    out[idx][0] *= -BIG_G / k_squared;
                    out[idx][1] *= -BIG_G / k_squared;
                } else {
                    out[idx][0] = 0.0;
                    out[idx][1] = 0.0;
                }
            }
        }

        fftw_execute(backward_plan);

        for (int i = 0; i < grid_size; ++i){
            for (int j = 0; j < grid_size; ++j){
                potential[i][j] = in[i * grid_size + j][0] / (grid_size * grid_size);
        }  
        } 

        for (size_t chunk_x = 0; chunk_x < num_threads; ++chunk_x) {
            for (size_t chunk_y = 0; chunk_y < num_threads; ++chunk_y) {
                threads.emplace_back([&, chunk_x, chunk_y]() {
                double x0 = -R + chunk_x * chunk_width;
                double x1 = x0 + chunk_width;
                double y0 = -R + chunk_y * chunk_width;
                double y1 = y0 + chunk_width;

               
                for (auto& body : universe.bodies) {

                    // Check if body is within this chunk
                    if (body.coordinates.data[0] < x0 || body.coordinates.data[0] >= x1 ||
                        body.coordinates.data[1] < y0 || body.coordinates.data[1] >= y1)
                        continue;

                    // NGP interpolation
                    int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                    int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);

                    if (grid_x >= 1 && grid_x < grid_size - 1 && grid_y >= 1 && grid_y < grid_size - 1) {
                        double fx = -(potential[grid_x + 1][grid_y] - potential[grid_x - 1][grid_y]) / (2.0 * cell_size);
                        double fy = -(potential[grid_x][grid_y + 1] - potential[grid_x][grid_y - 1]) / (2.0 * cell_size);

                        body.acceleration = Vector(fx, fy);
                        std::cout << "fx, fy: " << fx << " " << fy << " for Body " << body.title << "\n";
                    } else {
                        body.acceleration = Vector(0.0, 0.0);
                    }
                }
                
                });
            }
        }
        for (auto &thread : threads) thread.join();
        threads.clear();

        
        // === Save telemetry ===
        std::vector<Vector> positions;
        for (auto &body : universe.bodies){
             body.update(dt);
            positions.push_back(body.coordinates);
        }
        universe.telemetry.push_back(positions);

        if (DEBUG) {
            std::cout << "cell_size: " << cell_size << "\nStep: " << step << "\n";
            for (const auto &body : universe.bodies) {
                std::cout << "Body " << body.title << " Pos: (" << body.coordinates.data[0] << ", " << body.coordinates.data[1]
                          << ") Vel: (" << body.velocity.data[0] << ", " << body.velocity.data[1]
                          << ") Acc: (" << body.acceleration.data[0] << ", " << body.acceleration.data[1] << ")\n";
            }
        }



    }
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(in);
    fftw_free(out);


}
// parallelize over the number of bodies 
/*

void particle_mesh_simulation_parallel(System &universe, double dt, int grid_size) {
    universe.telemetry.clear();

    std::vector<Vector> initial_positions;
    for (const auto &body : universe.bodies) {
        initial_positions.push_back(body.coordinates);
    }
    universe.telemetry.push_back(initial_positions);

    Vector min_pos(-1.197985e+39, 0.0);
    Vector max_pos(4.843318e+38, 7.459491e+24);

    double cell_size = (max_pos.data[0] - min_pos.data[0]) / grid_size;

    std::vector<std::vector<double>> grid(grid_size, std::vector<double>(grid_size, 0.0));
    std::vector<std::vector<double>> potential(grid_size, std::vector<double>(grid_size, 0.0));

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_plan forward_plan = fftw_plan_dft_2d(grid_size, grid_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_2d(grid_size, grid_size, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int step = 0; step < STEP_COUNT; ++step) {
        // Zero grid
        for (auto &row : grid) {
            std::fill(row.begin(), row.end(), 0.0);
        }

        std::mutex grid_mutex;
        std::vector<std::thread> threads;

        // === 1. Parallel Mass Assignment ===
        for (const auto &body : universe.bodies) {
            threads.emplace_back([&grid, &body, &min_pos, &cell_size, &grid_size, step, &grid_mutex]() {
                int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);

                if (grid_x >= 0 && grid_x < grid_size && grid_y >= 0 && grid_y < grid_size) {
                    std::lock_guard<std::mutex> lock(grid_mutex);
                    grid[grid_x][grid_y] += body.m;
                } else if (DEBUG) {
                    std::lock_guard<std::mutex> lock(grid_mutex); // Ensure safe debug output
                    std::cerr << "Body " << body.title << " is out of bounds at step " << step << "\n";
                }
            });
        }
        for (auto &t : threads) t.join();
        threads.clear();

        // === 2. FFT to compute gravitational potential ===
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                in[i * grid_size + j][0] = grid[i][j];  // Real part
                in[i * grid_size + j][1] = 0.0;         // Imaginary part
            }
        }

        fftw_execute(forward_plan);

        for (int i = 0; i < grid_size; ++i) {
            int kx_index = (i <= grid_size / 2) ? i : i - grid_size;
            double kx = 2.0 * M_PI * kx_index / (grid_size * cell_size);

            for (int j = 0; j < grid_size; ++j) {
                int ky_index = (j <= grid_size / 2) ? j : j - grid_size;
                double ky = 2.0 * M_PI * ky_index / (grid_size * cell_size);

                double k_squared = kx * kx + ky * ky;
                int idx = i * grid_size + j;
                if (k_squared > 0.0) {
                    out[idx][0] *= -BIG_G / k_squared;
                    out[idx][1] *= -BIG_G / k_squared;
                } else {
                    out[idx][0] = 0.0;
                    out[idx][1] = 0.0;
                }
            }
        }

        fftw_execute(backward_plan);

        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                potential[i][j] = in[i * grid_size + j][0] / (grid_size * grid_size);
            }
        }

        // === 3. Parallel Acceleration Calculation ===
        for (auto &body : universe.bodies) {
            threads.emplace_back([&body, &potential, &min_pos, &cell_size, &grid_size]() {
                int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);

                if (grid_x >= 1 && grid_x < grid_size - 1 && grid_y >= 1 && grid_y < grid_size - 1) {
                    double fx = -(potential[grid_x + 1][grid_y] - potential[grid_x - 1][grid_y]) / (2.0 * cell_size);
                    double fy = -(potential[grid_x][grid_y + 1] - potential[grid_x][grid_y - 1]) / (2.0 * cell_size);
                    body.acceleration = Vector(fx, fy);
                } else {
                    body.acceleration = Vector(0.0, 0.0);
                }
            });
        }
        for (auto &t : threads) t.join();
        threads.clear();

        std::vector<Vector> positions(universe.bodies.size());

for (size_t i = 0; i < universe.bodies.size(); ++i) {
    threads.emplace_back([&, i]() {
        universe.bodies[i].update(dt);
        positions[i] = universe.bodies[i].coordinates;
    });
}

for (auto &t : threads) {
    t.join();
}
        universe.telemetry.push_back(positions);

        // === Optional Debug Info ===
        if (DEBUG) {
            std::cout << "Step: " << step << " | cell_size: " << cell_size << "\n";
            for (const auto &body : universe.bodies) {
                std::cout << "Body " << body.title << " Pos: (" << body.coordinates.data[0] << ", " << body.coordinates.data[1]
                          << ") Vel: (" << body.velocity.data[0] << ", " << body.velocity.data[1]
                          << ") Acc: (" << body.acceleration.data[0] << ", " << body.acceleration.data[1] << ")\n";
            }
        }
    }

    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
    fftw_free(in);
    fftw_free(out);
}
*/