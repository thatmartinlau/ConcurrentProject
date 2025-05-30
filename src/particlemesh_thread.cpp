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
#define DEBUG false
#define THREAD_CHECKING false


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

    size_t total_bodies = universe.bodies.size();

    for (int step = 0; step < STEP_COUNT; ++step) {
       
        for (auto &row : grid){
            std::fill(row.begin(), row.end(), 0.0);
        }

        std::vector<std::thread> threads;
        size_t bodies_per_thread = (universe.bodies.size() + num_threads) / num_threads;
       std::vector<std::mutex> grid_locks(grid_size * grid_size); // flat array of mutexes


        //if(THREAD_CHECKING){
        // Add a vector to track grids accessed by each thread
        std::mutex cout_mutex;  // for synchronized printing
        std::vector<std::vector<std::tuple<size_t, int, int>>> thread_grid_access(num_threads);
        //}
    
        for (size_t t = 0; t < num_threads; ++t) {
            threads.emplace_back([&, t=t]() {
                size_t start_idx = t * bodies_per_thread;
                size_t end_idx = std::min(start_idx + bodies_per_thread, total_bodies);
                
                for (size_t i = start_idx; i < end_idx; ++i) {
                    const auto& body = universe.bodies[i];
                    double x = body.coordinates.data[0];
                    double y = body.coordinates.data[1];

                    int gx = static_cast<int>((x - min_pos.data[0]) / cell_size);
                    int gy = static_cast<int>((y - min_pos.data[1]) / cell_size);

                    if (gx >= 0 && gx < grid_size && gy >= 0 && gy < grid_size) {
                        size_t lock_index = gx * grid_size + gy;
                        
                        if(THREAD_CHECKING){
                            thread_grid_access[t].emplace_back(i, gx, gy);  // Track body index and cell
                        }

                        {   
                            std::lock_guard<std::mutex> lock(grid_locks[lock_index]);
                            grid[gx][gy] += body.m;
                        }
                    }
                }
                
                if(THREAD_CHECKING){
                    {
                        std::lock_guard<std::mutex> lock(cout_mutex);
                        std::cout << "Thread " << t << " accessed grid cells:\n";
                        for (const auto& entry : thread_grid_access[t]) {
                            size_t body_index;
                            int gx, gy;
                            std::tie(body_index, gx, gy) = entry;
                            const auto& b = universe.bodies[body_index];
                            std::cout << "  Body " << body_index 
                                    << " (mass=" << b.m << ", pos=(" << b.coordinates.data[0] 
                                    << "," << b.coordinates.data[1] << ")) -> (" << gx << ", " << gy << ")\n";
                        }
                    }
                }
            });
        }

        for (auto& thread : threads) {
            thread.join();
        }
        threads.clear();


        
        //std::cout<<"managed first parallelization!"<< "\n";
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

        for (size_t t = 0; t < num_threads; ++t) {
            threads.emplace_back([&, t]() {
                size_t start_idx = t * bodies_per_thread;
                size_t end_idx = std::min(start_idx + bodies_per_thread, total_bodies);
                for (size_t i = start_idx; i < end_idx; ++i) {
                    auto& body = universe.bodies[i];
                    // NGP interpolation
                    int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
                    int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);

                    if (grid_x >= 1 && grid_x < grid_size - 1 && grid_y >= 1 && grid_y < grid_size - 1) {
                        double fx = -(potential[grid_x + 1][grid_y] - potential[grid_x - 1][grid_y]) / (2.0 * cell_size);
                        double fy = -(potential[grid_x][grid_y + 1] - potential[grid_x][grid_y - 1]) / (2.0 * cell_size);

                        body.acceleration = Vector(fx, fy);
                        //if(DEBUG){
                        //std::cout << "fx, fy: " << fx << " " << fy << " for Body " << body.title << "\n";
                        //}
                    } else {
                        body.acceleration = Vector(0.0, 0.0);
                    }
                }
            });
        }

        for (auto& thread : threads) thread.join();
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

