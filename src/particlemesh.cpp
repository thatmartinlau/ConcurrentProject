#include "core.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cmath> 
#include <fftw3.h>

#define M_PI 3.14159265358979323846
#define BIG_G 6.67e-11
#define DEBUG true

void particle_mesh_simulation(System &universe, double dt, int grid_size, double R) {
    universe.telemetry.clear();

    std::vector<Vector> initial_positions;
    for (const auto& body : universe.bodies) {
        initial_positions.push_back(body.coordinates);
    }
    universe.telemetry.push_back(initial_positions);

    //double R = 6.5e12;  // Slightly beyond Pluto
    Vector min_pos(-R, -R);
    Vector max_pos(R, R);
    //Vector min_pos(-1.197985e+39,0);  // Using the smallest x and y values
    //Vector max_pos(4.843318e+38, 7.459491e+24);  // Using the largest x and y values

    /*
    // Dynamic bounds based on initial body positions
    auto [min_pos, max_pos] = universe.exposeBounds();
    double padding = 0.1 * (max_pos - min_pos).norm();
    min_pos = min_pos - Vector(padding, padding);
    max_pos = max_pos + Vector(padding, padding);
    */

    double cell_size = (max_pos.data[0] - min_pos.data[0]) / grid_size;

    std::vector<std::vector<double>> grid(grid_size, std::vector<double>(grid_size, 0.0));
    std::vector<std::vector<double>> potential(grid_size, std::vector<double>(grid_size, 0.0));

    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * grid_size * grid_size);
    fftw_plan forward_plan = fftw_plan_dft_2d(grid_size, grid_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan backward_plan = fftw_plan_dft_2d(grid_size, grid_size, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int step = 0; step < STEP_COUNT; ++step) {
        // Clear grid mass density
        for (auto &row : grid){
        std::fill(row.begin(), row.end(), 0.0);
        }

        //for TSC method
        auto tsc_weight = [](double dx, double h) -> double {double r = std::abs(dx) / h;
        if (r < 0.5) {
                return 0.75 - r * r;
        } else if (r < 1.5) {
                return 0.5 * std::pow(1.5 - r, 2);
        } else {
                return 0.0;
            }
        };

        // Deposit mass to grid
        for (const auto &body : universe.bodies) {

            /*
            // NGP (Nearest Grid Point Method)
            int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
            int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);
            if (grid_x >= 0 && grid_x < grid_size && grid_y >= 0 && grid_y < grid_size) {
                grid[grid_x][grid_y] += body.m;
            } else if (DEBUG) {
                std::cerr << "Body " << body.title << " is out of bounds at step " << step << "\n";
            }
            */
            

            //CIC (Cloud in Cell)
            /*
            double fx = (body.coordinates.data[0] - min_pos.data[0]) / cell_size;
            double fy = (body.coordinates.data[1] - min_pos.data[1]) / cell_size;

            int i = static_cast<int>(fx);
            int j = static_cast<int>(fy);

            double dx = fx - i;
            double dy = fy - j;

            // Distribute mass to 4 surrounding cells with CIC weights
            if (i >= 0 && i < grid_size - 1 && j >= 0 && j < grid_size - 1) {
                grid[i][j]         += body.m * (1 - dx) * (1 - dy);
                grid[i+1][j]       += body.m * dx * (1 - dy);
                grid[i][j+1]       += body.m * (1 - dx) * dy;
                grid[i+1][j+1]     += body.m * dx * dy;
            }   else if (DEBUG) {
                std::cerr << "Body " << body.title << " is out of bounds at step " << step << "\n";
            }
            */
            

            //TSC method
            
            double x = (body.coordinates.data[0] - min_pos.data[0]) / cell_size;
            double y = (body.coordinates.data[1] - min_pos.data[1]) / cell_size;

            int i_center = static_cast<int>(std::round(x));
            int j_center = static_cast<int>(std::round(y));

            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                    int i = i_center + di;
                    int j = j_center + dj;

                    if (i >= 0 && i < grid_size && j >= 0 && j < grid_size) {
                        double dx = (i + 0.5) - x;
                        double dy = (j + 0.5) - y;

                        double weight_x = tsc_weight(dx, 1.0);  // normalized over cell size
                        double weight_y = tsc_weight(dy, 1.0);

                        grid[i][j] += body.m * weight_x * weight_y;
                    }
                }
            } 
            
        }
        // Copy grid mass into FFT input
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                in[i * grid_size + j][0] = grid[i][j];
                in[i *  grid_size + j][1] = 0.0;
            }
        }
        // Forward FFT
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
        // Inverse FFT to compute potential
        fftw_execute(backward_plan);
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                potential[i][j] =  in[i * grid_size + j][0] / (grid_size * grid_size);
            }
        }
        std::cout<< "cell_size" << cell_size <<"\n";

        // Compute acceleration from potential gradients
        for (auto &body : universe.bodies) {
            /*
            // NGP interpolation
            int grid_x = static_cast<int>((body.coordinates.data[0] - min_pos.data[0]) / cell_size);
            int grid_y = static_cast<int>((body.coordinates.data[1] - min_pos.data[1]) / cell_size);
            if (grid_x >= 1 && grid_x < grid_size - 1 && grid_y >= 1 && grid_y < grid_size - 1) {
                double fx = -(potential[grid_x + 1][grid_y] - potential[grid_x - 1][grid_y]) / (2.0 * cell_size);
                double fy = -(potential[grid_x][grid_y + 1] - potential[grid_x][grid_y - 1]) / (2.0 * cell_size);
                body.acceleration = Vector(fx, fy);
                std::cout <<"fx, fy: " << fx <<" "<< fy <<" for Body " << body.title << "\n";
            } else {
                body.acceleration = Vector(0.0, 0.0);
            }
            */
            //TSC method 
            double x = (body.coordinates.data[0] - min_pos.data[0]) / cell_size;
            double y = (body.coordinates.data[1] - min_pos.data[1]) / cell_size;
            //double x = (body.coordinates.data[0] ) / cell_size;
            //double y = (body.coordinates.data[1] ) / cell_size;
            std::cout <<"x, y: " << x <<" "<< y <<" for Body " << body.title << "\n";
            int i_center = static_cast<int>(std::round(x));
            int j_center = static_cast<int>(std::round(y));
            double fx = 0.0, fy = 0.0;
            for (int di = -1; di <= 1; ++di) {
                for (int dj = -1; dj <= 1; ++dj) {
                        int i = i_center + di;
                        int j = j_center + dj;

                    if (i >= 1 && i < grid_size - 1 && j >= 1 && j < grid_size - 1) {
                            double dx = (i + 0.5) - x;
                            double dy = (j + 0.5) - y;

                            double weight_x = tsc_weight(dx, 1.0);
                            double weight_y = tsc_weight(dy, 1.0);

                            // Central differences of potential for gradient
                        double dPhidx = (potential[i + 1][j] - potential[i - 1][j]) / (2.0 * cell_size);
                        double dPhidy = (potential[i][j + 1] - potential[i][j - 1]) / (2.0 * cell_size);

                        fx += -dPhidx * weight_x * weight_y;
                        fy += -dPhidy * weight_x * weight_y;
                        if(!DEBUG){
                            std::cout <<"x, y: " << x <<" "<< y <<" for Body " << body.title << "\n";
                            //std::cout <<"Forces DPhidx, DPhidy: " << dPhidx <<" "<< dPhidy <<" for Body " << body.title << "\n";
                            //std::cout <<"Forces f_x, f_y: " << fx <<" "<< fy <<" for Body " << body.title << ")\n";
                        }
                    }
                }
            }
        
            std::cout <<"fx, fy: " << fx <<" "<< fy <<" for Body " << body.title << "\n";
            body.acceleration = Vector(fx, fy);  
        }  
        std::vector<Vector> positions;
        for (auto &body : universe.bodies) {
            body.update(dt);
            positions.push_back(body.coordinates);
        }
        universe.telemetry.push_back(positions);

        // Debug output
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
