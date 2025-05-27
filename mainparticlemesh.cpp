#include "src/core.hpp"
#include "src/simplesimulation.hpp" 
#include "src/barneshutt.hpp"
#include "src/particlemesh.hpp"
#include "src/particlemesh_thread.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <chrono>

int main(int argc, char** argv) {
    // Default choice of method
    //std::string method = "particlemesh_thread";
   std::string method =  "particlemesh";
    // Parse args for different methods
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const std::string prefix = "-method=";
        if (arg.rfind(prefix, 0) == 0) {
            method = arg.substr(prefix.size());
        }
    }

    Magick::InitializeMagick(nullptr);
    MagickCore::SetMagickResourceLimit(MagickCore::AreaResource, 1024*1024*512);  // 256MB
    MagickCore::SetMagickResourceLimit(MagickCore::MemoryResource, 1024*1024*1024);  // 512MB
    MagickCore::SetMagickResourceLimit(MagickCore::MapResource, 1024*1024*1024);    // 512MB
    
    System universe;
    

    // Central mass (a star-like object)
    Body central(900.0,
        Vector(0.0, 0.0),  // Center
        Vector(0.0, 0.0),  // Stationary
        "yellow", 10, "Central");

    // Small orbiting body 1
    Body orbiter1(1.0,
        Vector(10.0, 0.0),   // 10 units away on x-axis
        Vector(0.0, 3.16),   // Roughly circular orbit
        "blue", 3, "Orbiter1");

    // Small orbiting body 2
    Body orbiter2(0.75,
        Vector(-15.0, 0.0),   // 15 units away on -x-axis
        Vector(0.0, -2.58),   // Opposite orbit
        "green", 3, "Orbiter2");

    Body orbiter3(0.5,
    Vector(5.0, 5.0),
    Vector(0, 1.4),
    "red", 3, "Orbiter3");

    universe.add(central);
    universe.add(orbiter1);
    universe.add(orbiter2);
    universe.add(orbiter3);

    double dt = 0.2; // fuck it, one hour
    universe.dt = dt;

    // Dispatch to the appropriate simulation method
    auto start = std::chrono::high_resolution_clock::now();
    if (method == "naive") {
        naive_simulation(universe);
    }
    else if (method == "barneshut") {
        // BarnesHut(universe, DT);
    }
    else if (method == "particlemesh"){
        int grid_size = 10; // added this 
        double R = 20000;
        particle_mesh_simulation(universe, dt,grid_size, R);
    }else if (method == "particlemesh_thread"){
        int grid_size =  10;
        double R = 20000;
        // number of cores of 8
        size_t num_threads = 5;
        particle_mesh_simulation_parallel(universe, dt, grid_size,num_threads, R);
    }
    #ifdef USE_CUDA
    else if (method == "gpu") {
        simulateBruteForceGPU(universe, DT);
    }
    #endif
    else {
        std::cerr << "Unknown method “" << method << "”\n";
        return 1;
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // For testing purposes
    bool print_telemetry = false;
    if (print_telemetry) {
        for (size_t i = 0; i < universe.telemetry.size(); i += 1) {
            std::cout << "Frame: " << i << ", Time: " << i*universe.dt << "s\n";
            for (size_t j = 0; j < universe.telemetry[i].size(); j++) {
                const Vector& pos = universe.telemetry[i][j];
                std::cout << "  Position " << universe.bodies[j].title << ": ("
                        << pos.data[0] << ", " << pos.data[1] << ")\n";
            }
        }
    }
    
    bool record_csv = true;
    if (record_csv) {
        // Export to CSV    
        std::ofstream file("telemetry.csv");
        if (!file.is_open()) {
            std::cerr << "Error: Could not open telemetry.csv for writing\n";
            return 1;
        }

        // Add CSV header
        file << "time";
        for (const auto& body : universe.bodies) {
            file << "," << body.title << "_x," << body.title << "_y";
        }
        file << "\n";

        // Write data
        for (size_t i = 0; i < universe.telemetry.size(); i++) {
            file << std::scientific << std::setprecision(6)  // Use scientific notation for large numbers
                << i * universe.dt; // time
            for (const auto& pos : universe.telemetry[i]) {
                file << "," << pos.data[0] << "," << pos.data[1];
            }
            file << "\n";
        }
        file.close();
    }

    std::cout << "Simulation done. Generating the visualization...\n";
    string out_name;
    if(method == "particlemesh"){
            out_name = "particle_mesh.gif";
    }
    else if(method == "particlemesh_thread"){
          out_name = "particle_mesh_thread.gif";
    }else{
         out_name = "rockyplanets.gif";
    }
 
    universe.visualize(out_name, true, true);
    std::cout << "Done! \n";

    std::cout << "\n\nSimulation time: " << time_taken.count() << " milliseconds.\n";
    return 0;
}