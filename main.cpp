#include "src/core.hpp"
#include "src/simplesimulation.hpp"
#include "src/barneshutt.hpp"
#include "src/particlemesh.hpp"

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
    std::string method = "naive";
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
    
    // Sun at approximate center
    Body sun(1.989e30,  // Mass of Sun in kg
        Vector(0, 0),   // Position at origin
        Vector(0, 0), "orange", 10, "Sun");  // Nearly stationary

    // Mercury
    Body mercury(3.285e23,    // Mass in kg
    Vector(57.9e9, 0),    // Position at 0.387 AU
    Vector(0, 43360),     // Orbital velocity
    "gray", 2, "Mercury");

    // Venus
    Body venus(4.867e24,      // Mass in kg
    Vector(108.2e9, 0),   // Position at 0.723 AU
    Vector(0, 35020),     // Orbital velocity
    "yellow", 3, "Venus");

    // Earth
    Body earth(5.972e24,      // Mass of Earth in kg
        Vector(149.6e9, 0),   // Position at 1 AU on x-axis
        Vector(0, 29780), "blue", 3, "Earth");    // Orbital velocity primarily in y direction

    // Mars
    Body mars(6.39e23,        // Mass of Mars in kg
        Vector(227.9e9, 0),   // Position at 1.52 AU on x-axis
        Vector(0, 24077), "red", 4, "Mars");    // Orbital velocity primarily in y direction

    // Jupiter
    Body jupiter(1.898e27,    // Mass in kg
    Vector(778.5e9, 0),   // Position at 5.203 AU
    Vector(0, 13070),     // Orbital velocity
    "brown", 7, "Jupiter");

    // Saturn
    Body saturn(5.683e26,     // Mass in kg
        Vector(1.434e12, 0),  // Position at 9.582 AU
        Vector(0, 9680),      // Orbital velocity
        "tan", 6, "Saturn");

    // Uranus
    Body uranus(8.681e25,     // Mass in kg
        Vector(2.871e12, 0),  // Position at 19.201 AU
        Vector(0, 6800),      // Orbital velocity
        "lightblue", 4, "Uranus");

    // Neptune
    Body neptune(1.024e26,    // Mass in kg
        Vector(4.495e12, 0),  // Position at 30.047 AU
        Vector(0, 5430),      // Orbital velocity
        "blue", 4, "Neptune");

    // Pluto (technically a dwarf planet now, but included for completeness)
    Body pluto(1.309e22,      // Mass in kg
        Vector(5.906e12, 0),  // Position at 39.482 AU
        Vector(0, 4670),      // Orbital velocity
        "gray", 1, "Pluto");

    

    universe.add(sun);
    universe.add(mercury);
    universe.add(venus);
    universe.add(earth);
    universe.add(mars);
    // universe.add(jupiter);
    // universe.add(neptune);
    // universe.add(saturn);
    // universe.add(uranus);
    // universe.add(pluto);

    const double AU = 149.6e9;           // 1 AU in meters
    const double MIN_DIST = 2.2 * AU;    // Inner asteroid belt (~2.2 AU)
    const double MAX_DIST = 3.2 * AU;    // Outer asteroid belt (~3.2 AU)
    const double AVG_VELOCITY = 17500;   // Average orbital velocity in m/s

    // Create 50 random asteroids
    for (int i = 0; i < 50; i++) {
        // Random mass between 1e13 and 1e17 kg (typical asteroid masses)
        double mass = (rand() % 10000) * 1e13;
        double distance = MIN_DIST + (rand() % (int)(MAX_DIST - MIN_DIST));
        double angle = (rand() % 360) * M_PI / 180.0;
        Vector position(distance * cos(angle), distance * sin(angle));
        double velocity_magnitude = AVG_VELOCITY * sqrt(2.2 * AU / distance);
        Vector velocity(-velocity_magnitude * sin(angle), velocity_magnitude * cos(angle));
        
        Body asteroid(mass, position, velocity, "green", 1, "");
        
        // universe.add(asteroid);
    }


    double dt = 3600; // fuck it, one hour
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
        int grid_size = 100; // added this 
        particle_mesh_simulation(universe, dt,grid_size);
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
    string out_name = "rockyplanets.gif";
    universe.visualize(out_name, true, true);
    std::cout << "Done! \n";

    std::cout << "\n\nSimulation time: " << time_taken.count() << " milliseconds.\n";
    return 0;
}