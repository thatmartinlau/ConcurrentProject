#include "src/core.hpp"
#include "src/simplesimulation.hpp"
#include "src/barneshutt.hpp"
#include "src/particlemesh.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

#ifndef VISUALIZE
#define VISUALIZE false
#endif 

#ifndef PRINT_TELEMETRY
#define PRINT_TELEMETRY false
#endif 

#ifndef EXPORT_CSV
#define EXPORT_CSV false
#endif 

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>
#include <chrono>

#define ASTEROIDS 200
#define DT 3600. // Not suggested to take a timestep larger than three hours - orbits start getting weird.

int main(int argc, char** argv) {
    // Default choice of method
    std::string method = "naive";
    bool useParallelBH = false;

    // Parse args for different methods
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        const std::string prefix = "-method=";
        if (arg.rfind(prefix, 0) == 0) {
            method = arg.substr(prefix.size());
        }
        else if (arg == "-parallel") {
            useParallelBH = true;
        }
    }

    Magick::InitializeMagick(nullptr);
    
    System universe;
    

    // Define planets
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

    // Not a good thing to add to the system, it strongly messes with the scaling inside the visualization. 
    // (Yes, a solution would be to change System::getBounds by taking the farthest distance point from the sun,
    // then projecting that in a ray across the sun, but I'm not paid enough to do that) - Martin
    // universe.add(jupiter);
    // universe.add(neptune);
    // universe.add(saturn);
    // universe.add(uranus);
    // universe.add(pluto);


    /* 
    Can't be asked to think of a better performance metric for testing the 
    efficiency of our programs, so we'll be testing increasing numbers of 
    asteroids to show how complexity increases super fast as we increase 
    number of bodies. - Martin
    */
    // Setup random number generation
    std::random_device rd;  // Used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine

    const double AU = 149.6e9;           // 1 AU in meters
    const double MIN_DIST = 2.2 * AU;    // Inner asteroid belt (~2.2 AU)
    const double MAX_DIST = 3.2 * AU;    // Outer asteroid belt (~3.2 AU)
    const double AVG_VELOCITY = 17500;   // Average orbital velocity in m/s

    // Create distributions for each random value
    std::uniform_real_distribution<> mass_dist(1e13, 1e17);
    std::uniform_real_distribution<> dist_dist(MIN_DIST, MAX_DIST);
    std::uniform_real_distribution<> angle_dist(0, 2 * M_PI);

    for (int i = 0; i < ASTEROIDS; i++) {
        double mass = mass_dist(gen);
        double distance = dist_dist(gen);
        double angle = angle_dist(gen);
        Vector position(distance * cos(angle), distance * sin(angle));
        double velocity_magnitude = sqrt((6.67430e-11 * 1.989e30) / distance);

        Vector velocity(-velocity_magnitude * sin(angle), velocity_magnitude * cos(angle));
        
        Body asteroid(mass, position, velocity, "green", 1, "");
        universe.add(asteroid);
    }

    universe.dt = DT;
    System universe_multithreaded = universe;
    System universe_multithreaded2 = universe;

    std::cout << "Starting the simulation: " << ASTEROIDS << " asteroids\n";

    // Dispatch to the appropriate simulation method: Do a simulation for unoptimized, then one for optimized.
    
    if (method == "naive") {
        auto start_seq = std::chrono::high_resolution_clock::now();
        naive_simulation(universe);
        auto end_seq = std::chrono::high_resolution_clock::now();
        auto time_taken_seq = std::chrono::duration_cast<std::chrono::milliseconds>(end_seq - start_seq);
        std::cout << "Simulation time sequential: " << time_taken_seq.count() << " milliseconds.\n";

        auto start_par = std::chrono::high_resolution_clock::now();
        optimized_simulation(universe_multithreaded);
        auto end_par = std::chrono::high_resolution_clock::now();
        auto time_taken_par = std::chrono::duration_cast<std::chrono::milliseconds>(end_par - start_par);
        std::cout << "Simulation time parallel: " << time_taken_par.count() << " milliseconds.\n";

        auto start_par2 = std::chrono::high_resolution_clock::now();
        optimized_simulationmk2(universe_multithreaded2);
        auto end_par2 = std::chrono::high_resolution_clock::now();
        auto time_taken_par2 = std::chrono::duration_cast<std::chrono::milliseconds>(end_par2 - start_par2);
        std::cout << "Simulation2 time parallel: " << time_taken_par2.count() << " milliseconds.\n";

    }
    else if (method == "barneshut") {
        auto start_bh = std::chrono::high_resolution_clock::now();

        universe.telemetry.clear();
        {
            std::vector<Vector> init;
            for (auto& b : universe.bodies)
                init.push_back(b.coordinates);
            universe.telemetry.push_back(init);
        }

        for (int step = 0; step < STEP_COUNT; ++step) {
            BarnesHutStep(universe.bodies, universe.dt, 0.5, useParallelBH);

            std::vector<Vector> frame;
            for (auto& b : universe.bodies) {
                frame.push_back(b.coordinates);
            }
            universe.telemetry.push_back(frame);
        }

        auto end_bh = std::chrono::high_resolution_clock::now();
        auto time_bh =
            std::chrono::duration_cast<std::chrono::milliseconds>(end_bh - start_bh);
        std::cout << "Barnes-Hut simulation time: "
                  << time_bh.count() << " milliseconds.\n";
    }
    else if (method == "particlemesh"){
        int grid_size = 100; // added this 
        // particle_mesh_simulation(universe, dt,grid_size);
    }
    else {
        std::cerr << "Unknown method “" << method << "”\n";
        return 1;
    }
    

    // For testing purposes
    bool print_telemetry = PRINT_TELEMETRY;
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
    
    bool record_csv = EXPORT_CSV;
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

    std::cout << "Simulations done.\n";

    bool visualization = VISUALIZE;
    if (visualization) {
        std::cout << "Generating the visualizations...\n";
        string out_name = "nbody_simulation";
        auto start2 = std::chrono::high_resolution_clock::now();
        universe.visualize2(out_name, false, false);
        auto end2 = std::chrono::high_resolution_clock::now();
        auto time_taken2 = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2);
        std::cout << "Visualization time:" << time_taken2.count() << " milliseconds.\n";

    } 

    return 0;
}
