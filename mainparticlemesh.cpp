#include "src/core.hpp"
#include "src/particlemesh.hpp"
#include "src/particlemesh.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

// Function to export telemetry data to a CSV file
void export_telemetry_to_csv(const System &universe, const std::string &filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open " << filename << " for writing.\n";
        return;
    }

    // Write header
    file << "step,body,x,y\n";

    for (size_t step = 0; step < universe.telemetry.size(); ++step) {
        const auto &positions = universe.telemetry[step];
        for (size_t i = 0; i < positions.size(); ++i) {
            file << step << "," << i << "," << positions[i].data[0] << "," << positions[i].data[1] << "\n";
        }
    }

    file.close();
    std::cout << "Telemetry data exported to " << filename << std::endl;
}


int main() {
    // Create a system
    System universe;

    // Add some bodies to the system
    universe.add(Body(1.0, Vector(0, 0), Vector(0, 0), "yellow", 10, "Body 1"));
    universe.add(Body(2.0, Vector(100, 0), Vector(0, 1), "blue", 10, "Body 2"));
    universe.add(Body(3.0, Vector(0, 100), Vector(-1, 0), "red", 10, "Body 3"));

    // Run the particle mesh simulation
    double dt = 0.1;
    int grid_size = 100;
    particle_mesh_simulation(universe, dt, grid_size);

    export_telemetry_to_csv(universe, "telemetry.csv");



    //std::cout << "Vizualization going to start" << std::endl;

    // Visualize the simulation and store it in a GIF file
    //universe.visualize("simulation_particle_mesh.gif", true, true);


    std::cout << "Simulation completed and stored in simulation_particle_mesh.gif" << std::endl;

    return 0;
}
