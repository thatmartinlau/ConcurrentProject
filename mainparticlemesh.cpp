
#include "src/core.hpp"
#include "src/particlemesh.hpp"
#include "src/partickemesh.hpp"
#include "src/particlemeshcuda.cpp"

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

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
    particle_mesh_simulation_CUDA(universe, dt, grid_size, 3);

    std::cout << "Vizualization going to start" << std::endl;

    // Visualize the simulation and store it in a GIF file
    //universe.visualize("simulation2.gif", true, true);

    std::cout << "Simulation completed and stored in simulation.gif" << std::endl;

    return 0;
}
