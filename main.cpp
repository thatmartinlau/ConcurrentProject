#include "src/core.hpp"
#include "src/simplesimulation.hpp"
#include "src/barneshutt.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

int main() {
    
    /* 
    System universe;
    
    // Sun at approximate center
    Body sun(1.989e30,  // Mass of Sun in kg
        Vector(0, 0),   // Position at origin
        Vector(0, 0));  // Nearly stationary

    // Earth
    Body earth(5.972e24,      // Mass of Earth in kg
        Vector(149.6e9, 0),   // Position at 1 AU on x-axis
        Vector(0, 29780));    // Orbital velocity primarily in y direction

    // Mars
    Body mars(6.39e23,        // Mass of Mars in kg
        Vector(227.9e9, 0),   // Position at 1.52 AU on x-axis
        Vector(0, 24077));    // Orbital velocity primarily in y direction

    universe.add(sun);
    universe.add(earth);
    universe.add(mars);

    naive_simulation(universe, 0.001);
    universe.visualize("naive_simulation.gif");
    */

    System universe;
    
    // Initial conditions
    double initial_velocity = 100.0;  // m/s
    double angle = 45.0;  // degrees
    double g = 9.81;  // m/s^2
    
    // Convert angle to radians
    double angle_rad = angle * M_PI / 180.0;
    
    // Initial velocity components
    double v0x = initial_velocity * cos(angle_rad);
    double v0y = initial_velocity * sin(angle_rad);
    
    // Create a body (mass doesn't matter for this simulation)
    Body projectile(1.0,  // mass (kg)
                   Vector(0, 0),  // starting at origin
                   Vector(v0x, v0y));  // initial velocity
    
    universe.add(projectile);
    
    // Simulation parameters
    double dt = 0.01;  // time step (s)
    int steps = 1000;  // number of steps
    
    // Clear any existing telemetry
    universe.telemetry.clear();
    
    // Manually compute trajectory
    for (int step = 0; step < steps; ++step) {
        double t = step * dt;
        
        // Compute position at time t using projectile motion equations
        double x = v0x * t;
        double y = v0y * t - 0.5 * g * t * t;
        
        // Store position in telemetry
        std::vector<Vector> positions;
        positions.push_back(Vector(x, y));
        universe.telemetry.push_back(positions);
        
        // Stop if we hit the ground
        if (y < 0) break;
    }

    for (size_t i = 0; i < universe.telemetry.size(); i += 50) {
        const Vector& pos = universe.telemetry[i][0];
        std::cout << "Time: " << i*dt << "s, Position: ("
                  << pos.data[0] << ", " << pos.data[1] << ")\n";
    }
    universe.visualize("out.gif");
    
    return 0;
}