#include "src/core.hpp"
#include "src/simplesimulation.hpp"
#include "src/particlemesh.hpp"
#include "src/barneshutt.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

int main() {
    Magick::InitializeMagick(nullptr);
    MagickCore::SetMagickResourceLimit(MagickCore::AreaResource, 1024*1024*256);  // 256MB
    MagickCore::SetMagickResourceLimit(MagickCore::MemoryResource, 1024*1024*512);  // 512MB
    MagickCore::SetMagickResourceLimit(MagickCore::MapResource, 1024*1024*512);    // 512MB
    
    System universe;
    
    // Initial conditions
    double initial_velocity = 100.0;  // m/s
    double initial_velocity2 = 200.0;  // m/s

    double angle = 45.0;  // degrees
    double angle2 = 60.0; 
    double g = 9.81;  // m/s^2
    
    // Convert angle to radians
    double angle_rad = angle * M_PI / 180.0;
    double angle_rad2 = angle2 * M_PI / 180.0;
    
    // Initial velocity components
    double v0x = initial_velocity * cos(angle_rad);
    double v0y = initial_velocity * sin(angle_rad);
    
    double v0x2 = initial_velocity2 * cos(angle_rad2);
    double v0y2 = initial_velocity2 * sin(angle_rad2);

    
    // Create two projectiles with different starting positions
    Body projectile1(1.0,  // mass (kg)
                    Vector(0, 0),  // starting at origin
                    Vector(v0x, v0y),  // initial velocity
                    "red", 4, "proj 1");  // color and size

    Body projectile2(2.0,  // mass (kg)
                    Vector(50, 0),  // starting 5m to the right
                    Vector(v0x2, v0y2),  // same initial velocity
                    "blue", 3, "proj 2");  // different color and size

    universe.add(projectile1);
    universe.add(projectile2);
    
    // Simulation parameters
    double dt = 0.1;  // time step (s)
    int steps = 10000;  // number of steps
    
    // Clear any existing telemetry
    universe.telemetry.clear();
    
    // Manually compute trajectories
    for (int step = 0; step < steps; ++step) {
        double t = step * dt;
        
        // Vector to store positions of all projectiles at this time step
        std::vector<Vector> positions;
        
        // Compute position for first projectile
        double x1 = v0x * t;
        double y1 = v0y * t - 0.5 * g * t * t;
        positions.push_back(Vector(x1, y1));
        
        // Compute position for second projectile
        double x2 = 5 + v0x2 * t;  // Start 5m to the right
        double y2 = v0y2 * t - 0.5 * g * t * t;
        positions.push_back(Vector(x2, y2));
        
        // Store positions in telemetry
        universe.telemetry.push_back(positions);
        
        // Stop if both projectiles hit the ground
        if (y1 < 0 || y2 < 0) break;
    }

    // Print positions for both projectiles
    for (size_t i = 0; i < universe.telemetry.size(); i += 1) {
        std::cout << "Time: " << i*dt << "s\n";
        for (size_t j = 0; j < universe.telemetry[i].size(); j++) {
            const Vector& pos = universe.telemetry[i][j];
            std::cout << "  Projectile " << j+1 << ": ("
                      << pos.data[0] << ", " << pos.data[1] << ")\n";
        }
    }
    
    universe.visualize("telemetry_test.gif", true, true);
    return 0;
}