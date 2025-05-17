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
    Magick::InitializeMagick(nullptr);
    MagickCore::SetMagickResourceLimit(MagickCore::AreaResource, 1024*1024*256);  // 256MB
    MagickCore::SetMagickResourceLimit(MagickCore::MemoryResource, 1024*1024*512);  // 512MB
    MagickCore::SetMagickResourceLimit(MagickCore::MapResource, 1024*1024*512);    // 512MB
    
    System universe;
    
    // Sun at approximate center
    Body sun(1.989e30,  // Mass of Sun in kg
        Vector(0, 0),   // Position at origin
        Vector(0, 0), "Orange", 10);  // Nearly stationary

    // Earth
    Body earth(5.972e24,      // Mass of Earth in kg
        Vector(149.6e9, 0),   // Position at 1 AU on x-axis
        Vector(0, 29780), "Blue", 3);    // Orbital velocity primarily in y direction

    // Mars
    Body mars(6.39e23,        // Mass of Mars in kg
        Vector(227.9e9, 0),   // Position at 1.52 AU on x-axis
        Vector(0, 24077), "Red", 4);    // Orbital velocity primarily in y direction

    universe.add(sun);
    universe.add(earth);
    universe.add(mars);

    naive_simulation(universe, 0.001);
    universe.visualize("naive_simulation.gif");
    
    return 0;
}