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

    // Dispatch to the appropriate simulation method
    if (method == "naive") {
        naive_simulation(universe, 0.001);
    }
    else if (method == "barneshut") {
        BarnesHut(universe, 0.001);
    }
#ifdef USE_CUDA
    else if (method == "gpu") {
        simulateBruteForceGPU(universe, 0.001);
    }
#endif
    else {
        std::cerr << "Unknown method “" << method << "”\n";
        return 1;
    }

    universe.visualize(method + "_simulation.gif");
    return 0;
}