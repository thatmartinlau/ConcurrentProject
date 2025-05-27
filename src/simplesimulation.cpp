#include "simplesimulation.hpp"
#include <iostream>

#define DEBUG false

void naive_simulation(System &universe) {
    universe.telemetry.clear();
    std::vector<Vector> initial_positions;
    for (const auto& body : universe.bodies) {
        initial_positions.push_back(body.coordinates);
    }
    universe.telemetry.push_back(initial_positions);

    for (int step = 0; step < STEP_COUNT; ++step) {
        for (auto& body : universe.bodies) {
            body.acceleration = Vector(0, 0);
        }

        for (size_t i = 0; i < universe.bodies.size(); ++i) {
            for (size_t j = i+1; j < universe.bodies.size(); ++j) {
                Body& body1 = universe.bodies[i];
                Body& body2 = universe.bodies[j];

                Vector f = force(universe.bodies[i], universe.bodies[j]);
                body1.acceleration += f / body1.m;
                body2.acceleration += f / (-body2.m);

            }
        }
        std::vector<Vector> positions;
        for (auto& body : universe.bodies) {
            body.update(universe.dt);
            positions.push_back(body.coordinates);
        }
        universe.telemetry.push_back(positions);
        if (DEBUG) {
            std::cout << "Step:" << step << "\n";
            for (const auto& body : universe.bodies) {
                std::cout << universe.dt << "    Body " << body.title 
                        << " velocity: (" << body.velocity.data[0] 
                        << ", " << body.velocity.data[1] << ")"
                        << " acceleration: (" << body.acceleration.data[0]
                        << ", " << body.acceleration.data[1] << ")\n";
            }
        }
        

    }
    std::cout << universe.telemetry.size() << "\n";
    
}