#include "simplesimulation.hpp"

void naive_simulation(System &universe, double dt) {
    universe.telemetry.clear();
    for (int step = 0; step < STEP_COUNT; ++step) {
        for (size_t i = 0; i < universe.bodies.size(); ++i) {
            for (size_t j = 0; j < universe.bodies.size(); ++j) {
                if (i != j) {
                    Vector f = force(universe.bodies[i], universe.bodies[j]);
                    universe.bodies[i].acceleration += f * (1.0/universe.bodies[i].m);
                }
            }
        }
        std::vector<Vector> positions;
        for (auto& body : universe.bodies) {
            body.update(dt);
            positions.push_back(body.coordinates);
        }
        universe.telemetry.push_back(positions);
    }
}