#include "simplesimulation.hpp"
#include <iostream>

#define DEBUG false

// Basic simulation algorithm for N bodies in a system.
void naive_simulation(System &universe) {
    // Clear the previous telemetries to make sure we're on a blank slate.
    // Record initial positions, so we have a starting frame for the video.
    universe.telemetry.clear();
    std::vector<Vector> initial_positions; 
    for (const auto& body : universe.bodies) {
        initial_positions.push_back(body.coordinates);
    }
    universe.telemetry.push_back(initial_positions);
    
    // Iterate over simulation steps
    for (int step = 0; step < STEP_COUNT; ++step) {
        // This part is O(n), I don't believe this part needs to be included in 
        // complexity reductions in the following.
        for (auto& body : universe.bodies) {
            body.acceleration = Vector(0, 0);
        }

        // The force computations could be done in parallel.
        for (size_t i = 0; i < universe.bodies.size(); ++i) {
            for (size_t j = i+1; j < universe.bodies.size(); ++j) {
                Body& body1 = universe.bodies[i];
                Body& body2 = universe.bodies[j];
                
                // Here we can already avoid computing the forces for two bodies twice this way
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
}

/*
Notes on optimizations
The bulk of the optimzation needs  to happen at the force calculation step. As we add more bodies, 
the complexity of computing forces grows at a rate of n^2. On the other hand, the body update step 
has a complexity at the level of O(n), so I don't see a need to really optimize this.

We should not overdo the parallelization with threading: Given how we need to open/join new threads 
for each step (there are 15000 steps), it would not be  a good idea to open threads for every small 
thing - we benchmark performance with parallelizing update step and without.


Proposed idea:
(determined constants: number of bodies N, number of threads to use)
On each simulation step: 
- create a copy of the system.
- Split as follows: -> Includes force calculation + update step.
    - split all the bodies on different number of threads
    - Using the original universe, compute the accumulation of accelerations for each body for that thread.
    - Once accumulated acceleration computation finished, update the original system coordinate with acceleration.
- Join vectors
*/

void optimized_aux(System& universe) {
    
}

void optimized_simulation(System &universe) {
    std::cout << "Not implemented yet.\n";
}