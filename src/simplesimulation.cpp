#include "simplesimulation.hpp"
#include <iostream>
#include <thread>
#include <vector>

#define DEBUG false

#define N_THREADS 20

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
            positions.push_back(body.coordinates); // This part is problematic for multithreading.
        }
        universe.telemetry.push_back(positions);   // This part is problematic for multithreading.
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
- create a copy of the system. (from 1 to N)
- Split as follows: -> Includes force calculation + update step. (from 1 to N/Thread_num)
    - split all the bodies on different number of threads
    - Using the original universe, compute the accumulation of accelerations for each body for that 
      thread.
    - Once accumulated acceleration computation finished, update the original system coordinate with 
      acceleration.
- Join vectors

!!! Need to separate auxiliary functions for force computation and updating positions:
    - Potential race condition for computing the forces if the bodies are updated too early.
      (realistically, in a planetary simulation, this hardly matters since the pull of other things
       is negligible compare to the pull of the sun, but it's the principle that counts.)
*/

void optimized_force_aux(const System& universe,  // Only read
                        System& universe_copy,    // Modify
                        const int start, const int end) {
    // Take copy of the universe, and compute all the forces for the concerned bodies.
    // We are only allowed to change the values inside universe_copy.
    
    // Reset the acceleration accumulation vector.
    for (int i = start; i < end; i++) {   
        universe_copy.bodies[i].acceleration = Vector(0,0);
    }

    // Compute all the forces, just like last time.
    for (int j = start; j < end; ++j) {
        for (size_t k = 0; k < universe_copy.bodies.size(); ++k) {
            if (j == k) {continue;}
            Body& body1 = universe_copy.bodies[j];
            Body& body2 = universe_copy.bodies[k];
            
            Vector f = force(universe_copy.bodies[j], universe_copy.bodies[k]);
            body1.acceleration += f / body1.m;   // Newton's 3rd law, duh
            // body2.acceleration += f / (-body2.m);

        }
    }
}

void optimized_update_aux(System& universe,          // Modify
                        const System& universe_copy, // Only read
                        const int start, const int end, std::vector<Vector>& positions,
                        const double dt) {
    /*     
        Calling Body.update will not cause any race conditions, as each thread is exclusively
        calling the  update function for each  body in our system, and the update method does
        not depend on anything outside of the Body class.
    */

    for (int i = start; i < end; i++) {
        universe.bodies[i].acceleration = universe_copy.bodies[i].acceleration;
        universe.bodies[i].update(dt);
        positions[i] = universe.bodies[i].coordinates;
    }
    // By this point, all bodies in the original universe would have their updated positions, 
    // Position vector would have everything in the correct position, as well.
}
 
void optimized_simulation(System &universe) {

    // Thread initialization
    const int length = universe.bodies.size();  // This value should be N. 
    const int block_size = length / N_THREADS;
    std::vector<std::thread> workers(N_THREADS-1);


    // Keep the same initialization as before.
    universe.telemetry.clear();
    universe.telemetry.resize(STEP_COUNT+1);

    std::vector<Vector> initial_positions(length); 
    for (int i = 0; i < length; ++i) {
        initial_positions[i] = universe.bodies[i].coordinates;
    }
    universe.telemetry[0] = initial_positions;


    // Loop through the steps here!
    for (int step=0; step < STEP_COUNT; ++step) {
        // Inside one step!
        // Make a copy of the given system
        System universe_copy = universe;
        std::vector<Vector> positions(length); // Pre initialize this this time, to conserve the order.

        int start_block = 0;
        for (int i = 0; i<N_THREADS-1; ++i) {
            int end_block = start_block + block_size;
            workers[i] = std::thread(optimized_force_aux, std::ref(universe), std::ref(universe_copy), 
                                    start_block, end_block);
            start_block = end_block;
        }
        optimized_force_aux(universe, universe_copy, start_block, length);
        for (int i=0; i < N_THREADS-1; ++i) {
            workers[i].join();
        }

        // By this point, universe_copy will contain all the updated positions of the bodies.
        start_block = 0;
        for (int i = 0; i<N_THREADS-1; ++i) {
            int end_block = start_block + block_size;
            workers[i] = std::thread(optimized_update_aux, std::ref(universe), std::ref(universe_copy),
                                    start_block, end_block, std::ref(positions), universe.dt);
            start_block = end_block;
        }
        optimized_update_aux(universe, universe_copy, start_block, length, positions, universe.dt);
        for (int i=0; i < N_THREADS-1; ++i) {
            workers[i].join();
        }

        universe.telemetry[step] = positions; // This doesn't cause race conditions.
        if (DEBUG) {
            std::cout << "Step: " << step << "\n";
            std::cout << "Telemetry size: " << universe.telemetry.size() << "\n";
            
            // Check if we have data for this step
            if (step < universe.telemetry.size()) {
                std::cout << "Bodies in this step: " << universe.telemetry[step].size() << "\n";
                
                // Print positions for each body
                for (size_t i = 0; i < universe.telemetry[step].size(); i++) {
                    const Vector& pos = universe.telemetry[step][i];
                    std::cout << "Body " << i << " position: ("
                            << pos.data[0] << ", " 
                            << pos.data[1] << ")\n";
                    
                    // Optional: Compare with actual body positions
                    std::cout << "  Current body position: ("
                            << universe.bodies[i].coordinates.data[0] << ", "
                            << universe.bodies[i].coordinates.data[1] << ")\n";
                }
            } else {
                std::cout << "ERROR: No telemetry data for step " << step << "\n";
            }
            
            // Add separator for readability
            std::cout << "----------------------------------------\n";
        }
    }
    // All steps done, Simulation done.
}