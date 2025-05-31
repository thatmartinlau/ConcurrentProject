#include "simplesimulation.hpp"
#include <iostream>
#include <thread>
#include <vector>

#define DEBUG false

#define N_THREADS 5

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

void optimized_force_aux(System& universe,  // Only read
                        const int start, const int end) {
    // Take copy of the universe, and compute all the forces for the concerned bodies.
    // We are only allowed to change the values inside universe_copy.
    
    // Reset the acceleration accumulation vector.
    for (int i = start; i < end; i++) {   
        universe.bodies[i].acceleration = Vector(0,0);
    }

    // Compute all the forces, just like last time.
    for (int j = start; j < end; ++j) {
        for (size_t k = 0; k < universe.bodies.size(); ++k) {
            if (j == k) {continue;}
            Body& body1 = universe.bodies[j];
            Body& body2 = universe.bodies[k];
            
            Vector f = force(universe.bodies[j], universe.bodies[k]);
            body1.acceleration += f / body1.m;   // Newton's 3rd law, duh
            // body2.acceleration += f / (-body2.m);

        }
    }
}

void optimized_update_aux(System& universe,          // Modify
                        const int start, const int end, std::vector<Vector>& positions,
                        const double dt) {
    /*     
        Calling Body.update will not cause any race conditions, as each thread is exclusively
        calling the  update function for each  body in our system, and the update method does
        not depend on anything outside of the Body class.
    */

    for (int i = start; i < end; i++) {
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
        std::vector<Vector> positions(length); // Pre initialize this this time, to conserve the order.

        int start_block = 0;
        for (int i = 0; i<N_THREADS-1; ++i) {
            int end_block = start_block + block_size;
            workers[i] = std::thread(optimized_force_aux, std::ref(universe),
                                    start_block, end_block);
            start_block = end_block;
        }
        optimized_force_aux(universe, start_block, length);
        for (int i=0; i < N_THREADS-1; ++i) {
            workers[i].join();
        }

        // By this point, universe_copy will contain all the updated positions of the bodies.
        start_block = 0;
        for (int i = 0; i<N_THREADS-1; ++i) {
            int end_block = start_block + block_size;
            workers[i] = std::thread(optimized_update_aux, std::ref(universe), 
                                    start_block, end_block, std::ref(positions), universe.dt);
            start_block = end_block;
        }
        optimized_update_aux(universe, start_block, length, positions, universe.dt);
        for (int i=0; i < N_THREADS-1; ++i) {
            workers[i].join();
        }

        universe.telemetry[step] = positions; // This doesn't cause race conditions.
    }
    // All steps done, Simulation done.
}


/*
    (Martin's home PC - performed on Core i5-8500 6 Core CPU)
    Simple sequential simulation
    10 asteroids: 63ms
    50 asteroids: 785ms
    100 asteroids: 2873ms
    150 asteroids: 6041ms
    200 asteroids: 10624ms
    500 asteroids: 65270ms
    1000 asteroids: 251897ms
    2000 asteroids: 1017985ms 
    3000 asteroids: 2295012ms (36 minutes LMAO)

    Simple parallelized simulation (5 threads)
    10 asteroids: 1732ms (27.5x slower)
    50 asteroids: 2087ms (2.65x slower)
    100 asteroids: 3047ms (6% slower)
    150 asteroids: 4213ms (43% faster)
    200 asteroids: 6194ms (71% faster)
    500 asteroids: 24611ms (2.49x faster)
    1000 asteroids: 90238ms (2.8x faster)
    2000 asteroids: 342844ms (3x faster)
    3000 asteroids: 778039ms (3x faster)
    
    Better parallelized simulation (5 threads)
    500 asteroids: 12772ms
    1000 asteroids: 57094ms
    2000 asteroids: 212045ms
    3000 asteroids: 480847ms



    Parameter optimisation: (on 500 asteroids)
    3 threads: 38119ms
    4 threads: 29924
    5 threads: 24359ms (<- actually the optimum, not at 6 threads, since we're opening 5 extra threads, and using the last)
    6 threads: 26222ms
    10 threads: 30261ms
    20 threads: 33278ms
    
*/


/*
    Proposed new algorithm:
    Split the simulation step in to 3 parts: filling out the force matrix, computing accelerations, updating step.
    Force matrix computation Example:
    N = 5
    total computations: N(N-1)/2 = 5*4/2 = 10 
    (we only concern ourselves with the lower triangle, since the matrix is anti symmetric)

       I
       0   1   2   3   4
    j  0  f01 f02 f03 f04
      f10  0  f12 f13 f14
      f20 f21  0  f23 f24
      f30 f31 f32  0  f34 
      f40 f41 f42 f43  0 
       -----------------
       a0  a1  a2  a3  a4

    acceleration
    a_i = sigma_0^4 f_ki (sum of columns)

    Splitting the work, round robin style.
    2 threads
    thread 0: handle 0, 2, 4
    thread 1: handle 1, 3
    equal workload.

    3 threads:
    thread 0: handle 0, 3
    thread 1: handle 1, 4
    thread 2: handle 2

    inside thread 1:
    loop int i=0, i < N, i+=thread_id
    skip 0.
    at 2: loop from 0 to 2-1: compute fi2 -> store at [0][2], [1][2]
    at 4: loop from 0 to 4-1: compute fi4 -> store at [0][4], [1][4], [2][4], [3][4]

*/

void better_force_aux(System& universe, 
                        const int thread_id) {
    // Better parallelized force computation, without doubled force computations.
    for (int i = thread_id; i < universe.bodies.size(); i += N_THREADS) {
        universe.force_matrix[i][i] = Vector(0.,0.);
        for (int j = 0; j < i ; j++) {
            Body& body1 = universe.bodies[i];
            Body& body2 = universe.bodies[j];
            Vector f = force(body1, body2);
            // universe.force_matrix[j][i] = f * (-1); // No need to store this
            universe.force_matrix[i][j] = f; 
        }
    }
}

void compute_accelerations_aux(System& universe, 
                        const int start, const int end) {
    for (int i =start; i< end; i++) {
        Vector acc(0, 0);
        // Add forces from upper triangle where i is row
        for (size_t j = 0; j < i; j++) {
            acc += universe.force_matrix[i][j];
        }
        // Add forces from lower triangle where i is column
        for (size_t j = i+1; j < universe.bodies.size(); j++) {
            acc -= universe.force_matrix[j][i];
        }
        universe.bodies[i].acceleration = acc/universe.bodies[i].m;
    }
}

void better_update_aux(System& universe,          // Modify
                        const int start, const int end, std::vector<Vector>& positions,
                        const double dt) {
    // Compute accelerations first
    // Normally, we don't have any race conditions here.
    // Then update everything
    for (int i = start; i < end; i++) {
        universe.bodies[i].update(dt);
        positions[i] = universe.bodies[i].coordinates;
    }
} 

void optimized_simulationmk2(System& universe) {
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

    // Initialize the force matrix;
    universe.force_matrix.resize(length);
    for(auto& row : universe.force_matrix) {
        row.resize(length, Vector(0,0));
    }

    // Loop through the steps here!
    for (int step=0; step < STEP_COUNT; ++step) {

        // Inside one step!
        std::vector<Vector> positions(length); // Pre initialize this this time, to conserve the order.

        // Split body calculation for all the threads, round robin style.
        for (int i = 0; i<N_THREADS-1; ++i) {
            workers[i] = std::thread(better_force_aux, std::ref(universe), i);
        }
        better_force_aux(universe, N_THREADS-1);
        for (int i=0; i < N_THREADS-1; ++i) {
            workers[i].join();
        }

        // By this point, universe.force_matrix should have all the values of forces stored. 
        // We just need to update accelerations

        int start_block = 0;
        for (int i = 0; i<N_THREADS-1; ++i) {
            int end_block = start_block + block_size;
            workers[i] = std::thread(compute_accelerations_aux, std::ref(universe),
                                    start_block, end_block);
            start_block = end_block;
        }
        compute_accelerations_aux(universe, start_block, length);
        for (int i=0; i<N_THREADS-1;++i) {
            workers[i].join();
        }

        start_block = 0;
        for (int i = 0; i<N_THREADS-1; ++i) {
            int end_block = start_block + block_size;
            workers[i] = std::thread(optimized_update_aux, std::ref(universe), 
                                    start_block, end_block, std::ref(positions), universe.dt);
            start_block = end_block;
        }
        optimized_update_aux(universe, start_block, length, positions, universe.dt);
        for (int i=0; i < N_THREADS-1; ++i) {
            workers[i].join();
        }

        universe.telemetry[step+1] = positions; // This doesn't cause race conditions.
    }
    // All steps done, Simulation done.

}