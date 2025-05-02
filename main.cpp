/*
    CSE305 PROJECT: N BODY SIMULATION
    Authors: Martin LAU, Ziyue QIU, Oscar PEYRON
*/

#include <stdio.h>
#include <math.h>

#include "string"
#include "iostream"
#include "algorithm"
#include "vector"

#define STEP_COUNT 500
#define G 6.67430e-11


class Vector {
    public:
        explicit Vector(double x = 0, double y=0) {
            data[0] = x;
            data[1] = y;
        }
        double data[2];

        Vector operator+(const Vector& other) const {
            return Vector(data[0] + other.data[0], data[1] + other.data[1]);
        }
        Vector operator*(double scalar) const {
            return Vector(data[0] * scalar, data[1] * scalar);
        }
        void operator+=(const Vector& other) {
            data[0] += other.data[0];
            data[1] += other.data[1];
        }
};

class Body {
    public:
    double m;
    Vector coordinates;
    Vector velocity;
    Vector acceleration;
    Body() : m(0), coordinates(), velocity(), acceleration() {}
    
    Body(double mass, const Vector& pos, const Vector& vel) 
        : m(mass), coordinates(pos), velocity(vel), acceleration() {}

    void update(double dt) {
        velocity += acceleration * dt;
        coordinates += velocity * dt;
        acceleration = Vector(0, 0);
    }
};

class System {
    public:
    std::vector<Body> bodies;
    std::vector<std::vector<Vector>> telemetry; 

    void add(Body body) {
        bodies.push_back(body);
    }

    void visualize(const std::string& name) {
        // Use telemetry here
        // Make the visualization based off the telemetry table of values.
    }
};

double sqr(double x) {
    return x*x;
}

Vector force(Body bi, Body bj) {
    /*
        Compute the gravitational force exerted by body B_i onto B_j
    */
    double x1 = bi.coordinates.data[0];
    double y1 = bi.coordinates.data[1];
    double x2 = bj.coordinates.data[0];
    double y2 = bj.coordinates.data[1]; 

    double dx = x2-x1;
    double dy = y2-y1;
    double magnitude = G * (bi.m * bj.m) / (sqr(dx) + sqr(dy));
    double dist = sqrt(sqr(dx) + sqr(dy));
    return Vector(magnitude*dx/dist, magnitude*dy/dist);
}

void naive_simulation(System &universe, double dt) {
    /*
        A straightforward, naive implementation of 
    */
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

void optimized_simulation(System &universe, double dt) {
    for (int i=0; i<STEP_COUNT; ++i) {
        // Build upon naive simulation, but make it better
        // Parallelization
        // Avoiding computing the forces twice
        // Parallelizing using atomic variables
        // Avoid mutexes and atomic variables as much as possible by dividing bodies in groups to be treated by different threads properly
    }
}

void BarnesHut(System &universe, double dt) {
    // Fill in later.
}



int main(){
    System universe;
    // Sun, Earth, Mars system (simplified 2D projection)
    // Note: Initial positions and velocities are simplified to 2D and approximate

    // These were provided by Claude
    // Sun at approximate center
    Body sun(1.989e30,  // Mass of Sun in kg
        Vector(0, 0),  // Position at origin
        Vector(0, 0)); // Nearly stationary

    // Earth
    // Average distance from Sun: 149.6e9 meters (1 AU)
    // Orbital velocity: ~29.78 km/s
    Body earth(5.972e24,  // Mass of Earth in kg
        Vector(149.6e9, 0),  // Position at 1 AU on x-axis
        Vector(0, 29780));   // Orbital velocity primarily in y direction

    // Mars
    // Average distance from Sun: 227.9e9 meters (1.52 AU)
    // Orbital velocity: ~24.077 km/s
    Body mars(6.39e23,    // Mass of Mars in kg
        Vector(227.9e9, 0),  // Position at 1.52 AU on x-axis
        Vector(0, 24077));   // Orbital velocity primarily in y direction

    universe.add(sun);
    universe.add(earth);
    universe.add(mars);

    naive_simulation(universe, 0.001);
    universe.visualize("naive_simulation.gif");
    // optimized_simulation(universe, 0.001);
    // universe.visualize("optimized_simulation.gif");
    // BarnesHut(universe, 0.001);
    
    return 0;
}