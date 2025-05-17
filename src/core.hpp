#ifndef CORE_H
#define CORE_H

#include <vector>
#include <string>
#include <Magick++.h>
using std::string; 

#define STEP_COUNT 500
#define G 6.67430e-11

// Vector class
class Vector {
    public:
        explicit Vector(double x = 0, double y = 0);
        double data[2];

        Vector operator+(const Vector& other) const;
        Vector operator*(double scalar) const;
        void operator+=(const Vector& other);
        double& operator[](int i) { return data[i]; }
        const double& operator[](int i) const { return data[i]; }
};

// Body class
class Body {
    public:
        double m;
        Vector coordinates;
        Vector velocity;
        Vector acceleration;

        // Aesthetic, for the visualizations
        string color;
        int size;
        string title;

        
        Body();
        Body(double mass, const Vector& pos, const Vector& vel, string color, int size, string title);
        void update(double dt);
};

// System class
class System {
    public:
        std::vector<Body> bodies;
        std::vector<std::vector<Vector>> telemetry;

        void add(Body body);
        void visualize(const std::string& name, bool time, bool axes);
    private:
        std::pair<Vector, Vector> getBounds() const;
        std::pair<int, int> worldToScreen(const Vector& pos, const Vector& min, const Vector& max, int width, int height) const;

};

// Utility functions
double sqr(double x);
Vector force(Body bi, Body bj);

#endif