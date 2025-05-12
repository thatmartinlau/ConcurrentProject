#include "core.hpp"
#include <cmath>
#include <fstream>
#include <Magick++.h>
using namespace Magick;

// Vector implementations
Vector::Vector(double x, double y) {
    data[0] = x;
    data[1] = y;
}

Vector Vector::operator+(const Vector& other) const {
    return Vector(data[0] + other.data[0], data[1] + other.data[1]);
}

Vector Vector::operator*(double scalar) const {
    return Vector(data[0] * scalar, data[1] * scalar);
}

void Vector::operator+=(const Vector& other) {
    data[0] += other.data[0];
    data[1] += other.data[1];
}

// Body implementations
Body::Body() : m(0), coordinates(), velocity(), acceleration() {}

Body::Body(double mass, const Vector& pos, const Vector& vel) 
    : m(mass), coordinates(pos), velocity(vel), acceleration() {}

void Body::update(double dt) {
    velocity += acceleration * dt;
    coordinates += velocity * dt;
    acceleration = Vector(0, 0);
}

// System implementations
void System::add(Body body) {
    bodies.push_back(body);
}

// Utility function implementations
double sqr(double x) {
    return x * x;
}

Vector force(Body bi, Body bj) {
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


std::pair<Vector, Vector> System::getBounds() const {
    if (telemetry.empty()) return {Vector(0,0), Vector(1,1)};
    
    Vector min_pos(INFINITY, INFINITY);
    Vector max_pos(-INFINITY, -INFINITY);
    
    for (const auto& frame : telemetry) {
        for (const auto& pos : frame) {
            min_pos.data[0] = std::min(min_pos.data[0], pos.data[0]);
            min_pos.data[1] = std::min(min_pos.data[1], pos.data[1]);
            max_pos.data[0] = std::max(max_pos.data[0], pos.data[0]);
            max_pos.data[1] = std::max(max_pos.data[1], pos.data[1]);
        }
    }
    
    // Add 10% padding
    Vector padding = (max_pos + min_pos * -1) * 0.1;
    return {min_pos + padding * -1, max_pos + padding};
}

void System::visualize(const std::string& name) {
    InitializeMagick(nullptr);
    
    const int width = 800;
    const int height = 600;
    std::vector<Image> frames;
    
    // Find bounds
    auto [min_pos, max_pos] = getBounds();
    double min_x = min_pos.data[0];
    double min_y = min_pos.data[1];
    double max_x = max_pos.data[0];
    double max_y = max_pos.data[1];
    
    // Create frames with memory management
    size_t batch_size = 50; // Process frames in batches
    for (size_t i = 0; i < telemetry.size(); i++) {
        Image image(Geometry(width, height), Color("white"));
        
        // Draw points
        for (const auto& pos : telemetry[i]) {
            int x = static_cast<int>((pos.data[0] - min_x) / (max_x - min_x) * (width - 20) + 10);
            int y = static_cast<int>(height - ((pos.data[1] - min_y) / (max_y - min_y) * (height - 20) + 10));
            
            // Draw a larger point
            image.fillColor("blue");  // Set fill color
            image.draw(DrawableCircle(x, y, x + 3, y + 3));  // Draw a small circle instead of a point
        }
        
        image.animationDelay(10);  // Slow down animation slightly
        frames.push_back(std::move(image));  // Use move semantics
        
        // Write in batches to manage memory
        if (frames.size() >= batch_size || i == telemetry.size() - 1) {
            writeImages(frames.begin(), frames.end(), name);
            frames.clear();
        }
    }
    
    // Resource cleanup
    frames.clear();

}
