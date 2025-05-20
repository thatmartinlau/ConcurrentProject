#include "core.hpp"
#include <cmath>
#include <fstream>
#include <Magick++.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
using std::string;
using namespace Magick;

// Vector implementations
Vector::Vector(double x, double y) {
    data[0] = x;
    data[1] = y;
}

    

Vector Vector::operator+(const Vector& other) const {
    return Vector(data[0] + other.data[0], data[1] + other.data[1]);
}

Vector Vector::operator-(const Vector& other) const {
    return Vector(data[0] - other.data[0], data[1] - other.data[1]);
}

Vector Vector::operator*(double scalar) const {
    return Vector(data[0] * scalar, data[1] * scalar);
}

Vector Vector::operator/(double scalar) const {
    return Vector(data[0]/scalar, data[1]/scalar);
}

Vector& Vector::operator=(const Vector& other) {
    data[0] = other.data[0];
    data[1] = other.data[1];
    return *this;
}


void Vector::operator+=(const Vector& other) {
    data[0] += other.data[0];
    data[1] += other.data[1];
}

// Body implementations
Body::Body() : m(0), coordinates(), velocity(), acceleration(), color(""), size(0), title("") {}

Body::Body(double mass, const Vector& pos, const Vector& vel, string color, int size, string title) 
    : m(mass), coordinates(pos), velocity(vel), acceleration(), color(color), size(size), title(title) {}

void Body::update(double dt) {
    // Update position first
    coordinates.data[0] += velocity.data[0] * dt + 0.5 * acceleration.data[0] * dt * dt;
    coordinates.data[1] += velocity.data[1] * dt + 0.5 * acceleration.data[1] * dt * dt;
    
    // Then update velocity
    velocity.data[0] += acceleration.data[0] * dt;
    velocity.data[1] += acceleration.data[1] * dt;

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

    double dx = x2 - x1;
    double dy = y2 - y1;
    double dist = sqrt(sqr(dx) + sqr(dy));  // Calculate distance once
    
    // Avoid division by zero
    if (dist < 1e-10) {
        return Vector(0, 0);
    }
    
    // F = G * m1 * m2 / r^2
    double magnitude = G * (bi.m * bj.m) / sqr(dist);  // Use dist² not (dx² + dy²)
    
    // Convert to vector components using direction cosines
    return Vector(magnitude * dx/dist, magnitude * dy/dist);
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

std::pair<Vector, Vector> System::exposeBounds() const{
    return getBounds();
} // to be able to get the bound

void System::visualize(const std::string& name, bool time=true, bool axes=true) {
    InitializeMagick(nullptr);
    const int width = 800;
    const int height = 600;
    const int padding = 50;
    std::vector<Image> frames;

    // Define different colors for different bodies
    std::vector<Color> colors = {
        Color("yellow"),
        Color("blue"),
        Color("red"),
        Color("green"),
        Color("purple"),
        Color("orange")
    };

    // Find bounds
    auto [min_pos, max_pos] = getBounds();

    // Create frames for each time step
    for (size_t i = 0; i < telemetry.size(); i++) {
        // Progress bar
        std::cout << i << " " << std::flush;

        Image image(Geometry(width, height), Color("black"));

        //just to debug
        std::cout << "telemetry.size(): " << telemetry.size() << ", i: " << i << std::endl;

        if (axes) {
            // Draw coordinate axes

            image.strokeColor("gray");
            image.draw(DrawableLine(padding, height-padding, width-padding, height-padding)); // X axis
            image.draw(DrawableLine(padding, height-padding, padding, padding)); // Y axis
            
            // Draw scale markers
            image.fillColor("white");
            for(int j = 0; j <= 10; j++) {
                int x = padding + j * (width - 2*padding) / 10;
                int y = height - padding + 20;
                std::string label = std::to_string(static_cast<int>(min_pos.data[0] + j * (max_pos.data[0] - min_pos.data[0]) / 10));
                image.draw(DrawableText(x, y, label));
            }
        }
        
        
        // For each body in the current frame
        for (size_t body_idx = 0; body_idx < telemetry[i].size(); body_idx++) {
            // Get color for this body
            Color bodyColor = colors[body_idx % colors.size()];
            //Color trailColor = bodyColor; //does not work with me
            ColorRGB trailColor(bodyColor); //For Oscar when running on his local computer
            trailColor.alpha(65535 * 0.3); // 30% opacity for trails
            
            // Draw trail for this body
            if (i > 0) {
                for (size_t j = 0; j < i; j++) {
                    const auto& trail_pos = telemetry[j][body_idx];
                    int trail_x = static_cast<int>((trail_pos.data[0] - min_pos.data[0]) / (max_pos.data[0] - min_pos.data[0]) 
                            * (width - 2*padding) + padding);
                    int trail_y = static_cast<int>(height - ((trail_pos.data[1] - min_pos.data[1]) / (max_pos.data[1] - min_pos.data[1]) 
                            * (height - 2*padding) + padding));
                    
                    image.fillColor(bodyColor);
                    image.draw(DrawableCircle(trail_x, trail_y, trail_x + 2, trail_y + 2));
                }
            }
            
            // Draw current position for this body
            const auto& pos = telemetry[i][body_idx];
            int x = static_cast<int>((pos.data[0] - min_pos.data[0]) / (max_pos.data[0] - min_pos.data[0]) 
                    * (width - 2*padding) + padding);
            int y = static_cast<int>(height - ((pos.data[1] - min_pos.data[1]) / (max_pos.data[1] - min_pos.data[1]) 
                    * (height - 2*padding) + padding));
            
            image.fillColor(bodyColor);
            image.draw(DrawableCircle(x, y, x + 5, y + 5));

            // Draw title for this body
                image.fillColor("white");
                // Position the title slightly above and to the right of the body
                image.draw(DrawableText(x + 10, y - 10, bodies[body_idx].title));
        }
        
        // Add time information
        if (time) {
            double current_time = i * dt;
            std::string timeInfo = "Time: " + std::to_string(current_time) + "s";
            image.fillColor("white");
            image.draw(DrawableText(padding, padding-20, timeInfo));
        }
        image.animationDelay(5);
        frames.push_back(std::move(image));
        
    // Image frame creation end
    }
    std::cout << "\nImages generated. Now writing to file...\n";

    // Write all frames at once
    writeImages(frames.begin(), frames.end(), name);
    frames.clear();
    std::cout << "File written.\n";
}