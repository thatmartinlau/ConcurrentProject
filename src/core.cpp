#include "core.hpp"
#include <cmath>
#include <fstream>

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

std::pair<int, int> System::worldToScreen(const Vector& pos, const Vector& min, const Vector& max, 
                                        int width, int height) const {
    double x_scale = width / (max.data[0] - min.data[0]);
    double y_scale = height / (max.data[1] - min.data[1]);
    
    int screen_x = static_cast<int>((pos.data[0] - min.data[0]) * x_scale);
    int screen_y = height - static_cast<int>((pos.data[1] - min.data[1]) * y_scale);
    
    return {screen_x, screen_y};
}

void System::writePPM(const std::string& filename, const std::vector<std::vector<uint8_t>>& image, 
                     int width, int height) {
    std::ofstream file(filename);
    file << "P3\n" << width << " " << height << "\n255\n";
    
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            uint8_t value = image[y][x];
            // Convert grayscale to RGB
            file << (int)value << " " << (int)value << " " << (int)value << " ";
        }
        file << "\n";
    }
}

void System::visualize(const std::string& name) {
    const int width = 800;
    const int height = 600;
    auto [min_pos, max_pos] = getBounds();
    
    // Create frames
    for (size_t frame = 0; frame < telemetry.size(); ++frame) {
        // Create a new white image
        std::vector<std::vector<uint8_t>> image(height, std::vector<uint8_t>(width, 255));
        
        // Draw trajectory trails (previous positions)
        for (size_t prev_frame = 0; prev_frame < frame; ++prev_frame) {
            for (const auto& pos : telemetry[prev_frame]) {
                auto [x, y] = worldToScreen(pos, min_pos, max_pos, width, height);
                if (x >= 0 && x < width && y >= 0 && y < height) {
                    image[y][x] = 200;  // Light gray for trails
                }
            }
        }
        
        // Draw current positions
        for (const auto& pos : telemetry[frame]) {
            auto [x, y] = worldToScreen(pos, min_pos, max_pos, width, height);
            if (x >= 0 && x < width && y >= 0 && y < height) {
                // Draw a larger point (3x3 pixels)
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dx = -1; dx <= 1; ++dx) {
                        int px = x + dx;
                        int py = y + dy;
                        if (px >= 0 && px < width && py >= 0 && py < height) {
                            image[py][px] = 0;  // Black for current position
                        }
                    }
                }
            }
        }
        
        // Save frame
        std::string frame_filename = name + "_frame_" + 
                                   std::to_string(frame) + ".ppm";
        writePPM(frame_filename, image, width, height);
    }
    
    // Create a simple shell script to view frames in sequence
    std::ofstream script("view_frames.sh");
    script << "#!/bin/bash\n";
    script << "for f in " << name << "_frame_*.ppm; do\n";
    script << "    cat $f\n";
    script << "    sleep 0.05\n";
    script << "done\n";
    script.close();
    
    // Make the script executable
    system("chmod +x view_frames.sh");
}
