#include "core.hpp"
#include <Magick++.h>

// #define STB_IMAGE_WRITE_IMPLEMENTATION
// #include "stb_image_write.h"

#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <cstdlib>
#include <atomic>
//#include <omp.h>
#include <thread>
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

void Vector::operator-=(const Vector& other) {
    data[0] -= other.data[0];
    data[1] -= other.data[1];
}

double Vector::norm() {
    return sqrt(data[0]*data[0] + data[1] + data[1]);
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
    double dist = sqrt(sqr(dx) + sqr(dy)); 

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
    const int width = 600;
    const int height = 400;
    const int padding = 25;
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
    for (size_t i = 0; i < telemetry.size(); i+=50) {
        // Progress bar; this is mainly for debug purposes
        std::cout << i << " " << std::flush;

        Image image(Geometry(width, height), Color("black"));
        if (axes) {
            // Draw coordinate axes
            image.strokeColor("gray");
            image.draw(DrawableLine(padding, height-padding, width-padding, height-padding)); // X axis
            image.draw(DrawableLine(padding, height-padding, padding, padding)); // Y axis
         }
         // For each body in the current frame, we draw its position on the image.
        for (size_t body_idx = 0; body_idx < telemetry[i].size(); body_idx++) {
            // Get color for this body
            Color bodyColor = colors[body_idx % colors.size()];
            //Color trailColor = bodyColor; //does not work with me
            ColorRGB trailColor(bodyColor); //For Oscar when running on his local computer
            trailColor.alpha(65535 * 0.3); // 30% opacity for trails
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
    auto start = std::chrono::high_resolution_clock::now();
    writeImages(frames.begin(), frames.end(), name);
    auto end = std::chrono::high_resolution_clock::now();
    auto time_taken = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Image Writing time: " << time_taken.count() << " milliseconds.";
    frames.clear();
    std::cout << "File written.\n";

///end marker of the vizualization function
}


// In this function, we allow ourselves to use pragma omp for the visualization purposes
// The visualization is not the most important part for multi threading, and we would just
// like something that works simply and easily for us. The other algorithm implementations
// will use proper multithreading, as seen in class.


void System::visualize2(const std::string& name, bool time=true, bool axes=true) {
    InitializeMagick(nullptr);
    const int width = 600;
    const int height = 600;
    const int padding = 25;

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
    
    //auto [min_pos, max_pos] = getBounds();
    auto bounds = getBounds();
    auto min_pos = bounds.first;
    auto max_pos = bounds.second;

    // Create directory if it doesn't exist
    std::string dir_name = name + "_frames";
    int a = system(("mkdir -p " + dir_name).c_str());

    // Calculate total frames
    size_t total_frames = FRAME_NUM;
    const size_t step_size = STEP_COUNT / FRAME_NUM;
    std::atomic<size_t> progress{0};

    // std::cout<< telemetry.size() << step_size << "\n";

    #pragma omp parallel
    {
        std::ostringstream progress_stream;

        #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < telemetry.size(); i+=step_size) {
            // Progress tracking
            size_t frame_num = i/step_size;
            progress_stream.str("");
            progress_stream << "Processing frame " << frame_num << "/" << total_frames << "\r";
            #pragma omp critical
            {
                std::cout << progress_stream.str() << std::flush;
            }

            Image image(Geometry(width, height), Color("black"));

            // Compression settings
            image.quantize(256);
            image.quality(50);
            image.compressType(LZWCompression);

            if (axes) {
                #pragma omp critical
                {
                    // Draw coordinate axes
                    image.strokeColor("gray");
                    image.draw(DrawableLine(padding, height-padding, width-padding, height-padding)); // X axis
                    image.draw(DrawableLine(padding, height-padding, padding, padding)); // Y axis
                }
            }
            
            // Draw each body
            for (size_t body_idx = 0; body_idx < telemetry[i].size(); body_idx++) {
                Color bodyColor = Color(bodies[body_idx].color);
                
                // Draw current position for this body
                const auto& pos = telemetry[i][body_idx];
                int x = static_cast<int>((pos.data[0] - min_pos.data[0]) / (max_pos.data[0] - min_pos.data[0]) 
                        * (width - 2*padding) + padding);
                int y = static_cast<int>(height - ((pos.data[1] - min_pos.data[1]) / (max_pos.data[1] - min_pos.data[1]) 
                        * (height - 2*padding) + padding));
                
                #pragma omp critical
                {
                    image.fillColor(bodyColor);
                    int radius = bodies[body_idx].size;
                    image.draw(DrawableCircle(x, y, x + radius, y + radius));

                    // Draw title
                    image.fillColor("white");
                    image.draw(DrawableText(x + 10, y - 10, bodies[body_idx].title));
                }
            }
            
            // Add time information
            if (time) {
                double current_time = i * dt;
                std::string timeInfo = "Time: " + std::to_string(current_time) + "s";
                #pragma omp critical
                {
                    image.fillColor("white");
                    image.draw(DrawableText(padding, padding-20, timeInfo));
                }
            }

            // Save individual frame
            std::ostringstream frame_name;
            frame_name << dir_name << "/frame_" << std::setw(6) << std::setfill('0') << frame_num << ".png";
            #pragma omp critical
            {
                image.write(frame_name.str());
            }
        }
    }

    std::cout << "\nAll frames generated.\n" << std::flush;
    
    // Create video using ffmpeg
    std::string ffmpeg_command = "ffmpeg -framerate 30 -i " + dir_name + "/frame_%06d.png -c:v libx264 -pix_fmt yuv420p " + name + ".mp4 2>/dev/null";
    int b = system(ffmpeg_command.c_str());

    // Clean up frames directory
    int c = system(("rm -rf " + dir_name).c_str());

    std::cout << "Video generated: " << name << ".mp4\n" << std::flush;
}
