#pragma once

#include "core.hpp"

struct Bounds {
    double minX, minY, maxX, maxY;
};

struct QuadNode {
    Bounds region;
    double totalMass;
    Vector centerOfMass;
    Body* singleBody;
    QuadNode* children[4];

    QuadNode(const Bounds& b);
};

void BarnesHutStep(std::vector<Body>& bodies,
                   double dt,
                   double theta,
                   bool useParallel);

void naive_with_record(System &universe,
                       std::vector<std::vector<Vector>> &groundPos);
