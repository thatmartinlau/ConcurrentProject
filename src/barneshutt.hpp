#ifndef BARNES_HUTT_HPP
#define BARNES_HUTT_HPP

#include <vector>
#include "core.hpp"

// Axis-aligned bounding box
struct Bounds {
    double minX, maxX;
    double minY, maxY;
};

// Quadtree node for Barnes–Hut
struct QuadNode {
    Bounds region;
    double totalMass;
    Vector centerOfMass;
    Body* singleBody;
    QuadNode* children[4];

    QuadNode(const Bounds& r);
    ~QuadNode();
};

// Compute simulation bounds
Bounds computeBounds(const System& universe);

// Build quadtree
QuadNode* createRootNode(const Bounds& bounds);
void insertBody(QuadNode* node, Body& body);

// Aggregate mass & center-of-mass
void computeMassDistribution(QuadNode* node);

// Compute force on a body via BH
Vector forceOnBody(const Body& bi, QuadNode* node, double theta);

// Parallel force computation
void computeForcesParallel(System& universe, QuadNode* root, double theta);
void computeForcesSerial(System& universe, QuadNode* root, double theta);

// Integrate positions & velocities
void updateBodies(System& universe, double dt);

// Cleanup the quadtree
void freeQuadTree(QuadNode* node);

// One step of Barnes–Hut
void BarnesHutStep(System& universe, double dt, double theta = 0.5, bool useParallel = false);

#endif // BARNES_HUTT_HPP
