#ifndef BARNES_HUTT_HPP
#define BARNES_HUTT_HPP

#include <vector>
#include <utility>
#include "core.hpp"    // Body, System, Vector

// Spatial bounds for quadtree region
struct Bounds {
    double minX, maxX;
    double minY, maxY;
};

// Quadtree node for Barnes–Hut algorithm
struct QuadNode {
    Bounds region;               // spatial extent of this node
    double totalMass;            // total mass in this node
    Vector centerOfMass;         // center-of-mass of bodies in this node
    Body* singleBody;            // non-nullptr if leaf containing exactly one body
    QuadNode* children[4];       // NW=0, NE=1, SW=2, SE=3

    QuadNode(const Bounds& r);
    ~QuadNode();
};

// Compute the axis-aligned bounding box for all bodies
std::pair<Bounds, Bounds> computeBounds(const System& universe);

// Quadtree construction
QuadNode* createRootNode(const Bounds& bounds);
void insertBody(QuadNode* node, Body& body);

// Compute mass distribution (post-order traversal)
void computeMassDistribution(QuadNode* node);

// Force calculation using Barnes–Hut opening angle theta
Vector forceOnBody(const Body& bi, QuadNode* node, double theta);

// Parallel force computation (shared-memory)
void computeForcesParallel(System& universe, QuadNode* root, double theta);

// Time integration (velocity & position updates)
void updateBodies(System& universe, double dt);

// Cleanup quadtree memory
void freeQuadTree(QuadNode* node);

// GPU brute-force PRAM illustration
#ifdef USE_CUDA
void simulateBruteForceGPU(System& universe, double dt);
#endif

#endif // BARNES_HUTT_HPP
