#include "core.hpp"
#include "barneshutt.hpp"

#ifdef G
#undef G
#endif

#include <vector>
#include <cmath>
#include <cstddef>
#include <omp.h>

static constexpr double BH_G = 6.67430e-11;

// Node pool stores all QuadNode objects for one step, so we don’t keep new/delete calls.
static std::vector<QuadNode> nodePool;
static std::size_t poolIndex = 0;

// QuadNode just sets its region; other fields will be filled in later.
QuadNode::QuadNode(const Bounds& b)
    : region(b), totalMass(0.0), centerOfMass(0.0, 0.0), singleBody(nullptr) {
    children[0] = children[1] = children[2] = children[3] = nullptr;
}

// Grab a node from the pool (reuse if possible, or add a new one).
QuadNode* allocateNode(const Bounds& region) {
    QuadNode* node;
    if (poolIndex < nodePool.size()) {
        node = &nodePool[poolIndex];
    } else {
        nodePool.emplace_back(region);
        node = &nodePool.back();
    }
    // Reset fields
    node->region = region;
    node->totalMass = 0.0;
    node->centerOfMass = Vector(0.0, 0.0);
    node->singleBody = nullptr;
    for (int i = 0; i < 4; ++i) {
        node->children[i] = nullptr;
    }
    poolIndex++;
    return node;
}

// Find the smallest square that contains all bodies (centered).
Bounds computeBounds(const std::vector<Body>& bodies) {
    Bounds b;
    if (bodies.empty()) {
        // If no bodies, just return a default square
        b.minX = b.minY = -1.0;
        b.maxX = b.maxY =  1.0;
        return b;
    }
    // Start with the first body’s position
    b.minX = b.maxX = bodies[0].coordinates.data[0];
    b.minY = b.maxY = bodies[0].coordinates.data[1];
    for (const auto& c : bodies) {
        double x = c.coordinates.data[0];
        double y = c.coordinates.data[1];
        if (x < b.minX) b.minX = x;
        if (x > b.maxX) b.maxX = x;
        if (y < b.minY) b.minY = y;
        if (y > b.maxY) b.maxY = y;
    }
    // Make it a square, keeping the same center
    double dx = b.maxX - b.minX;
    double dy = b.maxY - b.minY;
    double d = (dx > dy ? dx : dy);
    double cx = 0.5 * (b.minX + b.maxX);
    double cy = 0.5 * (b.minY + b.maxY);
    b.minX = cx - 0.5 * d;
    b.maxX = cx + 0.5 * d;
    b.minY = cy - 0.5 * d;
    b.maxY = cy + 0.5 * d;
    return b;
}

// Figure out which quadrant a point belongs to:
//  0 = NE, 1 = NW, 2 = SW, 3 = SE
int getQuadrant(const Bounds& region, const Vector& pos) {
    double midX = 0.5 * (region.minX + region.maxX);
    double midY = 0.5 * (region.minY + region.maxY);
    bool right = (pos.data[0] > midX);
    bool top   = (pos.data[1] > midY);
    if (right && top)        return 0; // NE
    if (!right && top)       return 1; // NW
    if (!right && !top)      return 2; // SW
    return 3;                        // SE
}

// Given a parent region, get the bounds for its child quadrant
Bounds childBounds(const Bounds& b, int quad) {
    double midX = 0.5 * (b.minX + b.maxX);
    double midY = 0.5 * (b.minY + b.maxY);
    switch (quad) {
        case 0: return Bounds{midX, midY, b.maxX, b.maxY}; // NE
        case 1: return Bounds{b.minX, midY, midX, b.maxY}; // NW
        case 2: return Bounds{b.minX, b.minY, midX, midY}; // SW
        case 3: return Bounds{midX, b.minY, b.maxX, midY}; // SE
    }
    return b; // should never reach here
}

// Insert a Body into the quadtree node.
// If the node is empty leaf, just store the body there.
// If it already has one, split into 4 children and re-insert.
void insertBody(QuadNode* node, Body& body) {
    if (node->singleBody == nullptr && node->children[0] == nullptr) {
        node->singleBody = &body;
        return;
    }
    if (node->children[0] == nullptr) {
        // Need to split this leaf because it already has one body
        Bounds b = node->region;
        Body* existing = node->singleBody;
        node->singleBody = nullptr;
        for (int i = 0; i < 4; ++i) {
            node->children[i] = allocateNode(childBounds(b, i));
        }
        int quadOld = getQuadrant(b, existing->coordinates);
        node->children[quadOld]->singleBody = existing;
    }
    // Not an empty leaf anymore, so recurse into the right child
    int quadNew = getQuadrant(node->region, body.coordinates);
    insertBody(node->children[quadNew], body);
}

// Compute totalMass and centerOfMass from the leaves up:
// If it’s a leaf with a body, just use that body.
// Otherwise combine children’s masses and centers.
void computeMassDistribution(QuadNode* node) {
    if (!node) return;
    if (node->children[0] == nullptr) {
        if (node->singleBody) {
            node->totalMass = node->singleBody->m;
            node->centerOfMass = node->singleBody->coordinates;
        }
        return;
    }
    // Internal node: first do children
    double msum = 0.0;
    Vector weightedSum(0.0, 0.0);
    for (int i = 0; i < 4; ++i) {
        computeMassDistribution(node->children[i]);
        if (node->children[i] && node->children[i]->totalMass > 0) {
            msum += node->children[i]->totalMass;
            weightedSum.data[0] += node->children[i]->centerOfMass.data[0] * node->children[i]->totalMass;
            weightedSum.data[1] += node->children[i]->centerOfMass.data[1] * node->children[i]->totalMass;
        }
    }
    node->totalMass = msum;
    if (msum > 0.0) {
        node->centerOfMass.data[0] = weightedSum.data[0] / msum;
        node->centerOfMass.data[1] = weightedSum.data[1] / msum;
    }
}

Vector computeForceIterative(const Body& bi, QuadNode* root, double theta) {
    Vector total(0.0, 0.0);
    if (!root || root->totalMass == 0.0) return total;
    std::vector<QuadNode*> stack;
    stack.reserve(64);
    stack.push_back(root);
    while (!stack.empty()) {
        QuadNode* node = stack.back();
        stack.pop_back();
        if (!node || node->totalMass == 0.0) continue;
        if (node->singleBody == &bi && node->children[0] == nullptr) continue;
        double dx = node->centerOfMass.data[0] - bi.coordinates.data[0];
        double dy = node->centerOfMass.data[1] - bi.coordinates.data[1];
        double r2 = dx*dx + dy*dy + 1e-12; // avoid divide-by-zero
        double r  = std::sqrt(r2);
        double size = node->region.maxX - node->region.minX;
        if (node->children[0] == nullptr || (size / r) < theta) {
            double forceMag = BH_G * bi.m * node->totalMass / r2;
            total.data[0] += dx / r * forceMag;
            total.data[1] += dy / r * forceMag;
        } else {
            for (int i = 0; i < 4; ++i) {
                if (node->children[i]) {
                    stack.push_back(node->children[i]);
                }
            }
        }
    }
    return total;
}

void BarnesHutStep(std::vector<Body>& bodies, double dt, double theta, bool useParallel) {
    // 1) Set up the pool
    poolIndex = 0;
    nodePool.clear();
    nodePool.reserve(bodies.size() * 4 + 1); // roughly 4N nodes

    // 2) Find bounding box
    Bounds b = computeBounds(bodies);

    // 3) Take one node from the pool as root
    QuadNode* root = allocateNode(b);

    // 4) Insert every body into the tree
    for (auto& c : bodies) {
        insertBody(root, c);
    }

    // 5) Compute masses and centers
    computeMassDistribution(root);

    int n = static_cast<int>(bodies.size());

    // 6) Compute accelerations (parallel if big enough)
    if (useParallel) { // remove the restriction of bodies > 2000
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i) {
            Vector f = computeForceIterative(bodies[i], root, theta);
            bodies[i].acceleration.data[0] = f.data[0] / bodies[i].m;
            bodies[i].acceleration.data[1] = f.data[1] / bodies[i].m;
        }
    } else {
        for (int i = 0; i < n; ++i) {
            Vector f = computeForceIterative(bodies[i], root, theta);
            bodies[i].acceleration.data[0] = f.data[0] / bodies[i].m;
            bodies[i].acceleration.data[1] = f.data[1] / bodies[i].m;
        }
    }

    // 7) Update velocities & positions
    for (auto& c : bodies) {
        c.velocity.data[0] += c.acceleration.data[0] * dt;
        c.velocity.data[1] += c.acceleration.data[1] * dt;
        c.coordinates.data[0] += c.velocity.data[0] * dt;
        c.coordinates.data[1] += c.velocity.data[1] * dt;
    }
}


// Additional function for test theta (this one is for computing the error between barneshut and naive_simulation)
void naive_with_record(System &universe,
                       std::vector<std::vector<Vector>> &groundPos) {
    std::cout << "Running naive baseline (recording every step)...\n";

    universe.telemetry.clear();

    int N = static_cast<int>(universe.bodies.size());

    for (int i = 0; i < N; ++i) {
        groundPos[0][i] = universe.bodies[i].coordinates;
    }

    for (int step = 0; step < STEP_COUNT; ++step) {
        for (auto &body : universe.bodies) {
            body.acceleration = Vector(0, 0);
        }
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                Body &b1 = universe.bodies[i];
                Body &b2 = universe.bodies[j];
                Vector f = force(b1, b2);
                b1.acceleration += f / b1.m;
                b2.acceleration += f / (-b2.m);
            }
        }
        
        for (int i = 0; i < N; ++i) {
            universe.bodies[i].update(universe.dt);
            groundPos[step+1][i] = universe.bodies[i].coordinates;
        }
    }
