#include "barneshutt.hpp"
#include <cmath>
#include <omp.h>

// QuadNode constructor: set bounds and init members
QuadNode::QuadNode(const Bounds& r)
    : region(r), totalMass(0), centerOfMass(0,0), singleBody(nullptr) {
    // children are set to nullptr
    for(int i=0; i<4; ++i) children[i] = nullptr;
}
// QuadNode destructor: nothing special
QuadNode::~QuadNode() {}

// 1. Compute bounds from system telemetry (just calls exposeBounds)
Bounds computeBounds(const System& universe) {
    // get min and max positions
    auto [minV, maxV] = universe.exposeBounds();
    // return as Bounds struct
    return { minV.data[0], maxV.data[0], minV.data[1], maxV.data[1] };
}

// Helper: decide which quadrant a vector v belongs to inside bounds b
static int quadrant(const Bounds& b, const Vector& v) {
    double mx = 0.5 * (b.minX + b.maxX); // mid x
    double my = 0.5 * (b.minY + b.maxY); // mid y
    bool east = v.data[0] > mx;   // is to the right?
    bool north = v.data[1] > my;  // is above?
    if(north) return east ? 1 : 0; // NE : NW
    else      return east ? 3 : 2; // SE : SW
}
// Helper: get bounds for a child quadrant index q
static Bounds childBounds(const Bounds& b, int q) {
    double mx = 0.5 * (b.minX + b.maxX);
    double my = 0.5 * (b.minY + b.maxY);
    switch(q) {
      case 0: return { b.minX, mx,    my,    b.maxY }; // NW
      case 1: return { mx,    b.maxX, my,    b.maxY }; // NE
      case 2: return { b.minX, mx,    b.minY, my   }; // SW
      case 3: return { mx,    b.maxX, b.minY, my   }; // SE
      default: return b; // fallback
    }
}

// 2a. Create a new root node for the quadtree
QuadNode* createRootNode(const Bounds& bounds) {
    return new QuadNode(bounds);
}
// 2b. Insert a body into the quadtree rooted at node
void insertBody(QuadNode* node, Body& body) {
    // if node is empty leaf (no body and no children), place body here
    if(!node->singleBody && !node->children[0]) {
        node->singleBody = &body;
        return;
    }
    // if leaf has one body, subdivide into 4 children
    if(node->singleBody) {
        Body* old = node->singleBody; // existing body
        node->singleBody = nullptr;   // clear this leaf
        // create children nodes
        for(int i=0; i<4; ++i)
            node->children[i] = new QuadNode(childBounds(node->region, i));
        // re-insert old body into appropriate child
        int qOld = quadrant(node->region, old->coordinates);
        insertBody(node->children[qOld], *old);
    }
    // now insert the new body into correct child
    int q = quadrant(node->region, body.coordinates);
    // if child does not exist, create it
    if(!node->children[q]) 
        node->children[q] = new QuadNode(childBounds(node->region, q));
    insertBody(node->children[q], body);
}

// 3. Compute mass and center-of-mass for each node (post-order)
void computeMassDistribution(QuadNode* node) {
    if(!node) return; // nothing to do
    bool leaf = !node->children[0]; // if no children, it's leaf
    if(leaf) {
        // if leaf has a body, set mass and center-of-mass
        if(node->singleBody) {
            node->totalMass    = node->singleBody->m;
            node->centerOfMass = node->singleBody->coordinates;
        }
        return;
    }
    // not a leaf: sum over children
    double msum = 0;
    Vector cm(0,0);
    for(int i=0; i<4; ++i) {
        if(node->children[i]) {
            computeMassDistribution(node->children[i]); // recurse
            double m = node->children[i]->totalMass;
            msum += m;
            cm += node->children[i]->centerOfMass * m; // weighted sum
        }
    }
    node->totalMass = msum;
    if(msum > 0)
        node->centerOfMass = cm / msum; // divide to get COM
}

// 4. Compute force on body bi from the quadtree node
Vector forceOnBody(const Body& bi, QuadNode* node, double theta) {
    Vector zero(0,0);
    if(!node || node->totalMass == 0) return zero; // no force
    // if this node is the same body and a leaf, skip
    if(node->singleBody == &bi && !node->children[0]) return zero;
    // vector from bi to node's center-of-mass
    Vector d = node->centerOfMass - bi.coordinates;
    double r2 = d.data[0]*d.data[0] + d.data[1]*d.data[1] + 1e-12;
    double r  = std::sqrt(r2);
    double size = node->region.maxX - node->region.minX;
    // if leaf or far enough (size/r < theta), treat as one body
    if(!node->children[0] || size/r < theta) {
        double f = G * bi.m * node->totalMass / r2; // magnitude
        // return force vector
        return Vector(d.data[0]/r * f, d.data[1]/r * f);
    }
    // otherwise, sum forces from children
    Vector sum(0,0);
    for(int i=0; i<4; ++i) {
        if(node->children[i])
            sum += forceOnBody(bi, node->children[i], theta);
    }
    return sum;
}

// 5. Compute forces for all bodies in parallel
void computeForcesParallel(System& universe, QuadNode* root, double theta) {
    int n = universe.bodies.size();
    std::vector<Vector> F(n);
    // for each body i, compute F[i] in parallel
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<n; ++i)
        F[i] = forceOnBody(universe.bodies[i], root, theta);
    // write back accelerations in parallel
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<n; ++i)
        universe.bodies[i].acceleration = F[i] / universe.bodies[i].m;
}

void computeForcesSerial(System& universe, QuadNode* root, double theta) {
    int n = universe.bodies.size();
    for (int i = 0; i < n; ++i) {
        Vector f = forceOnBody(universe.bodies[i], root, theta);
        universe.bodies[i].acceleration.data[0] = f.data[0] / universe.bodies[i].m;
        universe.bodies[i].acceleration.data[1] = f.data[1] / universe.bodies[i].m;
    }
}

// 6. Update all bodies' positions and velocities in parallel
void updateBodies(System& universe, double dt) {
    int n = universe.bodies.size();
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<n; ++i)
        universe.bodies[i].update(dt);
}

// 7. Recursively delete quadtree nodes
void freeQuadTree(QuadNode* node) {
    if(!node) return; // nothing
    for(int i=0; i<4; ++i)
        freeQuadTree(node->children[i]); // recurse
    delete node; // free this node
}

// Wrapper: do one full Barnes–Hut step (build, force, update, cleanup)
void BarnesHutStep(System& universe, double dt, double theta, bool useParallel) {
    Bounds b = computeBounds(universe);                 // get bounds
    QuadNode* root = createRootNode(b);                 // new root
    for(auto& c : universe.bodies)                      // insert each body
        insertBody(root, c);
    computeMassDistribution(root);                       // compute mass at each node
    if (useParallel) {
        computeForcesParallel(universe, root, theta);
    } else {
        computeForcesSerial(universe, root, theta);
    }
    updateBodies(universe, dt);                          // update positions
    freeQuadTree(root);                                  // delete tree
}
