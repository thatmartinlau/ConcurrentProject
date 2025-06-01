#include "barneshutt.hpp"
#include <cmath>
#include <omp.h>

// QuadNode implementation
QuadNode::QuadNode(const Bounds& r)
    : region(r), totalMass(0), centerOfMass(0,0), singleBody(nullptr) {
    for(int i=0;i<4;++i) children[i] = nullptr;
}
QuadNode::~QuadNode() {}

// 1. Bounds from telemetry
Bounds computeBounds(const System& universe) {
    auto [minV, maxV] = universe.exposeBounds();
    return { minV.data[0], maxV.data[0], minV.data[1], maxV.data[1] };
}

// Helpers
static int quadrant(const Bounds& b, const Vector& v) {
    double mx = 0.5*(b.minX + b.maxX);
    double my = 0.5*(b.minY + b.maxY);
    bool east = v.data[0] > mx;
    bool north= v.data[1] > my;
    if(north) return east?1:0; else return east?3:2;
}
static Bounds childBounds(const Bounds& b, int q) {
    double mx = 0.5*(b.minX + b.maxX);
    double my = 0.5*(b.minY + b.maxY);
    switch(q) {
      case 0: return {b.minX,mx, my,b.maxY}; // NW
      case 1: return {mx,b.maxX, my,b.maxY}; // NE
      case 2: return {b.minX,mx, b.minY,my}; // SW
      case 3: return {mx,b.maxX, b.minY,my}; // SE
      default: return b;
    }
}

// 2. Build
QuadNode* createRootNode(const Bounds& bounds) {
    return new QuadNode(bounds);
}
void insertBody(QuadNode* node, Body& body) {
    // empty leaf
    if(!node->singleBody && !node->children[0]){
        node->singleBody = &body;
        return;
    }
    // subdivide if leaf with one
    if(node->singleBody) {
        Body* old = node->singleBody;
        node->singleBody = nullptr;
        for(int i=0;i<4;++i)
            node->children[i] = new QuadNode(childBounds(node->region,i));
        insertBody(node->children[ quadrant(node->region,old->coordinates) ], *old);
    }
    // insert
    int q = quadrant(node->region, body.coordinates);
    if(!node->children[q])
        node->children[q] = new QuadNode(childBounds(node->region,q));
    insertBody(node->children[q], body);
}

// 3. Mass distribution
void computeMassDistribution(QuadNode* node) {
    if(!node) return;
    bool leaf = !node->children[0];
    if(leaf) {
        if(node->singleBody) {
            node->totalMass    = node->singleBody->m;
            node->centerOfMass = node->singleBody->coordinates;
        }
        return;
    }
    double msum=0; Vector cm(0,0);
    for(int i=0;i<4;++i) if(node->children[i]){
        computeMassDistribution(node->children[i]);
        double m = node->children[i]->totalMass;
        msum += m;
        cm += node->children[i]->centerOfMass * m;
    }
    node->totalMass = msum;
    if(msum>0) node->centerOfMass = cm / msum;
}

// 4. Force on body
Vector forceOnBody(const Body& bi, QuadNode* node, double theta) {
    Vector zero(0,0);
    if(!node || node->totalMass==0) return zero;
    if(node->singleBody==&bi && !node->children[0]) return zero;
    Vector d = node->centerOfMass - bi.coordinates;
    double r2 = d.data[0]*d.data[0] + d.data[1]*d.data[1] + 1e-12;
    double r  = std::sqrt(r2);
    double size = node->region.maxX - node->region.minX;
    if(!node->children[0] || size/r < theta) {
        double f = G * bi.m * node->totalMass / r2;
        return Vector(d.data[0]/r*f, d.data[1]/r*f);
    }
    Vector sum(0,0);
    for(int i=0;i<4;++i) if(node->children[i])
        sum += forceOnBody(bi, node->children[i], theta);
    return sum;
}

// 5. Parallel forces
void computeForcesParallel(System& universe, QuadNode* root, double theta) {
    int n = universe.bodies.size();
    std::vector<Vector> F(n);
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<n;++i)
        F[i] = forceOnBody(universe.bodies[i], root, theta);
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<n;++i)
        universe.bodies[i].acceleration = F[i] / universe.bodies[i].m;
}

// 6. Integrate
void updateBodies(System& universe, double dt) {
    int n = universe.bodies.size();
    #pragma omp parallel for schedule(dynamic)
    for(int i=0;i<n;++i)
        universe.bodies[i].update(dt);
}

// 7. Cleanup
void freeQuadTree(QuadNode* node) {
    if(!node) return;
    for(int i=0;i<4;++i) freeQuadTree(node->children[i]);
    delete node;
}

// Wrapper step
void BarnesHutStep(System& universe, double dt, double theta) {
    Bounds b = computeBounds(universe);
    QuadNode* root = createRootNode(b);
    for(auto& c:universe.bodies) insertBody(root,c);
    computeMassDistribution(root);
    computeForcesParallel(universe, root, theta);
    updateBodies(universe, dt);
    freeQuadTree(root);
}
