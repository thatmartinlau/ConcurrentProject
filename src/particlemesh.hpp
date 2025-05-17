#ifndef PARTICLE_MESH_H
#define PARTICLE_MESH_H

#include "core.hpp"
#include <vector>
#include <fftw3.h>


class ParticleMeshSystem : public System {
public:
    ParticleMeshSystem();
    ~ParticleMeshSystem();

    void simulate(double dt) override;

private:
    static constexpr int GRID_SIZE = 128;
    static constexpr double SOFTENING = 1e-2;

    std::vector<std::vector<double>> density;
    fftw_complex* potential_freq;
    fftw_complex* density_freq;
    double* potential;

    fftw_plan forward_plan;
    fftw_plan backward_plan;

    void computeDensity();
    void solvePotential();
    void computeForcesAndUpdate(double dt);
};

#endif 
