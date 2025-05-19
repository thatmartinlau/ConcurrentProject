#ifndef PARTICLE_MESH_CUDA_H
#define PARTICLE_MESH_CUDA_H

#include "core.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>


void particle_mesh_simulation_CUDA(System &universe, double dt, int grid_size, int num_bodies);

#endif