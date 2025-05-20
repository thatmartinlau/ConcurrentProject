#ifndef PARTICLE_MESH_THREAD_H
#define PARTICLE_MESH_THREAD_H



#include "core.hpp"
#include <vector>


void particle_mesh_simulation_parallel(System &universe, double dt, int grid_size, size_t num_threads);

#endif 
