#ifndef VANILLA_AVX_SIMULATION_H
#define VANILLA_AVX_SIMULATION_H
#include <stdio.h>
#include "preallocated_simulation.h"

typedef struct preallocated_simulation vanilla_avx_simulation;

vanilla_avx_simulation* new_vanilla_avx_simulation(
    size_t dimension, double size, double rho, double nu);

void advance_vanilla_avx_simulation(vanilla_avx_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_vanilla_avx_simulation(vanilla_avx_simulation* sim,
                               FILE* fp);

void destroy_vanilla_avx_simulation(vanilla_avx_simulation* sim);
#endif
