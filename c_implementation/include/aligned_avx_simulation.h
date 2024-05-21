#ifndef ALIGNED_AVX_SIMULATION_H
#define ALIGNED_AVX_SIMULATION_H
#include <stdio.h>
#include "faster_math_simulation.h"

typedef struct aligned_AVX_simulation{
    faster_math_simulation base;
    double *tmp;
} aligned_AVX_simulation;

aligned_AVX_simulation* new_aligned_AVX_simulation(
    size_t dimension, double size, double rho, double nu);

void init_aligned_AVX_simulation(aligned_AVX_simulation* sim,
    size_t dimension, double size, double rho, double nu);

void advance_aligned_AVX_simulation(aligned_AVX_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_aligned_AVX_simulation(aligned_AVX_simulation* sim,
                               FILE* fp);

void destroy_aligned_AVX_simulation(aligned_AVX_simulation* sim);
void deinit_aligned_AVX_simulation(aligned_AVX_simulation* sim);
#endif
