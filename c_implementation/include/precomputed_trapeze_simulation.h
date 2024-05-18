#ifndef PRECOMPUTED_TRAPEZE_SIMULATION_H
#define PRECOMPUTED_TRAPEZE_SIMULATION_H
#include <stdio.h>
#include "trapeze_simulation.h"

typedef struct trapeze_base_case_args{
  unsigned t;
  unsigned x0, x1;
  unsigned y0, y1;
} trapeze_base_case_args;

typedef struct precomputed_trapeze_simulation{
    trapeze_simulation base;
    size_t poisson_len, poisson_cap;
    trapeze_base_case_args* poisson_order;
} precomputed_trapeze_simulation;

precomputed_trapeze_simulation* new_precomputed_trapeze_simulation(
    size_t dimension, double size, double rho, double nu);

void advance_precomputed_trapeze_simulation(precomputed_trapeze_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_precomputed_trapeze_simulation(precomputed_trapeze_simulation* sim,
                               FILE* fp);

void destroy_precomputed_trapeze_simulation(precomputed_trapeze_simulation* sim);
#endif
