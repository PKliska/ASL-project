#ifndef CUTOFF_TRAPEZE_SIMULATION_H
#define CUTOFF_TRAPEZE_SIMULATION_H
#include <stdio.h>
#include "utils.h"
#include "trapeze_simulation.h"

typedef struct cutoff_trapeze_simulation{
    trapeze_simulation base;
    size_t poisson_len, poisson_cap;
    struct trapeze* poisson_order;
} cutoff_trapeze_simulation;

cutoff_trapeze_simulation* new_cutoff_trapeze_simulation(
    size_t dimension, double size, double rho, double nu);

void advance_cutoff_trapeze_simulation(cutoff_trapeze_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_cutoff_trapeze_simulation(cutoff_trapeze_simulation* sim,
                               FILE* fp);

void destroy_cutoff_trapeze_simulation(cutoff_trapeze_simulation* sim);
#endif
