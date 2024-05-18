#ifndef TRAPEZE_SIMULATION_H
#define TRAPEZE_SIMULATION_H
#include <stdio.h>
#include "preallocated_simulation.h"

typedef struct preallocated_simulation trapeze_simulation;

trapeze_simulation* new_trapeze_simulation(
    size_t dimension, double size, double rho, double nu);

void advance_trapeze_simulation(trapeze_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_trapeze_simulation(trapeze_simulation* sim,
                               FILE* fp);

void destroy_trapeze_simulation(trapeze_simulation* sim);
#endif
