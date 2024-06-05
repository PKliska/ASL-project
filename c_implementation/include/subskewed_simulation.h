#ifndef subskewed_SIMULATION_H
#define subskewed_SIMULATION_H
#include <stdio.h>
#include "faster_math_simulation.h"

typedef faster_math_simulation subskewed_simulation;

subskewed_simulation* new_subskewed_simulation(
    size_t dimension, double size, double rho, double nu);

void init_subskewed_simulation(subskewed_simulation* sim,
    size_t dimension, double size, double rho, double nu);

void advance_subskewed_simulation(subskewed_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_subskewed_simulation(subskewed_simulation* sim,
                               FILE* fp);

void destroy_subskewed_simulation(subskewed_simulation* sim);
void deinit_subskewed_simulation(subskewed_simulation* sim);
#endif
