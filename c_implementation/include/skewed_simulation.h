#ifndef SKEWED_SIMULATION_H
#define SKEWED_SIMULATION_H
#include <stdio.h>
#include "faster_math_simulation.h"

typedef faster_math_simulation skewed_simulation;

skewed_simulation* new_skewed_simulation(
    size_t dimension, double size, double rho, double nu);

void init_skewed_simulation(skewed_simulation* sim,
    size_t dimension, double size, double rho, double nu);

void advance_skewed_simulation(skewed_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_skewed_simulation(skewed_simulation* sim,
                               FILE* fp);

void destroy_skewed_simulation(skewed_simulation* sim);
void deinit_skewed_simulation(skewed_simulation* sim);
#endif
