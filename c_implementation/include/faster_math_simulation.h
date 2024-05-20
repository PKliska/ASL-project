#ifndef FASTER_MATH_SIMULATION_H
#define FASTER_MATH_SIMULATION_H
#include <stdio.h>
#include "preallocated_simulation.h"

typedef struct preallocated_simulation faster_math_simulation;

faster_math_simulation* new_faster_math_simulation(
    size_t dimension, double size, double rho, double nu);

void init_faster_math_simulation(faster_math_simulation* sim,
    size_t dimension, double size, double rho, double nu);

void faster_math_build_up_b(const faster_math_simulation* sim,
                            double dt);

void advance_faster_math_simulation(faster_math_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_faster_math_simulation(faster_math_simulation* sim,
                               FILE* fp);

void destroy_faster_math_simulation(faster_math_simulation* sim);
void deinit_faster_math_simulation(faster_math_simulation* sim);
#endif
