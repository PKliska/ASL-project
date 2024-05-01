#ifndef FASTER_MATH_SIMULATION_H
#define FASTER_MATH_SIMULATION_H
#include <stdio.h>
#include "preallocated_simulation.h"

typedef struct preallocated_simulation faster_math_simulation;

faster_math_simulation* new_faster_math_simulation(
    size_t nx, size_t ny, double rho, double nu);

void advance_faster_math_simulation(faster_math_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_faster_math_simulation(faster_math_simulation* sim,
                               FILE* fp);

void destroy_faster_math_simulation(faster_math_simulation* sim);
#endif
