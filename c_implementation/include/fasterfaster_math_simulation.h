#ifndef FASTERFASTER_MATH_SIMULATION_H
#define FASTERFASTER_MATH_SIMULATION_H
#include <stdio.h>
#include "preallocated_simulation.h"

typedef struct preallocated_simulation fasterfaster_math_simulation;

fasterfaster_math_simulation* new_fasterfaster_math_simulation(
    size_t dimension, double size, double rho, double nu);

void advance_fasterfaster_math_simulation(fasterfaster_math_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_fasterfaster_math_simulation(fasterfaster_math_simulation* sim,
                               FILE* fp);

void destroy_fasterfaster_math_simulation(fasterfaster_math_simulation* sim);
#endif
