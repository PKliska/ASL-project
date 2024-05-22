#ifndef BLOCKING_SIMULATION_H
#define BLOCKING_SIMULATION_H
#include <stdio.h>
#include "preallocated_simulation.h"

typedef struct preallocated_simulation blocking_simulation;

blocking_simulation *new_blocking_simulation(size_t dimension, double size, double rho, double nu);

void advance_blocking_simulation(blocking_simulation *sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_blocking_simulation(blocking_simulation *sim,
                               FILE *fp);

void destroy_blocking_simulation(blocking_simulation *sim);
#endif
