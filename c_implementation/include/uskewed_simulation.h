#ifndef USKEWED_SIMULATION_H
#define USKEWED_SIMULATION_H
#include <stdio.h>
#include "trapeze_macros.h"
#include "faster_math_simulation.h"

typedef faster_math_simulation uskewed_simulation;

uskewed_simulation* new_uskewed_simulation(
    size_t dimension, double size, double rho, double nu);

void init_uskewed_simulation(uskewed_simulation* sim,
    size_t dimension, double size, double rho, double nu);

void advance_uskewed_simulation(uskewed_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_uskewed_simulation(uskewed_simulation* sim,
                               FILE* fp);

void destroy_uskewed_simulation(uskewed_simulation* sim);
void deinit_uskewed_simulation(uskewed_simulation* sim);
#endif
