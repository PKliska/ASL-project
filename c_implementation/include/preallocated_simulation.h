#ifndef PREALLOCATED_SIMULATION_H
#define PREALLOCATED_SIMULATION_H
#include <stdio.h>
#include "simulation.h"

struct preallocated_simulation{
    struct simulation base;
    size_t d;
    double rho,nu,size;
    double* u;
    double* un;
    double* v;
    double* vn;
    double* p;
    double* pn;
    double* b;
};

struct preallocated_simulation* new_preallocated_simulation(
    size_t dimension, double size, double rho, double nu);

void advance_preallocated_simulation(struct preallocated_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_preallocated_simulation(struct preallocated_simulation* sim,
                               FILE* fp);

void destroy_preallocated_simulation(struct preallocated_simulation* sim);
#endif
