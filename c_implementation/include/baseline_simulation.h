#ifndef BASELINE_SIMULATION_H
#define BASELINE_SIMULATION_H
#include <stdio.h>
#include "simulation.h"

struct baseline_simulation{
    struct simulation base;
    size_t d;
    double rho,nu,size;
    double* u;
    double* v;
    double* p;
};

struct baseline_simulation* new_baseline_simulation(size_t dimension,
                                                    double size,
                                                    double rho,
                                                    double nu);

void advance_baseline_simulation(struct baseline_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_baseline_simulation(struct baseline_simulation* sim, FILE* fp);

void destroy_baseline_simulation(struct baseline_simulation* sim);
#endif
