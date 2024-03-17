#ifndef BASELINE_SIMULATION_H
#define BASELINE_SIMULATION_H
#include <stdio.h>

struct baseline_simulation{
    size_t nx, ny;
    double rho,nu;
    double* u;
    double* v;
    double* p;
};

struct baseline_simulation* new_baseline_simulation(size_t nx, size_t ny,
		                                    double rho, double nu);

void advance_baseline_simulation(struct baseline_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_baseline_simulation(struct baseline_simulation* sim, FILE* fp);

void destroy_baseline_simulation(struct baseline_simulation* sim);
#endif
