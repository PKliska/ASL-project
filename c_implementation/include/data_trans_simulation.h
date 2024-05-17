#ifndef DATA_TRANS_SIMULATION_H
#define DATA_TRANS_SIMULATION_H
#include <stdio.h>
#include "preallocated_simulation.h"

typedef struct preallocated_simulation data_trans_simulation;

data_trans_simulation* new_data_trans_simulation(
    size_t dimension, double size, double rho, double nu);

void advance_data_trans_simulation(data_trans_simulation* sim,
                                 unsigned int steps, unsigned int pit,
                                 double dt);

void write_data_trans_simulation(data_trans_simulation* sim,
                               FILE* fp);

void destroy_data_trans_simulation(data_trans_simulation* sim);
#endif
