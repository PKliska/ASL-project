#include <stdio.h>
#include <stdlib.h>
#include "baseline_simulation.h"
#include "tsc_x86.h"
#define NUM_RUNS 30
#define CYCLES_REQUIRED 1e8

void fill_matrix(double * A, int n, int m) {
    for(int i=0; i < n; i++) {
        for(int j=0; j < m; j++) {
            A[n*i+j] = (double) rand() / RAND_MAX;
        }
    }
}

/*
 * Timing function based on the TimeStep Counter of the CPU.
 * 
 * Arguments:
 * sim: simulation
 * b: intermediate matrix
 * dt: size of step
 * pit: number of steps for computation of pressure
 */
double rdtsc(struct baseline_simulation* sim,
				     double *b, unsigned int pit,
				     double dt) {
    int i, num_runs;
    myInt64 cycles = 0;
    myInt64 start;
    num_runs = NUM_RUNS;

    /*
     * The CPUID instruction serializes the pipeline.
     * Using it, we can create execution barriers around the code we want to time.
     * The calibrate section is used to make the computation large enough so as to
     * avoid measurements bias due to the timing overhead.
     */
    while(num_runs < (1 << 14)) {
        for (i = 0; i < num_runs; ++i) {
            start = start_tsc();
            step_baseline_simulation(sim, b, pit, dt);
            cycles += stop_tsc(start);
        }

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }

    cycles = 0;
    for (i = 0; i < num_runs; ++i) {
        start = start_tsc();
        step_baseline_simulation(sim, b, pit, dt);
        cycles += stop_tsc(start);
    }

    cycles = cycles/num_runs;
    return (double) cycles;
}

int main() {

    int n = 100; // array length
    int sizes[n];
     
    int start_value = 10;
    int step_size = (500 - 10) / (n - 1);
    for (int i = 0; i < n; i++) {
        sizes[i] = start_value + i * step_size;
    }

    for (int i = 0; i < n; i++) {
        struct baseline_simulation * sim = new_baseline_simulation(sizes[i], sizes[i],
                                    1, 0.1);
        double* b = (double *)malloc(sim->nx*sim->ny*sizeof(double));
        fill_matrix(b, sim->nx, sim->ny);
        double pit = 50;
        double dt = 0.001;
        double cycles = rdtsc(sim, b, pit, dt);
        printf("Size = %d, Cycles = %f \n", sizes[i], cycles);
    }

    return 0;
}