#include <stdio.h>
#include <stdlib.h>
#include <cargs.h>
#include <stdbool.h>
#include <string.h>
#include "simulation.h"
#include "baseline_simulation.h"
#include "tsc_x86.h"

#define NUM_RUNS 30
#define CYCLES_REQUIRED 1e8

struct implementation {
    const char* name;
    struct simulation* (*create)(size_t nx, size_t ny,
                                 double rho, double nu);
};

static const struct implementation IMPLEMENTATIONS[] = {
    {.name = "baseline",
     .create = new_baseline_simulation
    }
};

struct arguments {
    // Where we store the output of the simulation
    char *output_file;
    // Number of iterations of the simulation to run
    unsigned int num_iter;
    // Which implementation to use
    struct implementation implementation;
    // Should the simulation be timed?
    bool should_time;
};

static struct cag_option options[] = {
    {.identifier = 'o',
     .access_letters = "o",
     .access_name = "output_file",
     .value_name = "OUTPUT_FILE",
     .description = "File to which the result of the simulation"
		    "should be saved"
    },
    {.identifier = 'n',
     .access_letters = "n",
     .access_name = "num_iter",
     .value_name = "NUM_ITER",
     .description = "Number of simulation steps to simulate"
    },
    {.identifier = 'I',
     .access_letters = "I",
     .access_name = "implementetion",
     .value_name = "IMPLEMENTATION",
     .description = "Implementation of the simulation to use"
                    "(default: baseline)"
    },
    {.identifier = 't',
     .access_letters = "t",
     .access_name = "time",
     .description = "Time the execution of the simulation"
    }
};

static void parse_args(struct arguments* args, int argc, char* argv[]){
    cag_option_context context;
    cag_option_init(&context, options, CAG_ARRAY_SIZE(options),
		    argc, argv);
    while(cag_option_fetch(&context)){
        switch(cag_option_get_identifier(&context)){
            case 'o':
                args->output_file = cag_option_get_value(&context);
                break;
            case 'n':
                //todo: check that arg is really an unsigned int
                args->num_iter = atoi(cag_option_get_value(&context));
                break;
            case 'I':
                const char* impl_name = cag_option_get_value(&context);
                for(size_t i = 0;i<sizeof(IMPLEMENTATIONS)/sizeof(IMPLEMENTATIONS[0]);i++){
                    if(strcmp(IMPLEMENTATIONS[i].name, impl_name) == 0){
                        args->implementation = IMPLEMENTATIONS[i];
                    }
                }
                break;
            case 't':
                args->should_time = true;
                break;
            case '?':
                cag_option_print_error(&context, stderr);
                break;
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
double rdtsc(struct simulation* sim,
				     unsigned int steps, unsigned int pit,
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
            advance_simulation(sim, steps, pit, dt);
            cycles += stop_tsc(start);
        }

        if(cycles >= CYCLES_REQUIRED) break;

        num_runs *= 2;
    }

    cycles = 0;
    for (i = 0; i < num_runs; ++i) {
        start = start_tsc();
        advance_simulation(sim, steps, pit, dt);
        cycles += stop_tsc(start);
    }

    cycles = cycles/num_runs;
    return (double) cycles;
}


int main(int argc, char* argv[]){
    struct arguments arguments = {
        .output_file = "sim.csv",
        .num_iter = 100,
        .implementation = IMPLEMENTATIONS[0],
        .should_time = false
    };
    parse_args(&arguments, argc, argv);
    #ifdef DEBUG 
    fprintf(stderr, "output_file = %s\nnum_iter = %u\n",
        arguments.output_file, arguments.num_iter);
    #endif
    struct simulation * sim = arguments.implementation.create(41,41,
                                                               1, 0.1);
    advance_simulation(sim, arguments.num_iter, 50, 0.001);
    FILE* output_file = fopen(arguments.output_file, "w");
    if(!output_file){
        fprintf(stderr, "can't write to output file %s\n",
                arguments.output_file);
        return -1;
    }
    write_simulation(sim, output_file);
    if (arguments.should_time) {
        int n = 100; // array length
        int sizes[n];
     
        int start_value = 10;
        int step_size = (500 - 10) / (n - 1);
        for (int i = 0; i < n; i++) {
            sizes[i] = start_value + i * step_size;
        }

        for (int i = 0; i < n; i++) {
            struct simulation * sim = arguments.implementation.create(
                sizes[i], sizes[i], 1, 0.1
            );
            unsigned int steps = 10;
            unsigned int pit = 50;
            double dt = 0.001;
            double cycles = rdtsc(sim, steps, pit, dt);
            printf("Cycles = %f \n", cycles);
        }
    }


   destroy_simulation(sim);
}
