#include <stdio.h>
#include <stdlib.h>
#include <argp.h>
#include "baseline_simulation.h"

struct arguments {
    // Where we store the output of the simulation
    char *output_file;
    // Number of iterations of the simulation to run
    unsigned int num_iter;
};

static struct argp_option options[] = {
    {"output_file", 'o',
     "OUTPUT_FILE", 0,
     "Specify where the result of the simulation should be saved",
    },
    {"num_iter", 'n',
     "NUM_ITER", 0,
     "Specify the number of iterations to simulate",
    },
    {0}
};

static error_t parse_opt(int key, char *arg, struct argp_state *state){
    struct arguments *arguments = state->input;
    switch(key){
    case 'o':
	arguments->output_file = arg;
	break;
    case 'n':
	//todo: check that arg is really an unsigned int
	arguments->num_iter = atoi(arg);
	break;
    default:
	return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = {options, parse_opt, 0, 0};

int main(int argc, char* argv[]){
   struct arguments arguments;
   arguments.output_file = "sim.csv";
   arguments.num_iter = 100;
   argp_parse(&argp, argc, argv, 0, 0, &arguments);
   #ifdef DEBUG 
   fprintf(stderr, "output_file = %s\nnum_iter = %u\n",
	   arguments.output_file, arguments.num_iter);
   #endif
   struct baseline_simulation * sim = new_baseline_simulation(41,41,
							      1, 0.1);
   advance_baseline_simulation(sim, arguments.num_iter, 50, 0.001);
   FILE* output_file = fopen(arguments.output_file, "w");
   if(!output_file){
	fprintf(stderr, "can't write to output file %s\n",
			arguments.output_file);
	return -1;
   }
   write_baseline_simulation(sim, output_file);
   destroy_baseline_simulation(sim);
}
