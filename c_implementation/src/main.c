#include <stdio.h>
#include <stdlib.h>
#include <cargs.h>
#include "baseline_simulation.h"

struct arguments {
    // Where we store the output of the simulation
    char *output_file;
    // Number of iterations of the simulation to run
    unsigned int num_iter;
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
	case '?':
	    cag_option_print_error(&context, stderr);
	    break;
	}
    }
}


int main(int argc, char* argv[]){
   struct arguments arguments;
   arguments.output_file = "sim.csv";
   arguments.num_iter = 100;
   parse_args(&arguments, argc, argv);
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
