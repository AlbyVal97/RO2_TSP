#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>

#include "tsp.h"

// Command line parameter list examples:
// -f ../data/inst_max_100_nodes/berlin52.tsp -seed 123456 -model_type 1 -verbose 0 -time_limit 3600
// -m 0 -n_inst 20 -n_nodes 50 -f ../data/inst_max_100_nodes/berlin52.tsp -seed 123456 -model_type 3 -verbose 2 -time_limit 60

int main(int argc, char **argv) {

	if ( argc < 2 ) {
		printf("Usage: %s -help for help\n", argv[0]);
		exit(1);
	}

	instance inst;

	parse_command_line(argc, argv, &inst);

	if (inst.mode == CREATE_INSTANCES) {
		create_instances(&inst);
		return 0;
	}

	parse_input_file(&inst);

	double t1 = second();

	if (inst.verbose >= HIGH) { printf("Timestamp 1: %lf seconds.\n", t1); }

	if (TSPopt(&inst)) { print_error(" error within TSPopt()"); }

	double t2 = second();

	if (inst.verbose == TEST) {
		FILE* csv_file = fopen("../../outputs/test.csv", "a");
		if (csv_file == NULL) print_error("test.csv file not found inside \"outputs\" folder!");

		if (inst.model_type == MTZ_STATIC) {				// Print the instance name just for the first test execution (on the first model)
			fprintf(csv_file, "%s, %f, ", inst.inst_name, t2 - t1);
		}
		else if (inst.model_type == GG) {					// Go to next line only when the last test (on the last model) has been executed
			fprintf(csv_file, "%f\n", t2 - t1);
		}
		else {
			fprintf(csv_file, "%f, ", t2 - t1);
		}

		fclose(csv_file);
	}
    
	if ( inst.verbose >= TEST ) printf("TSP problem solved successfully in %lf seconds.\n\n", t2-t1);

	

	// Plot the resulting optimal tour using Gnuplot
	if (inst.verbose >= MEDIUM) system("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe ../outputs/gnuplot_commands.txt");

	free_instance(&inst);

	return 0; 

}