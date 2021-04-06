#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>

#include "tsp.h"

// Command line parameter list example: -f ../data/inst_max_100_nodes/berlin52.tsp -seed 123456 -model_type 1 -verbose 0 -time_limit 3600

int main(int argc, char **argv) {

	if ( argc < 2 ) {
		printf("Usage: %s -help for help\n", argv[0]);
		exit(1);
	}

	instance inst;

	parse_command_line(argc, argv, &inst);

	parse_input_file(&inst);

	double t1 = second();

	if (inst.verbose >= HIGH) { printf("Timestamp 1: %lf seconds.\n", t1); }

	if (TSPopt(&inst)) { print_error(" error within TSPopt()"); }

	double t2 = second();
    
	if ( inst.verbose >= LOW ) { printf("TSP problem solved successfully in %lf seconds.\n", t2-t1); }

	

	// Plot the resulting optimal tour using Gnuplot
	if (inst.verbose >= MEDIUM) { system("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe ../src/gnuplot_commands.txt"); }

	free_instance(&inst);

	return 0; 

}