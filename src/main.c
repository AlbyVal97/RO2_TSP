#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>

#include "tsp.h"
#include "tsp_utilities.h"


int main(int argc, char **argv) {

	if ( argc < 2 ) {
		printf("Usage: %s -help for help\n", argv[0]);
		exit(1);
	}

	double t1 = second();

	instance inst;

	parse_command_line(argc, argv, &inst);

	if (inst.verbose >= HIGH) { printf("Timestamp 1: %lf seconds.\n", t1); }
	
	//printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);

	parse_input_file(&inst);

	if ( TSPopt(&inst) ) print_error(" error within TSPopt()");

	double t2 = second();
    
	if ( inst.verbose >= LOW ) { printf("TSP problem solved successfully in %lf seconds.\n", t2-t1); }

	free_instance(&inst);

	// Plot the resulting optimal tour using Gnuplot
	system("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe ../src/gnuplot_commands.txt");

	return 0; 

}