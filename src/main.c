#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>

#include "tsp.h"
#include "tsp_utilities.h"

double second();
void parse_command_line(int argc, char** argv, instance *inst);
void read_input(instance *inst);



int main(int argc, char **argv) {

	if ( argc < 2 ) {
		printf("Usage: %s -help for help\n", argv[0]);
		exit(1);
	}

	//double t1 = second();
	//printf("Timestamp 1: %f", t1);

	instance inst;

	parse_command_line(argc, argv, &inst);     
	
	//printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);

	parse_input_file(&inst);

	/* 
	if ( TSPopt(&inst) ) print_error(" error within VRPopt()");
	double t2 = second(); 
    
	if ( VERBOSE >= 1 )   
	{
		printf("... VRP problem solved in %lf sec.s\n", t2-t1);  
	}
	*/

	free_instance(&inst);

	return 0; 

}