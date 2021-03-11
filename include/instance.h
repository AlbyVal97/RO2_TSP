#ifndef INSTANCE_H_

#define INSTANCE_H_

#define VERBOSE				    1000		// printing level  (debug_mode if VERBOSE >= 1000)

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ

//data structures
typedef struct {

	// input data
	int nnodes;
	double* xcoord;
	double* ycoord;
	int integer_costs;                      // = 1 for integer costs (rounded distances), 0 otherwise

	// parameters 
	int model_type;
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file

	// global data
	double zbest;							// value of the best sol. available

} instance;

#endif   /* INSTANCE_H_ */
