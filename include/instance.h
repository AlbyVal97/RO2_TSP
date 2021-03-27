#ifndef INSTANCE_H_

#define INSTANCE_H_

// Hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ


// Verbosity levels enumeration
typedef enum {
	LOW,
	MEDIUM,
	HIGH
} verbose;

// Model type enumeration
typedef enum {
	BASIC,
	MTZ_STATIC,
	MTZ_LAZY,
	MTZ_SUBTOUR_SIZE_2,
	GG
} model_type;

// Problem instance data structure
typedef struct {

	// input data
	int nnodes;
	char inst_name[100];
	double* xcoord;
	double* ycoord;
	int integer_costs;                      // = 1 for integer costs (rounded distances), 0 otherwise

	// parameters 
	int model_type;
	double timelimit;						// overall time limit, in seconds
	char input_file[1000];		  			// input file name
	int verbose;							// verbosity value

	// global data
	double zbest;							// value of the best sol. available

} instance;

#endif   /* INSTANCE_H_ */
