#ifndef INSTANCE_H_

#define INSTANCE_H_

// Hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ

#define N_MODELS 8
#define N_VERBOSITIES 4
#define N_MODES 3

static const char* models[] = { "BASIC", "MTZ_STATIC", "MTZ_LAZY", "MTZ_SEC2_STATIC", "MTZ_SEC2_LAZY", "GG", "BENDERS", "BRANCH_CUT" };
static const char* verbosities[] = { "TEST", "LOW", "MEDIUM", "HIGH" };
static const char* modes[] = { "DEFAULT", "CREATE_INSTANCES", "RUN_TEST" };

// Verbosity levels enumeration
typedef enum {
	TEST,
	LOW,
	MEDIUM,
	HIGH
} verbose;

// Model type enumeration
typedef enum {
	BASIC,
	MTZ_STATIC,
	MTZ_LAZY,
	MTZ_SEC2_STATIC,
	MTZ_SEC2_LAZY,
	GG,
	BENDERS,
	BRANCH_CUT
} model_type;

typedef enum {
	DEFAULT,
	CREATE_INSTANCES,
	RUN_TEST
} mode;

// Problem instance data structure
typedef struct {

	// Input data
	int nnodes;
	char inst_name[100];
	double* xcoord;
	double* ycoord;
	int integer_costs;                      // = 1 for integer costs (rounded distances), 0 otherwise

	// Command line parameters list
	int mode;								// working mode of the program

	// Command line arguments for CREATE_INSTANCES mode
	int n_instances;						// number of instances to generate
	int n_nodes_per_instance;				// number of nodes (fixed) for the instances to be generated
	
	char folder_istances[100];				// name of folder inside /data where to find the instances to be created in CREATE_INSTANCES mode or tested in RUN_TEST mode
	char instance_prefix_name[20];			// prefix of instance file name (used both in CREATE_INSTANCES and RUN_TEST modes)

	// Command line arguments for TEST mode
	char test_name[100];
	int n_models_test;
	int* models_to_test;

	// Command line arguments for DEFAULT mode
	int model_type;							// model number to be used to solve the instance
	double timelimit;						// overall time limit, in seconds
	char input_file[1000];		  			// input file name
	int verbose;							// verbosity value
	int seed;								// internal branching random seed used by Cplex. If fixed, leads to more consistent computational time
	int first_model;
	int last_model;

	// global data
	double zbest;							// value of the best sol. available

} instance;

#endif   /* INSTANCE_H_ */
