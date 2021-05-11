#ifndef INSTANCE_H_

#define INSTANCE_H_

// Hard-wired parameters
#define XSMALL		  		  1e-5 		// Tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// Very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// Cplex's ticks (depends on CPU)

#define N_MODELS 22						// Total number of TSP models/approaches implemented to date = |models[]|
#define N_VERBOSITIES 4					// Total number of verbosity levels considered = |verbosities[]|
#define N_MODES 3						// Total number of working modes of the program = |modes[]|


// List of the TSP models/approaches names implemented to date
static const char* models[] = { "BASIC", "MTZ_STATIC", "MTZ_LAZY", "MTZ_SEC2_STATIC", "MTZ_SEC2_LAZY", "GG",						// Compact models
								"BENDERS",																							// Benders-like solution scheme using SECs
								"BRANCH_CUT", "ADVBC_STD", "ADVBC_ROOT", "ADVBC_DEPTH_5", "ADVBC_PROB_50", "ADVBC_PROB_10",			// SEC separation using Cplex's callbacks
								"HEUR_HARD_FIX_50",	"HEUR_HARD_FIX_70",	"HEUR_HARD_FIX_90", "HEUR_HARD_FIX_VAR",					// Hard-fix heuristics variants
								"HEUR_SOFT_FIX_3", "HEUR_SOFT_FIX_5", "HEUR_SOFT_FIX_7", "HEUR_SOFT_FIX_9", "HEUR_SOFT_FIX_VAR" };	// Soft-fix heuristics variants

// List of verbosity levels considered, indicating which and how many details have to be shown on the console
static const char* verbosities[] = { "TEST",		// Prints only summary of input parameters and total time required;
													// Does not create any related file nor calls Gnuplot to avoid overhead during tests;
													// N.B. It can only be used in RUN_TEST mode, otherwise an error occurs.
									 "LOW",			// Prints only summary of input parameters and total time required;
													// Creates logfile, model_*.lp and model_*_edges.dat files;
													// Does not call Gnuplot to show the resulting graph.
									 "MEDIUM",		// Prints summary of input parameters and instance, details about final solution and total time required;
													// Creates logfile, model_*.lp and model_*_edges.dat files;
													// Calls Gnuplot to show the resulting graph.
									 "HIGH" };		// Prints intermediate solutions, loop states, concurrent nodes and cut values (only for debugging purpose);
													// More details are available in the specific logfile.

// List of working modes of the program: each of them requires different input parameters (see "instance" struct for more details)
static const char* modes[] = { "DEFAULT",			// Deafult mode: solves the provided instance with a single model/approach.
							   "CREATE_INSTANCES",	// Special mode: generates a number of random instances with a given number of nodes.
													//				 N.B. Uses a pre-determined seed for each instance of the sequence to allow for reproducibility.
							   "RUN_TEST" };		// Special mode: solves a list of instances using a list of models/approaches, so it consists of n_inst * n_models runs.
													//				 N.B. Uses TEST verbosity. Only generates a .bat script, runs it and returns a .csv file with the results.


// Verbosity levels enumeration
typedef enum {
	TEST, // 0
	LOW, // 1
	MEDIUM, // 2 (DEFAULT)
	HIGH // 3
} verbose;

// Models/approaches names enumeration
typedef enum {
	BASIC, // 0 (DEFAULT)
	MTZ_STATIC, // 1
	MTZ_LAZY, // 2
	MTZ_SEC2_STATIC, // 3
	MTZ_SEC2_LAZY, // 4
	GG, // 5
	BENDERS, // 6
	BRANCH_CUT, // 7
	ADVBC_STD, // 8
	ADVBC_ROOT, // 9
	ADVBC_DEPTH_5, // 10
	ADVBC_PROB_50, // 11
	ADVBC_PROB_10, // 12
	HEUR_HARD_FIX_50, // 13
	HEUR_HARD_FIX_70, // 14
	HEUR_HARD_FIX_90, // 15
	HEUR_HARD_FIX_VAR, // 16
	HEUR_SOFT_FIX_3, // 17
	HEUR_SOFT_FIX_5, // 18
	HEUR_SOFT_FIX_7, // 19
	HEUR_SOFT_FIX_9, // 20
	HEUR_SOFT_FIX_VAR // 21
} model_type;

// Working modes enumeration
typedef enum {
	DEFAULT, // 0 (DEFAULT)
	CREATE_INSTANCES, // 1
	RUN_TEST // 2
} mode;


// TSP instance data structure: includes both information about the instance[s] provided and about the model[s] chosen to solve it[them]
typedef struct {

	// Data from input instance
	int nnodes;								// Number of nodes of the provided instance
	int ncols;								// Number of columns/variables (depends on the TSP model used)
	char inst_name[100];					// Name of the provided instance (indicated inside the *.tsp file)
	double* xcoord;							// Array with X coordinates of the nodes
	double* ycoord;							// Array with Y coordinates of the nodes
	int integer_costs;                      // Flag to choose between using rounded distances (=1) or fractional distances (=0, default)
	int timelimit_exceeded;					// Flag to check if timelimit has been reached while solving the instance (=0, default)

	// Command line parameters list
	int mode;								// Working mode of the program

	// Command line arguments for CREATE_INSTANCES mode
	int n_instances;						// Number of new random instances to generate
	int n_nodes_per_instance;				// Number of nodes (fixed) of the instances to generate

	char folder_istances[100];				// Name of folder inside "/data" where to find the instances that are created in CREATE_INSTANCES mode or tested in RUN_TEST mode
	char instance_prefix_name[20];			// Prefix of each instance file name (used both in CREATE_INSTANCES and RUN_TEST modes)

	// Command line arguments for TEST mode
	char testname[100];
	int n_models_test;
	int* models_to_test;

	// Command line arguments for DEFAULT mode
	int model_type;							// model number to be used to solve the instance
	int tsp_solver;							// model number to be used as TSP solver inside Heuristic methods
	double timelimit;						// overall time limit, in seconds
	char input_file[100];		  			// input file name
	int verbose;							// verbosity value
	int seed;								// internal branching random seed used by Cplex. If fixed, leads to more consistent computational time
	int first_model;
	int last_model;

	//ADV_BRANCH_CUT_PARAM
	int adv_bc_param;

	// global data
	double z_best;							// Value/cost of the objective function to minimize
	double* best_sol;						// Array with the best current integer solution

} instance;

typedef struct {
	instance* inst;
	CPXCALLBACKCONTEXTptr context;
	int* index;
	double* value;
	int local;
} concorde_instance;

#endif   /* INSTANCE_H_ */
