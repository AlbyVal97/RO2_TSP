#ifndef INSTANCE_H_

#define INSTANCE_H_

// Hard-wired parameters
#define XSMALL				1e-5 		// Tolerance used to decide ingerality of 0-1 var.s
#define EPSILON				1e-9		// Very small numerical tolerance 
#define TICKS_PER_SECOND	1000.0  	// Cplex's ticks (depends on CPU)
#define SECONDS_PER_DAY		86400		// Number of seconds in a day

#define N_MODELS			32			// Total number of TSP models/approaches implemented to date = |models[]|
#define N_VERBOSITIES		4			// Total number of verbosity levels considered = |verbosities[]|
#define N_MODES				3			// Total number of working modes of the program = |modes[]|


// List of the TSP models/approaches names implemented to date
static const char* models[] = { "BASIC", "MTZ_STATIC", "MTZ_LAZY", "MTZ_SEC2_STATIC", "MTZ_SEC2_LAZY", "GG",						// Compact models
								"BENDERS",																							// Benders-like solution scheme using SECs
								"BRANCH_CUT", "BRANCH_CUT_2_OPT",																	// SEC separation using Cplex's callbacks
								"ADVBC_STD", "ADVBC_ROOT", "ADVBC_DEPTH_5", "ADVBC_PROB_50", "ADVBC_PROB_10",						// SEC separation using Cplex's callbacks + Concorde utility
								"HEUR_HARD_FIX_50",	"HEUR_HARD_FIX_70",	"HEUR_HARD_FIX_90", "HEUR_HARD_FIX_VAR",					// Hard-fix math heuristics variants
								"HEUR_SOFT_FIX_3", "HEUR_SOFT_FIX_5", "HEUR_SOFT_FIX_7", "HEUR_SOFT_FIX_9",							// Soft-fix (local branching) math heuristics variants
								"HEUR_GREEDY", "HEUR_GRASP_GREEDY", "HEUR_EXTRA_MILEAGE", "HEUR_GRASP_EXTRA_MILEAGE",				// Constructive heuristics
								"HEUR_2_OPT", "HEUR_MULTI_START", 													                // Refinement heuristics
								"HEUR_VNS", "HEUR_TABU", "HEUR_GENETIC", "HEUR_GENETIC_2_OPT" };									// Meta heuristics

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
	BRANCH_CUT_2_OPT, // 8
	ADVBC_STD, // 9
	ADVBC_ROOT, // 10
	ADVBC_DEPTH_5, // 11
	ADVBC_PROB_50, // 12
	ADVBC_PROB_10, // 13
	HEUR_HARD_FIX_50, // 14
	HEUR_HARD_FIX_70, // 15
	HEUR_HARD_FIX_90, // 16
	HEUR_HARD_FIX_VAR, // 17
	HEUR_SOFT_FIX_3, // 18
	HEUR_SOFT_FIX_5, // 19
	HEUR_SOFT_FIX_7, // 20
	HEUR_SOFT_FIX_9, // 21
	HEUR_GREEDY, // 22
	HEUR_GRASP_GREEDY, // 23
	HEUR_EXTRA_MILEAGE, // 24
	HEUR_GRASP_EXTRA_MILEAGE, // 25
	HEUR_2_OPT, // 26
	HEUR_MULTI_START, // 27
	HEUR_VNS, // 28
	HEUR_TABU, // 29
	HEUR_GENETIC, // 30
	HEUR_GENETIC_2_OPT // 31
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

	// Command line arguments for TEST mode
	char testname[100];						// Names of the .bat and .csv files of a single test session 
	int n_models_test;						// Number of models/approaches to test
	int* models_to_test;					// Array of integers (see "model_type" struct) specifying which models/approaches are tested
	int first_model;						// Integer associated to the first model/approach of the test session. The .bat script sets it automatically
	int last_model;							// Integer associated to the last model/approach of the test session. The .bat script sets it automatically

	// Command line arguments both for CREATE_INSTANCES and TEST mode
	char folder_instances[100];				// Name of folder inside "/data" where to find the instances that are created in CREATE_INSTANCES mode or tested in RUN_TEST mode
	char instance_prefix_name[20];			// Prefix of each instance file name. For example, if the instance file name is "rnd_200_1.tsp", then the prefix will be "rnd_200"

	// Command line arguments for DEFAULT mode
	int model_type;							// Integer associated to the model/approach (see "model_type" struct) that is used to solve the instance
	int tsp_solver;							// Integer associated to the model/approach to be used as main TSP solver inside Heuristic methods
	double timelimit;						// The overall time limit (in seconds), that should be set by the user
	char input_file[100];		  			// Path to the input instance file (for example: "../data/test_instances_200/rnd_200_1.tsp")
	int verbose;							// Verbosity level integer value (see "verbose" struct)
	int seed;								// Internal branching integer random seed used by Cplex. When fixed, leads to more consistent computational times
	int n_runs;								// Number of runs to execute: used in (meta)heuristics such as HEUR_GRASP. Of course it may affect the returned solution
	int use_2_opt;							// Flag to indicate if BRANCH_CUT has to make use of 2-opt refinement heuristics or not

	// Global data related to the instance
	double z_best;							// Value/cost of the objective function to minimize
	double* best_sol;						// Array with the best current integer (or almost integer) solution.
											// It can be seen as a list of 1s and 0s indicating if any potential edge of the graph is part of the solution or not, respectively.

} instance;

// Concorde struct used to pass information to the doit_fn_concorde callback (used for ADVBC methods)
typedef struct {
	instance* inst;							// Pointer to the "instance" struct
	CPXCALLBACKCONTEXTptr context;			// Pointer to Cplex callback context, defining the context in which a generic callback is invoked
	int* index;								// Array of indexes associated to the variables of a constraint/cut
	double* value;							// Array of coefficients associated to the variables of a constraint/cut
	int local;								// Flag to indicate if a cut is valid only in current node (=1) or is globally valid (=0, like SECs)
} concorde_instance;

#endif   /* INSTANCE_H_ */
