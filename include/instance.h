#ifndef INSTANCE_H_

#define INSTANCE_H_

// Hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ

#define N_MODELS 22
#define N_VERBOSITIES 4
#define N_MODES 3
#define N_ADV_BC_PARAM 5

static const char* models[] = { "BASIC", "MTZ_STATIC", "MTZ_LAZY", "MTZ_SEC2_STATIC", "MTZ_SEC2_LAZY", "GG",
								"BENDERS", "BRANCH_CUT", "ADVBC_STD", "ADVBC_ROOT", "ADVBC_DEPTH_5", "ADVBC_PROB_50", "ADVBC_PROB_10",
								"HEUR_HARD_FIX_50",	"HEUR_HARD_FIX_70",	"HEUR_HARD_FIX_90", "HEUR_HARD_FIX_VAR",
								"HEUR_SOFT_FIX_3", "HEUR_SOFT_FIX_5", "HEUR_SOFT_FIX_7", "HEUR_SOFT_FIX_9", "HEUR_SOFT_FIX_VAR" };
static const char* verbosities[] = { "TEST", "LOW", "MEDIUM", "HIGH" };
static const char* modes[] = { "DEFAULT", "CREATE_INSTANCES", "RUN_TEST" };
static const char* adv_bc_param[] = { "STANDARD", "ROOT_NODE_ONLY", "MAX_DEPTH_5", "PROB_50", "PROB_10" };

// Verbosity levels enumeration
typedef enum {
	TEST,
	LOW,
	MEDIUM,
	HIGH
} verbose;

// Model type enumeration
typedef enum {
	BASIC, // 0
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

typedef enum {
	DEFAULT,
	CREATE_INSTANCES,
	RUN_TEST
} mode;

typedef enum {
	STANDARD,
	ROOT_NODE_ONLY,
	MAX_DEPTH_5,
	PROB_50,
	PROB_10
} adv_branch_cut_param;

// Problem instance data structure
typedef struct {

	// Input data
	int nnodes;
	int ncols;
	char inst_name[100];
	double* xcoord;
	double* ycoord;
	int integer_costs;                      // = 1 for integer costs (rounded distances), 0 otherwise
	int timelimit_exceeded;

	// Command line parameters list
	int mode;								// working mode of the program

	// Command line arguments for CREATE_INSTANCES mode
	int n_instances;						// number of instances to generate
	int n_nodes_per_instance;				// number of nodes (fixed) for the instances to be generated
	
	char folder_istances[100];				// name of folder inside /data where to find the instances to be created in CREATE_INSTANCES mode or tested in RUN_TEST mode
	char instance_prefix_name[20];			// prefix of instance file name (used both in CREATE_INSTANCES and RUN_TEST modes)

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
