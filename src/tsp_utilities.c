#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include <time.h>

#include "tsp_utilities.h"


void parse_command_line(int argc, char** argv, instance* inst) {

	// Set default values for the instance fields
	inst->mode = DEFAULT;
	inst->model_type = BASIC;
	inst->verbose = MEDIUM;
	inst->tsp_solver = ADVBC_ROOT;
	inst->n_instances = 0;
	inst->n_nodes_per_instance = 0;
	strcpy(inst->input_file, "NULL\0");
	strcpy(inst->folder_instances, "NULL\0");
	strcpy(inst->instance_prefix_name, "NULL\0");
	strcpy(inst->testname, "NULL\0");
	inst->n_models_test = 0;
	inst->models_to_test = NULL;
	inst->timelimit = CPX_INFBOUND;
	inst->timelimit_exceeded = 0;
	inst->integer_costs = 0;
	inst->first_model = 0;
	inst->last_model = 0;
	inst->seed = 123456;
	inst->n_runs = 10000;
	inst->z_best = INFINITY;
	inst->use_2_opt = 0;

	// Chech if at least one command line argument is provided
	int help = 0;
	if (argc < 1) help = 1;

	// Set the instance fields values according to the command line arguments provided by the user
	// For further information about each argument, see the correspondent inst->* field description inside the definition of "instance" struct
	for (int i = 1; i < argc; i++) {

		if (strcmp(argv[i], "-end") == 0) break;
		if (strcmp(argv[i], "-m") == 0) {
			inst->mode = atoi(argv[++i]); 
			if (inst->mode < 0 || inst->mode >= N_MODES) print_error("Incorrect value of -m (working mode)!");
			continue;
		}							
		if (strcmp(argv[i], "-folder") == 0) { 
			strcpy(inst->folder_instances, argv[++i]);
			if (strlen(inst->folder_instances) > 100) print_error("Argument for -folder is too long!");
			continue;
		}
		if (strcmp(argv[i], "-prefix") == 0) {
			strcpy(inst->instance_prefix_name, argv[++i]);
			if (strlen(inst->instance_prefix_name) > 20) print_error("Argument for -prefix is too long!");
			continue;
		}
		if (strcmp(argv[i], "-test") == 0) {
			strcpy(inst->testname, argv[++i]);
			if (strlen(inst->testname) > 100) print_error("Argument for -test is too long!");
			continue;
		}
		if (strcmp(argv[i], "-n_inst") == 0) {
			inst->n_instances = atoi(argv[++i]);
			if (inst->n_instances <= 0 || inst->n_instances > 100) print_error("Values of -n_inst must an integer in the interval [1, 100]!");
			continue;
		}
		if (strcmp(argv[i], "-n_nodes") == 0) {
			inst->n_nodes_per_instance = atoi(argv[++i]);
			if (inst->n_nodes_per_instance > 10000) print_error("Values of -n_nodes greater than 10000 nodes are not supported!");
			continue;
		}
		if (strcmp(argv[i], "-f") == 0) {
			strcpy(inst->input_file, argv[++i]);
			if (strlen(inst->input_file) > 100) print_error("Argument for -f (path of input instance) is too long!");
			continue;
		}
		if (strcmp(argv[i], "-time_limit") == 0) {
			inst->timelimit = atof(argv[++i]);
			if (inst->timelimit > SECONDS_PER_DAY) print_error("Values of -time_limit greater than 24 hours are not supported!");
			continue;
		}
		if (strcmp(argv[i], "-model_type") == 0) {
			inst->model_type = atoi(argv[++i]);
			if (inst->model_type < 0 || inst->model_type >= N_MODELS) print_error("Incorrect value of -model_type!");
			continue;
		} 			
		if (strcmp(argv[i], "-integer_costs") == 0) { inst->integer_costs = 1; continue; }
		if (strcmp(argv[i], "-first") == 0) {
			inst->first_model = atoi(argv[++i]);
			if (inst->first_model < 0 || inst->first_model >= N_MODELS) print_error("Incorrect value of -first!");
			continue;
		}
		if (strcmp(argv[i], "-last") == 0) {
			inst->last_model = atoi(argv[++i]);
			if (inst->last_model < 0 || inst->last_model >= N_MODELS) print_error("Incorrect value of -last!");
			continue;
		}
		if (strcmp(argv[i], "-verbose") == 0) {
			inst->verbose = atoi(argv[++i]);
			if (inst->verbose < 0 || inst->verbose >= N_VERBOSITIES) print_error("Incorrect value of -verbose!");
			continue;
		} 				
		if (strcmp(argv[i], "-seed") == 0) {
			inst->seed = atoi(argv[++i]);
			if (inst->seed < 0 || inst->seed > CPX_BIGINT) print_error("Value of -seed (Cplex internal seed) must an integer in the interval [0, BIGINT]!");
			continue;
		}
		if (strcmp(argv[i], "-n_runs") == 0) {
			inst->n_runs = atoi(argv[++i]);
			if (inst->n_runs < 0 || inst->n_runs > 10000000) print_error("Value of -n_runs for (meta) heuristics must an integer in the interval [0, 10000000]!");
			continue;
		}
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; }
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; }

		if (strcmp(argv[i], "-n_models") == 0) {
			inst->n_models_test = atoi(argv[++i]);
			if (inst->n_models_test <= 0 || inst->n_models_test > N_MODELS) print_error("Incorrect value of -n_models!");
			inst->models_to_test = (int*)malloc(inst->n_models_test * sizeof(int));
			if (inst->models_to_test != NULL) {
				for (int j = 0; j < inst->n_models_test; j++) {
					int temp = atoi(argv[++i]);
					if (temp < 0 || temp >= N_MODELS) print_error("Incorrect value of at least one model in the list of -n_models!");
					inst->models_to_test[j] = temp;
				}
			}
			else {
				print_error("Error when allocating memory for models_to_test!");
			}
			continue;
		}
		help = 1;
	}

	if (inst->verbose >= HIGH) { printf("Running %s with %d parameters \n", argv[0], argc - 1); }

	// Print a summary of the command line parameters provided by the user before starting to solve the instance
	if (help || ((inst->verbose >= LOW) && inst->mode != RUN_TEST)) {

		printf("\nAvailable parameters (vers. 05-march-2021) --------------------------------------------------\n");
		printf("-m %s\n", modes[inst->mode]);
		printf("-folder %s\n", inst->folder_instances);
		printf("-prefix %s\n", inst->instance_prefix_name);
		printf("-n_inst %d\n", inst->n_instances);
		printf("-n_nodes %d\n", inst->n_nodes_per_instance);
		printf("-f %s\n", inst->input_file);
		printf("-time_limit %lf\n", inst->timelimit);
		printf("-model_type %s\n", models[inst->model_type]);
		printf("-integer_costs %d\n", inst->integer_costs);
		printf("-verbose %s\n", verbosities[inst->verbose]);
		printf("-seed %d\n", inst->seed);
		printf("-n_runs %d\n", inst->n_runs);
		printf("\nEnter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n");
	}
	else if (inst->verbose == TEST) {
		printf("-f %s\n", inst->input_file);
		printf("-time_limit %lf\n", inst->timelimit);
		printf("-model_type %s\n", models[inst->model_type]);
		printf("-test %s\n", inst->testname);
	}

	if (help) exit(1);
}


void parse_input_file(instance* inst) {

	FILE* fin = fopen(inst->input_file, "r");
	if (fin == NULL) {
		print_error("Input file not found!");
	}

	inst->nnodes = -1;
	inst->ncols = -1;

	char line[180];
	char* par_name;
	char* token1;
	char* token2;

	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 

	int debug_mode = (inst->verbose >= MEDIUM);

	while (fgets(line, sizeof(line), fin) != NULL) { // fgets returns NULL when error or EOF occurs

		// Skip empty lines
		if (strlen(line) <= 1) continue;

		// Gets the portion of line before separators " " and ":", that is the parameter name. 
		// Then replaces \0 to the first found separator
		par_name = strtok(line, " :\n");

		if (strncmp(par_name, "NAME", 4) == 0) {
			token1 = strtok(NULL, " :\n");
			strcpy(inst->inst_name, token1);
			if (debug_mode) { printf("Name of the instance: \"%s\" \n", token1); }
			continue;
		}

		if (strncmp(par_name, "COMMENT", 7) == 0) {
			token1 = strtok(NULL, ":\n");			// NULL is needed to parse the same line as before
			if (debug_mode) { printf("Solving instance \"%s\" with model \"%d\" \n", token1+1, inst->model_type); }
			continue;
		}

		if (strncmp(par_name, "TYPE", 4) == 0) {
			token1 = strtok(NULL, " :\n");
			if (strncmp(token1, "TSP", 3) != 0) { print_error("Format error: only TSP algorithm implemented yet."); }
			if (debug_mode) { printf("Problem type: \"%s\" \n", token1); }
			continue;
		}

		if (strncmp(par_name, "DIMENSION", 9) == 0) {
			if (inst->nnodes >= 0) { print_error("Number of nodes was alreay initialized before DIMENSION section in input file."); }
			token1 = strtok(NULL, " :\n");
			inst->nnodes = atoi(token1);
			if (debug_mode) { printf("Number of nodes: %d \n", inst->nnodes); }
			// Allocates memory for the arrays of coordinates.
			// Their size is equal to (number of nodes)*(number of bytes per node coordinate)
			inst->xcoord = (double*)calloc(inst->nnodes, sizeof(double));
			inst->ycoord = (double*)calloc(inst->nnodes, sizeof(double));
			continue;
		}

		if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0) {
			token1 = strtok(NULL, " :\n");
			if ((strncmp(token1, "EUC_2D", 6) != 0) && (strncmp(token1, "ATT", 3) != 0)) { print_error("Format error: only EUC_2D and ATT weight implemented yet."); }
			if (debug_mode) { printf("Edge weight type: \"%s\" \n", token1); }
			continue;
		}

		if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0) {
			if (inst->nnodes <= 0) { print_error("Number of nodes should have been already initialized before NODE_COORD_SECTION."); }
			active_section = 1;
			continue;
		}

		if (strncmp(par_name, "EOF", 3) == 0) {
			active_section = 0;
			break;
		}

		// To be executed only if we are inside NODE_COORD_SECTION
		if (active_section == 1) {
			int i = atoi(par_name) - 1;
			if (i < 0 || i >= inst->nnodes) { print_error("In NODE_COORD_SECTION a node has an illegal index."); }
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if (debug_mode) { printf("Node %4d coordinates: ( %15.7lf , %15.7lf )\n", i + 1, inst->xcoord[i], inst->ycoord[i]); }
			continue;
		}

		// If we are here, something has gone wrong with parsing
		print_error("Wrong format for this input file.");
	}

	fclose(fin);
}


void print_nodes_dat_file(const instance* inst) {

	FILE* nodes_file = fopen("nodes.dat", "w");

	if (nodes_file == NULL) {
		print_error("nodes.dat not found!");
	}

	for (int i = 0; i < inst->nnodes; ++i) {
		fprintf(nodes_file, "%f %f\n", inst->xcoord[i], inst->ycoord[i]);
	}

	fclose(nodes_file);

	return;
}


void create_instances(instance* inst) {

	char folder_path[50];
	sprintf(folder_path, "../data/%s", inst->folder_instances);
	if (mkdir(folder_path) == -1) printf("Folder \"%s\" already exists.\n", inst->folder_instances);
	else printf("Folder \"%s\" created for the first time!\n", inst->folder_instances);

	// Generate a list of seeds, one per instance
	int* seeds = (int*)malloc(sizeof(int) * inst->n_instances);
	int base_seed = 1236460 - 1236009;
	for (int i = 0; i < inst->n_instances; i++) {
		seeds[i] = base_seed + i * 9;
	}
	
	char inst_name[20];
	char inst_file_path[50];
	for (int i = 0; i < inst->n_instances; i++) {

		// Create the name for the current instance file
		sprintf(inst_name, "%s_%d", inst->instance_prefix_name, i+1);
		sprintf(inst_file_path, "../data/%s/%s.tsp", inst->folder_instances, inst_name);

		// Open the file for the current instance
		FILE* fin = fopen(inst_file_path, "w");
		if (fin == NULL) print_error("Input file not found!");

		fprintf(fin, "NAME : %s\n", inst_name);
		fprintf(fin, "COMMENT : Instance #%d out of %d with seed %d.\n", i+1, inst->n_instances, seeds[i]);
		fprintf(fin, "TYPE : TSP\n");
		fprintf(fin, "DIMENSION : %d\n", inst->n_nodes_per_instance);
		fprintf(fin, "EDGE_WEIGHT_TYPE : EUC_2D\n");
		fprintf(fin, "NODE_COORD_SECTION\n");

		// Set the seed for the current instance
		srand(seeds[i]);

		// Add the nodes coordinates (one node per row)
		for (int j = 0; j < inst->n_nodes_per_instance; j++) {
			double x_coord = (double)rand() / RAND_MAX;
			double y_coord = (double)rand() / RAND_MAX;

			fprintf(fin, "%d %f %f\n", j+1, x_coord, y_coord);
		}

		fprintf(fin, "EOF");

		fclose(fin);
	}

	return;
}


void print_error(const char* err) { 
	printf("\n\n ERROR: %s \n\n", err);
	fflush(NULL);											// If stream pointer is NULL, all open output streams are flushed
	exit(1);
}


double dist(int i, int j, instance* inst) {
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];
	if (!inst->integer_costs) return sqrt(dx * dx + dy * dy);
	
	// Only if integer_costs is 1, then we need to round to the nearest integer
	int dis = (int)(sqrt(dx * dx + dy * dy) + 0.499999999);
	return dis + 0.0;
}


void free_instance(instance* inst) {

	free(inst->xcoord);
	free(inst->ycoord);
	// N.B. only HEUR_HARD_FIX_* and HEUR_SOFT_FIX_K initialize "inst->model_type" !
	if (inst->model_type > ADVBC_PROB_10 && inst->model_type < HEUR_GREEDY) free(inst->best_sol);

	return;
}


double second() {
	return ((double)clock() / (double)CLK_TCK);
}