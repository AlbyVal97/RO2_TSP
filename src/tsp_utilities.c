#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include <time.h>

#include "tsp_utilities.h"


void parse_command_line(int argc, char** argv, instance* inst) {

	// Set default values for the instance
	inst->model_type = BASIC;
	strcpy(inst->input_file, "NULL");
	inst->timelimit = CPX_INFBOUND;
	inst->integer_costs = 0;
	inst->verbose = MEDIUM;
	inst->seed = 0;

	int help = 0;
	if (argc < 1) help = 1;
	for (int i = 1; i < argc; i++) {

		if (strcmp(argv[i], "-f") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 					// input file
		if (strcmp(argv[i], "-time_limit") == 0) { inst->timelimit = atof(argv[++i]); continue; }			// total time limit
		if (strcmp(argv[i], "-model_type") == 0) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		if (strcmp(argv[i], "-integer_costs") == 0) { inst->integer_costs = 1; continue; } 					// integer costs
		if (strcmp(argv[i], "-verbose") == 0) { inst->verbose = atoi(argv[++i]); continue; } 				// verbose
		if (strcmp(argv[i], "-seed") == 0) { inst->seed = atoi(argv[++i]); continue; } 				// random seed
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; } 											// help
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 										// help
		help = 1;
	}

	if (inst->verbose >= HIGH) { printf("Running %s with %d parameters \n", argv[0], argc - 1); }

	if (help || (inst->verbose >= MEDIUM)) {

		printf("\nAvailable parameters (vers. 05-march-2021) --------------------------------------------------\n");
		printf("-f %s\n", inst->input_file);
		printf("-time_limit %lf\n", inst->timelimit);
		printf("-model_type %d\n", inst->model_type);
		printf("-integer_costs %d\n", inst->integer_costs);
		printf("-verbose %d\n", inst->verbose);
		printf("-seed %d\n", inst->seed);
		printf("\nEnter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}

	if (help) exit(1);
}


void parse_input_file(instance* inst) {

	FILE* fin = fopen(inst->input_file, "r");
	if (fin == NULL) {
		print_error("Input file not found!");
	}

	inst->nnodes = -1;

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

		if (inst->verbose >= HIGH) { printf("Current parameter name: \"%s\" \n", par_name); fflush(NULL); }

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


void print_error(const char* err) { 
	printf("\n\n ERROR: %s \n\n", err);
	fflush(NULL);														// If stream pointer is NULL, all open output streams are flushed
	exit(1);
}


double dist(int i, int j, instance* inst) {
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];
	if (!inst->integer_costs) return sqrt(dx * dx + dy * dy);
	int dis = (int)(sqrt(dx * dx + dy * dy) + 0.499999999); 			// Round to the nearest integer
	return dis + 0.0;
}


void free_instance(instance* inst) { 
	free(inst->xcoord);
	free(inst->ycoord);
}


double second() {
	return ((double)clock() / (double)CLK_TCK);
}