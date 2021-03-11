#include "tsp.h"

double second();
void parse_command_line(int argc, char** argv, instance *inst);
void read_input(instance *inst);

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }
void free_instance(instance *inst) { free(inst->xcoord); free(inst->ycoord); }

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

	read_input(&inst);

	/* 
	if ( VRPopt(&inst) ) print_error(" error within VRPopt()");
	double t2 = second(); 
    
	if ( VERBOSE >= 1 )   
	{
		printf("... VRP problem solved in %lf sec.s\n", t2-t1);  
	}
	*/

	free_instance(&inst);

	return 0; 
}        

void parse_command_line(int argc, char** argv, instance *inst) { 
	
	if ( VERBOSE >= 100 ) {
		printf("Running %s with %d parameters \n", argv[0], argc - 1);
	}
		
	// default   
	//inst->model_type = 0;
	//inst->old_benders = 0;
	strcpy(inst->input_file, "NULL");
	strcpy(inst->node_file, "NULL");
	//inst->randomseed = 0; 
	//inst->num_threads = 0;
	inst->timelimit = 99999999;
	//inst->timelimit = CPX_INFBOUND;
	//inst->cutoff = CPX_INFBOUND; 
	//inst->integer_costs = 0;

	//inst->available_memory = 12000;   			// available memory, in MB, for Cplex execution (e.g., 12000)
	//inst->max_nodes = -1; 						// max n. of branching nodes in the final run (-1 unlimited)        

    int help = 0;
	if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) {

		//if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		//if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file, argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		//if ( strcmp(argv[i],"-model_type") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 	// model type
		//if ( strcmp(argv[i],"-model") == 0 ) { inst->model_type = atoi(argv[++i]); continue; } 			// model type
		//if ( strcmp(argv[i],"-old_benders") == 0 ) { inst->old_benders = atoi(argv[++i]); continue; } 	// old benders
		//if ( strcmp(argv[i],"-seed") == 0 ) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		//if ( strcmp(argv[i],"-threads") == 0 ) { inst->num_threads = atoi(argv[++i]); continue; } 		// n. threads
		//if ( strcmp(argv[i],"-memory") == 0 ) { inst->available_memory = atoi(argv[++i]); continue; }	// available memory (in MB)
		if ( strcmp(argv[i],"-node_file") == 0 ) { strcpy(inst->node_file,argv[++i]); continue; }		// cplex's node file
		//if ( strcmp(argv[i],"-max_nodes") == 0 ) { inst->max_nodes = atoi(argv[++i]); continue; } 		// max n. of nodes
		//if ( strcmp(argv[i],"-cutoff") == 0 ) { inst->cutoff = atof(argv[++i]); continue; }				// master cutoff
		//if ( strcmp(argv[i],"-int") == 0 ) { inst->integer_costs = 1; continue; } 						// integer costs
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }      

	if ( help || (VERBOSE >= 1000) )		// print current parameters
	{
		printf("\n\nAvailable parameters (vers. 05-march-2021) --------------------------------------------------\n");
		printf("-f %s\n", inst->input_file); 
		printf("-time_limit %lf\n", inst->timelimit); 
		//printf("-model_type %d\n", inst->model_type); 
		//printf("-old_benders %d\n", inst->old_benders); 
		//printf("-seed %d\n", inst->randomseed); 
		//printf("-threads %d\n", inst->num_threads);  
		//printf("-max_nodes %d\n", inst->max_nodes); 
		//printf("-memory %d\n", inst->available_memory); 
		//printf("-int %d\n", inst->integer_costs); 
		printf("-node_file %s\n", inst->node_file);
		//printf("-cutoff %lf\n", inst->cutoff); 
		printf("\nEnter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}

	if ( help ) exit(1);
}

void read_input(instance *inst) {
                            
	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) {
		print_error("Input file not found!");
	}

	inst->nnodes = -1;

	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	int debug_mode = ( VERBOSE >= 1000 );

	while ( fgets(line, sizeof(line), fin) != NULL ) { // fgets returns NULL when error or EOF occurs

		// Skip empty lines
		if ( strlen(line) <= 1 ) continue; 

		// Gets the portion of line before separators " " and ":", that is the parameter name. 
		// Then replaces \0 to the first found separator
	    par_name = strtok(line, " :\n");

		if ( debug_mode ) { printf("Current parameter name: \"%s\" \n", par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 ) { continue; }

		if ( strncmp(par_name, "COMMENT", 7) == 0 ) {
			// NULL is needed to parse the same line as before
			token1 = strtok(NULL, ":\n");
			if ( debug_mode ) { printf("Solving instance \"%s\" \n", token1+1); }
			//if ( VERBOSE >= 10 ) { printf("Solving instance %s with model %d\n\n", token1, inst->model_type); }
			continue;
		}   
		
		if ( strncmp(par_name, "TYPE", 4) == 0 ) {
			token1 = strtok(NULL, " :\n");  
			if ( strncmp(token1, "TSP", 3) != 0 ) { print_error("Format error: only TSP algorithm implemented yet."); }
			if ( debug_mode ) { printf("Problem type: \"%s\" \n", token1); }
			continue;
		}

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) {
			if ( inst->nnodes >= 0 ) { print_error("Number of nodes was alreay initialized before DIMENSION section in input file."); }
			token1 = strtok(NULL, " :\n");
			inst->nnodes = atoi(token1);
			if ( debug_mode ) { printf("Number of nodes: %d \n", inst->nnodes); }
			// Allocates memory for the arrays of coordinates.
			// Their size is equal to (number of nodes)*(number of bytes per node coordinate)
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double));
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
			continue;
		}

		if ( strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0 ) {
			token1 = strtok(NULL, " :\n");
			if ( strncmp(token1, "EUC_2D", 6) != 0 ) { print_error("Format error: only EUC_2D weight implemented yet."); }
			if ( debug_mode ) { printf("Edge weight type:: \"%s\" \n", token1); }
			continue;
		}            
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) {
			if ( inst->nnodes <= 0 ) { print_error("Number of nodes should have been already initialized before NODE_COORD_SECTION."); }
			active_section = 1;
			continue;
		}

		if ( strncmp(par_name, "EOF", 3) == 0 ) {
			active_section = 0;
			break;
		}
		
		// To be executed only if we are inside NODE_COORD_SECTION
		if ( active_section == 1 ) {
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) { print_error("In NODE_COORD_SECTION a node has an illegal index."); } 
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if ( debug_mode ) { printf("Node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]); }
			continue;
		}
		
		// If we are here, something has gone wrong with parsing
		print_error("Wrong format for this input file.");     
		    
	}                

	fclose(fin);    
	
}