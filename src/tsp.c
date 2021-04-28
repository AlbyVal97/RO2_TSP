#define _CRT_RAND_S
#define safe_rand(rval, value) \
	value = 0; \
	if (rand_s(&value)) \
		print_error("Error in random generation\n"); \
	*rval = ((double)value) / (UINT_MAX); 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
//#include <gnuplot_c.h>

#include "tsp.h"


int TSPopt(instance* inst) {

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	// Set the timelimit parameter according to user input or default value.
	if (CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit)) print_error("CPXsetdblparam() error in setting timelimit");

	// Set the seed parameter according to user input or default value.
	if (CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed)) print_error("CPXsetintparam() error in setting seed");

	char edges_file_path[100];
	char logfile_path[100];
	// Create "outputs" folder if needed
	printf("\n");
	if (mkdir("../outputs") == -1) {
		if (inst->verbose >= MEDIUM) printf("Folder \"outputs\" already exists.\n");
	}
	else {
		if (inst->verbose >= MEDIUM) printf("Folder \"outputs\" created for the first time!\n");
	}
	sprintf(edges_file_path, "../outputs/%s", inst->inst_name);
	sprintf(logfile_path, "../outputs/%s", inst->inst_name);
	if (inst->verbose != TEST) {
		if (mkdir(edges_file_path) == -1) {
			if (inst->verbose >= MEDIUM) printf("Folder for the current instance already exists.\n");
		}
		else {
			if (inst->verbose >= MEDIUM) printf("Folder for the current instance created for the first time!\n");
		}
	}

	switch (inst->model_type) {

		case BASIC:
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			printf("inst->ncols: %d\n", inst->ncols);

			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_BASIC.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			sprintf(edges_file_path, "%s/model_BASIC_edges.dat", edges_file_path);
			break;

		case MTZ_STATIC:
			build_model_MTZ_STATIC(inst, env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_MTZ_STATIC.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) { print_error("CPXsetdblparam() error in setting integer value tolerance"); }
			if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			sprintf(edges_file_path, "%s/model_MTZ_STATIC_edges.dat", edges_file_path);
			break;

		case MTZ_LAZY:
			build_model_MTZ_LAZY(inst, env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_MTZ_LAZY.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) { print_error("CPXsetdblparam() error in setting integer value tolerance"); }
			if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			sprintf(edges_file_path, "%s/model_MTZ_LAZY_edges.dat", edges_file_path);
			break;

		case MTZ_SEC2_STATIC:
			build_model_MTZ_SEC2_STATIC(inst, env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_MTZ_SEC2_STATIC.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) { print_error("CPXsetdblparam() error in setting integer value tolerance"); }
			if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			sprintf(edges_file_path, "%s/model_MTZ_SEC2_STATIC_edges.dat", edges_file_path);
			break;

		case MTZ_SEC2_LAZY:
			build_model_MTZ_SEC2_LAZY(inst, env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_MTZ_SEC2_LAZY.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) { print_error("CPXsetdblparam() error in setting integer value tolerance"); }
			if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			sprintf(edges_file_path, "%s/model_MTZ_SEC2_LAZY_edges.dat", edges_file_path);
			break;

		case GG:
			build_model_GG(inst, env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_GG.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			sprintf(edges_file_path, "%s/model_GG_edges.dat", edges_file_path);
			break;

		case BENDERS:
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_BENDERS.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			solve_benders(inst, env, lp);
			sprintf(edges_file_path, "%s/model_BENDERS_edges.dat", edges_file_path);
			break;

		case BRANCH_CUT:
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_BRANCH_CUT.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			solve_branch_cut(inst, env, lp);
			sprintf(edges_file_path, "%s/model_BRANCH_CUT_edges.dat", edges_file_path);
			break;
		
		case HEUR_HARD_FIX:
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_HEUR_HARD_FIX.txt", logfile_path);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			solve_heur_hard_fix(inst, env, lp);
			sprintf(edges_file_path, "%s/model_HEUR_HARD_FIX_edges.dat", edges_file_path);
			break;

		//case ADV_BRANCH_CUT:
		default:
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			if (inst->verbose >= LOW) {
				sprintf(logfile_path, "%s/logfile_%s.txt", logfile_path, models[inst->model_type]);
				if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");
			}
			solve_adv_branch_cut(inst, env, lp);
			sprintf(edges_file_path, "%s/model_%s_edges.dat", edges_file_path, models[inst->model_type]);
			break;

		/*default:
			print_error("Choose a correct value for the model to be used!");*/
	}

	if (inst->verbose >= MEDIUM) printf("Complete path for *.dat file: %s\n\n", edges_file_path);
	
	// Allocate memory for the optimal solution array
	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*)calloc(ncols, sizeof(double));

	// Copy the optimal solution from the Cplex environment to the new array "xstar"
	if (!inst->timelimit_exceeded && CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error");

	// Discern if a model is symmetric (solves the TSP for a directed or an undirected graph)
	int symmetric = -1;
	if (inst->model_type == BASIC || inst->model_type == BENDERS || inst->model_type == BRANCH_CUT || inst->model_type >= ADVBC_STD) symmetric = 0;
	else if (inst->model_type == MTZ_STATIC || inst->model_type == MTZ_LAZY ||inst->model_type == MTZ_SEC2_STATIC || inst->model_type == MTZ_SEC2_LAZY || inst->model_type == GG) symmetric = 1;
	
	// Fill the .dat file with the correctly formatted nodes of the found solution
	if (inst->verbose >= LOW) print_solution(inst, xstar, symmetric, edges_file_path);
	
	// Free allocated memory and close Cplex model
	free(xstar);

	if (CPXsetlogfilename(env, NULL, NULL)) { print_error("CPXsetlogfilename() error"); }
	if (CPXfreeprob(env, &lp)) { print_error("CPXfreeprob() error"); }
	if (CPXcloseCPLEX(&env)) { print_error("CPXcloseCPLEX() error"); }

	return 0;
}


void print_solution(instance* inst, double* xstar, int symmetric, char* edges_file_path) {

	if (symmetric != 0 && symmetric != 1) print_error("symmetric is a boolean variable\n");
	
	FILE* edges_plot_file_name = fopen(edges_file_path, "w");
	if (edges_plot_file_name == NULL) print_error("File *_edges.dat not found!");

	if (!symmetric) {
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = i + 1; j < inst->nnodes; j++) {

				if (xstar[xpos(i, j, inst)] > 0.5) {
				//if (xstar[xpos(i, j, inst)] > 0.0) {

					if (inst->verbose >= HIGH) {
						printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
					}
					fprintf(edges_plot_file_name, "%f %f\n%f %f\n\n", inst->xcoord[i], inst->ycoord[i], inst->xcoord[j], inst->ycoord[j]);
				}
			}
		}
	}
	else {
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = 0; j < inst->nnodes; j++) {
				if (i == j) continue;

				if (xstar[xpos_compact(i, j, inst)] > 0.5) {

					if (inst->verbose >= HIGH) {
						printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
					}
					fprintf(edges_plot_file_name, "%f %f %f %f\n", inst->xcoord[i], inst->ycoord[i], inst->xcoord[j] - inst->xcoord[i], inst->ycoord[j] - inst->ycoord[i]);
				}
			}
		}
	}

	FILE* gn_com = fopen("../outputs/gnuplot_commands.txt", "w");
	char* comand;
	if (gn_com == NULL) print_error("Error while opening gnuplot_commands.txt file\n");
	if (!symmetric) 
		comand = "set style line 1 \\\n\tlinecolor rgb '#FF0000' \\\n\tlinetype 1 linewidth 2 \\\n\tpointtype 7 pointsize 2 \\\n\nplot \"%s\" using 1:2 with linespoints linestyle 1\npause mouse close";
	else 
		comand = "set style line 1 \\\n\tlinecolor rgb '#FF0000' \\\n\tlinetype 1 linewidth 2 \\\n\tpointtype 7 pointsize 2 \\\n\nplot \"%s\" using 1:2:3:4 with vectors linestyle 1\npause mouse close";
	
	fprintf(gn_com, comand, edges_file_path);
	
	fclose(gn_com);
	fclose(edges_plot_file_name);
	
}


int xpos(int i, int j, instance* inst) {

	if (i == j) print_error(" i == j in xpos");
	if (i > j) return xpos(j, i, inst);
	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
	return pos;

	return 0;
}


int xpos_compact(int i, int j, instance* inst) {

	int pos = i * inst->nnodes + j;
	return pos;
}


int upos_compact(int i, instance* inst) {

	return (inst->nnodes) * (inst->nnodes) + i - 1;
}


int ypos_compact(int i, int j, instance* inst) {

	int pos = (inst->nnodes) * (inst->nnodes) + i * (inst->nnodes) + j;

	return pos;
}


int mip_solved_to_optimality(instance* inst, CPXENVptr env, CPXLPptr lp) {

	int lpstat = CPXgetstat(env, lp);
	switch (lpstat) {
		case CPXMIP_OPTIMAL:
			if (inst->verbose >= MEDIUM) printf("\nOptimal integer solution found.\n");
			break;
		case CPXMIP_OPTIMAL_INFEAS:
			if (inst->verbose >= MEDIUM) printf("\nProblem optimal with unscaled infeasibilities.\n");
			break;
		case CPXMIP_OPTIMAL_TOL:
			if (inst->verbose >= MEDIUM) printf("\nOptimal solution within epgap or epagap tolerance found.\n");
			break;
		case CPXMIP_TIME_LIM_FEAS:
			inst->timelimit_exceeded = 1;
			if (inst->verbose >= LOW) print_error("Time limit exceeded, integer solution exists.\n");
			break;
		case CPXMIP_TIME_LIM_INFEAS:
			inst->timelimit_exceeded = 1;
			if (inst->verbose >= LOW) print_error("Time limit exceeded, no integer solution.\n");
			break;
	}
	
	int solved = (lpstat == CPXMIP_OPTIMAL) ||
		(lpstat == CPXMIP_OPTIMAL_INFEAS) ||
		//( lpstat ==  CPXMIP_OPTIMAL_RELAXED ) ||
		(lpstat == CPXMIP_OPTIMAL_TOL);

	return solved;
}


void solve_benders(instance* inst, CPXENVptr env, CPXLPptr lp) {

	int n_comp = 2;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));	
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));
	double residual_timelimit = inst->timelimit;

	int ncols = inst->ncols;
	double* x = (double*)calloc(ncols, sizeof(double));
	
	int n_iter = 0;
	while (n_comp > 1) {									// Repeat iteratively until just one connected component is left

		if (n_iter > 0) {									// For the first iteration just solve the BASIC_MODEL
			// Add a new subtour elimination constraint for each connected components of the current solution
			update_benders_constraints(env, lp, inst, comp, n_comp, n_iter);
		}

		double t1 = second();
		// Optimize with the new constraint
		if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
		double t2 = second();
		
		// Check if CPXmipopt has ended correctly
		mip_solved_to_optimality(inst, env, lp);

		if (inst->verbose >= MEDIUM) printf("Time used for iteration number %d: %f\n", n_iter, t2 - t1);

		// Update the amount of time left before timelimit is reached and provide it to Cplex to check
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (inst->verbose >= MEDIUM) printf("New time limit: %f\n\n", residual_timelimit);
		if (CPXsetdblparam(env, CPX_PARAM_TILIM, residual_timelimit)) { print_error("CPXsetdblparam() error in setting timelimit"); }

		// Extract the new solution (but only if timelimit has not been reached)
		if (!inst->timelimit_exceeded && CPXgetx(env, lp, x, 0, ncols - 1)) { print_error("CPXgetx() error"); }

		// Update the number of connected components of the new graph
		update_connected_components(x, inst, succ, comp, &n_comp);
		if (inst->verbose >= MEDIUM) printf("Current number of connected components (n_comp): %d \n", n_comp);
		if (inst->verbose >= MEDIUM && n_comp == 1) printf("BENDERS method ended successfully!\n");

		// Increment current iteration number
		n_iter++;
	}

	if (inst->verbose >= LOW) create_lp_file(inst, env, lp, "model_BENDERS");

	free(succ);
	free(comp);
	free(x);
}


void update_benders_constraints(CPXCENVptr env, CPXLPptr lp, instance* inst, const int* comp, int n_comp, int n_iter) {
						
	int* index = (int*)calloc(inst->ncols, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->ncols, sizeof(double));			// Array of row variables coefficients
	int* comp_nodes = (int*)malloc(inst->nnodes * sizeof(int));				// Array of indexes of the nodes for one component

	for (int n = 1; n <= n_comp; n++) {					// n is the index of the current connected component

		// Scan connected components for its number of nodes -> value of comp_n_nodes
		int comp_n_nodes = 0;
		for (int i = 0; i < inst->nnodes; i++) {
			if (comp[i] == n) {
				comp_n_nodes++;
			}
		}

		// Scan connected components for its nodes indexes
		int j = 0;
		for (int i = 0; i < inst->nnodes; i++) {
			if (comp[i] == n) {
				comp_nodes[j++] = i;
			}
		}

		char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
		cname[0] = (char*)calloc(100, sizeof(char));
		sprintf(cname[0], "subtour_elimination_constraint(%d,%d)", n_iter + 1, n);
		double rhs = comp_n_nodes - 1.0;
		char sense = 'L';
		int i_zero = 0;
		// We should compute the number of edges connecting the nodes of the current connected component
		// as the binomial coefficient (n over 2) = n!/((n-2)!*(2!)) = n*(n-1)/2
		int n_edges_curr_comp = comp_n_nodes * (comp_n_nodes - 1) / 2;
		if (inst->verbose >= HIGH) printf("Number of edges in connected component #%d (n_edges_curr_comp): %d\n", n, n_edges_curr_comp);
		int nnz = n_edges_curr_comp;

		// Build (and add to the model) the subtour elimination constraint for the current connected component
		int k = 0;
		for (int i = 0; i < comp_n_nodes; i++) {
			for (int j = i + 1; j < comp_n_nodes; j++) {
				if (i == j) {
					continue;
				}
				index[k] = xpos(comp_nodes[i], comp_nodes[j], inst);
				value[k++] = 1.0;
			}
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for subtour elimination constraints!");

		free(cname[0]);
		free(cname);																	// Go to the next connected component to add more constraints
	}

	free(index);
	free(value);
	free(comp_nodes);
}


void solve_branch_cut(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// N.B. It is required to install a "lazy constraint" callback to cut infeasible integer solutions (es. found by heuristics) 
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
	if (CPXcallbacksetfunc(env, lp, contextid, branch_cut_callback, inst)) print_error("CPXcallbacksetfunc() error");

	if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

	// Check if CPXmipopt has ended correctly
	mip_solved_to_optimality(inst, env, lp);

	return;
}


static int CPXPUBLIC branch_cut_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {

	instance* inst = (instance*)userhandle;
	int ncols = inst->ncols;
	double* xstar = (double*)malloc(ncols * sizeof(double));
	double objval = CPX_INFBOUND;
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval)) print_error("CPXcallbackgetcandidatepoint error");

	// get some random information at the node (as an example for the students)
	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
	int mynode = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
	double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
	if (inst->verbose >= HIGH) {
		printf("mythread: %d\n", mythread);
		printf("mynode: %d\n", mynode);
		printf("incumbent: %f\n", incumbent);
	}

	int n_comp = 0;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));

	update_connected_components(xstar, inst, succ, comp, &n_comp);

	if (n_comp > 1) {								// means that the solution is infeasible and a violated cut has been found

		int* index = (int*)calloc(ncols, sizeof(int));					// Array of indexes associated to the row variables
		double* value = (double*)calloc(ncols, sizeof(double));			// Array of row variables coefficients
		int* comp_nodes = (int*)malloc(inst->nnodes * sizeof(int));		// Array of indexes of the nodes for one component
								
		for (int n = 1; n <= n_comp; n++) {			// n is the index of the current connected component

			// Scan connected components for its number of nodes -> value of comp_n_nodes
			int comp_n_nodes = 0;
			for (int i = 0; i < inst->nnodes; i++) {
				if (comp[i] == n) {
					comp_n_nodes++;
				}
			}

			// Scan connected components for its nodes indexes
			int j = 0;
			for (int i = 0; i < inst->nnodes; i++) {
				if (comp[i] == n) {
					comp_nodes[j++] = i;
				}
			}

			double rhs = comp_n_nodes - 1.0;
			char sense = 'L';
			int i_zero = 0;
			// We should compute the number of edges connecting the nodes of the current connected component
			// as the binomial coefficient (n over 2) = n!/((n-2)!*(2!)) = n*(n-1)/2
			int n_edges_curr_comp = comp_n_nodes * (comp_n_nodes - 1) / 2;
			int nnz = n_edges_curr_comp;

			// Build (and add to the model) the subtour elimination constraint for the current connected component
			int k = 0;
			for (int i = 0; i < comp_n_nodes; i++) {
				for (int j = i + 1; j < comp_n_nodes; j++) {
					if (i == j) {
						continue;
					}
					index[k] = xpos(comp_nodes[i], comp_nodes[j], inst);
					value[k++] = 1.0;
				}
			}
			// Rejects the current (infeasible) solution and adds one cut (a SEC)
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &i_zero, index, value)) print_error("CPXcallbackrejectcandidate() error");
		}

		free(index);
		free(value);
		free(comp_nodes);
	}

	free(xstar);
	free(succ);
	free(comp);

	return 0;
}


void solve_adv_branch_cut(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// N.B. It is required to install a "lazy constraint" callback to cut infeasible integer solutions (es. found by heuristics) 
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_RELAXATION;
	if (CPXcallbacksetfunc(env, lp, contextid, adv_branch_cut_callback_driver, inst)) print_error("CPXcallbacksetfunc() error");

	if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");

	// Check if CPXmipopt has ended correctly
	mip_solved_to_optimality(inst, env, lp);

	return;
}


static int CPXPUBLIC adv_branch_cut_callback_driver(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {

	instance* inst = (instance*)userhandle;
	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE) return branch_cut_callback(context, contextid, inst);
	int mynode = -1;
	int mydepth = -1;
	double rval;
	unsigned int value = 0;
	switch (inst->model_type) {
		
		case ADVBC_ROOT:
			CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
			if (mynode != 0) return 0;
			break;
		case ADVBC_DEPTH_5:
			CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODEDEPTH, &mydepth);
			if (mydepth > 5) return 0;
			break;
		case ADVBC_PROB_50:
			safe_rand(&rval, value);
			if (rval > 0.5) return 0;
			break;
		case ADVBC_PROB_10:
			safe_rand(&rval, value);
			if (rval > 0.1) return 0;
			break;
	}
	if (contextid == CPX_CALLBACKCONTEXT_RELAXATION) return adv_branch_cut_callback(context, contextid, inst);
	print_error("Unknown \"contextid\" in adv_branch_cut_callback_driver");
	return 1;
}


static int CPXPUBLIC adv_branch_cut_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {

	instance* inst = (instance*)userhandle;
	int ncols = inst->ncols;
	double* xstar = (double*)malloc(inst->ncols * sizeof(double));
	double objval = CPX_INFBOUND;
	if (CPXcallbackgetrelaxationpoint(context, xstar, 0, ncols - 1, &objval)) print_error("CPXcallbackgetrelaxationpoint error");

	int ncount = inst->nnodes;
	int ecount = inst->nnodes * (inst->nnodes - 1) / 2;
	int* elist = malloc(2 * ecount * sizeof(int));
	int k = 0;
	for (int i = 0; i < ncount; i++) {
		for (int j = i + 1; j < ncount; j++) {
			elist[k++] = i;
			elist[k++] = j;
		}
	}
	
	int ncomp = 0;
	int* compscount;
	int* comps;
	
	if (CCcut_connect_components(ncount, ecount, elist, xstar, &ncomp, &compscount, &comps)) print_error("Error during concorde connect comps algorithm!");
	if (inst->verbose >= HIGH) printf("#connected components: %d\n", ncomp);

	int* index = (int*)calloc(ncols, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(ncols, sizeof(double));			// Array of row variables coefficients


	if (ncomp > 1) {

		for (int n = 0; n < ncomp; n++) {			// n is the index of the current connected component

			double rhs = compscount[n] - 1.0;
			char sense = 'L';
			int i_zero = 0;
			int purgeable = CPX_USECUT_FILTER;
			int local = 0;
			int nnz = 0;
			int pos = 0;
			for (int i = pos; i < pos + compscount[n]; i++) {
				for (int j = i + 1; j < pos + compscount[n]; j++) {
					if (i == j) {
						continue;
					}
					index[nnz] = xpos(comps[i], comps[j], inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
			if (CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &i_zero, index, value, &purgeable, &local)) print_error("CPXcallbackaddusercuts() error");
		}
		
	}

	else if (ncomp == 1) {
		
		concorde_instance cc;
		cc.inst = inst;
		cc.context = context;
		cc.local = 0;
		cc.index = malloc(inst->ncols * sizeof(int));
		cc.value = calloc(inst->ncols, sizeof(double));
		if (CCcut_violated_cuts(ncount, ecount, elist, xstar, 2 - 0.1, doit_fn_concorde, &cc)) print_error("Error during concorde violated cuts algorithm!");
		free(cc.index);
		free(cc.value);
	}

	free(elist);
	free(index);
	free(value);
	free(xstar);
	free(compscount);
	free(comps);
	
	return 0;
}


int doit_fn_concorde(double cutval, int cutcount, int* cut, void* in_param) {
	concorde_instance* cc = (concorde_instance*)in_param;
	instance* inst = cc->inst;
	int* index = cc->index;
	double* value = cc->value;
	if (cutcount > 2 && cutcount < inst->nnodes) {
		int i_zero = 0;
		char sense = 'G';
		double rhs = 2.0;
		int purgeable = CPX_USECUT_FILTER;
		int nnz = 0;

		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = 0; j < cutcount; j++) {
				if (i < cut[j]) {
					index[nnz] = xpos(i, cut[j], inst);
					value[nnz] = 1.0;
					nnz++;
				}
			}
		}
		if (CPXcallbackaddusercuts(cc->context, 1, nnz, &rhs, &sense, &i_zero, index, value, &purgeable, &cc->local)) print_error("CPXcallbackaddusercuts() error");
	}

	if (inst->verbose >= HIGH) printf("cutval %f\n", cutval);

	return 0;
}


void solve_heur_hard_fix(instance* inst, CPXENVptr env, CPXLPptr lp) {

	double residual_timelimit = inst->timelimit;

	// Build an array of indices of variables/columns x(i,j) of the model
	int* indices = (int*)calloc(inst->ncols, sizeof(int));
	int k = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			indices[k++] = xpos(i, j, inst);
		}
	}

	if (CPXsetintparam(env, CPXPARAM_Advance, 1)) print_error("CPXsetintparam() error in setting CPXPARAM_Advance");

	// Set to solve just the root node: by doing that we will get a feasible solution, but not the optimal one => good as a starting point for heuristics
	if (CPXsetintparam(env, CPX_PARAM_NODELIM, 0)) print_error("CPXsetintparam() error in setting seed");
	inst->model_type = ADVBC_ROOT;

	int n_iter = 0;
	char lu = 'L';
	int beg = 0;
	double temp_obj_val = INFINITY;
	double low_bound_value = 1.0;
	double* curr_best_sol = (double*)calloc(inst->ncols, sizeof(double));
	if (CPXsetdblparam(env, CPX_PARAM_TILIM, CPX_INFBOUND)) { print_error("CPXsetdblparam() error in setting timelimit"); }
	while (1) {

		double t1 = second();
		solve_adv_branch_cut(inst, env, lp);
		double t2 = second();
		if (CPXgetx(env, lp, curr_best_sol, 0, inst->ncols - 1)) print_error("CPXgetx() error");

		if (inst->verbose >= MEDIUM) printf("Time used for iteration number %d: %f\n", n_iter, t2 - t1);

		// Update the amount of time left before timelimit is reached and provide it to Cplex to check
		residual_timelimit = residual_timelimit - (t2 - t1);

		// If the timelimit is reached for the current iteration => exit the loop
		if (residual_timelimit <= 0) break;

		CPXgetobjval(env, lp, &temp_obj_val);
		if (inst->z_best > temp_obj_val) {
			inst->z_best = temp_obj_val;
			inst->best_sol = curr_best_sol;

			// Set the feasible solution (not optimal) from which to start with the heuristics
			//if (CPXcopy(env, lp, inst->ncols, indices, curr_best_sol)) print_error("CPXcopymipstart() error in setting known solution");
			if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, indices, inst->best_sol, CPX_MIPSTART_AUTO, NULL)) print_error("CPXaddmipstarts() error in setting known solution");
		}
		if (inst->verbose >= HIGH) printf("inst->z_best: %f\n", inst->z_best);

		int n_fixed_edges = 0;
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = i + 1; j < inst->nnodes; j++) {
				if (curr_best_sol[xpos(i, j, inst)] > 0.5 && ((double)rand() / RAND_MAX) <= 0.5) {
					int var_index = xpos(i, j, inst);
					if (CPXchgbds(env, lp, 1, &var_index, &lu, &low_bound_value)) print_error("CPXchgbds() error in setting edge lower bound");
					n_fixed_edges++;
				}
			}
		}
		if (inst->verbose >= HIGH) printf("n_fixed_edges: %d\n", n_fixed_edges);

		// Increment current iteration number
		n_iter++;
	}

	free(indices);
	free(curr_best_sol);

	return;
}


void update_connected_components(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp) {

	// Start with 0 connected components
	*ncomp = 0;
	// Initialize all cells of "successors" and "components" arrays to -1 (meaning that the node has not been visited yet)
	for (int i = 0; i < inst->nnodes; i++) {
		succ[i] = -1;
		comp[i] = -1;
	}

	// Start visiting from node 0, then from the first node not yet visited, and so on...
	for (int start = 0; start < inst->nnodes; start++) {
		if (comp[start] >= 0) continue;								// Node "start" was already visited, just skip it

		// If we arrive here, it means that a new component is found (count it)
		(*ncomp)++;
		int i = start;
		int done = 0;
		while (!done) {												// Visit the current component
			comp[i] = *ncomp;
			done = 1;
			// Scan for feasible nodes to draw an edge and close the component
			for (int j = 0; j < inst->nnodes; j++) {
				// The edge [i,j] is selected in xstar and j was not visited before
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && comp[j] == -1) {
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}
		succ[i] = start;											// Last edge to close the cycle

		// Now go to the next component, if there is one...
	}
}


void build_model_BASIC(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// N.B. This model lacks subtour elimination constraints!

	double zero = 0.0;
	char binary = 'B';

	char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary variables x(i,j) for i < j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);		// Set a name for the new column/varible
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) { print_error("Wrong CPXnewcols on x variables"); }
			// Verify if xpos returns the right position inside the tableu (xpos starts from 0)
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst)) { print_error("Wrong position for x variables"); }
		}
	}

	// Add the degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
	int i_zero = 0;
	int* index = (int*)calloc(inst->nnodes, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->nnodes, sizeof(double));			// Array of row variables coefficients
	double rhs = 2.0;														// "rhs" = right-hand side of the degree constraints
	char sense = 'E';														// 'L' for less-or-equal constraint
	int nnz = inst->nnodes;													// number of cells != 0 in the row
	for (int h = 0; h < inst->nnodes; h++) {

		sprintf(cname[0], "degree(%d)", h + 1);								// Set a name for the new row/constraint

		int j = 0;
		for (int i = 0; i < inst->nnodes; i++) {

			if (i == h ) {
				continue;
			}

			index[j] = xpos(i, h, inst);
			value[j] = 1.0;
			j++;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for degree constraints!");
	}

	// Outputs to file "basic_model.lp" the built model
	if (inst->verbose >= LOW) create_lp_file(inst, env, lp, "model_BASIC");

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


void build_model_MTZ_STATIC(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// Model MTZ with static subtour elimination constraints.

	double zero = 0.0;
	char binary = 'B';
	char integer = 'I';

	char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary variables x(i,j) for any i,j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {

			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);		// Set a name for the new column/varible
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (i == j) { ub = 0.0; }
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) { print_error("Wrong CPXnewcols on x variables"); }

			//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
			// Verify if xpos returns the right position inside the tableu (xpos_compact starts from 0)
			if (CPXgetnumcols(env, lp) - 1 != xpos_compact(i, j, inst)) {
				printf("Curr interation: %d %d", i, j);
				print_error("Wrong position for x variables");
			}
		}
	}

	// Add u variables for each node that is not the first one
	for (int i = 1; i < inst->nnodes; i++) {

		sprintf(cname[0], "u(%d)", i + 1);		// Set a name for the new column/varible
		double obj = 0.0;
		double lb = 0.0;
		double ub = inst->nnodes - 2.0;
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) { print_error("Wrong CPXnewcols on x variables"); }

		//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
		if (CPXgetnumcols(env, lp) - 1 != upos_compact(i, inst)) {
			printf("Curr interation: %d", i);
			print_error("Wrong position for x variables");
		}
	}

	// Add the inner and outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
	int i_zero = 0;
	int* index = (int*)calloc(inst->nnodes, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->nnodes, sizeof(double));			// Array of row variables coefficients
	double rhs = 1.0;														// "rhs" = right-hand side of the degree constraints
	char sense = 'E';														// 'L' for less-or-equal constraint
	int nnz = inst->nnodes;													// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "inner_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(j, i, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "outer_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

				index[j] = xpos_compact(i, j, inst);
				value[j] = 1.0;
				if (j == i)
					value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}

	// Add the u consistency constraints for each edge (i,j)
	double M = inst->nnodes - 1.0;													// Smallest M value for big M trick
	rhs = M - 1;																	// "rhs" = right-hand side of the degree constraints
	sense = 'L';																	// 'L' for less-or-equal constraint
	nnz = 3;																		// number of cells != 0 in the row
	for (int i = 1; i < inst->nnodes; i++) {
		for (int j = 1; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "u_consistency_for_arc(%d,%d)", i + 1, j + 1);
			index[0] = upos_compact(i, inst);										// +1.0 * Ui
			value[0] = 1.0;
			index[1] = upos_compact(j, inst);										// -1.0 * Uj
			value[1] = -1.0;
			index[2] = xpos_compact(i, j, inst);									// +M * Xij
			value[2] = M;
			if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for u-consistency");
		}
	}

	// Outputs to file "model_MTZ_static.lp" the built model
	if (inst->verbose >= LOW) create_lp_file(inst, env, lp, "model_MTZ_STATIC");

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


void build_model_MTZ_LAZY(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// Model MTZ with lazy subtour elimination constraints.

	double zero = 0.0;
	char binary = 'B';
	char integer = 'I';

	char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary variables x(i,j) for any i,j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {

			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);		// Set a name for the new column/varible
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (i == j) { ub = 0.0; }
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) { print_error("Wrong CPXnewcols on x variables"); }

			//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
			// Verify if xpos returns the right position inside the tableu (xpos_compact starts from 0)
			if (CPXgetnumcols(env, lp) - 1 != xpos_compact(i, j, inst)) {
				printf("Curr interation: %d %d", i, j);
				print_error("Wrong position for x variables");
			}
		}
	}

	// Add u variables for each node that is not the first one
	for (int i = 1; i < inst->nnodes; i++) {

		sprintf(cname[0], "u(%d)", i + 1);		// Set a name for the new column/varible
		double obj = 0.0;
		double lb = 0.0;
		double ub = inst->nnodes - 2.0;
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) { print_error("Wrong CPXnewcols on x variables"); }

		//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
		if (CPXgetnumcols(env, lp) - 1 != upos_compact(i, inst)) {
			printf("Curr interation: %d", i);
			print_error("Wrong position for x variables");
		}
	}

	// Add the inner and outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
	int i_zero = 0;
	int* index = (int*)calloc(inst->nnodes, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->nnodes, sizeof(double));			// Array of row variables coefficients
	double rhs = 1.0;														// "rhs" = right-hand side of the degree constraints
	char sense = 'E';														// 'L' for less-or-equal constraint
	int nnz = inst->nnodes;													// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "inner_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(j, i, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "outer_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(i, j, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}

	// Add the u consistency constraints for each edge (i,j)
	double M = inst->nnodes - 1.0;								// Smallest M value for big M trick
	rhs = M - 1;											// "rhs" = right-hand side of the degree constraints
	sense = 'L';											// 'L' for less-or-equal constraint
	nnz = 3;												// number of cells != 0 in the row
	for (int i = 1; i < inst->nnodes; i++) {
		for (int j = 1; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "u_consistency_for_arc(%d,%d)", i + 1, j + 1);
			index[0] = upos_compact(i, inst);										// +1.0 * Ui
			value[0] = 1.0;
			index[1] = upos_compact(j, inst);										// -1.0 * Uj
			value[1] = -1.0;
			index[2] = xpos_compact(i, j, inst);									// +M * Xij
			value[2] = M;
			if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &i_zero, index, value, cname)) print_error("wrong CPXlazyconstraints() for u-consistency");
		}
	}

	// Outputs to file "model_MTZ_lazy_u_consistency.lp" the built model
	if (inst->verbose >= LOW) create_lp_file(inst, env, lp, "model_MTZ_LAZY");

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


void build_model_MTZ_SEC2_STATIC(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// Model MTZ with static (added directly to the TSP model) subtour elimination constraints of dimension 2.

	double zero = 0.0;
	char binary = 'B';
	char integer = 'I';

	char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary variables x(i,j) for any i,j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {

			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);		// Set a name for the new column/varible
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (i == j) { ub = 0.0; }
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) { print_error("Wrong CPXnewcols on x variables"); }

			//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
			// Verify if xpos returns the right position inside the tableu (xpos_compact starts from 0)
			if (CPXgetnumcols(env, lp) - 1 != xpos_compact(i, j, inst)) {
				printf("Curr interation: %d %d", i, j);
				print_error("Wrong position for x variables");
			}
		}
	}

	// Add u variables for each node that is not the first one
	for (int i = 1; i < inst->nnodes; i++) {

		sprintf(cname[0], "u(%d)", i + 1);		// Set a name for the new column/varible
		double obj = 0.0;
		double lb = 0.0;
		double ub = inst->nnodes - 2.0;
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) { print_error("Wrong CPXnewcols on x variables"); }

		//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
		if (CPXgetnumcols(env, lp) - 1 != upos_compact(i, inst)) {
			printf("Curr interation: %d", i);
			print_error("Wrong position for x variables");
		}
	}

	// Add the inner and outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
	int i_zero = 0;
	int* index = (int*)calloc(inst->nnodes, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->nnodes, sizeof(double));			// Array of row variables coefficients
	double rhs = 1.0;														// "rhs" = right-hand side of the degree constraints
	char sense = 'E';														// 'L' for less-or-equal constraint
	int nnz = inst->nnodes;													// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "inner_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(j, i, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "outer_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(i, j, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}

	// Add the u consistency constraints for each edge (i,j)
	double M = inst->nnodes - 1.0;													// Smallest M value for big M trick
	rhs = M - 1;																	// "rhs" = right-hand side of the degree constraints
	sense = 'L';																	// 'L' for less-or-equal constraint
	nnz = 3;																		// number of cells != 0 in the row
	for (int i = 1; i < inst->nnodes; i++) {
		for (int j = 1; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "u_consistency_for_arc(%d,%d)", i + 1, j + 1);
			index[0] = upos_compact(i, inst);										// +1.0 * Ui
			value[0] = 1.0;
			index[1] = upos_compact(j, inst);										// -1.0 * Uj
			value[1] = -1.0;
			index[2] = xpos_compact(i, j, inst);									// +M * Xij
			value[2] = M;
			if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for u-consistency!");
		}
	}

	// Add the subtour elimination constraints of dimension 2 (2-node SECs): Xij + Xji <= 1, for any i < j
	rhs = 1.0;													// "rhs" = right-hand side of the degree constraints
	sense = 'L';												// 'L' for less-or-equal constraint
	nnz = 2;													// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "_2_node_SEC(%d,%d)", i + 1, j + 1);
			index[0] = xpos_compact(i, j, inst);										// +1.0 * Xij
			value[0] = 1.0;
			index[1] = xpos_compact(j, i, inst);										// +1.0 * Xji
			value[1] = 1.0;
			if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for 2-node SECs!");
		}
	}

	// Outputs to file "model_MTZ_SEC2_STATIC.lp" the built model
	if (inst->verbose >= LOW) create_lp_file(inst, env, lp, "model_MTZ_SEC2_STATIC");

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


void build_model_MTZ_SEC2_LAZY(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// Model MTZ with lazy subtour elimination constraints of dimension 2.

	double zero = 0.0;
	char binary = 'B';
	char integer = 'I';

	char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary variables x(i,j) for any i,j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {

			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);		// Set a name for the new column/varible
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (i == j) { ub = 0.0; }
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) { print_error("Wrong CPXnewcols on x variables"); }

			//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
			// Verify if xpos returns the right position inside the tableu (xpos_compact starts from 0)
			if (CPXgetnumcols(env, lp) - 1 != xpos_compact(i, j, inst)) {
				printf("Curr interation: %d %d", i, j);
				print_error("Wrong position for x variables");
			}
		}
	}

	// Add u variables for each node that is not the first one
	for (int i = 1; i < inst->nnodes; i++) {

		sprintf(cname[0], "u(%d)", i + 1);		// Set a name for the new column/varible
		double obj = 0.0;
		double lb = 0.0;
		double ub = inst->nnodes - 2.0;
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) { print_error("Wrong CPXnewcols on x variables"); }

		//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
		if (CPXgetnumcols(env, lp) - 1 != upos_compact(i, inst)) {
			printf("Curr interation: %d", i);
			print_error("Wrong position for x variables");
		}
	}

	// Add the inner and outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
	int i_zero = 0;
	int* index = (int*)calloc(inst->nnodes, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->nnodes, sizeof(double));			// Array of row variables coefficients
	double rhs = 1.0;														// "rhs" = right-hand side of the degree constraints
	char sense = 'E';														// 'L' for less-or-equal constraint
	int nnz = inst->nnodes;													// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "inner_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(j, i, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "outer_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(i, j, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}

	// Add the u consistency constraints for each edge (i,j)
	double M = inst->nnodes - 1.0;													// Smallest M value for big M trick
	rhs = M - 1;																	// "rhs" = right-hand side of the degree constraints
	sense = 'L';																	// 'L' for less-or-equal constraint
	nnz = 3;																		// number of cells != 0 in the row
	for (int i = 1; i < inst->nnodes; i++) {
		for (int j = 1; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "u_consistency_for_arc(%d,%d)", i + 1, j + 1);
			index[0] = upos_compact(i, inst);										// +1.0 * Ui
			value[0] = 1.0;
			index[1] = upos_compact(j, inst);										// -1.0 * Uj
			value[1] = -1.0;
			index[2] = xpos_compact(i, j, inst);									// +M * Xij
			value[2] = M;
			if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &i_zero, index, value, cname)) print_error("wrong CPXlazyconstraints() for u-consistency");
		}
	}

	// Add the subtour elimination constraints of dimension 2 (2-node SECs): Xij + Xji <= 1, for any i < j
	rhs = 1.0;													// "rhs" = right-hand side of the degree constraints
	sense = 'L';												// 'L' for less-or-equal constraint
	nnz = 2;													// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "_2_node_SEC(%d,%d)", i + 1, j + 1);
			index[0] = xpos_compact(i, j, inst);										// +1.0 * Xij
			value[0] = 1.0;
			index[1] = xpos_compact(j, i, inst);										// +1.0 * Xji
			value[1] = 1.0;
			if (CPXaddlazyconstraints(env, lp, 1, nnz, &rhs, &sense, &i_zero, index, value, cname)) print_error("wrong CPXlazyconstraints() for 2-node SECs");
		}
	}

	// Outputs to file "model_MTZ_SEC2_LAZY.lp" the built model
	if (inst->verbose >= LOW) create_lp_file(inst, env, lp, "model_MTZ_SEC2_LAZY");

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


void build_model_GG(instance* inst, CPXENVptr env, CPXLPptr lp) {

	// Model GG with flow constraints

	double zero = 0.0;
	char binary = 'B';
	char integer = 'I';

	char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	// Add binary variables x(i,j) for any i,j
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {

			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);		// Set a name for the new column/varible
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (i == j) { ub = 0.0; }
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) { print_error("Wrong CPXnewcols on x variables"); }

			//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
			// Verify if xpos returns the right position inside the tableu (xpos_compact starts from 0)
			if (CPXgetnumcols(env, lp) - 1 != xpos_compact(i, j, inst)) {
				printf("Curr interation: %d %d", i, j);
				print_error("Wrong position for x variables");
			}
		}
	}

	// Add y variables for each edge
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {

			sprintf(cname[0], "y(%d,%d)", i + 1, j + 1);		// Set a name for the new column/varible
			double obj = 0.0;
			double lb = 0.0;
			double ub = inst->nnodes - 1.0;
			if (i == j || j == 0) ub = 0.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) { print_error("Wrong CPXnewcols on y variables"); }

			//printf("Curr num_cols: %d\n", CPXgetnumcols(env, lp));
			if (CPXgetnumcols(env, lp) - 1 != ypos_compact(i, j, inst)) {
				printf("Curr interation: %d", i);
				print_error("Wrong position for x variables");
			}
		}
	}

	// Add the inner and outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
	int i_zero = 0;
	int* index = (int*)calloc(inst->nnodes, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->nnodes, sizeof(double));			// Array of row variables coefficients
	double rhs = 1.0;														// "rhs" = right-hand side of the degree constraints
	char sense = 'E';														// 'L' for less-or-equal constraint
	int nnz = inst->nnodes;													// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "inner_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(j, i, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}
	for (int i = 0; i < inst->nnodes; i++) {

		sprintf(cname[0], "outer_degree(%d)", i + 1);								// Set a name for the new row/constraint

		for (int j = 0; j < inst->nnodes; j++) {

			index[j] = xpos_compact(i, j, inst);
			value[j] = 1.0;
			if (j == i)
				value[j] = 0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for inner degree constraints!");
	}

	// Add the linking constraints for the variable y
	// Add the u consistency constraints for each edge (i,j)
	rhs = 0.0;																		// "rhs" = right-hand side of the degree constraints
	sense = 'L';																	// 'L' for less-or-equal constraint
	nnz = 2;																		// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "linking_constr(%d,%d)", i + 1, j + 1);
			index[0] = xpos_compact(i, j, inst);									// +(2.0 - N) * Xij
			value[0] = 1.0 - inst->nnodes;
			index[1] = ypos_compact(i, j, inst);									// +1.0 * Yij
			value[1] = 1.0;
			if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for linking_constr");
		}
	}

	// Add the starting flow constraint (the starting flow value from node 1 is n - 1)
	rhs = inst->nnodes - 1.0;														// "rhs" = right-hand side of the degree constraints
	sense = 'E';																	// 'E' for equality constraint 
	nnz = inst->nnodes;
	index[0] = ypos_compact(0, 0, inst);
	value[0] = 0.0;
	sprintf(cname[0], "starting_flow_constr");
	for (int j = 1; j < inst->nnodes; j++) {
		index[j] = ypos_compact(0, j, inst);
		value[j] = 1.0;
	}
	if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for linking_constr");

	// Add the flow balance constraints
	rhs = 1.0;
	sense = 'E';
	nnz = inst->nnodes * 2;
	index = realloc(index, nnz * sizeof(int));
	value = realloc(value, nnz * sizeof(double));
	for (int h = 1; h < inst->nnodes; h++) {		

		sprintf(cname[0], "flow_balance(%d)", h + 1);					// Set a name for the new row/constraint

		for (int i = 0; i < inst->nnodes; i++) {
			if (i == h) {
				index[i] = ypos_compact(i, h, inst);
				value[i] = 0.0;
				index[i + inst->nnodes] = ypos_compact(h, i, inst);
				value[i + inst->nnodes] = 0.0;
				continue;
			}
			index[i] = ypos_compact(i, h, inst);
			value[i] = 1.0;
			index[i + inst->nnodes] = ypos_compact(h, i, inst);
			value[i + inst->nnodes] = -1.0;
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for linking_constr");
	}

	// Outputs to file "model_GG.lp" the built model
	if (inst->verbose >= LOW) create_lp_file(inst, env, lp, "model_GG");

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
}


void create_lp_file(instance* inst, CPXENVptr env, CPXLPptr lp, const char* model_name) {

	char model_file_path[100];
	sprintf(model_file_path, "../outputs/%s/%s.lp", inst->inst_name, model_name);
	if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
	CPXwriteprob(env, lp, model_file_path, NULL);
}

