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
#include "convex_hull.h"


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

	if (inst->verbose >= LOW) {
		sprintf(logfile_path, "%s/logfile_%s.txt", logfile_path, models[inst->model_type]);
		if (CPXsetlogfilename(env, logfile_path, "w")) print_error("CPXsetlogfilename() error in setting logfile name");

		sprintf(edges_file_path, "%s/model_%s_edges.dat", edges_file_path, models[inst->model_type]);
	}

	// Discern if a model is symmetric (solves the TSP for a directed or an undirected graph)
	int symmetric = -1;

	switch (inst->model_type) {

		case BASIC:
			symmetric = 0;
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);

			if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			break;

		case MTZ_STATIC:
			symmetric = 1;
			build_model_MTZ_STATIC(inst, env, lp);
		
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) print_error("CPXsetdblparam() error in setting integer value tolerance");
			if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			break;

		case MTZ_LAZY:
			symmetric = 1;
			build_model_MTZ_LAZY(inst, env, lp);
		
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) print_error("CPXsetdblparam() error in setting integer value tolerance");
			if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			break;

		case MTZ_SEC2_STATIC:
			symmetric = 1;
			build_model_MTZ_SEC2_STATIC(inst, env, lp);
			
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) print_error("CPXsetdblparam() error in setting integer value tolerance");
			if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			break;

		case MTZ_SEC2_LAZY:
			symmetric = 1;
			build_model_MTZ_SEC2_LAZY(inst, env, lp);
		
			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) print_error("CPXsetdblparam() error in setting integer value tolerance");
			if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			break;

		case GG:
			symmetric = 1;
			build_model_GG(inst, env, lp);
		
			if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
			mip_solved_to_optimality(inst, env, lp);												// Check if CPXmipopt has ended correctly
			break;

		case BENDERS:
			symmetric = 0;
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			
			solve_benders(inst, env, lp);

			break;

		case BRANCH_CUT:
			symmetric = 0;
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			inst->use_2_opt = 0;						// Solve without 2-opt refinement heuristics

			solve_branch_cut(inst, env, lp);

			break;

		case BRANCH_CUT_2_OPT:
			symmetric = 0;
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			inst->use_2_opt = 1;						// Solve using 2-opt refinement heuristics

			solve_branch_cut(inst, env, lp);

			break;
		
		case ADVBC_STD:
		case ADVBC_ROOT:
		case ADVBC_DEPTH_5:
		case ADVBC_PROB_50:
		case ADVBC_PROB_10:
			symmetric = 0;
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
		
			solve_adv_branch_cut(inst, env, lp);

			break;

		case HEUR_HARD_FIX_50:
		case HEUR_HARD_FIX_70:
		case HEUR_HARD_FIX_90:
		case HEUR_HARD_FIX_VAR:
			symmetric = 0;
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);
			
			solve_heur_hard_fix(inst, env, lp);

			break;

		case HEUR_SOFT_FIX_3:
		case HEUR_SOFT_FIX_5:
		case HEUR_SOFT_FIX_7:
		case HEUR_SOFT_FIX_9:
			symmetric = 0;
			build_model_BASIC(inst, env, lp);
			inst->ncols = CPXgetnumcols(env, lp);

			solve_heur_soft_fix(inst, env, lp);

			break;

		case HEUR_GREEDY:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_greedy = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_greedy(inst, x_greedy);

			if (inst->verbose >= LOW) print_solution(inst, x_greedy, symmetric, edges_file_path);
			free(x_greedy);
			break;

		case HEUR_GRASP_GREEDY:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_grasp_greedy = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_grasp_greedy(inst, x_grasp_greedy, inst->timelimit);

			if (inst->verbose >= LOW) print_solution(inst, x_grasp_greedy, symmetric, edges_file_path);
			free(x_grasp_greedy);
			break;

		case HEUR_EXTRA_MILEAGE:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_extra_mileage = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_extra_mileage(inst, x_extra_mileage);

			if (inst->verbose >= LOW) print_solution(inst, x_extra_mileage, symmetric, edges_file_path);
			free(x_extra_mileage);
			break;

		case HEUR_GRASP_EXTRA_MILEAGE:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_grasp_extra_mileage = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_grasp_extra_mileage(inst, x_grasp_extra_mileage, inst->timelimit);

			if (inst->verbose >= LOW) print_solution(inst, x_grasp_extra_mileage, symmetric, edges_file_path);
			free(x_grasp_extra_mileage);
			break;

		case HEUR_2_OPT:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_2_opt = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_grasp_greedy(inst, x_2_opt, inst->timelimit / 10);
			solve_heur_2_opt(inst, x_2_opt, NULL, (inst->timelimit / 10) * 9);

			if (inst->verbose >= LOW) print_solution(inst, x_2_opt, symmetric, edges_file_path);
			free(x_2_opt);
			break;

		case HEUR_MULTI_START:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_multi_start = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_multi_start(inst, x_multi_start);

			if (inst->verbose >= LOW) print_solution(inst, x_multi_start, symmetric, edges_file_path);
			free(x_multi_start);
			break;

		case HEUR_VNS:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_vns = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_vns(inst, x_vns);

			if (inst->verbose >= LOW) print_solution(inst, x_vns, symmetric, edges_file_path);
			free(x_vns);
			break;

		case HEUR_TABU:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_tabu = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_tabu(inst, x_tabu);

			if (inst->verbose >= LOW) print_solution(inst, x_tabu, symmetric, edges_file_path);
			free(x_tabu);
			break;

		case HEUR_GENETIC:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_genetic = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_genetic(inst, x_genetic, 1000, 0.0);

			if (inst->verbose >= LOW) print_solution(inst, x_genetic, symmetric, edges_file_path);
			free(x_genetic);
			break;

		case HEUR_GENETIC_2_OPT:
			symmetric = 0;
			inst->ncols = (inst->nnodes * (inst->nnodes - 1)) / 2;
			double* x_genetic_2_opt = (double*)calloc(inst->ncols, sizeof(double));

			solve_heur_genetic(inst, x_genetic_2_opt, 100, 0.1);

			if (inst->verbose >= LOW) print_solution(inst, x_genetic_2_opt, symmetric, edges_file_path);
			free(x_genetic_2_opt);
			break;

		default:
			print_error("Choose a correct value for the model to be used!");
	}

	if (inst->verbose >= MEDIUM) printf("Complete path for *.dat file: %s\n\n", edges_file_path);

	if (inst->model_type < HEUR_GREEDY) {													// Let's manage how models/approaches using Cplex print their solutions

		// Allocate memory for the optimal solution array
		int ncols = CPXgetnumcols(env, lp);
		double* xstar = (double*)calloc(ncols, sizeof(double));

		// Copy the optimal solution from the Cplex environment to the new array "xstar" (does not work for HEUR_SOFT_FIX_K because it edits the model at the end)
		if ((inst->model_type < HEUR_SOFT_FIX_3 || inst->model_type > HEUR_SOFT_FIX_9) && !inst->timelimit_exceeded && CPXgetx(env, lp, xstar, 0, ncols - 1)) print_error("CPXgetx() error");

		// Fill the .dat file with the correctly formatted nodes of the found solution
		if (inst->verbose >= LOW && (inst->model_type < HEUR_SOFT_FIX_3 || inst->model_type > HEUR_SOFT_FIX_9)) print_solution(inst, xstar, symmetric, edges_file_path);
		else if (inst->verbose >= LOW) print_solution(inst, inst->best_sol, symmetric, edges_file_path);

		// Free allocated memory and close Cplex model
		free(xstar);
	}

	if (CPXsetlogfilename(env, NULL, NULL)) print_error("CPXsetlogfilename() error");
	if (CPXfreeprob(env, &lp)) print_error("CPXfreeprob() error");
	if (CPXcloseCPLEX(&env)) print_error("CPXcloseCPLEX() error");

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
	char* command;
	if (gn_com == NULL) print_error("Error while opening gnuplot_commands.txt file\n");
	if (!symmetric) 
		command = "set style line 1 \\\n\tlinecolor rgb '#FF0000' \\\n\tlinetype 1 linewidth 2 \\\n\tpointtype 7 pointsize 2 \\\n\nplot \"%s\" using 1:2 with linespoints linestyle 1\npause mouse close";
	else 
		command = "set style line 1 \\\n\tlinecolor rgb '#FF0000' \\\n\tlinetype 1 linewidth 2 \\\n\tpointtype 7 pointsize 2 \\\n\nplot \"%s\" using 1:2:3:4 with vectors linestyle 1\npause mouse close";
	
	fprintf(gn_com, command, edges_file_path);
	
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

	return i * inst->nnodes + j;
}


int upos_compact(int i, instance* inst) {

	return (inst->nnodes) * (inst->nnodes) + i - 1;
}


int ypos_compact(int i, int j, instance* inst) {

	return (inst->nnodes) * (inst->nnodes) + i * (inst->nnodes) + j;
}


double extra_mileage(int i, int j, int new_node, instance* inst) {

	return dist(i, new_node, inst) + dist(new_node, j, inst) - dist(i, j, inst);
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
			if (inst->model_type >= HEUR_HARD_FIX_50) return 0;
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
			// Add a new subtour elimination constraint for each connected component of the current solution
			update_benders_constraints(env, lp, inst, comp, n_comp, n_iter);
		}

		double t1 = second();
		// Optimize with the new constraint
		if (CPXmipopt(env, lp)) print_error("CPXmipopt() error");
		double t2 = second();
		
		// Check if CPXmipopt has ended correctly
		mip_solved_to_optimality(inst, env, lp);

		if (inst->verbose >= MEDIUM) printf("Time used for iteration number %d: %f\n", n_iter, t2 - t1);

		// Update the amount of time left before timelimit is reached and provide it to Cplex to check
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (inst->verbose >= MEDIUM) printf("New time limit: %f\n\n", residual_timelimit);
		if (CPXsetdblparam(env, CPX_PARAM_TILIM, residual_timelimit)) print_error("CPXsetdblparam() error in setting timelimit");

		// Extract the new solution (but only if timelimit has not been reached)
		if (!inst->timelimit_exceeded && CPXgetx(env, lp, x, 0, ncols - 1)) print_error("CPXgetx() error");

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
	int use_2_opt = inst->use_2_opt;
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
	else if (n_comp == 1 && use_2_opt) {									// Means that the solution is feasible (without cycles => candidate solution which can be improved with 2-opt)

		compute_succ(inst, xstar, succ);

		// Compute from scratch the cost of the solution
		double curr_sol_cost = 0.0;
		for (int i = 0; i < inst->nnodes; i++) curr_sol_cost += dist(i, succ[i], inst);
		inst->z_best = curr_sol_cost;

		// Call the 2-opt refinement heuristics
		solve_heur_2_opt(inst, xstar, succ, INFINITY);

		// Build the array of indexes associated to solution variables
		int* index = (int*)calloc(ncols, sizeof(int));						
		for (int k = 0; k < inst->ncols; k++) index[k] = k;

		// Ask Cplex to use our handcrafted improved solution (which is always complete, so no solution repair strategy is provided)
		CPXcallbackpostheursoln(context, inst->ncols, index, xstar, inst->z_best, CPXCALLBACKSOLUTION_NOCHECK);

		free(index);
	}

	// De-comment below to get a sequence of plots showing how the solution evolves
	/*
	if (mythread == 0 && incumbent < 10000) {
		printf("n_comp: %d\n\n", n_comp);
		print_solution(inst, xstar, 0, "../outputs/B&C_plot.dat");
		system("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe ../outputs/gnuplot_commands.txt");
		system("PAUSE");
	}
	*/

	free(xstar);
	free(succ);
	free(comp);

	return 0;
}


void solve_adv_branch_cut(instance* inst, CPXENVptr env, CPXLPptr lp) {

	if (inst->model_type >= HEUR_HARD_FIX_50) inst->tsp_solver = ADVBC_ROOT;
	else inst->tsp_solver = inst->model_type;

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
	switch (inst->tsp_solver) {
		
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
	
	// Returns the connected components of the graph given by the edgeset.
	// More details (look for "CCcut_connect_components") can be found here: http://www.math.uwaterloo.ca/tsp/concorde/DOC/cut.html
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
		// Computes the global minimum cut, but calls the doit_fn_concorde function for any cut the algorithm encounters that has capacity at most cutoff.
		// More details (look for "CCcut_violated_cuts") can be found here: http://www.math.uwaterloo.ca/tsp/concorde/DOC/cut.html
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



/*CODICE PER GRAFICO
FILE* hard_best_solution_csv_file;
char hard_best_solution_csv_path[120];
sprintf(hard_best_solution_csv_path, "../outputs/%s_best_solution_cost.csv", models[inst->model_type]);
hard_best_solution_csv_file = fopen(hard_best_solution_csv_path, "a");
fprintf(hard_best_solution_csv_file, "%s", inst->inst_name);
fprintf(hard_best_solution_csv_file, "%f %f %f\n", inst->timelimit - residual_timelimit, curr_sol_cost, min_sol_cost);
fclose(hard_best_solution_csv_file);
*/



void solve_heur_hard_fix(instance* inst, CPXENVptr env, CPXLPptr lp) {

	double param[] = { 0.9, 0.7, 0.5 };
	int param_id = 0;
	double prob_param = param[param_id];
	double temp = prob_param;
	switch (inst->model_type) {

		case HEUR_HARD_FIX_50:
			prob_param = param[2];
			break;
		case HEUR_HARD_FIX_70:
			prob_param = param[1];
			break;
		case HEUR_HARD_FIX_90:
			prob_param = param[0];
			break;
	}

	double residual_timelimit = inst->timelimit;
	double small_timelimit = 60.0;

	// Build an array of indices of variables/columns x(i,j) of the model
	int* indices = (int*)calloc(inst->ncols, sizeof(int));
	int k = 0;
	for (int i = 0; i < inst->ncols; i++) {
		indices[i] = i;
	}

	if (CPXsetintparam(env, CPXPARAM_Advance, 1)) print_error("CPXsetintparam() error in setting CPXPARAM_Advance");

	// Set to solve just the root node: by doing that we will get a feasible solution, but not the optimal one => good as a starting point for heuristics
	if (CPXsetintparam(env, CPX_PARAM_NODELIM, 0)) print_error("CPXsetintparam() error in setting seed");

	// Run TSP solver (which is the branch and cut with fractional subtour elimination constraints only applied on the root node)
	double t1 = second();
	inst->use_2_opt = 0;
	solve_branch_cut(inst, env, lp);				// solve the first time with nodelimit = 0 and using Branch & Cut with 2-opt refinement as TSP solver
	double t2 = second();

	int n_iter = 0;
	char lu = 'L';
	int beg = 0;
	double temp_obj_val = INFINITY;
	double low_bound_value = 1.0;
	double low_bound_reset = 0.0;
	double* curr_best_sol = (double*)calloc(inst->ncols, sizeof(double));

	if (CPXsetintparam(env, CPX_PARAM_NODELIM, 2100000000)) print_error("CPXsetintparam() error in setting seed");
	
	while (1) {
		
		int status = 0;
		if (status = CPXgetx(env, lp, curr_best_sol, 0, inst->ncols - 1)) {
			printf("Status: %d", status);
			print_error("CPXgetx() error in hard_fix");
		}

		if (inst->verbose >= MEDIUM) printf("Time used for iteration number %d: %f\n", n_iter, t2 - t1);

		CPXgetobjval(env, lp, &temp_obj_val);
		if (inst->z_best > temp_obj_val) {
			inst->z_best = temp_obj_val;

			// Set the feasible solution (not optimal) from which to start with the heuristics
			if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, indices, curr_best_sol, CPX_MIPSTART_AUTO, NULL)) print_error("CPXaddmipstarts() error in setting known solution");
		}
		else if (inst->model_type == HEUR_HARD_FIX_VAR) {
			// If there has been no improvement in the last iteration and we are using HARD_FIX heuristics with variable prob_param, then decrease prob_param value
			if (param_id < 2)
				param_id++;
			prob_param = param[param_id];
		}

		if (inst->verbose >= HIGH) printf("inst->z_best: %f\n", inst->z_best);
		if (inst->verbose == MEDIUM) printf("Iteration: %d, param: %f, best_sol: %f\n", n_iter, temp, inst->z_best);
		temp = prob_param;

		// Update the amount of time left before timelimit is reached and provide it to Cplex to check
		residual_timelimit = residual_timelimit - (t2 - t1);
		// If the timelimit is reached for the current iteration => exit the loop
		if (residual_timelimit <= 0) {
			if (inst->verbose >= LOW) printf("TOTAL time limit reached\n");
			break;
		}

		if (residual_timelimit <= small_timelimit) {
			if (CPXsetdblparam(env, CPX_PARAM_TILIM, residual_timelimit)) print_error("CPXsetdblparam() error in setting timelimit");
		}
		else { 
			if (CPXsetdblparam(env, CPX_PARAM_TILIM, small_timelimit)) print_error("CPXsetdblparam() error in setting timelimit");
		}

		for (int i = 0; i < inst->ncols; i++) {
			if (CPXchgbds(env, lp, 1, &i, &lu, &low_bound_reset)) print_error("CPXchgbds() error in setting edge lower bound");
		}

		int n_fixed_edges = 0;
		for (int i = 0; i < inst->ncols; i++) {
			if (curr_best_sol[i] > 0.5 && ((double)rand() / RAND_MAX) <= prob_param) {
				if (CPXchgbds(env, lp, 1, &i, &lu, &low_bound_value)) print_error("CPXchgbds() error in setting edge lower bound");
				n_fixed_edges++;
			}
		}
		if (inst->verbose >= HIGH) printf("n_fixed_edges: %d\n", n_fixed_edges);

		// Run TSP solver (which is the branch and cut with fractional subtour elimination constraints only applied on the root node)
		t1 = second();
		inst->use_2_opt = 0;
		solve_branch_cut(inst, env, lp);
		t2 = second();

		// Increment current iteration number
		n_iter++;
	}

	free(indices);
	free(curr_best_sol);

	return;
}


void solve_heur_soft_fix(instance* inst, CPXENVptr env, CPXLPptr lp) {

	int param[] = { 3, 5, 7, 9 };
	int K = param[0];
	int temp = K;
	switch (inst->model_type) {

		case HEUR_SOFT_FIX_3:
			K = param[0];
			temp = K;
			break;
		case HEUR_SOFT_FIX_5:
			K = param[1];
			temp = K;
			break;
		case HEUR_SOFT_FIX_7:
			K = param[2];
			temp = K;
			break;
		case HEUR_SOFT_FIX_9:
			K = param[3];
			temp = K;
			break;
	}

	double residual_timelimit = inst->timelimit;
	double small_timelimit = 60.0;

	// Build an array of indices of variables/columns x(i,j) of the model
	int* indices = (int*)calloc(inst->ncols, sizeof(int));
	for (int i = 0; i < inst->ncols; i++) indices[i] = i;

	// Ask CPLEX to continue with a partially explored MIP tree, if one is available.
	if (CPXsetintparam(env, CPXPARAM_Advance, 1)) print_error("CPXsetintparam() error in setting CPXPARAM_Advance");

	// Set to solve just the root node: by doing that we will get a feasible solution, but not the optimal one => good as a starting point for heuristics
	//if (CPXsetintparam(env, CPX_PARAM_NODELIM, 0)) print_error("CPXsetintparam() error in setting node limit");
	if (CPXsetdblparam(env, CPX_PARAM_TILIM, residual_timelimit / 3)) print_error("CPXsetdblparam() error in setting timelimit");

	// Run TSP solver (which is the branch and cut with fractional subtour elimination constraints only applied on the root node)
	double t1 = second();
	inst->use_2_opt = 1;
	solve_branch_cut(inst, env, lp);						// solve the first time with nodelimit = 0 and using Branch & Cut with 2-opt refinement as TSP solver
	double t2 = second();
	if (inst->verbose >= MEDIUM) printf("Time used for iteration number 0: %f\n", t2 - t1);
	
	residual_timelimit = residual_timelimit - (t2 - t1);
	t1 = second();

	int n_iter = 0;
	int beg = 0;
	double temp_obj_val = INFINITY;
	double* curr_best_sol = (double*)calloc(inst->ncols, sizeof(double));
	inst->best_sol = (double*)calloc(inst->ncols, sizeof(double));

	int* index = (int*)calloc(inst->nnodes, sizeof(int));							// Array of indexes associated to the row variables
	double* value = (double*)calloc(inst->nnodes, sizeof(double));					// Array of row variables coefficients
	char** cname = (char**)calloc(1, sizeof(char*));								// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	//if (CPXsetintparam(env, CPX_PARAM_NODELIM, 2100000000)) print_error("CPXsetintparam() error in setting node limit");

	while (1) {

		int status = 0;
		if (status = CPXgetx(env, lp, curr_best_sol, 0, inst->ncols - 1)) {
			printf("Status: %d", status);
			print_error("CPXgetx() error in soft_fix");
		}

		CPXgetobjval(env, lp, &temp_obj_val);
		if (inst->z_best > temp_obj_val) {
			inst->z_best = temp_obj_val;
			
			// Set the feasible solution (not optimal) from which to start with the heuristics
			if (CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, indices, inst->best_sol, CPX_MIPSTART_AUTO, NULL)) print_error("CPXaddmipstarts() error in setting known solution");
		}
		else {
			// If there has been no improvement in the last iteration and we are using SOFT_FIX heuristics with variable K, then increase K value up to 10
			if (K < 10) K++;
		}

		if (inst->verbose >= HIGH) printf("inst->z_best: %f\n", inst->z_best);
		if (inst->verbose == MEDIUM) printf("Iteration: %d, K: %d, best_sol: %f\n", n_iter, temp, inst->z_best);
		temp = K;

		t2 = second();

		if (inst->verbose >= MEDIUM && n_iter >= 1) printf("Time used for iteration number %d: %f\n", n_iter, t2 - t1);

		// After the first iteration, remove the previous local branching constraint
		if (n_iter >= 1) {
			int last;
			if ((last = CPXgetnumrows(env, lp)) == 0) print_error("CPXgetnumrows() error!");
			if (CPXdelrows(env, lp, last - 1, last - 1)) print_error("CPXdelrows() error in deleting local branching constraints!");
		}

		// Update the amount of time left before timelimit is reached and provide it to Cplex to check
		residual_timelimit = residual_timelimit - (t2 - t1);
		// If the timelimit is reached for the current iteration => exit the loop
		if (residual_timelimit <= 0) {
			if (inst->verbose >= LOW) printf("TOTAL time limit reached\n");
			memcpy(inst->best_sol, curr_best_sol, sizeof(double) * inst->ncols);
			break;
		}

		t1 = second();

		if (residual_timelimit <= small_timelimit) {
			if (CPXsetdblparam(env, CPX_PARAM_TILIM, residual_timelimit)) print_error("CPXsetdblparam() error in setting timelimit");
		}
		else {
			if (CPXsetdblparam(env, CPX_PARAM_TILIM, small_timelimit)) print_error("CPXsetdblparam() error in setting timelimit");
		}

		// Add the new local branching constraint
		int i_zero = 0;
		double rhs = inst->nnodes - K;														// "rhs" = right-hand side of the degree constraints
		char sense = 'G';																	// 'G' for greater-or-equal constraint
		int nnz = inst->nnodes;																// number of cells != 0 in the row
		sprintf(cname[0], "local_branching_constraint_(%d)", n_iter);						// Set a name for the new row/constraint
		int j = 0;
		for (int i = 0; i < inst->ncols; i++) {
			if (curr_best_sol[i] > 0.5) {
				index[j] = i;
				value[j] = 1.0;
				j++;
			}
		}
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("CPXaddrows() error in adding local branching constraints!");

		// Uncomment the following to debug adding/removing local branching constraints:
		/*char lp_temp_name[50];
		sprintf(lp_temp_name, "DEBUG_SOFT_FIX_%d", n_iter);
		create_lp_file(inst, env, lp, lp_temp_name);*/

		// Run the previously used TSP solver again
		inst->use_2_opt = 1;
		solve_branch_cut(inst, env, lp);

		n_iter++;
	}

	free(index);
	free(value);
	free(cname[0]);
	free(cname);
	free(curr_best_sol);
	free(indices);

	return;
}


void solve_heur_greedy(instance* inst, double* x) {

	// Array with all precomputed edges length = distances between all pairs of nodes
	double* edges_length = (double*)malloc(inst->ncols * sizeof(double));			
	int k = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			edges_length[k++] = dist(i, j, inst);
		}
	}

	// Array to record which nodes have already been visited
	int* nodes_visited = (int*)calloc(inst->nnodes, sizeof(int));

	// Array to memorize temporarily the new solution before comparing it to the best one (collected in array x)
	double* temp_x = (double*)calloc(inst->ncols, sizeof(double));

	int best_start_node = -1;
	double min_sol_cost = INFINITY;
	for (int start_node = 0; start_node < inst->nnodes; start_node++) {			// Try multiple runs starting from all possible nodes (affordable because it is very fast)

		double curr_sol_cost = 0.0;												// Restore the cost of the current solution before starting the new run
		for (int j = 0; j < inst->nnodes; j++) nodes_visited[j] = 0;			// Restore the array of visited nodes before starting the new run
		for (int j = 0; j < inst->ncols; j++) temp_x[j] = 0.0;					// Restore the array that will contain the new solution before starting the new run

		// Get the greedy solution (collected in array temp_x) associated to the current start_node
		int i = start_node;
		nodes_visited[start_node] = 1;											// Set the starting node as visited => it will not be reached automatically by the algorithm
		int next_node_index = -1;
		int curr_edge_index = -1;
		for (int n = 0; n < inst->nnodes - 1; n++) {							// Repeat until all edges apart from the "closing loop" one are found

			double min_dist = INFINITY;											// Reset the minimum distance as the highest possibile value
			for (int j = 0; j < inst->nnodes; j++) {							// Try all possible "landing nodes" from node i
				if (i != j && nodes_visited[j] == 0) {							// Avoid loop edges and check that the new node has not already been visited
					double curr_dist = edges_length[xpos(i, j, inst)];			// Retrieve the distance between the two current nodes from the array computed before
					if (curr_dist < min_dist) {
						min_dist = curr_dist;
						next_node_index = j;									// After the last iteration, next_node_index will be the best candidate node to move on
					}
				}
			}
			if (inst->verbose >= HIGH) printf("Edge #%d : [ %d -> %d ]\n", n, i + 1, next_node_index + 1);

			nodes_visited[next_node_index] = 1;									// Set the node on which we just landed as visited
			curr_edge_index = xpos(i, next_node_index, inst);
			curr_sol_cost += edges_length[curr_edge_index];						// Add up the cost of the edge just found
			temp_x[curr_edge_index] = 1.0;										// Set the edge just found as part of the solution

			i = next_node_index;												// Move to the next node and start searching again from it
		}

		int last_edge_index = xpos(i, start_node, inst);
		curr_sol_cost += edges_length[last_edge_index];							// Add up the cost of the last "closing loop" edge
		temp_x[last_edge_index] = 1.0;											// Add the last "closing loop" edge manually

		if (inst->verbose >= HIGH) printf("Cost of the solution with starting node #%d : %f\n", start_node + 1, curr_sol_cost);

		// Compare the new solution cost with the (currently) best one, then memorize the best solution and the start_node that lead to it
		if (curr_sol_cost < min_sol_cost) {
			min_sol_cost = curr_sol_cost;
			best_start_node = start_node;
			for (int j = 0; j < inst->ncols; j++) x[j] = temp_x[j];
		}
	}

	if (inst->verbose >= MEDIUM) printf("\nCost of the best solution: %f\nStarting node of the best solution: %d\n\n", min_sol_cost, best_start_node + 1);

	free(temp_x);
	free(nodes_visited);
	free(edges_length);

	return;
}


void solve_heur_grasp_greedy(instance* inst, double* x, double max_time) {

	double residual_timelimit = max_time;

	// Array with all precomputed edges length = distances between all pairs of nodes
	double* edges_length = (double*)malloc(inst->ncols * sizeof(double));
	int k = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			edges_length[k++] = dist(i, j, inst);
		}
	}

	// Array to record which nodes have already been visited
	int* nodes_visited = (int*)calloc(inst->nnodes, sizeof(int));

	// Array to memorize temporarily the new solution before comparing it to the best one (collected in array x)
	double* temp_x = (double*)calloc(inst->ncols, sizeof(double));

	int best_start_node = -1;
	double min_sol_cost = INFINITY;
	while(1) {																	// Try as many runs as made possible by max_time timelimit

		double t1 = second();

		int start_node = rand() % inst->nnodes;									// Choose randomly a starting node

		double curr_sol_cost = 0.0;												// Restore the cost of the current solution before starting the new run
		for (int j = 0; j < inst->nnodes; j++) nodes_visited[j] = 0;			// Restore the array of visited nodes before starting the new run
		for (int j = 0; j < inst->ncols; j++) temp_x[j] = 0.0;					// Restore the array that will contain the new solution before starting the new run

		// Get a greedy solution (collected in array temp_x) associated to the current start_node and to the random decisions taken during the building of the solution
		int i = start_node;
		nodes_visited[start_node] = 1;											// Set the starting node as visited => it will not be reached automatically by the algorithm
		int next_node_index = -1;
		int curr_edge_index = -1;
		for (int n = 0; n < inst->nnodes - 1; n++) {							// Repeat until all edges apart from the "closing loop" one are found

			int first_refused_node = -1;
			int second_refused_node = -1;
			int node_accepted = 0;												// Flag to indicate if the candidate node found hat to be accepted or not
			while (!node_accepted) {											// Repeat until the candidate node (either the 1st, 2nd, 3rd best one) is accepted	

				double min_dist = INFINITY;										// Reset the minimum distance as the highest possibile value
				for (int j = 0; j < inst->nnodes; j++) {						// Try all possible "landing nodes" from node i
					if (i != j &&
						first_refused_node != j &&
						second_refused_node != j &&
						nodes_visited[j] == 0) {								// Avoid loop edges and check that the new node has not already been visited

						double curr_dist = edges_length[xpos(i, j, inst)];		// Retrieve the distance between the two current nodes from the array computed before
						if (curr_dist < min_dist) {
							min_dist = curr_dist;
							next_node_index = j;								// After the last iteration, next_node_index will be the best candidate node to move on
						}
					}
				}

				// Once the best candidate node to move on is found, choose randomly to accept it or to choose the second or third best candidate
				if (first_refused_node == -1) {									// If the first candidate node has not been refused, refuse it with 10% probability
					if (((double)rand() / RAND_MAX) <= 0.1) {
						first_refused_node = next_node_index;
						node_accepted = 0;
					}
					else {														// With 90% probability just accept the first candidate node
						node_accepted = 1;
						if (inst->verbose >= HIGH) printf("First candidate node has been accepted!\n");
					}
				}
				else {															// If the first candidate node has already been refused...
					if (second_refused_node == -1) {							// ...then if the second candidate node has not been refused, refuse it with 10% probability
						if (((double)rand() / RAND_MAX) <= 0.1) {
							second_refused_node = next_node_index;
							node_accepted = 0;
						}
						else {													// With 90% probability just accept the second candidate node
							node_accepted = 1;
							if (inst->verbose >= HIGH) printf("Second candidate node has been accepted!\n");
						}
					}
					else {														// If even the second candidate node has already been refused...
						node_accepted = 1;										// ...then just accept the third candidate node
						if (inst->verbose >= HIGH) printf("Third candidate node has been accepted!\n");
					}
				}
			}

			if (inst->verbose >= HIGH) printf("Edge #%d : [ %d -> %d ]\n", n, i + 1, next_node_index + 1);

			nodes_visited[next_node_index] = 1;									// Set the node on which we just landed as visited
			curr_edge_index = xpos(i, next_node_index, inst);
			curr_sol_cost += edges_length[curr_edge_index];						// Add up the cost of the edge just found
			temp_x[curr_edge_index] = 1.0;										// Set the edge just found as part of the solution

			i = next_node_index;												// Move to the next node and start searching again from it
		}

		int last_edge_index = xpos(i, start_node, inst);
		curr_sol_cost += edges_length[last_edge_index];							// Add up the cost of the last "closing loop" edge
		temp_x[last_edge_index] = 1.0;											// Add the last "closing loop" edge manually

		// Compare the new solution cost with the (currently) best one, then memorize the best solution
		if (curr_sol_cost < min_sol_cost) {
			min_sol_cost = curr_sol_cost;
			for (int j = 0; j < inst->ncols; j++) x[j] = temp_x[j];
		}

		double t2 = second();

		// Update the residual timelimit and check if it has been reached (if so, stop generating new solutions and keep the best one found up to now)
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (residual_timelimit <= 0) break;
	}

	if (inst->verbose >= MEDIUM) printf("\nHEUR_GRASP -> Cost of the best solution: %f\n\n", min_sol_cost);
	inst->z_best = min_sol_cost;

	free(temp_x);
	free(nodes_visited);
	free(edges_length);

	return;
}


void solve_heur_extra_mileage(instance* inst, double* x) {

	// Fill the following data structure with the coordinates of nodes from the instance
	Point* instance_nodes = (Point*)calloc(inst->nnodes, sizeof(Point));
	for (int i = 0; i < inst->nnodes; i++) {
		instance_nodes[i].id = i;
		instance_nodes[i].x = inst->xcoord[i];
		instance_nodes[i].y = inst->ycoord[i];
	}
	
	// Compute the convex hull of the set of nodes
	int hull_size;
	Point* hull = convexHull(instance_nodes, inst->nnodes, &hull_size);

	if (inst->verbose >= HIGH) {
		printf("Nodes of the convex hull: ");
		printPoints(hull, hull_size);
		printf("\n");
	}

	/*
	// Since the convex hull is returned as a list of consecutive nodes arranged counter-clockwise,
	// we just needed to connect all of them to get the first part of the solution.
	for (int i = 0; i <= hull_size - 2; i++) {
		x[xpos(hull[i].id, hull[i + 1].id, inst)] = 1.0;				// Add all the edges connecting the nodes of the hull, apart from the one from last to first node
	}
	x[xpos(hull[hull_size - 1].id, hull[0].id, inst)] = 1.0;			// Then add the last edge from the last to the first node
	*/
	
	// Build the list of "successors" of each node. If a node is not part of the solution yet, then its successor will be -1
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++) succ[i] = -1;
	for (int i = 0; i <= hull_size - 2; i++) succ[hull[i].id] = hull[i + 1].id;
	succ[hull[hull_size - 1].id] = hull[0].id;

	if (inst->verbose >= HIGH) {
		printf("List of successors: [ ");
		for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, succ[i]);
		printf("%d->%d ]\n", inst->nnodes - 1, succ[inst->nnodes - 1]);
	}

	// Now we start inserting to the partial solution, at each iteration, the node which lead to the lowest "extra mileage" (as a locally optimal choice)
	int n_nodes_to_insert = inst->nnodes - hull_size;
	for (int n_iter = 1; n_iter <= n_nodes_to_insert; n_iter++) {

		int curr_node_index = -1;
		int next_node_index = -1;
		for (int i = 0; i < inst->nnodes; i++) {
			if (succ[i] != -1) {										// The node i has already been inserted => edge (i, succ[i]) good candidate to be removed

				double min_x_mil = INFINITY;
				for (int j = 0; j < inst->nnodes; j++) {
					if (succ[j] == -1) {								// The node j has not been inserted yet => good candidate node to be inserted

						double x_mil = extra_mileage(i, succ[i], j, inst);
						if (x_mil < min_x_mil) {
							min_x_mil = x_mil;
							curr_node_index = i;
							next_node_index = j;						// After the last iteration, next_node_index will be the best candidate node to insert
						}
					}
				}
			}
		}

		// Update the successors list accordingly
		succ[next_node_index] = succ[curr_node_index];
		succ[curr_node_index] = next_node_index;
	}

	if (inst->verbose >= HIGH) {
		printf("FINAL list of successors: [ ");
		for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, succ[i]);
		printf("%d->%d ]\n", inst->nnodes - 1, succ[inst->nnodes - 1]);
	}

	// Fill the solution array "x" from successors list
	for (int i = 0; i <= inst->nnodes - 1; i++) {
		x[xpos(i, succ[i], inst)] = 1.0;
	}

	// Compute from scratch the cost of the solution
	double curr_sol_cost = 0.0;
	for (int i = 0; i < inst->nnodes; i++) curr_sol_cost += dist(i, succ[i], inst);

	if (inst->verbose >= MEDIUM) printf("\nHEUR_EXTRA_MILEAGE -> Cost of the best solution: %f\n\n", curr_sol_cost);
	inst->z_best = curr_sol_cost;

	free(succ);
	free(hull);
	free(instance_nodes);

	return;
}


void solve_heur_grasp_extra_mileage(instance* inst, double* x, double max_time) {

	double residual_timelimit = max_time;

	// Fill the following data structure with the coordinates of nodes from the instance
	Point* instance_nodes = (Point*)calloc(inst->nnodes, sizeof(Point));
	for (int i = 0; i < inst->nnodes; i++) {
		instance_nodes[i].id = i;
		instance_nodes[i].x = inst->xcoord[i];
		instance_nodes[i].y = inst->ycoord[i];
	}

	// Compute the convex hull of the set of nodes
	int hull_size;
	Point* hull = convexHull(instance_nodes, inst->nnodes, &hull_size);

	// Build the list of "successors" of each node. If a node is not part of the solution yet, then its successor will be -1
	int* hull_succ = (int*)malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++) hull_succ[i] = -1;
	for (int i = 0; i <= hull_size - 2; i++) hull_succ[hull[i].id] = hull[i + 1].id;
	hull_succ[hull[hull_size - 1].id] = hull[0].id;

	if (inst->verbose >= HIGH) {
		printf("Convex hull list of successors: [ ");
		for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, hull_succ[i]);
		printf("%d->%d ]\n", inst->nnodes - 1, hull_succ[inst->nnodes - 1]);
	}

	// Now we start inserting to the partial solution, at each iteration, the node which lead to the lowest "extra mileage" (as a locally optimal choice)
	int n_nodes_to_insert = inst->nnodes - hull_size;
	int* temp_succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* final_succ = (int*)malloc(inst->nnodes * sizeof(int));

	double min_sol_cost = INFINITY;
	int n_runs = 1;
	while (1) {																		// Try as many runs as made possible by max_time timelimit

		if (inst->verbose >= HIGH) printf("Restarting from convex hull: run #%d\n", n_runs);

		double t1 = second();

		temp_succ = memcpy(temp_succ, hull_succ, sizeof(int) * inst->nnodes);		// Start from the original succ data structure with the same convex hull

		for (int n_iter = 1; n_iter <= n_nodes_to_insert; n_iter++) {				// Repeat until a complete solution is built

			int curr_node_index = -1;
			int next_node_index = -1;

			// Choose randomly to get the pair of nodes that lead to minimum "extra mileage" or to get them randomly (with 1% probability)
			if (((double)rand() / RAND_MAX) <= 0.01) {
				do {																// We  need that the chosen starting node has already been inserted (es. part of the hull)
					curr_node_index = rand() % inst->nnodes;
				} while (temp_succ[curr_node_index] == -1);
				do {																// We need that the chosen landing node has not already been inserted
					next_node_index = rand() % inst->nnodes;
				} while (temp_succ[next_node_index] != -1);
			}
			else {
				for (int i = 0; i < inst->nnodes; i++) {
					if (temp_succ[i] != -1) {										// The node i has already been inserted => edge (i, succ[i]) good candidate to be removed

						double min_x_mil = INFINITY;
						for (int j = 0; j < inst->nnodes; j++) {
							if (temp_succ[j] == -1) {								// The node j has not been inserted yet nor already discarded => good candidate node to be inserted

								double x_mil = extra_mileage(i, temp_succ[i], j, inst);
								if (x_mil < min_x_mil) {
									min_x_mil = x_mil;
									curr_node_index = i;
									next_node_index = j;							// After the last iteration, next_node_index will be the best candidate node to insert
								}
							}
						}
					}
				}
			}

			// Update the successors list accordingly
			temp_succ[next_node_index] = temp_succ[curr_node_index];
			temp_succ[curr_node_index] = next_node_index;
		}

		// Compute from scratch the cost of the new solution
		double curr_sol_cost = 0.0;
		for (int i = 0; i < inst->nnodes; i++) curr_sol_cost += dist(i, temp_succ[i], inst);
		if (inst->verbose >= HIGH) printf("New solution cost: %f\n", curr_sol_cost);

		// If the newly built solution is better than all the previous ones, then memorize it
		if (curr_sol_cost < min_sol_cost) {
			min_sol_cost = curr_sol_cost;
			final_succ = memcpy(final_succ, temp_succ, sizeof(int) * inst->nnodes);

			if (inst->verbose >= MEDIUM) printf("Found an improved solution of cost: %f at run #%d\n", curr_sol_cost, n_runs);
		}

		if (inst->verbose >= HIGH) {
			printf("New solution list of successors: [ ");
			for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, temp_succ[i]);
			printf("%d->%d ]\n\n", inst->nnodes - 1, temp_succ[inst->nnodes - 1]);
		}

		double t2 = second();

		// Update the residual timelimit and check if it has been reached (if so, stop generating new solutions and keep the best one found up to now)
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (residual_timelimit <= 0) break;

		n_runs++;
	}

	// Fill the solution array "x" from successors list
	for (int i = 0; i <= inst->nnodes - 1; i++) {
		x[xpos(i, final_succ[i], inst)] = 1.0;
	}
	
	if (inst->verbose >= MEDIUM) printf("\nHEUR_GRASP_EXTRA_MILEAGE -> Cost of the best solution: %f\n\n", min_sol_cost);
	inst->z_best = min_sol_cost;
	
	free(final_succ);
	free(temp_succ);
	free(hull_succ);
	free(hull);

	return;
}


void compute_succ(instance* inst, double* x, int* succ) {

	for (int i = 0; i < inst->nnodes; i++) succ[i] = -1;
	int i = 0;
	for (int k = 0; k < inst->nnodes - 1; k++) {
		for (int j = 1; j < inst->nnodes; j++) {
			if (i != j && succ[j] == -1 && x[xpos(i, j, inst)] > 0.5) {
				succ[i] = j;
				i = j;
				break;
			}
		}
	}
	succ[i] = 0;
}


void compute_best_node(instance* inst, int* succ, int* best_a, int* best_b, double* min_delta_cost) {

	for (int a = 0; a < inst->nnodes; a++) {
		for (int b = 0; b < inst->nnodes; b++) {

			// The two selected nodes must be non-consecutive
			if (a != b && b != succ[a] && a != succ[b]) {
				double curr_delta_cost = (dist(a, b, inst) + dist(succ[a], succ[b], inst)) - (dist(a, succ[a], inst) + dist(b, succ[b], inst));
				if (curr_delta_cost < *min_delta_cost) {
					*min_delta_cost = curr_delta_cost;
					*best_a = a;
					*best_b = b;
				}
			}
			
		}
	}
}


void _2opt_move(instance* inst, int a, int b, int* succ) {
	
	// Since the successors of best_a and best_b will soon be updated, memorize the original ones for later use
	int old_succ_a = succ[a];
	int old_succ_b = succ[b];
	// Change verse to all edges between nodes old_succ_a and b
	int temp_curr = old_succ_a;										// Start from node old_succ_a
	int temp_succ = -1;
	while (temp_curr != b) {										// Continue until node best_b is reached
		if (temp_succ == -1) temp_succ = succ[temp_curr];			// N.B. both successor (temp_succ) and successor of successor (temp_succ_succ)
		int temp_succ_succ = succ[temp_succ];						// of temp_curr have to be memorized for the next iteration!

		succ[temp_succ] = temp_curr;								// Change the verse of the edge (temp_curr -> temp_succ)

		temp_curr = temp_succ;										// Shift one position towards node best_b
		temp_succ = temp_succ_succ;

	}

	// Rearrange the connections between nodes a and b
	succ[old_succ_a] = old_succ_b;
	succ[a] = b;
	
}


void solve_heur_2_opt(instance* inst, double* x, int* succ, double timelimit) {

	if (succ == NULL && x == NULL) print_error("solve_heur_2_opt() cannot be called with both x and succ values set to NULL!");

	// If not provided, compute the list of successors from the solution provided by GRASP constructive heuristic
	int is_succ_provided = 1;
	if (succ == NULL) {
		is_succ_provided = 0;
		succ = (int*)malloc(inst->nnodes * sizeof(int));
		compute_succ(inst, x, succ);
	}

	// Retrieve the cost of the solution provided by GRASP or by an iteration of VNS
	double curr_sol_cost = inst->z_best;
	double residual_timelimit = timelimit;

	if (inst->verbose >= HIGH) {
		printf("STARTING list of successors: [ ");
		for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, succ[i]);
		printf("%d->%d ]\n", inst->nnodes - 1, succ[inst->nnodes - 1]);
	}
	
	// Find the pair of edges (associated to their two starting nodes a and b) whose substitution leads to 
	// the maximum "cut" in current solution cost, which corresponds to the most negative value of curr_delta_cost
	int n_iter = 1;
	while (1) {

		double t1 = second();

		if (inst->verbose >= HIGH) printf("n_iter: %d\n", n_iter);

		int best_a = -1;
		int best_b = -1;
		double min_delta_cost = INFINITY;

		compute_best_node(inst, succ, &best_a, &best_b, &min_delta_cost);

		// When the new solution cost is no longer able to become lower => local optimal solution achieved!
		if (min_delta_cost >= 0) break;

		// If min_delta_cost is negative, lower the current solution cost by its value
		curr_sol_cost += min_delta_cost;

		if (inst->verbose >= HIGH) {
			printf("best_a: %d\n", best_a);
			printf("best_b: %d\n", best_b);
		}

		_2opt_move(inst, best_a, best_b, succ);
		
		if (inst->verbose >= HIGH) {
			printf("IMPROVED list of successors: [ ");
			for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, succ[i]);
			printf("%d->%d ]\n", inst->nnodes - 1, succ[inst->nnodes - 1]);
		}
		
		double t2 = second();

		// Update the timelimit and check if it has been reached for this call of 2-opt
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (residual_timelimit <= 0) {
			if (inst->verbose >= MEDIUM) printf("Timelimit reached after %f secs!\n", timelimit);
			break;
		}

		n_iter++;
	}
	
	if (inst->verbose >= HIGH) {
		printf("FINAL list of successors: [ ");
		for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, succ[i]);
		printf("%d->%d ]\n", inst->nnodes - 1, succ[inst->nnodes - 1]);
	}

	// Remove the values of the solution provided by GRASP and insert the 2-opt improved one
	if (x != NULL) {
		for (int i = 0; i < inst->ncols; i++) x[i] = 0.0;
		for (int i = 0; i <= inst->nnodes - 1; i++) {
			x[xpos(i, succ[i], inst)] = 1.0;
		}
	}

	if (inst->verbose >= MEDIUM) printf("HEUR_2_OPT -> Cost of the best solution (after 2-opt refinement): %f\n\n", curr_sol_cost);
	inst->z_best = curr_sol_cost;

	if (!is_succ_provided) free(succ);

	return;
}


void solve_heur_multi_start(instance* inst, double* x) {

	double residual_timelimit = inst->timelimit;
	double best_sol_cost = INFINITY;

	// Array to memorize temporarily the new solution before comparing it to the best one (collected in array x)
	double* temp_x = (double*)calloc(inst->ncols, sizeof(double));

	int n_iter = 1;
	while (1) {

		double grasp_timelimit = (double)(((double)rand() / RAND_MAX) * 10.0);

		double t1 = second();
		solve_heur_grasp_greedy(inst, temp_x, grasp_timelimit);					// First let GRASP find an always different reference solution (no need to reset solution array temp_x)
		solve_heur_2_opt(inst, temp_x, NULL, INFINITY);							// Then use 2-opt refinement heuristic to improve the reference solution
		double t2 = second();

		// Keep the best solution found up to now
		if (inst->z_best < best_sol_cost) {
			best_sol_cost = inst->z_best;
			for (int j = 0; j < inst->ncols; j++) x[j] = temp_x[j];
		}

		if (inst->verbose >= MEDIUM) printf("Time used for iteration number %d (HEUR_GRASP + HEUR_2_OPT): %f\n", n_iter, t2 - t1);

		// Update the amount of time left before timelimit is reached and provide it to Cplex to check
		residual_timelimit = residual_timelimit - (t2 - t1);
		// If the timelimit is reached for the current iteration => exit the loop
		if (residual_timelimit <= 0) {
			if (inst->verbose >= LOW) printf("TOTAL time limit reached\n");
			break;
		}

		n_iter++;
	}

	if (inst->verbose >= MEDIUM) printf("HEUR_MULTI_START -> Cost of the best solution: %f\n\n", best_sol_cost);
	inst->z_best = best_sol_cost;

	free(temp_x);

	return;
}


void _apply_rnd_3_opt_move(instance* inst, int* succ) {

	// Select 3 random edges (corresponding to as many starting nodes) to be removed
	int a, b, c = -1;
	a = rand() % inst->nnodes;
	while ((b = rand() % inst->nnodes) == a);
	while ((c = rand() % inst->nnodes) == a || c == b);

	if (inst->verbose >= HIGH) printf("a = %d, b = %d, c = %d\n", a, b, c);
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d\n", succ[a], succ[b], succ[c]);

	// Make sure that nodes (a, b, c) are chosen in this order inside the "succ" data structure
	int n_nodes_chosen = 0;
	int temp_a = a;
	int temp_b = b;
	int temp_c = c;
	int curr_node = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		int is_curr_node_chosen = (curr_node == temp_a || curr_node == temp_b || curr_node == temp_c);
		if (n_nodes_chosen == 0 && is_curr_node_chosen) {
			a = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 1 && is_curr_node_chosen) {
			b = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 2 && is_curr_node_chosen) {
			c = curr_node;
		}
		curr_node = succ[curr_node];
	}
	if (inst->verbose >= HIGH) printf("a = %d, b = %d, c = %d\n", a, b, c);
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d\n", succ[a], succ[b], succ[c]);

	// Perform a 3-opt worsening move ("kick" the solution out of the local optimum)
	int temp_succ_a = succ[a];
	succ[a] = succ[b];
	int temp_succ_c = succ[c];
	succ[c] = temp_succ_a;
	succ[b] = temp_succ_c;
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d\n", succ[a], succ[b], succ[c]);

	return;
}


void _apply_rnd_5_opt_move(instance* inst, int* succ) {

	// Select 5 random edges (corresponding to as many starting nodes) to be removed
	int a, b, c, d, e = -1;
	a = rand() % inst->nnodes;
	while ((b = rand() % inst->nnodes) == a);
	while ((c = rand() % inst->nnodes) == a || c == b);
	while ((d = rand() % inst->nnodes) == a || d == b || d == c);
	while ((e = rand() % inst->nnodes) == a || e == b || e == c || e == d);

	if (inst->verbose >= HIGH) printf("a = %d, b = %d, c = %d, d = %d, e = %d\n", a, b, c, d, e);
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d, succ_d = %d, succ_e = %d\n", succ[a], succ[b], succ[c], succ[d], succ[e]);

	// Make sure that nodes (a, b, c, d, e) are chosen in this order inside the "succ" data structure
	int n_nodes_chosen = 0;
	int temp_a = a;
	int temp_b = b;
	int temp_c = c;
	int temp_d = d;
	int temp_e = e;
	int curr_node = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		int is_curr_node_chosen = (curr_node == temp_a || curr_node == temp_b || curr_node == temp_c || curr_node == temp_d || curr_node == temp_e);
		if (n_nodes_chosen == 0 && is_curr_node_chosen) {
			a = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 1 && is_curr_node_chosen) {
			b = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 2 && is_curr_node_chosen) {
			c = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 3 && is_curr_node_chosen) {
			d = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 4 && is_curr_node_chosen) {
			e = curr_node;
		}
		curr_node = succ[curr_node];
	}
	if (inst->verbose >= HIGH) printf("a = %d, b = %d, c = %d, d = %d, e = %d\n", a, b, c, d, e);
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d, succ_d = %d, succ_e = %d\n", succ[a], succ[b], succ[c], succ[d], succ[e]);

	// Perform a 5-opt worsening move ("kick" the solution out of the local optimum)
	int temp_succ_a = succ[a];
	succ[a] = succ[c];
	int temp_succ_d = succ[d];
	succ[d] = temp_succ_a;
	int temp_succ_b = succ[b];
	succ[b] = temp_succ_d;
	int temp_succ_e = succ[e];
	succ[e] = temp_succ_b;
	succ[c] = temp_succ_e;
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d, succ_d = %d, succ_e = %d\n", succ[a], succ[b], succ[c], succ[d], succ[e]);

	return;
}


void _apply_rnd_7_opt_move(instance* inst, int* succ) {

	// Select 7 random edges (corresponding to as many starting nodes) to be removed
	int a, b, c, d, e, f, g = -1;
	a = rand() % inst->nnodes;
	while ((b = rand() % inst->nnodes) == a);
	while ((c = rand() % inst->nnodes) == a || c == b);
	while ((d = rand() % inst->nnodes) == a || d == b || d == c);
	while ((e = rand() % inst->nnodes) == a || e == b || e == c || e == d);
	while ((f = rand() % inst->nnodes) == a || f == b || f == c || f == d || f == e);
	while ((g = rand() % inst->nnodes) == a || g == b || g == c || g == d || g == e || g == f);

	if (inst->verbose >= HIGH) printf("a = %d, b = %d, c = %d, d = %d, e = %d, f = %d, g = %d\n", a, b, c, d, e, f, g);
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d, succ_d = %d, succ_e = %d, succ_f = %d, succ_g = %d\n", succ[a], succ[b], succ[c], succ[d], succ[e], succ[f], succ[g]);

	// Make sure that nodes (a, b, c, d, e, f, g) are chosen in this order inside the "succ" data structure
	int n_nodes_chosen = 0;
	int temp_a = a;
	int temp_b = b;
	int temp_c = c;
	int temp_d = d;
	int temp_e = e;
	int temp_f = f;
	int temp_g = g;
	int curr_node = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		int is_curr_node_chosen = (curr_node == temp_a || curr_node == temp_b || curr_node == temp_c || curr_node == temp_d || curr_node == temp_e || curr_node == temp_f || curr_node == temp_g);
		if (n_nodes_chosen == 0 && is_curr_node_chosen) {
			a = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 1 && is_curr_node_chosen) {
			b = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 2 && is_curr_node_chosen) {
			c = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 3 && is_curr_node_chosen) {
			d = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 4 && is_curr_node_chosen) {
			e = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 5 && is_curr_node_chosen) {
			f = curr_node;
			n_nodes_chosen++;
		}
		else if (n_nodes_chosen == 6 && is_curr_node_chosen) {
			g = curr_node;
		}
		curr_node = succ[curr_node];
	}
	if (inst->verbose >= HIGH) printf("a = %d, b = %d, c = %d, d = %d, e = %d, f = %d, g = %d\n", a, b, c, d, e, f, g);
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d, succ_d = %d, succ_e = %d, succ_f = %d, succ_g = %d\n", succ[a], succ[b], succ[c], succ[d], succ[e], succ[f], succ[g]);

	// Perform a 7-opt worsening move ("kick" the solution out of the local optimum)
	int temp_succ_a = succ[a];
	succ[a] = succ[d];
	int temp_succ_e = succ[e];
	succ[e] = temp_succ_a;
	int temp_succ_b = succ[b];
	succ[b] = temp_succ_e;
	int temp_succ_f = succ[f];
	succ[f] = temp_succ_b;
	int temp_succ_c = succ[c];
	succ[c] = temp_succ_f;
	int temp_succ_g = succ[g];
	succ[g] = temp_succ_c;
	succ[d] = temp_succ_g;
	if (inst->verbose >= HIGH) printf("succ_a = %d, succ_b = %d, succ_c = %d, succ_d = %d, succ_e = %d, succ_f = %d, succ_g = %d\n", succ[a], succ[b], succ[c], succ[d], succ[e], succ[f], succ[g]);

	return;
}

void solve_heur_vns(instance* inst, double* x) {

	double residual_timelimit = inst->timelimit;
	double min_sol_cost = INFINITY;


	// Use GRASP to generate the reference solution, allowing it to use up to 1/10 of the total timelimit
	double t1 = second();
	solve_heur_grasp_greedy(inst, x, residual_timelimit / 10);
	double t2 = second();

	// Update remaining timelimit, but don't check it has been reached because it should just be at about 1/10
	residual_timelimit = residual_timelimit - (t2 - t1);

	if (inst->verbose >= MEDIUM) printf("Starting incumbent cost (obtained from GRASP): %f\n\n", inst->z_best);

	int* succ = (int*)malloc(inst->nnodes * sizeof(int));

	t1 = second();

	while (1) {

		// Apply the 2-opt refinement heuristics to the current solution to get to a local optimum
		solve_heur_2_opt(inst, x, NULL, INFINITY);

		t2 = second();

		double curr_sol_cost = inst->z_best;
		//if (inst->verbose >= MEDIUM) printf("Local Optimum reached: %f\n", inst->z_best);
		
		// Compare the new solution cost with the (currently) best one, then memorize the best solution
		if (curr_sol_cost < min_sol_cost) {
			min_sol_cost = curr_sol_cost;
		}


		// Check if the timelimit has been reached: if so => exit the loop
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (residual_timelimit <= 0) break;
		
		t1 = second();

		// Compute the list of successors of the local optimum solution 
		compute_succ(inst, x, succ);

		if (inst->verbose >= HIGH) {
			printf("STARTING list of successors: [ ");
			for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, succ[i]);
			printf("%d->%d ]\n", inst->nnodes - 1, succ[inst->nnodes - 1]);
		}

		// Apply a worsening 3-opt move with 50% prob, a 5-opt move with 30% prob and a 7-opt move with 20% prob
		int n_move_chosen = -1;
		double rnd_value = ((double)rand() / RAND_MAX);
		if (rnd_value < 0.5) {
			_apply_rnd_3_opt_move(inst, succ);
			n_move_chosen = 3;
		}
		else if (rnd_value >= 0.5 && rnd_value < 0.8) {
			_apply_rnd_5_opt_move(inst, succ);
			n_move_chosen = 5;
		}
		else {
			_apply_rnd_7_opt_move(inst, succ);
			n_move_chosen = 7;
		}

		if (inst->verbose >= HIGH) {
			printf("POST %d-OPT MOVE list of successors: [ ", n_move_chosen);
			for (int i = 0; i <= inst->nnodes - 2; i++) printf("%d->%d, ", i, succ[i]);
			printf("%d->%d ]\n", inst->nnodes - 1, succ[inst->nnodes - 1]);
		}

		// Update the current solution x with the (3-opt, 5-opt, 7-opt) worsening move
		for (int i = 0; i < inst->ncols; i++) x[i] = 0.0;
		for (int i = 0; i < inst->nnodes; i++) x[xpos(i, succ[i], inst)] = 1.0;

		// Compute from scratch the cost of the new (worsened) solution
		curr_sol_cost = 0.0;
		for (int i = 0; i < inst->nnodes; i++) curr_sol_cost += dist(i, succ[i], inst);

		if (inst->verbose >= MEDIUM) printf("COST after a %d-opt kick: %f\n\n", n_move_chosen, curr_sol_cost);
		inst->z_best = curr_sol_cost;
	}

	if (inst->verbose >= MEDIUM) printf("\nHEUR_VNS -> Cost of the best solution: %f\n\n", min_sol_cost);
	inst->z_best = min_sol_cost;

	free(succ);
}


void solve_heur_tabu(instance* inst, double* x) {

	double residual_timelimit = inst->timelimit;
	double min_sol_cost = INFINITY;

	// Set a default value for the "tenure", which indicates for how many iterations a tabu rule is applied (then it is removed)
	int tenure = 20;									

	// Use GRASP to generate the reference solution, allowing it to use up to 1/10 of the total timelimit
	double t1 = second();
	solve_heur_grasp_greedy(inst, x, residual_timelimit / 10);
	double t2 = second();

	// Update remaining timelimit, but don't check it has been reached because it should just be at about 1/10
	residual_timelimit = residual_timelimit - (t2 - t1);

	if (inst->verbose >= MEDIUM) printf("Starting incumbent cost: %f\n\n", inst->z_best);

	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int* final_succ = (int*)malloc(inst->nnodes * sizeof(int));
	compute_succ(inst, x, succ);

	// Initialize the "tabu list", which memorizes at which iteration a node has been declared tabu/untouchable
	int* tabu = malloc(inst->nnodes * sizeof(int));
	for (int i = 0; i < inst->nnodes; i++) tabu[i] = -1;

	int n_iter = 0;
	int tenure_change_flag = 0;
	int improving_phase = 1;
	double curr_cost = inst->z_best;

	t1 = second();

	while (1) {

		while (1) {

			int a, b;
			double min_delta_cost = INFINITY;


			// Every 500 iterations change the tenure value (es. for 1000 nodes, it alternates between 20 and 100)
			if (n_iter % 500 == 0) {
				if (tenure_change_flag == 0) {
					tenure = inst->nnodes / 10;
					tenure_change_flag = 1;
				}
				else {
					tenure = inst->nnodes / 50;
					tenure_change_flag = 0;
				}
				if (inst->verbose >= MEDIUM) printf("Switching to TENURE = %d\n", tenure);
			}
			
			// Find the pair of edges (associated to their two starting nodes a and b) whose substitution leads to 
			// the maximum "cut" in current solution cost, which corresponds to the most negative value of curr_delta_cost
			compute_best_node(inst, succ, &a, &b, &min_delta_cost);

			// If one of the nodes a,b involved in the worsening move are tabu, then exit the inner loop
			if (improving_phase == 0 && tabu[a] != -1 && n_iter - tabu[a] <= tenure) {
				if (inst->verbose >= HIGH) printf("Node a = %d is TABU!\n", a);
				if (curr_cost < min_sol_cost) {
					min_sol_cost = curr_cost;
					final_succ = memcpy(final_succ, succ, sizeof(int) * inst->nnodes);
				}
				break;
			}
			if (improving_phase == 0 && tabu[b] != -1 && n_iter - tabu[b] <= tenure) {
				if (inst->verbose >= HIGH) printf("Node b = %d is TABU!\n", b);
				if (curr_cost < min_sol_cost) {
					min_sol_cost = curr_cost;
					final_succ = memcpy(final_succ, succ, sizeof(int) * inst->nnodes);
				}
				break;
			}

			// If min_delta_cost is a positive value, it means that there was no pair of edges that could lead to a cut in current solution cost
			// => after a certain number of (improving) 2-opt moves, a local optimum has been reached!
			if (min_delta_cost >= 0) {
				
				if (inst->verbose >= MEDIUM) printf("LOCAL optimum reached: %f\n", curr_cost);

				// If the new locally optimal solution is better than the global one, then replace it inside "final_succ"
				if (curr_cost < min_sol_cost) {
					min_sol_cost = curr_cost;
					final_succ = memcpy(final_succ, succ, sizeof(int) * inst->nnodes);
				}

				// Switch to "worsening phase", where we worsen the current solution with 2-opt random kicks, and exit the inner loop
				improving_phase = 0;
				break;
			}

			// Apply an improving 2-opt move and update the cost accordingly
			improving_phase = 1;
			curr_cost += min_delta_cost;
			_2opt_move(inst, a, b, succ);
			
			n_iter++;
		}

		t2 = second();

		// Check if timelimit has been reached: if so, exit the loop
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (residual_timelimit <= 0) break;

		t1 = second();

		// Perform a 2-opt random kick to avoid loop
		int node_a, node_b;
		
		node_a = rand() % inst->nnodes;
		while ((node_b = rand() % inst->nnodes) == node_a || node_b == succ[node_a]);
			
		// Apply a worsening move
		double worse_cost = (dist(node_a, node_b, inst) + dist(succ[node_a], succ[node_b], inst)) - (dist(node_a, succ[node_a], inst) + dist(node_b, succ[node_b], inst));
		curr_cost = curr_cost + worse_cost;
		_2opt_move(inst, node_a, node_b, succ);
		if (inst->verbose >= HIGH) printf("Cost after 2-OPT kick: %f\n", curr_cost);

		// Update the "tabu list" of nodes node_a and node_b, so that they will be tabu for "tenure" iterations from now
		tabu[node_a] = n_iter;
		tabu[node_b] = n_iter;

		n_iter++;
	}

	if (inst->verbose >= MEDIUM) printf("\nHEUR_TABU -> Cost of the best solution: %f\n\n", min_sol_cost);
	inst->z_best = min_sol_cost;
	
	// Write the solution memorized in final_succ inside x
	for (int i = 0; i < inst->ncols; i++) x[i] = 0.0;
	for (int i = 0; i <= inst->nnodes - 1; i++) {
		x[xpos(i, final_succ[i], inst)] = 1.0;
	}

	free(succ);
	free(final_succ);
	free(tabu);
}


void generate_random_solution(instance* inst, double* x) {

	int* nodes_list = (int*)malloc(inst->nnodes * sizeof(int));				// List of node indexes from which we draw randomly the nodes
	for (int i = 0; i < inst->nnodes; i++) nodes_list[i] = i;

	int max_index = inst->nnodes - 1;										// max_index is the max value of node index we want to draw at each iteration
	int start_node_index = rand() % max_index;								// Chooose a random starting node
	int start_node = nodes_list[start_node_index];

	int temp_node = nodes_list[max_index];									// Every time we draw an index from nodes_list, we swap the node value with the last node
	nodes_list[max_index] = start_node;										// of the list and we lower max_index, so that the same node cannot be drawn multiple times
	nodes_list[start_node_index] = temp_node;

	max_index--;															

	//for (int i = 0; i < inst->nnodes; i++) printf("nodes_list[%d]: %d\n", i, nodes_list[i]);

	int last_node = start_node;
	int curr_edge_index = -1;
	int next_node_index = -1;
	int next_node = -1;
	double curr_sol_cost = 0.0;
	int n_edges = 0;
	while (max_index >= 0) {												// Repeat until all edges apart from the "closing loop" one are found

		if (max_index != 0) next_node_index = rand() % max_index;			// Generate a random node index between 0 and max_index. If max_index is 0, just take node at index 0
		else next_node_index = 0;
		next_node = nodes_list[next_node_index];

		curr_edge_index = xpos(last_node, next_node, inst);
		curr_sol_cost += dist(last_node, next_node, inst);					// Add up the cost of the edge just found
		x[curr_edge_index] = 1.0;											// Set the edge just found as part of the solution

		//if (inst->verbose >= HIGH) printf("Edge #%d : [ %d -> %d ]\n", n_edges, last_node, next_node);
		
		temp_node = nodes_list[max_index];									// Every time we draw an index from nodes_list, we swap the node value with the last node
		nodes_list[max_index] = next_node;									// of the list and we lower max_index, so that the same node cannot be drawn multiple times
		nodes_list[next_node_index] = temp_node;

		last_node = next_node;												// Move to the next node
		max_index--;
		n_edges++;
	}												

	int last_edge_index = xpos(last_node, start_node, inst);				// The node last_node is the last picked from the list
	curr_sol_cost += dist(last_node, start_node, inst);						// Add up the cost of the last "closing loop" edge
	x[last_edge_index] = 1.0;

	//if (inst->verbose >= HIGH) printf("Edge #%d : [ %d -> %d ]\n", n_edges, last_node, start_node);

	inst->z_best = curr_sol_cost;											// Memorize the cost of the generated solution

	if (inst->verbose >= MEDIUM) printf("Cost of the random solution (before 2-opt refinement): %f\n", curr_sol_cost);

	return;
}


void _generate_feasible_nodes_list(instance* inst, int* nodes_list) {

	int* all_nodes_list = (int*)malloc(inst->nnodes * sizeof(int));			// List of node indexes from which we draw randomly the nodes
	for (int i = 0; i < inst->nnodes; i++) all_nodes_list[i] = i;

	int max_index = inst->nnodes - 1;										// max_index is the max value of node index we want to draw at each iteration
	int start_node_index = rand() % max_index;								// Chooose a random starting node
	int start_node = all_nodes_list[start_node_index];

	int last_node_index = 0;
	nodes_list[last_node_index] = start_node;								// Insert the first node into the final solution nodes list
	//printf("nodes_list[%d]: %d\n", last_node_index, nodes_list[last_node_index]);
	last_node_index++;

	int temp_node = all_nodes_list[max_index];								// Every time we draw an index from nodes_list, we swap the node value with the last node
	all_nodes_list[max_index] = start_node;									// of the list and we lower max_index, so that the same node cannot be drawn multiple times
	all_nodes_list[start_node_index] = temp_node;

	max_index--;

	int last_node = start_node;
	int next_node_index = -1;
	int next_node = -1;
	double curr_sol_cost = 0.0;
	while (max_index >= 0) {												// Repeat until all edges apart from the "closing loop" one are found

		if (max_index != 0) next_node_index = rand() % max_index;			// Generate a random node index between 0 and max_index. If max_index is 0, just take node at index 0
		else next_node_index = 0;
		next_node = all_nodes_list[next_node_index];

		nodes_list[last_node_index] = next_node;							// Insert the new node into the final solution nodes list
		//printf("nodes_list[%d]: %d\n", last_node_index, nodes_list[last_node_index]);
		last_node_index++;

		curr_sol_cost += dist(last_node, next_node, inst);					// Add up the cost of the edge just found

		temp_node = all_nodes_list[max_index];								// Every time we draw an index from nodes_list, we swap the node value with the last node
		all_nodes_list[max_index] = next_node;								// of the list and we lower max_index, so that the same node cannot be drawn multiple times
		all_nodes_list[next_node_index] = temp_node;

		last_node = next_node;												// Move to the next node
		max_index--;
	}

	curr_sol_cost += dist(last_node, start_node, inst);						// Add up the cost of the last "closing loop" edge

	inst->z_best = curr_sol_cost;											// Memorize the cost of the generated solution

	return;
}


void get_solution_from_nodes_list(instance* inst, int* nodes_list, double* x) {

	// Clean the array of edges where the solution has to be returned
	for (int i = 0; i < inst->ncols; i++) x[i] = 0.0;
	
	int start_node = nodes_list[0];
	int last_node = start_node;
	int next_node = -1;
	for (int j = 1; j < inst->nnodes; j++) {
		next_node = nodes_list[j];
		x[xpos(last_node, next_node, inst)] = 1.0;
		last_node = next_node;
	}
	x[xpos(last_node, start_node, inst)] = 1.0;

	return;
}


void get_succ_from_nodes_list(instance* inst, int* nodes_list, int* succ, int* start_node) {

	*start_node = nodes_list[0];

	for (int j = 0; j < inst->nnodes - 1; j++) {
		succ[nodes_list[j]] = nodes_list[j + 1];
	}
	succ[nodes_list[inst->nnodes - 1]] = nodes_list[0];

	return;
}


void get_nodes_list_from_succ(instance* inst, int* nodes_list, int* succ, int start_node) {

	nodes_list[0] = start_node;

	for (int j = 1; j < inst->nnodes; j++) {
		nodes_list[j] = succ[nodes_list[j - 1]];
	}

	return;
}

void solve_heur_genetic(instance* inst, double* x, int pop_size, double ratio_2_opt) {
	 
	double residual_timelimit = inst->timelimit;
	double min_sol_cost = INFINITY;
	double fitness_spread_at_first_epoch = 0.0;
	double t1, t2 = 0.0;
	int champion_index = -1;

	// Generate the starting population of pop_size (ex. 1000) random solutions (already with 2-opt refinement) and keep their costs (called "fitness")
	int** population = (int**)malloc(pop_size * sizeof(int*));
	if (population == NULL) print_error("Unable to allocate memory for pop_size solutions in HEUR_GENETIC.");
	double* fitness = (double*)malloc(pop_size * sizeof(double));

	t1 = second();

	int start_node = -1;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));								// Shared successors list data structure
	for (int i = 0; i < pop_size; i++) {
		int* nodes_list = (int*)malloc(inst->nnodes * sizeof(int));						// Allocate memory for a single list of nodes

		_generate_feasible_nodes_list(inst, nodes_list);								// Generate it randomly

		// Only if ratio_2_opt parameter is > 0, then apply 2-opt refinement to this fraction of starting solutions
		if (i < (pop_size * ratio_2_opt)) {
			get_succ_from_nodes_list(inst, nodes_list, succ, &start_node);				// Get the successors list from the nodes list and the start node
			solve_heur_2_opt(inst, NULL, succ, 0.5);									// Refine it with 2-opt for a small period of time
			get_nodes_list_from_succ(inst, nodes_list, succ, start_node);				// Get back the nodes list given the successors list and the start node
		}
		
		population[i] = nodes_list;														// Add it to population list
		fitness[i] = inst->z_best;														// Memorize also the solution cost
		if (inst->verbose >= HIGH) printf("Random solution #%d generated with fitness: %f\n", i, fitness[i]);
	}

	t2 = second();

	// Lower the residual timelimit, but wait at least for the first epoch to be concluded to check if timelimit has been reached
	residual_timelimit = residual_timelimit - (t2 - t1);

	// Compute the worst fitness (highest cost) of the first epoch. Then it will be updated after the offspring generation
	double worst_fitness = -INFINITY;
	for (int i = 0; i < pop_size; i++) {
		if (fitness[i] > worst_fitness) worst_fitness = fitness[i];
	}
	if (inst->verbose >= HIGH) printf("Worst fitness of the starting population: %f\n\n", worst_fitness);

	int n_epoch = 1;
	while (1) {

		t1 = second();

		if (inst->verbose >= HIGH) printf("\nEpoch #%d\n", n_epoch);

		// Generate the new (pop_size / 10) offsprings merging the chromosomes and kill among the worst (pop_size / 10) solutions
		int offspring_size = pop_size / 10;
		for (int n = 0; n < offspring_size; n++) {

			// Choose the pair of population members to merge together in a probabilistic way, so that the solutions with best fitness are mostly chosen
			double norm_fitness = 1.0;
			double sol_thresh = 1.0;
			int first_parent_index = -1;
			while (1) {
				first_parent_index = rand() % pop_size;											// Draw a random solution
				norm_fitness = fitness[first_parent_index] / worst_fitness;						// Normalize its fitness w.r.t. the worst fitness value of the current epoch
				sol_thresh = 1.0 - norm_fitness * norm_fitness;									// Compute a threshold that favours solutions with fitness far from the worst one

				if (((double)rand() / RAND_MAX) <= sol_thresh + EPSILON) break;
			}
			if (inst->verbose >= HIGH) printf("first_parent_index: %d with fitness: %f\n", first_parent_index, fitness[first_parent_index]);

			int second_parent_index = -1;
			while (1) {
				second_parent_index = rand() % pop_size;										// Draw another random solution, making sure that it is different from the first parent
				while (second_parent_index == first_parent_index) {
					second_parent_index = rand() % pop_size;
				}

				norm_fitness = fitness[second_parent_index] / worst_fitness;					// Normalize its fitness w.r.t. the worst fitness value of the current epoch
				sol_thresh = 1.0 - norm_fitness * norm_fitness;									// Compute a threshold that favours solutions with fitness far from the worst one

				if (((double)rand() / RAND_MAX) <= sol_thresh + EPSILON) break;
			}
			if (inst->verbose >= HIGH) printf("second_parent_index: %d with fitness: %f\n", second_parent_index, fitness[second_parent_index]);

			// Choose the population member that will die in this epoch, which will be replaced by the new solution (the offspring).
			// N.B. We assume that the Champion of the last epoch cannot die (because it's the most likely to survive the environment).
			int offspring_index = -1;
			while (1) {
				offspring_index = rand() % pop_size;											// Draw another random solution, making sure that it is different
				while (	offspring_index == first_parent_index ||								// from both parents and from the Champion of the last epoch
						offspring_index == second_parent_index || 
						offspring_index == champion_index)
				{
					offspring_index = rand() % pop_size;
				}

				norm_fitness = fitness[offspring_index] / worst_fitness;						// Normalize its fitness w.r.t. the worst fitness value of the current epoch
				sol_thresh = norm_fitness * norm_fitness;										// Compute a threshold that favours solutions with fitness close to the worst one

				if (((double)rand() / RAND_MAX) < sol_thresh + EPSILON) break;
			}
			if (inst->verbose >= HIGH) printf("Index of the killed member: %d with fitness: %f\n", offspring_index, fitness[offspring_index]);

			// N.B. Notice that the killed member will (rarely) be a previous offspring of the same epoch, but this is how evolution works (simulates infant mortality)

			// Get the chromosomes (= list of ordered nodes corresponding to each solution) of these two population members
			// N.B. The chromosomes are basically: population[first_parent_index] and population[second_parent_index]. 
			if (inst->verbose >= HIGH) {
				printf("1st parent nodes list: ");
				for (int i = 0; i < inst->nnodes; i++) printf("%d ", population[first_parent_index][i]);
				printf("\n2nd parent nodes list: ");
				for (int i = 0; i < inst->nnodes; i++) printf("%d ", population[second_parent_index][i]);
			}
			
			// Merge the chromosomes of the two parents such that the offspring solution is feasible:
			int chromosome_cutting_index = rand() % (inst->nnodes / 2 - 1) + 1;								// Cut the chromosomes at index between 1 and inst->nnodes / 2 - 1
			// 1) Copy the first half from the first parent chromosome
			for (int i = 0; i < chromosome_cutting_index; i++) {
				population[offspring_index][i] = population[first_parent_index][i];
			}	
			// 2) Copy the second half from the second parent chromosome such that in the offspring solution no node is repeated and all are included
			int next_index_to_replace = chromosome_cutting_index;
			int n_repeated_nodes = 0;
			int is_node_repeated = 0;
			for (int i = chromosome_cutting_index; i < inst->nnodes; i++) {								// Visit the second parent nodes list only until its end
				for (int j = 0; j < chromosome_cutting_index; j++) {
					if (population[second_parent_index][i] == population[first_parent_index][j]) {
						is_node_repeated = 1;
						break;
					}
				}
				// If we arrive here without having spotted a repeated node => we can add the candidate node to the new solution
				if (!is_node_repeated) {
					next_index_to_replace = i - n_repeated_nodes;
					population[offspring_index][next_index_to_replace] = population[second_parent_index][i];
				}
				else {
					n_repeated_nodes++;
				}
				is_node_repeated = 0;
			}
			// 3) If at least one node of the second half of the second parent is repeated => start visiting the second parent chromosome from the beginning
			int i = 0;
			while (i < n_repeated_nodes) {																// Keep on visiting the second parent nodes list until no node is repeated
				for (int j = 0; j < chromosome_cutting_index; j++) {
					if (population[second_parent_index][i] == population[first_parent_index][j]) {
						is_node_repeated = 1;
						break;
					}
				}
				// If we arrive here without having spotted a repeated node => we can add the candidate node to the new solution
				if (!is_node_repeated) {
					next_index_to_replace++;
					population[offspring_index][next_index_to_replace] = population[second_parent_index][i];
				}
				else {
					n_repeated_nodes++;
				}
				is_node_repeated = 0;
				i++;
			}
			// N.B. The previous chromosomes merging procedure works only under the assumption that both parents' nodes lists are of fesible solutions

			// Compute the new solution fitness
			double new_sol_cost = 0.0;
			int first_node = population[offspring_index][0];
			int prev_node = first_node;
			int next_node = -1;
			for (int i = 1; i < inst->nnodes; i++) {
				next_node = population[offspring_index][i];
				new_sol_cost += dist(prev_node, next_node, inst);
				prev_node = next_node;
			}
			new_sol_cost += dist(next_node, first_node, inst);

			// Make sure that, if this option is enabled (ratio_2_opt > 0), then 2-opt refinement is applied to about ratio_2_opt of the offspring of the current epoch
			// => "ratio_2_opt" parameter controls both the ratio of starting and new (offspring) solutions to which 2-opt refinement is applied
			if ((ratio_2_opt > 0) && ((double)rand() / RAND_MAX) < ratio_2_opt) {
				inst->z_best = new_sol_cost;
				get_succ_from_nodes_list(inst, population[offspring_index], succ, &start_node);				// Get the successors list from the nodes list and the start node
				solve_heur_2_opt(inst, NULL, succ, 0.5);													// Refine it with 2-opt for a small period of time
				get_nodes_list_from_succ(inst, population[offspring_index], succ, start_node);				// Get back the nodes list given the successors list and the start node
				new_sol_cost = inst->z_best;
			}

			// Finally update the fitness data structure
			fitness[offspring_index] = new_sol_cost;

			if (inst->verbose >= HIGH) {
				printf("\n Offspring nodes list: ");
				for (int i = 0; i < inst->nnodes; i++) printf("%d ", population[offspring_index][i]);
				printf("\n\n");
			}
		}

		// Compute the average fitness of the current epoch
		double fitness_sum = 0.0;
		for (int i = 0; i < pop_size; i++) fitness_sum += fitness[i];
		double avg_fitness = fitness_sum / (double)pop_size;
		if (inst->verbose >= MEDIUM) printf("Average fitness of epoch #%d: %f\n", n_epoch, avg_fitness);
		
		// Find the solution with best (= smallest) fitness value (called "Champion") of the current epoch
		double best_fitness = INFINITY;
		for (int i = 0; i < pop_size; i++) {
			if (fitness[i] < best_fitness) {
				best_fitness = fitness[i];
				inst->z_best = best_fitness;
				champion_index = i;
			}
		}
		if (inst->verbose >= MEDIUM) printf("Champion fitness of epoch #%d: %f\n", n_epoch, best_fitness);

		// Also update the worst fitness after generating the offspring
		worst_fitness = -INFINITY;
		for (int i = 0; i < pop_size; i++) {
			if (fitness[i] > worst_fitness) worst_fitness = fitness[i];
		}
		if (inst->verbose >= HIGH) printf("Worst fitness of epoch #%d (after offspring generation): %f\n", n_epoch, worst_fitness);

		t2 = second();

		// Update the residual timelimit and check if the global timelimit has been reached. If so, return the champion solution.
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (residual_timelimit <= 0) {
			if (inst->verbose >= MEDIUM) printf("Timelimit reached! The Champion of the last epoch will be returned.\n");
			get_solution_from_nodes_list(inst, population[champion_index], x);
			break;
		}

		// If the champion (best) fitness turns out to be equal to the worst one, it means that all population members have become equal
		// => stop the algorithm, since it can't find any better solution and will loop until timelimit is reached.
		// N.B. This should never happen!
		if (best_fitness >= worst_fitness - EPSILON && best_fitness <= worst_fitness + EPSILON) {
			if (inst->verbose >= MEDIUM) printf("Population has become completely homogeneous! No further improvements are expected!\n");
			get_solution_from_nodes_list(inst, population[champion_index], x);
			break;
		}

		t1 = second();

		// Before starting the new epoch, let's introduce some random mutations, where a mutation consistes in swapping 2 nodes of a solution.
		// N.B. We want to apply a bunch of mutations when the spread between the worst and best solution is too low compared to that of the first epoch (ex. < 10%)
		double fitness_spread = (worst_fitness - best_fitness);
		if (n_epoch == 1) fitness_spread_at_first_epoch = fitness_spread;					// Memorize the fitness spread of the first epoch as a reference to keep "close" to
		double relative_fitness_spread = fitness_spread / fitness_spread_at_first_epoch;
		if (inst->verbose >= MEDIUM) printf("Relative worst-best fitness spread: %f\n\n", relative_fitness_spread);

		if (relative_fitness_spread < 0.1) {

			int n_mutations = rand() % (inst->nnodes) + 1;									// Select a number of mutations between 1 and inst->nnodes
			for (int i = 0; i < n_mutations; i++) {

				int mutant_member_index = rand() % pop_size;								// Choose a random solution to mutate, avoiding the current "Champion" one
				while (mutant_member_index == champion_index) {
					mutant_member_index = rand() % pop_size;
				}

				int first_node_to_swap_index = rand() % inst->nnodes;						// Choose two random (different) nodes of the selected solution to swap
				int second_node_to_swap_index = rand() % inst->nnodes;
				while (second_node_to_swap_index == first_node_to_swap_index) {
					second_node_to_swap_index = rand() % inst->nnodes;
				}

				int temp_node = population[mutant_member_index][first_node_to_swap_index];
				population[mutant_member_index][first_node_to_swap_index] = population[mutant_member_index][second_node_to_swap_index];
				population[mutant_member_index][second_node_to_swap_index] = temp_node;

				// Of course, now we need to re-compute the fitness value of the mutant solution
				double mutant_sol_cost = 0.0;
				int start_node = population[mutant_member_index][0];
				int last_node = start_node;
				int next_node = -1;
				for (int j = 1; j < inst->nnodes; j++) {
					next_node = population[mutant_member_index][j];
					mutant_sol_cost += dist(last_node, next_node, inst);
					last_node = next_node;
				}
				mutant_sol_cost += dist(last_node, start_node, inst);
				fitness[mutant_member_index] = mutant_sol_cost;
			}
			if (inst->verbose >= MEDIUM) printf("\nNumber of random mutations performed: %d\n\n", n_mutations);
		}

		t2 = second();

		// Update the residual timelimit and check if the global timelimit has been reached. If so, return the champion solution.
		residual_timelimit = residual_timelimit - (t2 - t1);
		if (residual_timelimit <= 0) {
			if (inst->verbose >= MEDIUM) printf("Timelimit reached! The Champion of the last epoch will be returned.\n");
			get_solution_from_nodes_list(inst, population[champion_index], x);
			break;
		}

		n_epoch++;
	}
	
	// Free all the memory allocated for the population solutions
	for (int i = 0; i < pop_size; i++) free(population[i]);
	free(population);
	free(fitness);
	free(succ);

	return;
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
			double obj = dist(i, j, inst);						// cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error("Wrong CPXnewcols on x variables");
			// Verify if xpos returns the right position inside the tableu (xpos starts from 0)
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst)) print_error("Wrong position for x variables");
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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error("Wrong CPXnewcols on x variables");

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
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) print_error("Wrong CPXnewcols on x variables");

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

	// Model MTZ with lazy consistency constraints.

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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error("Wrong CPXnewcols on x variables");

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
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) print_error("Wrong CPXnewcols on x variables");

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
	rhs = M - 1;												// "rhs" = right-hand side of the degree constraints
	sense = 'L';												// 'L' for less-or-equal constraint
	nnz = 3;													// number of cells != 0 in the row
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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error("Wrong CPXnewcols on x variables");

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
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) print_error("Wrong CPXnewcols on x variables");

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

	// Model MTZ with both lazy consistency constraints and lazy subtour elimination constraints of dimension 2.

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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error("Wrong CPXnewcols on x variables");

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
		if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) print_error("Wrong CPXnewcols on x variables");

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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error("Wrong CPXnewcols on x variables");

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
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &integer, cname)) print_error("Wrong CPXnewcols on y variables");

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
	rhs = 0.0;																		// "rhs" = right-hand side of the degree constraints
	sense = 'L';																	// 'L' for less-or-equal constraint
	nnz = 2;																		// number of cells != 0 in the row
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = 0; j < inst->nnodes; j++) {
			if (i == j) continue;

			sprintf(cname[0], "linking_constr(%d,%d)", i + 1, j + 1);
			index[0] = xpos_compact(i, j, inst);									// +(1.0 - N) * Xij
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

