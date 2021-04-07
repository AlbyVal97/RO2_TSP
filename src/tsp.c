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
	if (CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit)) { print_error("CPXsetdblparam() error in setting timelimit"); }

	// Set the seed parameter according to user input or default value.
	if (CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed)) { print_error("CPXsetdblparam() error in setting seed"); }

	char edges_file_path[100];
	//create outputs folder if needed
	if (mkdir("../outputs") == -1) printf("Folder outputs already exists\n");
	else printf("Folder outputs created for the first time\n");
	sprintf(edges_file_path, "../outputs/%s", inst->inst_name);
	if (mkdir(edges_file_path) == -1) printf("Folder for the current instance already exists\n");
	else printf("Folder for the current instance created for the first time\n");

	// Build the Cplex model according to model_type chosen
	// -> See "model_type" enum in instance.h for the list of models available
	// Then executes the actual optimization procedure
	if (inst->model_type != BENDERS) {							// Model creation + optimization for all methods other than	BENDERS	
		build_model(inst, env, lp);
		if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
	}
	else {														// Model creation + optimization specifically for BENDERS method	
		int n_comp = 0;			
		int* succ = (int*)malloc(inst->nnodes * sizeof(int));	// Array of the successors for each node
		int* comp = (int*)malloc(inst->nnodes * sizeof(int));	// Array of the connected component index for each node
		
		// For BENDERS model, build the model at least once, then start to solve
		build_model(inst, env, lp);

		if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }
			
		int ncols = CPXgetnumcols(env, lp);
		double* x = (double*)calloc(ncols, sizeof(double));
		if (CPXgetx(env, lp, x, 0, ncols - 1)) { print_error("CPXgetx() error"); }

		update_components(x, inst, succ, comp, &n_comp);
		if (inst->verbose >= MEDIUM) printf("Current n_comp: %d \n", n_comp);

		int n_iter = 0;
		while (n_comp > 1) {									// Repeat iteratively until just one connected component is left

			// Add a new subtour elimination constraint for the first connected component of the current solution
			update_benders_constraints(env, lp, inst, comp, n_iter);

			// Optimize with the new constraint
			if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }

			// Extract the new solution
			ncols = CPXgetnumcols(env, lp);
			if (CPXgetx(env, lp, x, 0, ncols - 1)) { print_error("CPXgetx() error"); }

			// Update the numbero of connected components of the new graph
			update_components(x, inst, succ, comp, &n_comp);
			if (inst->verbose >= MEDIUM) printf("Current n_comp: %d \n", n_comp);

			// Increment current iteration number
			n_iter++;
		}

		// Outputs only the final model to file "benders_model.lp"
		char model_file_path[100];
		sprintf(model_file_path, "../outputs/%s/benders_model.lp", inst->inst_name);
		if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
		CPXwriteprob(env, lp, model_file_path, NULL);

		free(succ);
		free(comp);
		free(x);

	}
	

	// use the optimal solution found by CPLEX

	// Create the edges.dat file associated to the correct model type
	//FILE* edges_plot_file_name = NULL;
	
	if (inst->verbose >= MEDIUM) {

		switch (inst->model_type) {

			case BASIC:
				sprintf(edges_file_path, "%s/basic_model_edges.dat", edges_file_path);
				break;

			case MTZ_STATIC:
				sprintf(edges_file_path, "%s/model_MTZ_static_edges.dat", edges_file_path);
				break;

			case MTZ_LAZY:
				sprintf(edges_file_path, "%s/model_MTZ_lazy_u_consistency_edges.dat", edges_file_path);
				break;

			case MTZ_SUBTOUR_SIZE_2:
				sprintf(edges_file_path, "%s/model_MTZ_lazy_2_node_SECs_edges.dat", edges_file_path);
				break;

			case GG:
				sprintf(edges_file_path, "%s/model_GG_edges.dat", edges_file_path);
				break;

			case BENDERS:
				sprintf(edges_file_path, "%s/model_BENDERS_edges.dat", edges_file_path);
				break;
		}
		printf("Complete path for *.dat file: %s\n\n", edges_file_path);

		//edges_plot_file_name = fopen(edges_file_path, "w");
		//if (edges_plot_file_name == NULL) print_error("File edges.dat not found!");
	}
	

	int ncols = CPXgetnumcols(env, lp);
	// Allocate memory for the optimal solution array
	double* xstar = (double*)calloc(ncols, sizeof(double));
	// Copy the optimal solution from the Cplex environment to the new array "xstar"
	if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { print_error("CPXgetx() error"); }

	int symmetric = -1;
	if (inst->model_type == BASIC || inst->model_type == BENDERS) symmetric = 0;
	else if (inst->model_type == MTZ_STATIC || inst->model_type == MTZ_LAZY || inst->model_type == MTZ_SUBTOUR_SIZE_2 || inst->model_type == GG) symmetric = 1;
	
	if (inst->verbose >= MEDIUM) {
		print_solution(inst, xstar, symmetric, edges_file_path);
	}

	// Scan all legal edges and print the ones involved (with x ~ 1) in the optimal tour
	/*if (inst->model_type == BASIC || inst->model_type == BENDERS) {
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = i + 1; j < inst->nnodes; j++) {

				if (xstar[xpos(i, j, inst)] > 0.5) {

					if (inst->verbose >= MEDIUM) {
						printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
						fprintf(edges_plot_file_name, "%f %f\n%f %f\n\n", inst->xcoord[i], inst->ycoord[i], inst->xcoord[j], inst->ycoord[j]);
					}
				}
			}
		}
	}
	else if (inst->model_type == MTZ_STATIC || inst->model_type == MTZ_LAZY || inst->model_type == MTZ_SUBTOUR_SIZE_2 || inst->model_type == GG) {
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = 0; j < inst->nnodes; j++) {
				if (i == j) continue;

				if (xstar[xpos_compact(i, j, inst)] > 0.5) {

					if (inst->verbose >= MEDIUM) {
						printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
						fprintf(edges_plot_file_name, "%f %f %f %f\n", inst->xcoord[i], inst->ycoord[i], inst->xcoord[j] - inst->xcoord[i], inst->ycoord[j] - inst->ycoord[i]);
					}
				}
			}
		}
	}*/

	//if (inst->verbose >= MEDIUM) {
		//fclose(edges_plot_file_name);
	//}

	// Free allocated memory and close Cplex model
	free(xstar);

	if (CPXfreeprob(env, &lp)) { print_error("CPXfreeprob() error"); }
	if (CPXcloseCPLEX(&env)) { print_error("CPXcloseCPLEX() error"); }

	return 0; // or an appropriate nonzero error code
}

void print_solution(instance* inst, double* xstar, int symmetric, char* edges_file_path) {
	if (symmetric != 0 && symmetric != 1) print_error("symmetric is a boolean variable\n");
	
	FILE* edges_plot_file_name = fopen(edges_file_path, "w");
	if (edges_plot_file_name == NULL) print_error("File edges.dat not found!");

	if (!symmetric) {
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = i + 1; j < inst->nnodes; j++) {

				if (xstar[xpos(i, j, inst)] > 0.5) {

					if (inst->verbose >= MEDIUM) {
						printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
						fprintf(edges_plot_file_name, "%f %f\n%f %f\n\n", inst->xcoord[i], inst->ycoord[i], inst->xcoord[j], inst->ycoord[j]);
					}
				}
			}
		}
	}
	else {
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = 0; j < inst->nnodes; j++) {
				if (i == j) continue;

				if (xstar[xpos_compact(i, j, inst)] > 0.5) {

					if (inst->verbose >= MEDIUM) {
						printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
						fprintf(edges_plot_file_name, "%f %f %f %f\n", inst->xcoord[i], inst->ycoord[i], inst->xcoord[j] - inst->xcoord[i], inst->ycoord[j] - inst->ycoord[i]);
					}
				}
			}
		}
	}

	FILE* gn_com = fopen("../src/gnuplot_commands.txt", "w");
	char* comand;
	if (gn_com == NULL) print_error("Erroe while opening gnuplot_commands file\n");
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


void update_benders_constraints(CPXCENVptr env, CPXLPptr lp, instance* inst, const int* comp, int n_iter) {
	// Scan first connencted component for its number of nodes -> value of first_comp_n_nodes
	int first_comp_n_nodes = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		if (comp[i] == 1) {
			first_comp_n_nodes++;
		}
	}
	
	// Scan first connected component for its nodes indexes
	int* first_comp_nodes = (int*)malloc(first_comp_n_nodes * sizeof(int));
	int j = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		if (comp[i] == 1) {
			first_comp_nodes[j++] = i;
		}
	}

	char** cname = (char**)calloc(1, sizeof(char*));			// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));
	sprintf(cname[0], "subtour_elimination_constraint(%d)", n_iter);
	double rhs = first_comp_n_nodes - 1.0;
	char sense = 'L';
	int i_zero = 0;
	// We should compute the number of edges connecting the nodes of the current connected component
	// as the binomial coefficient (n over 2) = n!/((n-2)!*(2!)), but it's not feasible for "big" values of n
	// => we just allocate the square of the number of nodes, which is not that much bigger than needed:
	int n_edges_curr_comp = first_comp_n_nodes * first_comp_n_nodes;
	if (inst->verbose >= MEDIUM) printf("n_edges_curr_comp: %d\n", n_edges_curr_comp);
	int* index = (int*)calloc(n_edges_curr_comp, sizeof(int));					// Array of indexes associated to the row variables
	double* value = (double*)calloc(n_edges_curr_comp, sizeof(double));			// Array of row variables coefficients
	int nnz = n_edges_curr_comp;

	int k = 0;
	for (int i = 0; i < first_comp_n_nodes; i++) {
		for (int j = i + 1; j < first_comp_n_nodes; j++) {
			if (i == j) {
				continue;
			}
			index[k] = xpos(first_comp_nodes[i], first_comp_nodes[j], inst);
			value[k++] = 1.0;
		}

	}
	if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for subtour elimination constraints!");

	free(cname[0]);
	free(cname);
	free(index);
	free(value);
}


void update_components(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp) {

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


void build_model(instance* inst, CPXENVptr env, CPXLPptr lp) {

	char model_file_path[100];

	switch (inst->model_type) {

		case BASIC:
		{
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
			for (int h = 0; h < inst->nnodes; h++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 2.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "degree(%d)", h + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int i = 0; i < inst->nnodes; i++) {
					if (i == h) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Outputs to file "basic_model.lp" the built model
			sprintf(model_file_path, "../outputs/%s/basic_model.lp", inst->inst_name);
			if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
			CPXwriteprob(env, lp, model_file_path, NULL);

			free(cname[0]);
			free(cname);

			break;
		}

		case MTZ_STATIC:
		{
			// Model MTZ with static subtour elimination constraints.

			double zero = 0.0;
			char binary = 'B';
			char integer = 'I';

			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) { print_error("CPXsetdblparam() error in setting integer value tolerance"); }

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

			// Add the inner degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "inner_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(j, i, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "outer_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the u consistency constraints for each edge (i,j)
			int i_zero = 0;
			int index[3];
			double value[3];
			double M = inst->nnodes - 1.0;								// Smallest M value for big M trick
			double rhs = M - 1;											// "rhs" = right-hand side of the degree constraints
			char sense = 'L';											// 'L' for less-or-equal constraint
			int nnz = 3;												// number of cells != 0 in the row
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
			sprintf(model_file_path, "../outputs/%s/model_MTZ_static.lp", inst->inst_name);
			if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
			CPXwriteprob(env, lp, model_file_path, NULL);

			free(cname[0]);
			free(cname);

			break;
		}

		case MTZ_LAZY:
		{
			// Model MTZ with lazy subtour elimination constraints.

			double zero = 0.0;
			char binary = 'B';
			char integer = 'I';

			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) { print_error("CPXsetdblparam() error in setting integer value tolerance"); }

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

			// Add the inner degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "inner_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(j, i, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "outer_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the u consistency constraints for each edge (i,j)
			int i_zero = 0;
			int index[3];
			double value[3];
			double M = inst->nnodes - 1.0;								// Smallest M value for big M trick
			double rhs = M - 1;											// "rhs" = right-hand side of the degree constraints
			char sense = 'L';											// 'L' for less-or-equal constraint
			int nnz = 3;												// number of cells != 0 in the row
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
			sprintf(model_file_path, "../outputs/%s/model_MTZ_lazy_u_consistency.lp", inst->inst_name);
			if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
			CPXwriteprob(env, lp, model_file_path, NULL);

			free(cname[0]);
			free(cname);

			break;
		}

		case MTZ_SUBTOUR_SIZE_2:
		{
			// Model MTZ with lazy subtour elimination constraints of dimension 2.

			double zero = 0.0;
			char binary = 'B';
			char integer = 'I';

			if (CPXsetdblparam(env, CPX_PARAM_EPINT, 0.0)) { print_error("CPXsetdblparam() error in setting integer value tolerance"); }

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

			// Add the inner degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "inner_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(j, i, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "outer_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the u consistency constraints for each edge (i,j)
			int i_zero = 0;
			int index[3];
			double value[3];
			double M = inst->nnodes - 1.0;								// Smallest M value for big M trick
			double rhs = M - 1;											// "rhs" = right-hand side of the degree constraints
			char sense = 'L';											// 'L' for less-or-equal constraint
			int nnz = 3;												// number of cells != 0 in the row
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

			// Outputs to file "model_MTZ_lazy_2_node_SECs.lp" the built model
			sprintf(model_file_path, "../outputs/%s/model_MTZ_lazy_2_node_SECs.lp", inst->inst_name);
			if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
			CPXwriteprob(env, lp, model_file_path, NULL);

			free(cname[0]);
			free(cname);

			break;
		}

		case GG:
		{
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

			// Add the inner degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "inner_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(j, i, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the outer degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int i = 0; i < inst->nnodes; i++) {					// degree constraints
				int lastrow = CPXgetnumrows(env, lp);
				double rhs = 1.0;										// "rhs" = right-hand side of the degree constraints
				char sense = 'E';										// 'E' for equality constraint 
				sprintf(cname[0], "outer_degree(%d)", i + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int j = 0; j < inst->nnodes; j++) {
					if (j == i) continue;
					if (CPXchgcoef(env, lp, lastrow, xpos_compact(i, j, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Add the linking constraints for the variable y
			// Add the u consistency constraints for each edge (i,j)
			int i_zero = 0;
			int index[2];
			double value[2];
			double rhs = 0.0;											// "rhs" = right-hand side of the degree constraints
			char sense = 'L';											// 'L' for less-or-equal constraint
			int nnz = 2;												// number of cells != 0 in the row
			for (int i = 0; i < inst->nnodes; i++) {
				for (int j = 0; j < inst->nnodes; j++) {
					if (i == j) continue;

					sprintf(cname[0], "linking_constr(%d,%d)", i + 1, j + 1);
					index[0] = xpos_compact(i, j, inst);										// +(2.0 - N) * Xij
					value[0] = 1.0 - inst->nnodes;
					index[1] = ypos_compact(i, j, inst);										// +1.0 * Yij
					value[1] = 1.0;
					if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for linking_constr");
				}
			}

			// Add the starting flow constraint (the starting flow value from node 1 is n - 1)
			int lastrow = CPXgetnumrows(env, lp);
			rhs = inst->nnodes - 1.0;												// "rhs" = right-hand side of the degree constraints
			sense = 'E';															// 'E' for equality constraint 
			sprintf(cname[0], "starting_flow_constr");								// Set a name for the new row/constraint
			if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

			// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
			// Set coefficient 1 to any variable associated to edges linked to node h
			for (int j = 1; j < inst->nnodes; j++) {
				if (CPXchgcoef(env, lp, lastrow, ypos_compact(0, j, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
			}

			// Add the flow balance constraints
			// Add the degree constraints for all nodes (because in TSP model the final tour is hamiltonian)
			for (int h = 1; h < inst->nnodes; h++) {							// degree constraints
				lastrow = CPXgetnumrows(env, lp);
				rhs = 1.0;														// "rhs" = right-hand side of the degree constraints
				sense = 'E';													// 'E' for equality constraint 
				sprintf(cname[0], "flow_balance(%d)", h + 1);					// Set a name for the new row/constraint
				if (CPXnewrows(env, lp, 1, &rhs, &sense, NULL, cname)) { print_error(" wrong CPXnewrows [degree]"); }

				// TODO: put all coefficients in an array, not creating empty constraints and then changing coefficients one by one
				// Set coefficient 1 to any variable associated to edges linked to node h
				for (int i = 0; i < inst->nnodes; i++) {
					if (i == h) continue;
					if (CPXchgcoef(env, lp, lastrow, ypos_compact(i, h, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
					if (CPXchgcoef(env, lp, lastrow, ypos_compact(h, i, inst), -1.0)) print_error(" wrong CPXchgcoef [degree]");
				}
			}

			// Outputs to file "model_GG.lp" the built model
			sprintf(model_file_path, "../outputs/%s/model_GG.lp", inst->inst_name);
			if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
			CPXwriteprob(env, lp, model_file_path, NULL);

			free(cname[0]);
			free(cname);

			break;
		}

		case BENDERS:
		{
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

				for (int i = 0; i < inst->nnodes; i++) {
					
					if (i != h) {
						index[i] = xpos(i, h, inst);
						value[i] = 1.0;
					}
					else {
						continue;
					}
				}
				if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &i_zero, index, value, NULL, cname)) print_error("wrong CPXaddrows() for degree constraints!");
			}

			// Outputs to file "basic_model.lp" the built model
			sprintf(model_file_path, "../outputs/%s/benders_model.lp", inst->inst_name);
			if (inst->verbose >= MEDIUM) printf("\nComplete path for *.lp file: %s\n", model_file_path);
			CPXwriteprob(env, lp, model_file_path, NULL);

			free(cname[0]);
			free(cname);
			free(index);
			free(value);

			

			break;
		}

	}

	

	return;
}
