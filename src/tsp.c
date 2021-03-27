#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>
#include <gnuplot_c.h>

#include "tsp.h"


int TSPopt(instance* inst) {

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	// Set the time_limit parameter according to user input or default value.
	if (CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit)) { print_error("CPXsetdblparam() error"); }

	build_model(inst, env, lp);

	// Cplex's parameter setting
	// ...

	// Executes the actual optimization procedure
	if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }

	// use the optimal solution found by CPLEX

	// Create the edges.dat file associated to the correct model type
	FILE* edges_plot_file_name = NULL;
	if (inst->verbose >= MEDIUM) {

		char edges_file_path[100];
		sprintf(edges_file_path, "../outputs/%s", inst->inst_name);

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
		}
		printf("Complete path for *.dat file: %s\n\n", edges_file_path);

		edges_plot_file_name = fopen(edges_file_path, "w");
		if (edges_plot_file_name == NULL) print_error("File edges.dat not found!");
	}
	

	int ncols = CPXgetnumcols(env, lp);
	// Allocate memory for the optimal solution array
	double* xstar = (double*)calloc(ncols, sizeof(double));
	// Copy the optimal solution from the Cplex environment to the new array "xstar"
	if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { print_error("CPXgetx() error"); }
	// Scan all legal edges and print the ones involved (with x ~ 1) in the optimal tour
	if (inst->model_type == BASIC) {
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
	}

	if (inst->verbose >= MEDIUM) {
		fclose(edges_plot_file_name);
	}

	// Free allocated memory and close Cplex model
	free(xstar);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
	
	return 0; // or an appropriate nonzero error code
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

	}

	

	return;
}
