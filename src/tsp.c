#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>

#include "tsp.h"


int TSPopt(instance* inst) {

	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	build_model(inst, env, lp);

	// Cplex's parameter setting
	// ...

	if (CPXmipopt(env, lp)) { print_error("CPXmipopt() error"); }

	// use the optimal solution found by CPLEX

	int ncols = CPXgetnumcols(env, lp);
	// Allocate memory for the optimal solution array
	double* xstar = (double*)calloc(ncols, sizeof(double));
	// Copy the optimal solution from the Cplex environment to the new array "xstar"
	if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { print_error("CPXgetx() error"); }
	// Scan all legal edges and print the ones involved (with x ~ 1) in the optimal tour
	for (int i = 0; i < inst->nnodes; i++) {
		for (int j = i + 1; j < inst->nnodes; j++) {
			if (xstar[xpos(i, j, inst)] > 0.5) printf("x(%3d,%3d) = 1\n", i + 1, j + 1);
		}
	}
	free(xstar);

	// free and close cplex model   
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


void build_model(instance* inst, CPXENVptr env, CPXLPptr lp) {

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
		// Set coefficient 1 to any variable associated to edges linked to node h
		for (int i = 0; i < inst->nnodes; i++) {
			if (i == h) continue;
			if (CPXchgcoef(env, lp, lastrow, xpos(i, h, inst), 1.0)) print_error(" wrong CPXchgcoef [degree]");
		}
	}

	// Outputs to file "model.lp" the built model
	if (inst->verbose >= MEDIUM) { CPXwriteprob(env, lp, "model.lp", NULL); }

	free(cname[0]);
	free(cname);

	return;
}