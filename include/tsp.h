#ifndef TSP_H_  

#define TSP_H_

#include "tsp_utilities.h"

// Optimizes the provided instance of the TSP
int TSPopt(instance* inst);

// Computes the position of edge [i,j] inside the Cplex tableau
int xpos(int i, int j, instance* inst);

// Computes the position of edge (i,j) inside the Cplex tableau
int xpos_compact(int i, int j, instance* inst);

// Computes the position of variable u associated to each node inside the Cplex tableau
int upos_compact(int i, instance* inst);

// Computes the position of variable y associated to edge (i, j) inside the Cplex tableau
int ypos_compact(int i, int j, instance* inst);

// Returns if the results are from an optimal solution or not (es. timelimit has been reached)
int mip_solved_to_optimality(instance* inst, CPXENVptr env, CPXLPptr lp);

// Add a new subtour elimination constraint for the first connected component of the current solution
void update_benders_constraints(CPXCENVptr env, CPXLPptr lp, instance* inst, const int* comp, int n_iter);

// Finds the number of connected components in the graph defined inside xstar
void update_connected_components(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);

// Apply the loop method (called Benders) to solve the tsp problem without a compact model
void solve_benders(instance* inst, CPXENVptr env, CPXLPptr lp);

// Apply the Branch & Cut method to solve the tsp problem without a compact model
void solve_branch_cut(instance* inst, CPXENVptr env, CPXLPptr lp);

// Callback function to be used internally by Cplex only in BRANCH_CUT method
static int CPXPUBLIC branch_cut_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle);

//void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);

// Builds the tableau (variables and constraints)
void build_model_BASIC(instance* inst, CPXENVptr env, CPXLPptr lp);
void build_model_MTZ_STATIC(instance* inst, CPXENVptr env, CPXLPptr lp);
void build_model_MTZ_LAZY(instance* inst, CPXENVptr env, CPXLPptr lp);
void build_model_MTZ_SEC2_STATIC(instance* inst, CPXENVptr env, CPXLPptr lp);
void build_model_MTZ_SEC2_LAZY(instance* inst, CPXENVptr env, CPXLPptr lp);
void build_model_GG(instance* inst, CPXENVptr env, CPXLPptr lp);

// Create and fill the .lp file containing the complete model definition used by Cplex
void create_lp_file(instance* inst, CPXENVptr env, CPXLPptr lp, const char* model_name);

// Prints the .dat file with the list of nodes from the optimized solution, with formatting depending on graph symmetry (directed/undirected)
void print_solution(instance* inst, double* xstar, int symmetric, char* edges_plot_file_name);

#endif   /* TSP_H_ */ 