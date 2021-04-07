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

// Add a new subtour elimination constraint for the first connected component of the current solution
void update_benders_constraints(CPXCENVptr env, CPXLPptr lp, instance* inst, const int* comp, int n_iter);

// Finds the number of connected components in the graph defined inside xstar
void update_components(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);

// Builds the tableau (variables and constraints)
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);

// Prints the .dat file with the list of nodes from the optimized solution, with formatting depending on graph symmetry (directed/undirected)
void print_solution(instance* inst, double* xstar, int symmetric, char* edges_plot_file_name);

extern int mkdir(const char* dir_name);

#endif   /* TSP_H_ */ 