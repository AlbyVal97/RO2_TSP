#ifndef TSP_H_  

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <cplex.h> 

#define VERBOSE				    1000		// printing level  (debug_mode if VERBOSE >= 1000)

//hard-wired parameters
#define XSMALL		  		  1e-5 		// 1e-4*	// tolerance used to decide ingerality of 0-1 var.s
#define EPSILON		  		  1e-9		// 1e-9		// very small numerical tolerance 
#define TICKS_PER_SECOND 	  1000.0  	// cplex's ticks on Intel Core i7 quadcore @2.3GHZ

//data structures  

typedef struct {

    //input data
	int nnodes;   
	double *xcoord;
	double *ycoord;
	int integer_costs;                   // = 1 for integer costs (rounded distances), 0 otherwise

	// parameters 
	int model_type; 
	double timelimit;						// overall time limit, in sec.s
	char input_file[1000];		  			// input file
							
	double zbest;							// value of the best sol. available

} instance;

//inline
inline int imax(int i1, int i2) { return ( i1 > i2 ) ? i1 : i2; } 
inline double dmin(double d1, double d2) { return ( d1 < d2 ) ? d1 : d2; } 
inline double dmax(double d1, double d2) { return ( d1 > d2 ) ? d1 : d2; }

// Prints in a standard way the provided type of error
void print_error(const char* err);

// Computes specifically the euclidean distance between nodes i and j
double dist(int i, int j, instance* inst);

// Optimizes the provided instance of the TSP
int TSPopt(instance* inst);

// Computes the position of edge [i,j] inside the Cplex tableau
int xpos(int i, int j, instance* inst);

// Builds the tableau (variables and constraints)
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);

#endif   /* TSP_H_ */ 