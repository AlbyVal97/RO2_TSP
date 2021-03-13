#ifndef TSP_H_  

#define TSP_H_

#include "tsp_utilities.h"

// Optimizes the provided instance of the TSP
int TSPopt(instance* inst);

// Computes the position of edge [i,j] inside the Cplex tableau
int xpos(int i, int j, instance* inst);

// Builds the tableau (variables and constraints)
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);

#endif   /* TSP_H_ */ 