#ifndef TSP_UTILITIES_H_

#define TSP_UTILITIES_H_

#include "instance.h"

// Parses all command line arguments passed by the user and sets the instance parameters
void parse_command_line(int argc, char** argv, instance* inst);

// Opens and parses input_file and fills the instance with data
void parse_input_file(instance* inst);

// Creates and fills a "nodes.dat" file to allow for graph plotting
void print_nodes_dat_file(const instance* inst);

// Create the given number of random (but deterministic) instances with a fixed number of nodes (all parameters are inside inst)
void create_instances(instance* inst);

// Prints in a standard way the provided type of error
void print_error(const char* err);

// Computes specifically the euclidean distance between nodes i and j
double dist(int i, int j, instance * inst);

// Releases allocated memory after the istance is solved
void free_instance(instance* inst);

// Returns the processor clock time used since the beginning of the program (in seconds)
double second();

extern int mkdir(const char* dir_name);

// Inline useful functions
inline int imax(int i1, int i2) { return (i1 > i2) ? i1 : i2; }
inline double dmin(double d1, double d2) { return (d1 < d2) ? d1 : d2; }
inline double dmax(double d1, double d2) { return (d1 > d2) ? d1 : d2; }

#endif   /* TSP_UTILITIES_H_ */ 