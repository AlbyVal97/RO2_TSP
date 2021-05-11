#ifndef TSP_UTILITIES_H_

#define TSP_UTILITIES_H_

#include "instance.h"

/**
 * Parses all command line arguments passed by the user and sets the instance fields
 * @param argc Number of command line arguments passed by the user
 * @param argv Array of command line arguments passed by the user
 * @param inst Instance struct whose fields are set by this function
 */
void parse_command_line(int argc, char** argv, instance* inst);

/**
 * Opens and parses inst->input_file and fills the instance with data retrieved from it
 * @param inst Instance struct whose fields are set by this function
 */
void parse_input_file(instance* inst);

/**
 * Creates and fills a "nodes.dat" file that could be useful for graph plotting
 * @param inst Instance struct from which nodes coordinates are retrieved
 */
void print_nodes_dat_file(const instance* inst);

/**
 * Creates the given number of random (but deterministic) instances with a fixed number of nodes
 * @param inst Instance struct from which number of instances to create and number of nodes per instance are retrieved
 */
void create_instances(instance* inst);

/**
 * Prints in a standard way the provided type of error and immediately ends the program
 * @param err Custom string explaining what happened or what caused the error
 */
void print_error(const char* err);

/**
 * Computes specifically the Euclidean distance between nodes i and j
 * If inst->integer_costs is set, then this distance is rounded to the nearest integer value
 * @param i Index of the first node
 * @param j Index of the second node
 * @param inst Instance struct from which nodes coordinates are retrieved
 * @return The Euclidean distance between nodes i and j
 */
double dist(int i, int j, instance* inst);

/**
 * Releases allocated memory after the instance is solved
 * @param inst Instance struct whose dynamically allocated memory is released
 */
void free_instance(instance* inst);

/**
 * @return The processor clock time used since the beginning of the program (in seconds)
 */
double second();

/**
 * Creates a directory with at the given path
 * @param dir_name The path, including the name, of the directory to create
 * @return 0 if the directory is successfully created, -1 if an error occurs
 */
extern int mkdir(const char* dir_name);

#endif   /* TSP_UTILITIES_H_ */ 