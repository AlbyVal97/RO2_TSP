#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <cplex.h>

#include "tsp.h"

// Command line parameter list examples:
// -f ../data/inst_max_100_nodes/berlin52.tsp -seed 123456 -model_type 1 -verbose 1 -time_limit 3600
// -m 0 -n_inst 20 -n_nodes 50 -f ../data/inst_max_100_nodes/berlin52.tsp -seed 123456 -model_type 3 -verbose 2 -time_limit 60
// -f ../data/test_instances_50/rnd_50_1.tsp -seed 123456 -model_type 7 -verbose 2 -time_limit 60

// -m 1 -folder test_instances -prefix rnd -n_inst 5 -n_nodes 15
// -m 1 -folder test_instances_10 -prefix rnd_10 -n_inst 20 -n_nodes 10
// -m 1 -folder test_instances_50 -prefix rnd_50 -n_inst 20 -n_nodes 50
// -m 1 -folder test_instances_200 -prefix rnd_200 -n_inst 20 -n_nodes 200
// -m 1 -folder test_instances_300 -prefix rnd_300 -n_inst 40 -n_nodes 300
// -m 1 -folder test_instances_500 -prefix rnd_500 -n_inst 20 -n_nodes 500
// -m 1 -folder test_instances_750 -prefix rnd_750 -n_inst 40 -n_nodes 750
// -m 1 -folder test_instances_1000 -prefix rnd_1000 -n_inst 40 -n_nodes 1000
// -m 1 -folder test_instances_2000 -prefix rnd_2000 -n_inst 40 -n_nodes 2000
// -m 1 -folder test_instances_5000 -prefix rnd_5000 -n_inst 20 -n_nodes 5000


// Command line arguments list for compact models official tests:
// -m 2 -test test_compact_models_50 -folder test_instances_50 -n_inst 20 -n_models 5 1 2 3 4 5 -prefix rnd_50 -time_limit 1800

// Command line arguments list for tuning hyperparameters of advanced branch and cut method:
// -m 2 -test test_tuning_advbc -folder test_instances_300_tuning -n_inst 20 -n_models 5 8 9 10 11 12 -prefix rnd_300 -time_limit 1800

// Command line arguments for Benders + B&C + B&C_2_opt + advB&C methods official tests:
// -m 2 -test test_benders_bc_models_300 -folder test_instances_300 -n_inst 20 -n_models 4 6 7 8 10 -prefix rnd_300 -time_limit 1800
// -m 2 -test test_benders_bc_models_500 -folder test_instances_500 -n_inst 20 -n_models 4 6 7 8 10 -prefix rnd_500 -time_limit 1800

// Command line arguments list for tuning hyperparameters of hard fix heuristics (only on 10 instances because each run takes about 15 minutes):
// -m 2 -test test_tuning_hard_fix -folder test_instances_750_tuning -n_inst 10 -n_models 4 14 15 16 17 -prefix rnd_750 -time_limit 900
// -m 2 -test test_tuning_1000_hard_fix -folder test_instances_1000_tuning -n_inst 10 -n_models 4 14 15 16 17 -prefix rnd_1000 -time_limit 900

// Command line arguments list for tuning hyperparameters of soft fix heuristics (only on 10 instances because each run takes about 15 minutes):
// -m 2 -test test_tuning_soft_fix -folder test_instances_750_tuning -n_inst 10 -n_models 4 18 19 20 21 -prefix rnd_750 -time_limit 900
// -m 2 -test test_tuning_1000_soft_fix -folder test_instances_1000_tuning -n_inst 10 -n_models 4 14 15 16 17 -prefix rnd_1000 -time_limit 900

// Command line arguments lists for hard fix and soft fix (best) methods official tests (split in two parts because each run takes about 30 minutes):
// -m 2 -test test_math_heuristics_part_1_bc_2opt -folder test_instances_750_part_1 -n_inst 10 -n_models 2 17 21 -prefix rnd_750 -time_limit 1800
// -m 2 -test test_math_heuristics_part_2_bc_2opt -folder test_instances_750_part_2 -n_inst 10 -n_models 2 17 21 -prefix rnd_750 -time_limit 1800
// -m 2 -test test_math_heuristics_part_1_bc_2opt -folder test_instances_1000_part_1 -n_inst 10 -n_models 2 17 21 -prefix rnd_1000 -time_limit 1800
// -m 2 -test test_math_heuristics_part_2_bc_2opt -folder test_instances_1000_part_2 -n_inst 10 -n_models 2 17 21 -prefix rnd_1000 -time_limit 1800

// Command line arguments list for constructive heuristics official tests (each run takes about 10 minutes):
// -m 2 -test test_constructive_heuristics -folder test_instances_1000 -n_inst 20 -n_models 3 23 25 27 -prefix rnd_1000 -time_limit 600

// Command line arguments lists for metaheuristics official tests (each run takes about 30 minutes):
// -m 2 -test test_metaheuristics_part_1 -folder test_instances_2000_part_1 -n_inst 5 -n_models 4 28 29 30 31 -prefix rnd_2000 -time_limit 1800
// -m 2 -test test_metaheuristics_part_2 -folder test_instances_2000_part_2 -n_inst 5 -n_models 4 28 29 30 31 -prefix rnd_2000 -time_limit 1800
// -m 2 -test test_metaheuristics_part_3 -folder test_instances_2000_part_3 -n_inst 5 -n_models 4 28 29 30 31 -prefix rnd_2000 -time_limit 1800
// -m 2 -test test_metaheuristics_part_4 -folder test_instances_2000_part_4 -n_inst 5 -n_models 4 28 29 30 31 -prefix rnd_2000 -time_limit 1800

// Command line arguments lists for heuristics final test: VNS vs HEUR_SOFT_FIX_9 (each run takes about 30 minutes):
// -m 2 -test test_soft_fixing_2000_6_minutes_first_timelimit -folder test_instances_2000_final_test -n_inst 20 -n_models 1 21 -prefix rnd_2000 -time_limit 1800
// -m 2 -test test_soft_fixing_2000_10_minutes_first_timelimit -folder test_instances_2000_final_test -n_inst 20 -n_models 1 21 -prefix rnd_2000 -time_limit 1800
// -m 2 -test test_vns_2000 -folder test_instances_2000_final_test -n_inst 20 -n_models 1 28 -prefix rnd_2000 -time_limit 1800


void update_csvfile(instance* inst, int first_model, int last_model, double time, double z_best);
void run_test(instance* inst);

int main(int argc, char **argv) {

	if ( argc < 2 ) {
		printf("Usage: %s -help for help\n", argv[0]);
		exit(1);
	}

	instance inst;

	parse_command_line(argc, argv, &inst);

	if (inst.mode == CREATE_INSTANCES) {
		create_instances(&inst);
		return 0;
	}

	if (inst.mode == RUN_TEST) {
		run_test(&inst);
		return 0;
	}

	parse_input_file(&inst);

	double t1 = second();

	if (inst.verbose >= HIGH) printf("Timestamp 1: %lf seconds.\n", t1);

	if (TSPopt(&inst)) print_error("Error within TSPopt()");

	double t2 = second();

	if (!(inst.timelimit_exceeded)) printf("TSP problem solved successfully in %lf seconds.\n\n", t2 - t1);
	else if (inst.verbose == TEST && inst.timelimit_exceeded) printf("TSP problem ended with timelimit of %lf seconds.\n\n", inst.timelimit);

	if (inst.verbose == TEST) {	
		if (inst.timelimit_exceeded == 1) {
			update_csvfile(&inst, inst.first_model, inst.last_model, 10.0 * inst.timelimit, inst.z_best);
		}
		else {
			update_csvfile(&inst, inst.first_model, inst.last_model, t2 - t1, inst.z_best);
		}
	}

	// Plot the resulting optimal tour using Gnuplot
	if (inst.verbose >= MEDIUM) system("C:/\"Program Files\"/gnuplot/bin/gnuplot.exe ../outputs/gnuplot_commands.txt");

	free_instance(&inst);

	return 0;
}


void run_test(instance* inst) {
	printf("*** RUN %s ***\n\n", inst->testname);
	printf("Istances from folder %s\n", inst->folder_instances);
	printf("Models used:\n");
	for (int j = 0; j < inst->n_models_test; j++)
		printf("%s\n", models[inst->models_to_test[j]]);

	if (mkdir("../outputs") == -1) {
		if (inst->verbose >= MEDIUM) printf("Folder \"outputs\" already exists.\n");
	}
	else {
		if (inst->verbose >= MEDIUM) printf("Folder \"outputs\" created for the first time!\n");
	}

	if (inst->n_instances <= 0) print_error("Number of instances to be tested must be greater than zero");

	char csv_pathname[100];
	sprintf(csv_pathname, "../outputs/%s.csv", inst->testname);
	FILE* csv_file = fopen(csv_pathname, "w");
	if (csv_file == NULL) print_error("test.csv file not found inside \"outputs\" folder!");

	fprintf(csv_file, "%d, ", inst->n_models_test);
	for (int j = 0; j < inst->n_models_test - 1; j++)
		fprintf(csv_file, "%s, ", models[inst->models_to_test[j]]);
	fprintf(csv_file, "%s\n", models[inst->models_to_test[inst->n_models_test-1]]);

	fclose(csv_file);

	char bat_pathname[100];
	sprintf(bat_pathname, "../outputs/%s.bat", inst->testname);
	FILE* bat_file = fopen(bat_pathname, "w");
	if (bat_file == NULL) print_error("bat file not found inside \"outputs\" folder!");

	for (int i = 0; i < inst->n_instances; i++) {
		for (int j = 0; j < inst->n_models_test; j++) {

			fprintf(bat_file, ".\\Release\\tsp -f ../data/%s/%s_%d.tsp -test %s -model_type %d -first %d -last %d -verbose 0 -seed 123456 -time_limit %f -end\n", inst->folder_instances, inst->instance_prefix_name, i+1, inst->testname, inst->models_to_test[j], inst->models_to_test[0], inst->models_to_test[inst->n_models_test-1], inst->timelimit);
		}
	}
	fprintf(bat_file, "@echo *** TEST ENDED ***\n");
	fprintf(bat_file, "pause\n");

	fclose(bat_file);

	free(inst->models_to_test);

	sprintf(bat_pathname, "..\\outputs\\%s.bat", inst->testname);
	system(bat_pathname);

	return;
}


void update_csvfile(instance* inst, int first_model, int last_model, double time, double z_best) {
	char csv_path[120];
	if (strcmp(inst->testname, "NULL\0") == 0) print_error("Using TEST verbosity not in RUN_TEST mode: no testname is provided!");
	sprintf(csv_path, "../outputs/%s.csv", inst->testname);
	FILE* csv_file = fopen(csv_path, "a");
	if (csv_file == NULL) print_error("csv file not found inside \"outputs\" folder!");

	// If the method to test is a heuristic (inst->model_type >= HEUR_HARD_FIX_50), then the .csv file must report the final value of the objective function
	// instead of the computational time, which is fixed since they all run until timelimit is reached
	if (inst->model_type <= ADVBC_PROB_10) {
		if (inst->model_type == first_model) {						// Print the instance name only before the test on the first model
			fprintf(csv_file, "%s, %f, ", inst->inst_name, time);
		}
		else if (inst->model_type == last_model) {					// Go to next line only after the test on the last model has been executed
			fprintf(csv_file, "%f\n", time);
		}
		else {
			fprintf(csv_file, "%f, ", time);
		}
	}
	else {
		if (inst->model_type == first_model) {						// Print the instance name only before the test on the first model
			fprintf(csv_file, "%s, %f, ", inst->inst_name, z_best);
		}
		else if (inst->model_type == last_model) {					// Go to next line only after the test on the last model has been executed
			fprintf(csv_file, "%f\n", z_best);
		}
		else {
			fprintf(csv_file, "%f, ", z_best);
		}
	}

	fclose(csv_file);
	
	return;
}