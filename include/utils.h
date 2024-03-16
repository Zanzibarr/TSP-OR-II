#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

// DEBUGGING
#define TSP_VERBOSE 10
/**
 * <0 for quiet                                 (nothing)
 * [0, 10[ for normal                           (basic info for final user)
 * == 5 for thread info                         (multithreading)
 * >=10 partial solutions                       (basic info for debugging)
 * >=100 for integrity checks                   (integrity checks enabled)      <--     suggested while in development
 * >=500 to see the path in the solution        (advanced debugging)
 * >=1000 for super-verbose                     (full verbose)
*/

// MULTITHREADING
#define N_THREADS 12

// PARSING
#define TSP_FILE_P "-file"
#define TSP_TIME_LIMIT "-tl"
#define TSP_SEED "-seed"
#define TSP_NNODES "-nodes"
#define TSP_HELP "-help"
#define TSP_ALGORITHM "-alg"
#define TSP_MT "-mt"

// DEFAULTS
#define TSP_DEF_TL 3.6e+6  //number of ms in an hour
#define TSP_DEF_NNODES 300  //default number of nodes
#define TSP_GRID_SIZE 10000
#define TSP_EDGE_W_TYPE "ATT"

// FILE NAMES
#define TSP_SOL_FOLDER "solutions"
#define TSP_INST_FOLDER "instances"
#define TSP_PLOT_FILE "solution_plot.png"
#define TSP_SOLUTION_FILE "solution_file.txt"
#define TSP_COORDS_FILE "coords_file.txt"
#define TSP_COMMAND_FILE "command_file.txt"

// USEFUL NUMBERS
#define TSP_TABU_TENUE 100
#define TSP_EPSILON 1e-8    //to round double values

extern pthread_t tsp_threads[N_THREADS];
extern int tsp_available_threads[N_THREADS];

extern int tsp_mt_choice;
extern int tsp_stoplight_update_sol;

typedef struct {
    int counter;
    int size;
    int* table;
} tabu;

typedef struct {    //temp struct used to sort nodes by cost
    int key;
    double value;
} tsp_entry;

typedef struct {    //node expressed as coordinates
    double x, y;
} tsp_pair;

typedef struct {    //instance

    int nnodes;             //number of nodes
    tsp_pair* coords;       //list of nodes

    double* costs;         //cost "matrix"
    int* sort_edges;       //min edges "matrix": row i contains a permutation of the nodes, ordered by increasing distance from node i

    int* best_solution; //store the best solution found so far
    double best_cost;   //store the best cost found so far
    double best_time;   //store the time of the best solution found so far (in seconds)

} tsp_instance;

typedef struct {    //struct used to pass parameters to functions in threads

    int t_index;

    tsp_instance* inst;
    int s_node;
    char alg;
    int* path;
    double* cost;
    int  (*swap_function)(const tsp_instance*, int*, double*);

} tsp_mt_parameters;

// TIME MANAGEMENT
extern double tsp_initial_time;
extern double tsp_total_time;
extern int tsp_over_time;
extern double tsp_time_limit;

// GLOBAL VARIABLES
extern uint64_t tsp_seed;
extern tabu tsp_tabu_table[N_THREADS];
extern char tsp_alg_type[20];
extern char tsp_file_name[100];

// USEFUL METHODS
void tsp_init_rand();
double tsp_rnd_coord();
double tsp_time_elapsed();

#endif