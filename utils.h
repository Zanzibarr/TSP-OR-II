#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

typedef struct {
    int key;
    double value;
} tsp_entry;

typedef struct {    //node expressed as coordinates
    double x, y;
} tsp_pair;

// PARSING
#define TSP_FILE_P "-file"
#define TSP_TIME_LIMIT "-tl"
#define TSP_SEED "-seed"
#define TSP_NNODES "-nodes"
#define TSP_HELP "-help"
#define TSP_ALGORITHM "-alg"
#define TSP_QUIET "-quiet"
#define TSP_VERBOSE "-verbose"

// TIME MANAGEMENT
#define TSP_DEF_TL 3.6e+6  //number of ms in a day
clock_t tsp_initial_time;
double tsp_total_time;

// INSTANCE INFO
char tsp_file_name[1000];
#define TSP_DEF_NNODES 300  //default number of nodes
#define TSP_GRID_SIZE 10000
uint64_t tsp_seed;
time_t tsp_time_limit;
char tsp_edge_weight_type[10];
char tsp_alg_type[10];  //store the type of algorithm using

// NUMBERS
#define TSP_EPSYLON 1e-9    //to round double values

// LOGGING
int tsp_verbose; // -1 for quiet, 0 for normal, 1 for verbose

void tsp_init_rand() { for (int i = 0; i < 100; i++) rand(); }  //fixing first random values being small
double tsp_rnd_coord() { return (double)rand()/RAND_MAX*TSP_GRID_SIZE; }  //generate a random number between 0 and GRID_SIZE
double tsp_time_elapsed() { return ((double)clock() - tsp_initial_time)/(CLOCKS_PER_SEC); }

#endif