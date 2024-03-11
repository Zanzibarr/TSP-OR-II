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

// DEBUGGING
#define TSP_VERBOSE 0 // <0 for quiet, 0 for normal, >0 for verbose

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

// TIME MANAGEMENT
clock_t tsp_initial_time;

// USEFUL NUMBERS
#define TSP_EPSYLON 1e-9    //to round double values


// USEFUL METHODS

void tsp_init_rand() { for (int i = 0; i < 100; i++) rand(); }  //fixing first random values being small

double tsp_rnd_coord() { return (double)rand()/RAND_MAX*TSP_GRID_SIZE; }  //generate a random number between 0 and GRID_SIZE

double tsp_time_elapsed() { return ((double)clock() - tsp_initial_time)/(CLOCKS_PER_SEC); }

#endif