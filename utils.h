#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <ctype.h>

// PARSING
#define FILE_P "-file"
#define TIME_LIMIT "-tl"
#define SEED "-seed"
#define NNODES "-nodes"
#define HELP "-help"
#define QUIET "-quiet"
#define VERBOSE "-verbose"

// TIME MANAGEMENT
#define MS_SEC 1e+3 //number of ms in s
#define DEF_TL 8.64e+7  //number of ms in a day

// INSTANCE INFO
char file_name[1000];
#define DEF_NNODES 300
#define GRID_SIZE 10000
uint64_t seed;
uint64_t time_limit;
char edge_weight_type[10];

// NUMBERS
#define EPSYLON 1e-9 // to round double values

// LOGGING
int verbose; // -1 for quiet, 0 for normal, 1 for verbose

double rnd() { return (double)rand()/RAND_MAX*GRID_SIZE; }  //generate a random number between 0 and GRID_SIZE

#endif