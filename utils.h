#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
//#include <gnuplot_c.h>

// PARSING
#define FILE "-file"
#define TIME_LIMIT "-tl"
#define SEED "-seed"
#define NNODES "-nodes"
#define HELP "-help"

// TIME MANAGEMENT
#define MS_SEC 1e+3
#define DEF_TL 8.64e+7

// GRID
#define DEF_NNODES 300
#define GRID_SIZE 10000

// NUMBERS
#define EPSYLON 1e-9

uint64_t time_limit;
uint64_t seed;

double rnd() { return (double)rand()/RAND_MAX*GRID_SIZE; }

#endif