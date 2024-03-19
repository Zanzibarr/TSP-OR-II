#ifndef _ALG_H
#define _ALG_H

#include "tsp.h"
#include "threads.h"

// 2OPT
int tsp_find_2opt_swap(const tsp_instance* inst, int* path, double* cost);
int tsp_find_2opt_best_swap(const tsp_instance* inst, int* path, double* cost);

// GREEDY / G2OPT
int     tsp_solve_greedy(tsp_instance* inst, int (*swap_function)(const tsp_instance*, int*, double*));
int     tsp_solve_greedy_st(tsp_instance* inst, int (*swap_function)(const tsp_instance*, int*, double*));

// TABU
int     tsp_solve_tabu(tsp_instance* inst);

#endif