#ifndef _ALG_H
#define _ALG_H

#include "tsp.h"

// GREEDY
int     tsp_solve_greedy(tsp_instance* inst, const char g2opt);
int     tsp_solve_greedy_mt(tsp_instance* inst, const char g2opt);

double  tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node);
void*   tsp_greedy_from_node_mt(void* params);

// 2OPT
void    tsp_2opt(const tsp_instance* inst, int* path, double* cost, double (*swap_function)(const tsp_instance*, int*, double*));
double  tsp_find_2opt_swap(const tsp_instance* inst, int* path, double* cost);
double  tsp_find_2opt_best_swap(const tsp_instance* inst, int* path, double* cost);

#endif