#ifndef _ALG_H
#define _ALG_H

#include "tsp.h"

// TABU
int     tsp_solve_tabu(tsp_instance* inst);
void    tsp_find_2opt_best_swap_tabu(tsp_instance* inst, int* path, double* cost);

// GREEDY
int     tsp_solve_greedy(tsp_instance* inst, const char g2opt);
int     tsp_solve_greedy_mt(tsp_instance* inst, const char g2opt);

double  tsp_greedy_from_node(const tsp_instance* inst, int* path, int start_node);
void*   tsp_greedy_from_node_mt(void* params);

// 2OPT
void    tsp_2opt(const tsp_instance* inst, int* path, double* cost, int (*swap_function)(const tsp_instance*, int*, double*));
int  tsp_find_2opt_swap(const tsp_instance* inst, int* path, double* cost);
int  tsp_find_2opt_best_swap(const tsp_instance* inst, int* path, double* cost);

#endif