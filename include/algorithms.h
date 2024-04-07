#ifndef _ALG_H
#define _ALG_H

#include "tsp.h"
#include "threads.h"

// 2OPT

/**
 * @brief Looks for a swap improving the current cost for the current path.
 *
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 *
 * @return 1 : found and applied a good swap, -1 : didn't found a swap (path and cost unchanged)
 */
int tsp_find_2opt_swap(int *path, double *cost);

/**
 * @brief Looks for a the swap improving the most the current cost for the current path
 * 
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 *
 * @return 1 : found and applied the best swap, -1 : didn't found a swap (path and cost unchanged)
 */
int tsp_find_2opt_best_swap(int *path, double *cost);

// GREEDY / G2OPT

/**
 * @brief (MULTITHREAD) Finds the best greedy + (eventually) 2opt optimization among the starting nodes
 * 
 * @param swap_function Specify the swap function to use
 *
 * @return 1 if the algorithm ended before the time limit, -1 otherwise
 */
int tsp_solve_greedy(int (*swap_function)(int *, double *));

/**
 * @brief (UNUSED) Finds the best greedy + (eventually) 2opt optimization among the starting nodes
 * 
 * @param swap_function Specify the swap function to use
 *
 * @return 1 if the algorithm ended before the time limit, -1 otherwise
 */
int tsp_solve_greedy_st(int (*swap_function)(int *, double *));

// TABU

/**
 * @brief (MULTITHREAD) Execute the tabu algorithm on multiple instances at the same time
 *
 * @return -1 (Reached the time limit)
 */
int tsp_solve_tabu();

// VNS

/**
 * @brief (MULTITHREAD) Execute the vns algorithm on multiple instances at the same time
 * 
 * @return int -1 (Reached the time limit)
 */
int tsp_solve_vns();

/**
 * @brief (MULTITHREAD) Execute the f2opt algorithm
 * 
 * @return int -1 (Reached the time limit)
*/
int tsp_solve_fvns();


// CPLEX

/**
 * @brief Execute the cplex algorithm
 * 
 * @return int 1 if the algorithm ended before the time limit, -1 otherwise
*/
int tsp_cplex_solve();

/**
 * @brief get cplex to solve the currently stored model and store solution and its cost
 */
void tsp_cplex_solve_model();

void tsp_cplex_benders_loop();

#endif