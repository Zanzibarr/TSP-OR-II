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


// GREEDY

/**
 * @brief (MULTITHREAD) Execute the greedy algorithm on multiple starting nodes at the same time
 */
void tsp_solve_greedy();


// G2OPT

/**
 * @brief (MULTITHREAD) Execute the g2opt (+ greedy) algorithm on multiple starting nodes at the same time
*/
void tsp_solve_g2opt();


// TABU

/**
 * @brief (MULTITHREAD) Execute the tabu algorithm on multiple instances at the same time
 */
void tsp_solve_tabu();


// VNS

/**
 * @brief (MULTITHREAD) Execute the vns algorithm with parallel kicks
 */
void tsp_solve_vns();


// CPLEX

/**
 * @brief (MULTITHREAD) Execute the cplex algorithm
*/
void tsp_solve_cplex();

#endif