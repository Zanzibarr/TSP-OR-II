#ifndef _ALG_H
#define _ALG_H

#include "tsp.h"
#include "exact.h"
#include "threads.h"


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