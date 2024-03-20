#ifndef _TSP_H
#define _TSP_H

#include "utils.h"
#include "threads.h"

// PRECOMPUTING

/**
 * Precomputes the sort_edges list
 */
void tsp_precompute_sort_edges();

/**
 * Precomputes the costs list
 */
void tsp_precompute_costs();

// ALGORITHMS TOOLS

/**
 * (THREAD SAFE) Checks and updates the incumbent of the instance
 * @param path Candidate solution for the update
 * @param cost Cost of the candidate solution
 * @param time Time at which the candidate solution was found
 */
void tsp_check_best_sol(const int *path, const double cost, const double time);

/**
 * Reverse a list from start to end
 *
 * @param path The list to reverse
 * @param start Beginning index of the reverse (included)
 * @param end Ending index of the reverse (included)
 */
void tsp_reverse(int *path, int start, int end);

/**
 * (THREAD SPECIFIC) Checks if an edge is in the tabu list
 *
 * @param t_index The index of the thread (used to determine which tabu list to look into)
 * @param from Node 1 (order doesn't matter)
 * @param to Node 2 (order doesn't matter)
 *
 * @returns 1 if the move is a tabu in the specified table, 0 otherwise
 */
int tsp_check_tabu(const int t_index, const int from, const int to);

/**
 * (THREAD SPECIFIC) Adds an edge to the tabu list
 *
 * @param t_index The index of the thread (used to determine which tabu list to look into)
 * @param from Node 1 (order doesn't matter)
 * @param to Node 2 (order doesn't matter)
 */
void tsp_add_tabu(const int t_index, const int from, const int to);

// INITIALIZATIONS

/**
 * Initialize default variables (sort of a constructor for the problem)
 */
void tsp_init_defs();

/**
 * Initialize the incumbent
 */
void tsp_init_solution();

// SAVING FILES

/**
 * Prints to stdout the best solution found so far
 */
void tsp_print_solution();

/**
 * Save to file the best solution found so far
 */
void tsp_save_solution();

/**
 * Plot the best solution found so far
 */
void tsp_plot_solution();

// DEBUGGING TOOLS

/**
 * Prints to stdout the instance parameters
 */
void tsp_instance_info();

/**
 * Checks the correctness of an intermediate solution
 *
 * @param path The intermediate solution
 * @param cost The cost of the intermediate solution
 */
void tsp_check_integrity(const int *path, const double cost);

// MEMORY MANAGEMENT

/**
 * Dinamically allocate the coords list
 */
void tsp_allocate_coords_space();

/**
 * Frees dinamically allocated memory created by tsp_allocate methods and concludes eventual other finishing operations
 */
void tsp_free_instance();

#endif