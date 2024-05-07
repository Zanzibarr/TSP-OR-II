#ifndef _TSP_H
#define _TSP_H

#include "utils.h"
#include "threads.h"

#pragma region PRECOMPUTING

/**
 * @brief Precomputes the sort_edges list
 */
void tsp_precompute_sort_edges();

/**
 * @brief Precomputes the costs list
 */
void tsp_precompute_costs();

#pragma endregion


#pragma region ALGORITHMS TOOLS

/**
 * @brief Returns and index corresponding to the algorithm chosen
*/
int tsp_find_alg();

/**
 * @brief Comparator used by the qsort method to sort the sort_edges list
 * 
 * @param arg1 The first element to compare (casted as tsp_entry)
 * @param arg2 The second element to compare (casted as tsp_entry)
 * 
 * @return -1 if arg1.value < arg2.value, 0 if arg1.value == arg2.value, 1 if arg1.value > arg2.value
*/
int tsp_compare_entries(const void* arg1, const void* arg2);

/**
 * @brief (THREAD SAFE) Checks and updates the incumbent of the instance
 * 
 * @param path Candidate solution for the update
 * @param cost Cost of the candidate solution
 * @param time Time at which the candidate solution was found
 */
void tsp_check_best_sol(const int *path, const double cost, const double time);

/**
 * @brief Reverse a list from start to end
 *
 * @param path The list to reverse
 * @param start Beginning index of the reverse (included)
 * @param end Ending index of the reverse (included)
 */
void tsp_reverse(int *path, int start, int end);

/**
 * @brief (THREAD SPECIFIC) Checks if an edge is in the tabu list
 *
 * @param t_index The index of the thread (used to determine which tabu list to look into)
 * @param from Node 1 (order doesn't matter)
 * @param to Node 2 (order doesn't matter)
 *
 * @return 1 if the move is a tabu in the specified table, 0 otherwise
 */
int tsp_check_tabu(const int t_index, const int from, const int to);

/**
 * @brief (THREAD SPECIFIC) Adds an edge to the tabu list
 *
 * @param t_index The index of the thread (used to determine which tabu list to look into)
 * @param from Node 1 (order doesn't matter)
 * @param to Node 2 (order doesn't matter)
 */
void tsp_add_tabu(const int t_index, const int from, const int to);

/**
 * @brief Computes the cost of a path
 * 
 * @param path The path whose cost is calculated
 * 
 * @return The cost of the path
*/
double tsp_compute_path_cost(const int* path);

/**
 * @brief return cost of edge (i,j)
 * 
 * @param i first edge node
 * @param j second edge node
 * @return double cost of edge (i,j) as stored in tsp_inst.costs
 */
double tsp_get_edge_cost(const int i, const int j);

/**
 * @brief convert a succ type solution to a path type solution
 * 
 * @param succ the starting list with the succ type solution
 * @param path list where to store the converted solution
 * 
 * @return The cost of the path calculated
*/
double tsp_succ_to_path(const int* succ, int* path);

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

/**
 * @brief Applies a swap policy till the solution cannot be improved further
 * 
 * @param path Path considered (will be changed if it finds a swap)
 * @param cost Cost of the current path (will be changed if it finds a swap)
 * @param swap_function The swap function to use
*/
void tsp_2opt(int* path, double* cost, const int (*swap_function)(int*, double*));

#pragma endregion


#pragma region INITIALIZATIONS

/**
 * @brief Initialize default variables (sort of a constructor for the problem)
 */
void tsp_init_defs();

/**
 * @brief Initialize the incumbent
 */
void tsp_init_solution();

#pragma endregion


#pragma region SAVING FILES

/**
 * @brief Save to file the best solution found so far
 */
int tsp_save_solution();

/**
 * @brief Plot the best solution found so far
 */
void tsp_plot_solution(const int unique);

#pragma endregion


#pragma region DEBUGGING TOOLS

/**
 * @brief Prints to stdout the instance parameters
 */
void tsp_instance_info();

/**
 * @brief Prints to stdout the best solution found so far
 */
void tsp_print_solution();

/**
 * @brief Checks the correctness of an intermediate solution
 *
 * @param path The intermediate solution
 * @param cost The cost of the intermediate solution
 * @param message The message to print in case of a violation
 */
void tsp_check_integrity(const int *path, const double cost, const char* message);

#pragma endregion


#pragma region MEMORY MANAGEMENT

/**
 * @brief Dinamically allocate the coords list
 */
void tsp_allocate_coords_space();

/**
 * @brief Dinamically allocate the space for the tabu list
*/
void tsp_allocate_tabu_space();

/**
 * @brief Frees dinamically allocated memory created by tsp_allocate methods and concludes eventual other finishing operations
 */
void tsp_free_all();

#pragma endregion

#endif