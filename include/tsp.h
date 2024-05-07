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
 * @return -1 if arg1<arg2, 0 if arg1 == arg2, 1 if arg1 > arg2
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
double tsp_get_edge_cost(int i, int j);

/**
 * @brief convert a succ type solution to a path type solution
 * 
 * @param succ the starting list with the succ type solution
 * @param path list where to store the converted solution
 * 
 * @return The cost of the path calculated
*/
double tsp_succ_to_path(const int* succ, int* path);

// CPLEX

/**
 * @brief Initialize the cplex model and lp
 * 
 * @param env cplex environment
 * @param lp cplex lp
 * @param error cplex error code
*/
void tsp_cplex_init(CPXENVptr* env, CPXLPptr* lp, int* error);

/**
 * @brief Computes the cost of the cplex solution
 */
void tsp_cplex_compute_xstar_cost(double* xstar, double* cost);

/**
 * @brief add a SEC to the cplex model
 * 
 * @param ncomp number of connected components in support graph for current solution
 * @param comp vector storing the connected component for each node of support graph for current solution
 */
void tsp_cplex_add_sec(CPXENVptr env, CPXLPptr lp, int* ncomp, int* comp, int* succ);

void tsp_cplex_patching(int* ncomp, int* comp, int* succ);

/**
 * @brief Store the solution found by cplex inside the instance ONLY IF this solution is the best found so far.
 * MEANT ONLY FOR SOLUTIONS WITHOUT CYCLES (not handled if it's not)
 * 
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
*/
void tsp_cplex_check_best_sol(const int ncomp, const int* comp, const int* succ);

/**
 * @brief Convert a path type solution to a cplex type solution
 * 
 * @param ncols Number of columns (for cplex)
 * @param path path type solution to convert
 * @param indexes indexes type solution (for cplex)
 * @param values values type solution (for cplex)
*/
void tsp_cplex_path_to_xstar(const int ncols, const int* path, int* indexes, double* values);

/**
 * @brief Decompose xstar into comp and succ
 * 
 * @param xstar cplex type solution
 * @param comp list containing the component index of each node
 * @param succ list containing successor of each node
 * @param ncomp number of connected components found
*/
void tsp_cplex_decompose_xstar(const double* xstar, int* comp, int* succ, int* ncomp);

/**
 * @brief cplex callback for candidate solution
 * 
 * @param context cplex context
 * @param ncols number of columns cplex uses
 * 
 * @return cplex error code
*/
int tsp_cplex_callback_candidate(CPXCALLBACKCONTEXTptr context, const int ncols);

/**
 * @brief cplex callback for relaxation solution
 * 
 * @param context cplex context
 * @param ncols number of columns cplex uses
 * 
 * @return cplex error code
*/
int tsp_cplex_callback_relaxation(CPXCALLBACKCONTEXTptr context, const int ncols);

/**
 * @brief Closes the env and lp and frees any intermediate variables
 * 
 * @param env cplex env
 * @param lp cplex lp
 * @param comp list containing the component index of each node
 * @param succ list containing successor of each node
*/
void tsp_cplex_close(CPXENVptr env, CPXLPptr lp, int* comp, int* succ);
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
 * @brief Prints to stdout the best solution found so far
 */
void tsp_print_solution();

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
 * @brief Checks the correctness of an intermediate solution
 *
 * @param path The intermediate solution
 * @param cost The cost of the intermediate solution
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