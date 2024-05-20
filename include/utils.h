#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <ctype.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>
#include <signal.h>
#include <cplex.h>

/**
 * @brief Number of threads
 */
#define N_THREADS 16


// PARSING CLI ARGUMENTS

#define TSP_PARSING_FILE        "-file"                 //parsing cli argument to select the file
#define TSP_PARSING_TIME_LIMIT  "-tl"                   //parsing cli argument to set the time limit
#define TSP_PARSING_SEED        "-seed"                 //parsing cli argument to set the seed
#define TSP_PARSING_NNODES      "-nodes"                //parsing cli argument to set the number of nodes
#define TSP_PARSING_ALGORITHM   "-alg"                  //parsing cli argument to set the algorithm to use
#define TSP_PARSING_VERBOSE     "-verbose"              //parsing cli argument to set verbose parameter
#define TSP_PARSING_NOPLOT      "-noplot"               //parsing cli argument to avoid creating the plot for the solution

#define TSP_PARSING_GREEDY      "greedy"                //parsing cli argument to use the greedy algorithm
#define TSP_PARSING_G2OPT       "g2opt"                 //parsing cli argument to use the g2opt algorithm
#define TSP_PARSING_TABU        "tabu"                  //parsing cli argument to use the tabu algorithm
#define TSP_PARSING_VNS         "vns"                   //parsing cli argument to use the vns algorithm
#define TSP_PARSING_CPLEX       "cplex"                 //parsing cli argument to use the cplex algorithm

#define TSP_PARSING_BEST_SWAP               "-bs"               //parsing cli argument to set first swap as swapping policy in g2opt
#define TSP_PARSING_F2OPT                   "-f2opt"            //parsing cli argument to use the f2opt algorithm
#define TSP_PARSING_TENURE                  "-tenure"           //parsing cli argument to set the tenure to use
#define TSP_PARSING_TENURE_A                "-tenure-a"         //parsing cli argument to set the amplitude parameter for the dinamic tenure
#define TSP_PARSING_TENURE_F                "-tenure-f"         //parsing cli argument to set the frequency parameter for the dinamic tenure
#define TSP_PARSING_FVNS                    "-fvns"             //parsing cli argument to use the fast vns algorithm
#define TSP_PARSING_MIPSTART                "-mipstart"         //parsing cli argument to use a mipstart in cplex
#define TSP_PARSING_CPLEX_BENDERS           "-benders"          //parsing cli argument to use benders loop with cplex
#define TSP_PARSING_CPLEX_PATCHING          "-patching"         //parsing cli argument to use patching with cplex
#define TSP_PARSING_CPLEX_CANDIDATE         "-cb-comps"         //parsing cli argument to use the candidate callback in cplex
#define TSP_PARSING_RELAX_CALLBACK          "-cb-fract"         //parsing cli argument to use the relaxation callback in cplex
#define TSP_PARSING_CPLEX_CB_PATCHING       "-cb-patching"      //parsing cli argument to use greedy patching with cplex
#define TSP_PARSING_CPLEX_HARD_FIXING       "-hard-fixing"      //parsing cli argument to use hard fixing matheuristic
#define TSP_PARSING_CPLEX_LOCAL_BRANCHING   "-local-branching"  //parsing cli argument to use local branching matheuristic

#define TSP_PARSING_TMP_CHOICE              "-tmp"              //parsing cli argument to use the temp choice


// FILE NAMES

#define TSP_SOL_FOLDER          "solutions"             // path to the solutions folder
#define TSP_INST_FOLDER         "instances"             // path to the instances folder
#define TSP_SOLUTION_FILE       "solution.txt"          // suffix for the solutions files
#define TSP_CPLEX_LP_FOLDER     "cplex_outputs/lp"      // folder for cplex lp files
#define TSP_CPLEX_LOG_FOLDER    "cplex_outputs/logs"    // folder for cplex logs


// USEFUL NUMBERS

#define TSP_EPSILON                 1e-7        // to round double values
#define TSP_CPLEX_ZERO_THRESHOLD    0.5         // threshold used by exact algorithms to determine 0/1 values
#define TSP_GRID_SIZE               10000       // grid size
#define TSP_EDGE_W_TYPE             "EUC_2D"    // default edge weight type

#define TSP_DEF_TL                  3.6e+6      // number of ms in an hour
#define TSP_DEF_NNODES              300         // default number of nodes
#define TSP_DEF_VERBOSE             100         // default verbose value

#define TSP_F2OPT_MAX_DEPTH         6           // maximum depth for the f2opt algorithm
#define TSP_DEF_TABU_TENURE         80          // tenure base size
#define TSP_CBREL_PERCENTAGE        20          // default value to use the relaxation callback (1/20 = 5%)


#pragma region USEFUL METHODS

/**
 * @brief Safely deallocates the pointer memory
 * 
 * @param ptr pointer to the memory to deallocate
*/
void safe_free(void* ptr);

/**
 * @brief Initialize the rand() function to avoid having small random numbers at the beginning
*/
void init_rand();

/**
 * @brief Get a random value in the [0, TSP_GRID_SIZE] interval
 * 
 * @return The random value in the [0, TSP_GRID_SIZE] interval
*/
double rnd_coord();

/**
 * @brief Get the time elapsed since the beginning of the execution of the code
 * 
 * @return The time (in seconds) elapsed from the beginning of the execution of the code
*/
double time_elapsed();

/**
 * @brief Prints the info passed as parameters
*/
void print_info(const char* str, ...);

/**
 * @brief Prints the warning passed as parameters
*/
void print_warn(const char* str, ...);

/**
 * @brief Prints the error passed as parameters and terminates the execution of the code
*/
void raise_error(const char* str, ...);

#pragma endregion

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


#pragma region CONVERSIONS

/**
 * @brief Converts a pair of coordinates to the edge index (in cplex notation)
 * 
 * @param i first node coordinate
 * @param j second node coordinate
 * 
 * @return edge's index
*/
int tsp_convert_coord_to_xpos(const int i, const int j);

/**
 * @brief Converts a list from a permutation to indexes and values (cplex format)
 * 
 * @param ncols number of columns the list has
 * @param path permutation type list
 * @param ind indexes type list
 * @param val values type list
*/
void tsp_convert_path_to_indval(const int ncols, const int* path, int* ind, double* val);

/**
 * @brief Converts a specific component to indexes and values (cplex format)
 * 
 * @param kcomp the index of the component
 * @param ncomps the number of components
 * @param ncols number of columns the list has
 * @param comp comp type list
 * @param ind indexes type list
 * @param val values type list
 * @param nnz number of non-zeros
 * @param rhs right-hand side
*/
void tsp_convert_comp_to_cutindval(const int kcomp, const int ncomps, const int ncols, const int* comp, int* ind, double* val, int* nnz, double* rhs);

/**
 * @brief Converts a solution from succ to indexes and values
 * 
 * @param succ succ type list
 * @param ncols number of columns
 * @param ind indexes type list
 * @param val values type list
*/
void tsp_convert_succ_to_solindval(const int* succ, const int ncols, int* ind, double* val);

/**
 * @brief Converts a list from xstar to comp, ncomp and succ
 * 
 * @param xstar xstar type list
 * @param comp comp type list
 * @param ncomp number of components the list has
 * @param succ succ type list
*/
void tsp_convert_xstar_to_compsucc(const double* xstar, int* comp, int* ncomp, int* succ);

/**
 * @brief Converts a list from succ to a permutation
 * 
 * @param succ succ type list
 * @param ncomp number of components the list has
 * @param path permutation type list
*/
void tsp_convert_succ_to_path(const int* succ, const double ncomp, int* path);

/**
 * @brief Converts a list from a permutation so succ
 * 
 * @param path permutation type list
 * @param succ succ type list
*/
void tsp_convert_path_to_succ(const int* path, int* succ);

/**
 * @brief Converts a list from xstar to nxstar and builds elist
 * 
 * @param xstar xstar type list
 * @param nnodes number of nodes the problem has
 * @param elist elist type list (for concorde)
 * @param nxstar new xstar type list
 * @param nedges number of edges that are non-zero
*/
void tsp_convert_xstar_to_elistnxstar(const double* xstar, const int nnodes, int* elist, double* nxstar, int* nedges);

/**
 * @brief Converts a list from xstar to indexes and values (cplex format)
 * 
 * @param xstar xstar type list
 * @param ncols number of columns the list has
 * @param ind indexes type list
 * @param val values type list
*/
//FIXME: Never used
void tsp_convert_xstar_to_indval(const double* xstar, const int ncols, int* ind, double* val);

/**
 * @brief Converts a list from cut_index (concorde) to index and value (cplex)
 * 
 * @param cut_index cut_index type list
 * @param cut_nnodes number of nodes in the cut
 * @param index indexes type list
 * @param value values type list
 * @param nnz number of non zero edges
*/
void tsp_convert_cutindex_to_indval(const int* cut_index, const int cut_nnodes, int* index, double* value, int* nnz);

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
 * @param path Candidate permutation solution for the update (pass NULL if you wish to use succ and ncomp)
 * @param succ Candidate successors solution for the update (pass NULL with ncomp if you wish to use path)
 * @param ncomp Number of components to specify if you use succ (pass NULL with succ if you wish to use path)
 * @param cost Cost of the candidate solution (pass NULL if you don't have it at hand)
 * @param time Time at which the candidate solution was found
 */
void tsp_check_best_sol(const int* path, const int* succ, const int* ncomp, const double* cost, const double time);

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
 * @brief Computes the cost of the cplex solution
 * 
 * @param xstar xstar cplex's solution
 * 
 * @return the computed cost
 */
double tsp_compute_xstar_cost(const double* xstar);

/**
 * @brief Computes the cost of a path (succ type list)
 * 
 * @param path The succ list whose cost is calculated
 * 
 * @return The cost of the path
*/
double tsp_compute_succ_cost(const int* succ);

/**
 * @brief return cost of edge (i,j)
 * 
 * @param i first edge node
 * @param j second edge node
 * 
 * @return double cost of edge (i,j) as stored in tsp_inst.costs
 */
double tsp_get_edge_cost(const int i, const int j);

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
void tsp_2opt(int* path, double* cost, int (*swap_function)(int*, double*));

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
void tsp_save_solution();

/**
 * @brief Plot the best solution found so far
 * 
 * @param filename the solution file to plot
 */
void tsp_plot_solution();

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
