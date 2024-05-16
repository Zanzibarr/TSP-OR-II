#ifndef _TSP_H
#define _TSP_H

#include "utils.h"
#include "threads.h"


/**
 * @brief Tabu entry used to save the edge (in pair with the tsp_tabu.list)
 */
typedef struct {
    int node_1, counter_1;
    int node_2, counter_2;
} tsp_tabu_entry;

/**
 * @brief Tabu list
 */
typedef struct {
    int counter;
    tsp_tabu_entry *list;
} tsp_tabu;

/**
 * @brief Entry (int key, double value)
 */
typedef struct {
    int key;
    double value;
} tsp_entry;

/**
 * @brief Coordinates of the node
 */
typedef struct {
    double x, y;
} tsp_pair;

/**
 * @brief Problem instance
 */
typedef struct {

    int         nnodes;             // number of nodes
    tsp_pair*   coords;             // list of nodes

    double*     costs;              // cost "matrix"
    int*        sort_edges;         // min edges "matrix":
                                    // row i contains a permutation of the nodes, ordered by increasing distance from node i

    int*        solution_succ;      // best solution found so far (successors type list)
    int         ncomp;              // number of connected components in support graph

    double      best_cost;          // cost of the best solution found so far
    double      best_time;          // time of the best solution found so far (in seconds)

} tsp_instance;

/**
 * @brief Problem environment
*/
typedef struct {

    /**
     * @brief problem's status code
     *  0 : OK
     *  1 : time limit, found solution
     *  2 : time limit, didn't find a solution
     *  3 : infeasible
     *  4 : terminated by user, found solution
     *  5 : terminated by user, didn't find a solution
     *  6 : cplex didn't even start
    */
    int         status;

    /**
     * @brief problem's effort level (was previously called "verbose")
     * <0 for quiet                                         (prints nothing)
     * [0, 10[ for normal                                   (basic info for final user)
     * >=1 to view best sol updates                         (basic info for final user)
     * ==5 for thread info                                  (multithreading info)
     * >=10 for cplex choices info                          (debugging)
     * >=100 for integrity checks                           (integrity checks enabled)      <--     suggested while in development
     * >=200 to see who finds new solutions                 (advanced debugging)
     * >=500 to see the path in the solution                (advanced debugging)
     * >=1000 for super-verbose                             (full verbose)
    */
    int         effort_level;

    char        file_name[100];             // name of the file where to read the instance (if not random)
    uint64_t    seed;                       // seed used for random algorithms
    char        alg_type[20];               // name of the algorithm using
    double      time_limit;                 // time limit
    int         cplex_terminate;            // variable used to manually terminate cplex

    char        solution_file[500];

    double      time_start;                 // initial time           
    double      time_total;                 // total execution time

    tsp_tabu    tabu_tables[N_THREADS];     // list of tabu tables needed to solve the tabu algorithm

    int         noplot;
    int         tmp_choice;                 //variable used for temporary implementation choices

    int         g2opt_swap_pol;             // swap policy for the g2opt algorithm
    int         g2opt_f2opt;                // choice for using f2opt algorithm
    int         tabu_tenure;                // tenure for the tabu algorithm
    int         tabu_tenure_a;              // tenure variability for the tabu algorithm
    double      tabu_tenure_f;              // tenure frequency for the tabu algorithm
    int         vns_fvns;                   // choice for fast/normal vns algorithm
    int         cplex_mipstart;             // choice for mipstart in cplex
    int         cplex_benders;              // choice for benders loop in cplex
    int         cplex_patching;             // choice for patching the solution in cplex
    int         cplex_can_cb;               // choice for using the candidate callback in cplex
    int         cplex_rel_cb;               // choice for using the relaxation callback in cplex
    int         cplex_cb_patching;          // choice for using patching in the callback functions

} tsp_environment;

/**
 * @brief Problem statistics
*/
typedef struct {

    int         n_solutions_found;              // store the number of feasible solutions found
    int         n_solutions_improv;             // store the number of feasible solutions that improve the best solution

    double      time_for_conversions;           // store the time lost in conversions
    double      time_for_candidate_callback;    // store the time used in the candidate callback
    double      time_for_relaxation_callback;   // store the time used in the relaxation callback

} tsp_statistics;


// SOLVING STRUCTs

extern tsp_environment  tsp_env;        // environment for the problem
extern tsp_instance     tsp_inst;       // instance of the problem
extern tsp_statistics   tsp_stat;       // statistics for the problem

#endif