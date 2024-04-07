#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>
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
 * @brief Debugging options
 *
 * <0 for quiet                                 (nothing)
 * [0, 10[ for normal                           (basic info for final user)
 * == 5 for thread info                         (multithreading)
 * >=10 for new best solutions                  (visual info)
 * >=50 to plot intermediate costs              (plotting)
 * >=100 for integrity checks                   (integrity checks enabled)      <--     suggested while in development
 * >=500 to see the path in the solution        (advanced debugging)
 * >=1000 for super-verbose                     (full verbose)
 */
#define TSP_VERBOSE 100

/**
 * @brief Number of threads
 */
#define N_THREADS 16

// PARSING CLI ARGUMENTS

#define TSP_PARSING_FILE        "-file"        //parsing cli argument to select the file
#define TSP_PARSING_TIME_LIMIT  "-tl"    //parsing cli argument to set the time limit
#define TSP_PARSING_SEED        "-seed"        //parsing cli argument to set the seed
#define TSP_PARSING_NNODES      "-nodes"     //parsing cli argument to set the number of nodes
#define TSP_PARSING_HELP        "-help"        //parsing cli argument to ask for cli help
#define TSP_PARSING_ALGORITHM   "-alg"    //parsing cli argument to set the algorithm to use
#define TSP_PARSING_TENURE      "-tenure"        //parsing cli argument to set the tenure to use
#define TSP_PARSING_TENURE_A    "-tenure-a"    //parsing cli argument to set the amplitude parameter for the dinamic tenure
#define TSP_PARSING_TENURE_F    "-tenure-f"    //parsing cli argument to set the frequency parameter for the dinamic tenure
//#define TSP_TEST_RUN "-test"

// DEFAULTS VALUES

#define TSP_DEF_TL      3.6e+6  // number of ms in an hour
#define TSP_DEF_NNODES  300     // default number of nodes
#define TSP_GRID_SIZE   10000   // grid size
#define TSP_EDGE_W_TYPE "ATT"   // default edge weight type
#define TSP_ALG_NUMBER  8       // number of available algorithms

// FILE NAMES

#define TSP_SOL_FOLDER      "solutions"                      // path to the solutions folder
#define TSP_INST_FOLDER     "instances"                     // path to the instances folder
#define TSP_PLOT_FOLDER     "plotting"                      // path to the plotting folder
#define TSP_PLOT_FILE       "solution_plot.png"               // suffix for the plots
#define TSP_SOLUTION_FILE   "solution_file.txt"           // suffix for the solutions files
#define TSP_COORDS_FILE     "coords_file.txt"               // temporary suffix for the plotting
#define TSP_COMMAND_FILE    "command_file.txt"             // temporary suffix for the plotting

/**
// TEST RUN EXCLUSIVE PARAMETERS

#define TSP_TEST_NNODES             5000                     // number of nodes for test instances
#define TSP_TEST_TL                 120                      // timelimit for algorithms during test runs
#define TSP_TEST_NINST              10                       // number of test instances to generate
#define TSP_TEST_RUN_FILES_FOLDER   "perf_prof/test_spec"   // path to the test run files folder
#define TSP_TEST_RUN_RESULTS_FOLDER "perf_prof/results"     // path to the test run results folder
*/

// USEFUL NUMBERS

#define TSP_F2OPT_MAX_DEPTH         6
#define TSP_DEF_TABU_TENURE         80      // tenure base size
#define TSP_EPSILON                 1e-7    // to round double values
#define TSP_CPLEX_ZERO_THRESHOLD    0.5     // threshold used by exact algorithms to determine 0/1 values

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
/**
 * @brief Hyperparameters needed for an algorithm used in a test run (all possible hyperparameters for all possible algorithms containeed; the values that are not needed are to be put to 0)
 
typedef struct {
    int tabu_tenure;
    int tabu_tenure_a;
    double tabu_tenure_f;
} tsp_test_hyperparameters;
*/

/**
 * @brief Problem instance
 */
typedef struct {

    int         nnodes;             // number of nodes
    tsp_pair   *coords;             // list of nodes

    double     *costs;              // cost "matrix"
    int        *sort_edges;         // min edges "matrix":
                                    // row i contains a permutation of the nodes, ordered by increasing distance from node i

    int        *best_solution;      // best solution found so far
                                    // solution is written and read differently depending on whether it was found through heuristics or cplex algorithms
    double      best_cost;          // cost of the best solution found so far
    double      best_time;          // time of the best solution found so far (in seconds)                         

} tsp_instance;

// TIME MANAGEMENT

extern double tsp_initial_time; // "time" at which the algorithm has tarted 
extern double tsp_total_time;   // time in seconds that the algorithm took to conclude
extern double tsp_time_limit;   // time limit for the algorithm

extern int tsp_over_time;           // flag to see whether the algorithm has exceeded the time limit 
extern int tsp_forced_termination;  // flag to see whether the algorithm has been stopped by the user

// SOLVING PARAMETERS

extern char         tsp_algorithms[TSP_ALG_NUMBER][50];  // list of available algorithms

extern uint64_t tsp_seed;                   // seed used for random algorithms
extern tsp_tabu tsp_tabu_tables[N_THREADS]; // list of tabu tables needed to solve the tabu algorithm
extern char     tsp_alg_type[20];           // name of the algorithm using
extern char     tsp_file_name[100];         // name of the file where to read the instance (if not random)

extern int      tsp_tabu_tenure;    // tenure for the tabu algorithm
extern int      tsp_tabu_tenure_a;  // tenure variability for the tabu algorithm
extern double   tsp_tabu_tenure_f;  // tenure frequency for the tabu algorithm

extern char tsp_intermediate_costs_files[N_THREADS][30];    //save the intermediate costs file names

extern tsp_instance tsp_inst;   // Problem instance

// CPLEX VARIABLES

extern CPXENVptr    tsp_cplex_env;              // environment variable for cplex
extern CPXLPptr     tsp_cplex_lp;               // lp variable for cplex
extern double*      tsp_cplex_solution;         // last solution found by cplex (before being attached to tsp_inst)
extern double       tsp_cplex_solution_cost;    // cost of the last solution found by cplex

/**
// TEST RUN EXCLUSIVE PARAMETERS

extern char tsp_test_flag;          // whether the program needs to solve an instance or perform a test run
extern char tsp_test_run_file[100]; // textual file to use for the test run
*/

// USEFUL METHODS


/**
 * @brief Initialize the rand() function to avoid having small random numbers at the beginning
*/
void tsp_init_rand();

/**
 * @brief Get a random value in the [0, TSP_GRID_SIZE] interval
 * 
 * @return The random value in the [0, TSP_GRID_SIZE] interval
*/
double tsp_rnd_coord();

/**
 * @brief Set the initial time of the solving algorithm
 
void tsp_set_initial_time();
*/

/**
 * @brief Get the time elapsed since the beginning of the execution of the code
 * 
 * @return The time (in seconds) elapsed from the beginning of the execution of the code
*/
double tsp_time_elapsed();

#endif
