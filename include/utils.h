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
#define TSP_PARSING_CPLEX_GREEDY_PATCHING   "-patching-greedy"  //parsing cli argument to use greedy patching with cplex
#define TSP_PARSING_CPLEX_CANDIDATE         "-cb-comps"         //parsing cli argument to use the candidate callback in cplex
#define TSP_PARSING_RELAX_CALLBACK          "-cb-fract"         //parsing cli argument to use the relaxation callback in cplex

#define TSP_PARSING_TMP_CHOICE  "-tmp"                  //parsing cli argument to use the temp choice


// FILE NAMES

#define TSP_SOL_FOLDER          "solutions"             // path to the solutions folder
#define TSP_INST_FOLDER         "instances"             // path to the instances folder
#define TSP_SOLUTION_FILE       "solution.txt"     // suffix for the solutions files
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


// STRUCTs

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
     * 0 : OK
     * 1 : time limit reached, found solution
     * 2 : terminated by the user
     * 3 : time limit reached, solution not found
     * 4 : infeasible
    */
    int         status;

    char        file_name[100];             // name of the file where to read the instance (if not random)
    uint64_t    seed;                       // seed used for random algorithms
    char        alg_type[20];               // name of the algorithm using
    double      time_limit;                 // time limit

    char        solution_file[500];

    double      time_start;                 // initial time           
    double      time_total;                 // total execution time

    double      time_for_conversions;       // store the time lost in conversions

    tsp_tabu    tabu_tables[N_THREADS];     // list of tabu tables needed to solve the tabu algorithm
    
    int         tmp_choice;                 //variable used for temporary implementation choices

    int         g2opt_swap_pol;             // swap policy for the g2opt algorithm
    int         g2opt_f2opt;                      // choice for using f2opt algorithm
    int         tabu_tenure;                // tenure for the tabu algorithm
    int         tabu_tenure_a;              // tenure variability for the tabu algorithm
    double      tabu_tenure_f;              // tenure frequency for the tabu algorithm
    int         vns_fvns;                   // choice for fast/normal vns algorithm
    int         cplex_mipstart;             // choice for mipstart in cplex
    int         cplex_benders;              // choice for benders loop in cplex
    int         cplex_patching;             // choice for patching the solution in cplex
    int         cplex_can_cb;               // choice for using the candidate callback in cplex
    int         cplex_rel_cb;               // choice for using the relaxation callback in cplex

} tsp_environment;


// VERBOSE

/**
 * @brief Debugging options
 *
 * <0 for quiet                                 (prints nothing)
 * [0, 10[ for normal                           (basic info for final user)
 * == 5 for thread info                         (multithreading)
 * >=10 for new best solutions                  (visual info)
 * >=50 to plot intermediate costs              (plotting)
 * >=100 for integrity checks                   (integrity checks enabled)      <--     suggested while in development
 * >=500 to see the path in the solution        (advanced debugging)
 * >=1000 for super-verbose                     (full verbose)
 */
extern int tsp_verbose;


// SOLVING STRUCTs

extern tsp_environment  tsp_env;
extern tsp_instance     tsp_inst;

extern int tsp_cplex_terminate;


// USEFUL METHODS

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

#endif
