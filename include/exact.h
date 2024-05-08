#ifndef _EXACT_H
#define _EXACT_H

#include "tsp.h"
#include "mincut.h"

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
 * 
 * @param xstar xstar cplex's solution
 * 
 * @return the computed cost
 */
double tsp_cplex_compute_xstar_cost(const double* xstar);

/**
 * @brief Convert a path type solution to a cplex type solution
 * 
 * @param ncols Number of columns (for cplex)
 * @param path path type solution to convert
 * @param indexes indexes type solution (for cplex)
 * @param values values type solution (for cplex)
*/
void tsp_cplex_path_to_ind_val(const int ncols, const int* path, int* indexes, double* values);

/**
 * @brief (THREAD SAFE) Checks and updates the incumbent of the instance
 * Stores also a multitour solution
 * 
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 * @param cost the cost of the solution to check
*/
void tsp_cplex_check_best_sol(const int ncomp, const int* comp, const int* succ, const double cost);

/**
 * @brief Decompose xstar into comp and succ
 * 
 * @param xstar xstar cplex's solution
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 * @param ncomp number of components
*/
void tsp_cplex_decompose_xstar(const double* xstar, int* comp, int* succ, int* ncomp);

/**
 * @brief add a SEC to the cplex model
 * 
 * @param env cplex environment
 * @param lp cplex lp
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 */
void tsp_cplex_add_sec(CPXENVptr env, CPXLPptr lp, const int* ncomp, const int* comp, const int* succ);

/**
 * @brief apply patching to the cplex solution
 * 
 * @param xstar xstar cplex's solution
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 * @param cost cost of the solution
*/
void tsp_cplex_patching(const double* xstar, int* ncomp, int* comp, int* succ, double* cost);

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
 * @param xstar xstar cplex's solution
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
*/
void tsp_cplex_close(CPXENVptr env, CPXLPptr lp, double* xstar, int* comp, int* succ);

#endif