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
 * @param type type of patching to use
 * @param xstar xstar cplex's solution
 * @param ncomp number of components
 * @param comp list containing the component index of each node
 * @param succ successors type list containing the solution found by cplex
 * @param cost cost of the solution (pass a pointer to -1 if you don't have it at hand)
*/
void tsp_cplex_patching(const int type, const double* xstar, int* ncomp, int* comp, int* succ, double* cost);

/**
 * @brief cplex callback for candidate solution
 * 
 * @param context cplex context
 * @param userhandle userhandle (NULL)
 * 
 * @return cplex error code
*/
int tsp_cplex_callback_candidate(CPXCALLBACKCONTEXTptr context, const void* userhandle);

/**
 * @brief cplex callback for relaxation solution
 * 
 * @param context cplex context
 * @param userhandle userhandle (NULL)
 * 
 * @return cplex error code
*/
int tsp_cplex_callback_relaxation(CPXCALLBACKCONTEXTptr context, const void* userhandle);

/**
 * @brief fix edges in hard fixing matheuristic (only completely random choice for fixing now)
 * 
 * @param env cplex env
 * @param lp cplex lp
 * @param fix_size number of edges to be fixed
 * @param fixed_edges empty array to be filled with indices of fixed edges
 */
void tsp_cplex_dive_fix(CPXENVptr env, CPXLPptr lp, const int fix_size, int* fixed_edges);

/**
 * @brief unfix edges in hard fixing matheuristic (only completely random choice for fixing now)
 * 
 * @param env cplex env
 * @param lp cplex lp
 * @param fix_size number of edges to be fixed
 * @param fixed_edges empty array to be filled with indices of fixed edges
 */
void tsp_cplex_dive_unfix(CPXENVptr env, CPXLPptr lp, const int fix_size, int* fixed_edges);

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