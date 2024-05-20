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
 * @brief set warm start for cplex solver
 * 
 * @param env cplex env
 * @param lp cplex lp
 * @param path path of solution to be set as warm start (expressed as permutation of nodes)
 * @return int code returned by CPXaddmipstarts
 */
int tsp_cplex_set_mipstarts(CPXENVptr env, CPXLPptr lp, const int* path);

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
 * @brief fix or unfix edges in hard fixing matheuristic (only completely random choice for fixing now)
 * 
 * @param env cplex env
 * @param lp cplex lp
 * @param fix_size number of edges to be fixed
 * @param fixed_edges empty array to be filled with indices of fixed edges
 * @param fix 0 to unfix edges, 1 to fix edges
 */
void tsp_cplex_hard_fixing_manage_edges(CPXENVptr env, CPXLPptr lp, const int fix_size, int* fixed_edges, char fix);

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