#ifndef _TSP_H
#define _TSP_H

#include "utils.h"
#include "threads.h"

// PRECOMPUTING
void    tsp_precompute_sort_edges(tsp_instance* inst);
void    tsp_precompute_costs(tsp_instance* inst);

// ALGORITHMS TOOLS
void    tsp_check_best_sol(tsp_instance* inst, int* path, double cost, double time);
void    tsp_reverse(int* path, int start, int end);
int     tsp_check_tabu(int t_index, int from, int to);
void    tsp_add_tabu(int t_index, int from, int to);

// INITIALIZATIONS
void    tsp_init_defs(tsp_instance* inst);
void    tsp_init_solution(tsp_instance* inst);

// SAVING FILES
void    tsp_print_solution(const tsp_instance* inst);
void    tsp_save_solution(const tsp_instance* inst);
void    tsp_plot_solution(const tsp_instance* inst);

// DEBUGGING TOOLS
void    tsp_instance_info(const tsp_instance* inst);
void    tsp_check_integrity(const tsp_instance* inst, double cost, int* path);

// MEMORY MANAGEMENT
void    tsp_allocate_coords_space(tsp_instance* inst);

void    tsp_free_instance(tsp_instance* inst);

#endif