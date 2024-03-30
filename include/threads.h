#ifndef _THREADS_H
#define _THREADS_H

#include "utils.h"

extern pthread_t tsp_threads[N_THREADS];
extern int tsp_available_threads[N_THREADS];

extern pthread_mutex_t tsp_mutex_update_sol;

/**
 * @brief Struct used to pass parameters to functions in threads
 */
typedef struct {

    int t_index;

    int s_node, e_node;
    int *path;
    double *cost;
    int (*swap_function)(int *, double *);

} tsp_mt_parameters;

/**
 * @brief Initialize thread related variables
 */
void tsp_init_threads();

/**
 * @brief Looks for a free thread and waits if none is free
 *
 * @return The index of the thread (referring to tsp_threads[...]), whenever it can find a free thread
 */
int tsp_wait_for_thread();

/**
 * @brief Frees a thread so that other methods can use it
 *
 * @param index The index of the thread to free
 */
void tsp_free_thread(const int t_index);

/**
 * @brief Waits for all thread to be free before continuing
 */
void tsp_wait_all_threads();

#endif