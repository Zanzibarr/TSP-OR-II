#ifndef _THREADS_H
#define _THREADS_H

#include "utils.h"

extern pthread_t tsp_threads[N_THREADS];
extern int tsp_available_threads[N_THREADS];

extern pthread_mutex_t tsp_mutex_update_sol;

/**
 * Struct used to pass parameters to functions in threads
 */
typedef struct {

    int t_index;

    tsp_instance *inst;
    int s_node;
    int *path;
    double *cost;
    int (*swap_function)(const tsp_instance *, int *, double *);

} tsp_mt_parameters;

/**
 * Initialize thread related variables
 */
void tsp_init_threads();

/**
 * Looks for a free thread and waits if none is free
 *
 * @returns The index of the thread (referring to tsp_threads[...]), whenever it can find a free thread
 */
int tsp_wait_for_thread();

/**
 * Frees a thread so that other methods can use it
 *
 * @param index The index of the thread to free
 */
void tsp_free_thread(const int index);

/**
 * Waits for all thread to be free before continuing
 */
void tsp_wait_all_threads();

#endif