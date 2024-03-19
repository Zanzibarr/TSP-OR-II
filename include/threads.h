#ifndef _THREADS_H
#define _THREADS_H

#include "utils.h"

extern pthread_t tsp_threads[N_THREADS];
extern int tsp_available_threads[N_THREADS];

extern pthread_mutex_t tsp_mutex_update_sol;

typedef struct {    //struct used to pass parameters to functions in threads

    int t_index;

    tsp_instance* inst;
    int s_node;
    int* path;
    double* cost;
    int  (*swap_function)(const tsp_instance*, int*, double*);

} tsp_mt_parameters;

void tsp_init_threads();
int tsp_wait_for_thread();
void tsp_free_thread(const int index);
void tsp_wait_all_threads();

#endif