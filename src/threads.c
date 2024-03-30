#include "../include/threads.h"

pthread_t tsp_threads[N_THREADS];
int tsp_available_threads[N_THREADS];
pthread_mutex_t tsp_mutex_update_sol;

void tsp_init_threads() {

    #if TSP_VERBOSE == 5
    printf("Initializing threads.\n");
    #endif

    for (int i = 0; i < N_THREADS; i++) tsp_available_threads[i] = 1;

    pthread_mutex_init(&tsp_mutex_update_sol, NULL);

}

int tsp_wait_for_thread() {

    #if TSP_VERBOSE == 5
    int rnd_index = (int)tsp_rnd_coord();
    printf("- %4d - Waiting for thread.\n", rnd_index);
    #endif

    while (1)
        for (int i = 0; i < N_THREADS; i++)
            if (tsp_available_threads[i]) {
                #if TSP_VERBOSE == 5
                printf("- %4d - Thread %d available.\n", rnd_index, i);
                #endif
                tsp_available_threads[i] = 0;
                return i;
            }

}

void tsp_free_thread(const int t_index) {

    #if TSP_VERBOSE == 5
    printf("Freeing thread %d.\n", t_index);
    #endif

    pthread_join(tsp_threads[t_index], NULL);
    tsp_available_threads[t_index] = 1;

}

void tsp_wait_all_threads() {

    #if TSP_VERBOSE == 5
    printf("Waiting for all threads to finish.\n");
    #endif

    int free = 0;
    while (!free) {
        free = 1;
        for (int i = 0; i < N_THREADS; i++)
            if (!tsp_available_threads[i]) free = 0;
    }

    #if TSP_VERBOSE == 5
    printf("All threads finished.\n");
    #endif

}