#include "../include/threads.h"

pthread_t tsp_threads[N_THREADS];
int tsp_available_threads[N_THREADS];
pthread_mutex_t tsp_mutex_update_sol;
pthread_mutex_t tsp_mutex_update_info;

void tsp_init_threads() {

    if (tsp_verbose == 5) print_info("Initializing threads.\n");

    for (int i = 0; i < N_THREADS; i++) tsp_available_threads[i] = 1;

    pthread_mutex_init(&tsp_mutex_update_sol, NULL);
    pthread_mutex_init(&tsp_mutex_update_info, NULL);

}

int tsp_wait_for_thread() {

    int rnd_index = 0;
    if (tsp_verbose == 5) {
        rnd_index = (int)rnd_coord();
        print_info("- %4d - Waiting for thread.\n", rnd_index);
    }

    while (1)
        for (int i = 0; i < N_THREADS; i++)
            if (tsp_available_threads[i]) {
                if (tsp_verbose == 5) print_info("- %4d - Thread %d available.\n", rnd_index, i);
                tsp_available_threads[i] = 0;
                return i;
            }

}

void tsp_free_thread(const int t_index) {

    if (tsp_verbose == 5) print_info("Freeing thread %d.\n", t_index);

    pthread_join(tsp_threads[t_index], NULL);
    tsp_available_threads[t_index] = 1;

}

void tsp_wait_all_threads() {

    if (tsp_verbose == 5) print_info("Waiting for all threads to finish.\n");

    int free = 0;
    while (!free) {
        free = 1;
        for (int i = 0; i < N_THREADS; i++)
            if (!tsp_available_threads[i]) free = 0;
    }

    if (tsp_verbose == 5) print_info("All threads finished.\n");

}